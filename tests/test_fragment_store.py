"""The fragment-energy store: keys, insert-or-ignore with the duplicate log,
verified vs unverified records, best-effort reads, and harvesting a job's
CBH:<SMILES> systems."""

import json
import sqlite3

import pytest

from seamm_thermochemistry import FragmentRecord, FragmentStore, harvest_job
from seamm_thermochemistry import fragment_store as F


@pytest.fixture()
def store(tmp_path):
    with FragmentStore(tmp_path / "fragments.db", create=True) as s:
        yield s


def _record(species="CC", E=-100.0, fp="orca:abc123", **kw):
    kw.setdefault("model_chemistry", "ORCA:DFT@B3LYP/bse:def2-SVP")
    return FragmentRecord.make(species, E=E, energy_fingerprint=fp, **kw)


def test_insert_duplicate_and_log(store, caplog):
    r = _record(ZPE=10.0, H=-90.0, frequency_fingerprint="orca:freq1", job="1")
    assert store.put(r) == "inserted"
    assert len(store) == 1
    # Same key, same values (within 0.1): ignored
    r2 = _record(E=-100.05, ZPE=10.0, H=-90.0, frequency_fingerprint="orca:freq1")
    assert store.put(r2) == "duplicate"
    # Same key, a different value: logged, stored value kept
    r3 = _record(
        E=-101.0, ZPE=10.0, H=-90.0, frequency_fingerprint="orca:freq1", job="2"
    )
    assert store.put(r3) == "logged"
    assert len(store) == 1
    got = store.get("CC", energy_fingerprint="orca:abc123")
    assert got["E"] == -100.0 and got["verified"] and got["job"] == "1"
    dups = store.duplicates()
    assert len(dups) == 1 and dups[0]["difference"] == pytest.approx(1.0)
    assert dups[0]["job"] == "2"
    assert "differing by 1.000" in caplog.text


def test_key_includes_charge_multiplicity_and_frequency_level(store):
    store.put(_record("[CH2]C", multiplicity=2, E=1.0))
    store.put(_record("[CH2]C", multiplicity=2, charge=1, E=2.0))
    store.put(
        _record("[CH2]C", multiplicity=2, E=3.0, ZPE=5.0, frequency_fingerprint="f")
    )
    assert len(store) == 3
    assert (
        store.get("[CH2]C", multiplicity=2, energy_fingerprint="orca:abc123")["E"]
        == 1.0
    )
    assert (
        store.get("[CH2]C", multiplicity=2, charge=1, energy_fingerprint="orca:abc123")[
            "E"
        ]
        == 2.0
    )
    r = store.get(
        "[CH2]C",
        multiplicity=2,
        energy_fingerprint="orca:abc123",
        frequency_fingerprint="f",
    )
    assert r["E"] == 3.0 and r["ZPE"] == 5.0
    assert store.get("[CH2]C", multiplicity=1, energy_fingerprint="orca:abc123") is None


def test_unverified_fallback(store):
    mc = "MOPAC:SQM@PM7"
    u = FragmentRecord.make("C", E=-70.0, model_chemistry=mc, version="23.0.1")
    assert not u["verified"] and u["energy_fingerprint"].startswith("unverified:")
    store.put(u)
    # No fingerprint given: the unverified record by model chemistry
    r = store.get("C", model_chemistry=mc)
    assert r is not None and not r["verified"] and r["E"] == -70.0
    assert store.get("C", model_chemistry=mc, version="22.0") is None
    # A fingerprint given but nothing verified: still the unverified one
    r = store.get("C", energy_fingerprint="mopac:xyz", model_chemistry=mc)
    assert not r["verified"]
    # Once a verified record exists it wins
    store.put(_record("C", E=-71.0, fp="mopac:xyz", model_chemistry=mc))
    r = store.get("C", energy_fingerprint="mopac:xyz", model_chemistry=mc)
    assert r["verified"] and r["E"] == -71.0
    with pytest.raises(ValueError, match="fingerprint or a model chemistry"):
        FragmentRecord.make("C", E=1.0)


def test_best_effort_reads(tmp_path, caplog):
    # No file: an empty store, no error
    s = FragmentStore(tmp_path / "none.db")
    assert len(s) == 0 and s.get("C", model_chemistry="x") is None
    assert s.records() == [] and s.duplicates() == []
    with pytest.raises(PermissionError):
        s.put(_record())
    # A corrupt file: None and a warning
    bad = tmp_path / "bad.db"
    bad.write_text("not a database")
    s = FragmentStore(bad)
    assert s.get("C", model_chemistry="x") is None
    assert "could not be read" in caplog.text


def test_connection_pragmas(store):
    db = store._connect()
    assert db.execute("PRAGMA mmap_size").fetchone()[0] == 0
    assert db.execute("PRAGMA journal_mode").fetchone()[0].lower() == "delete"
    assert db.execute("PRAGMA busy_timeout").fetchone()[0] == 30000


def test_readonly_connection_cannot_write(store, tmp_path):
    store.put(_record())
    store.close()
    reader = FragmentStore(store.path)
    assert len(reader) == 1
    with pytest.raises(sqlite3.OperationalError):
        reader._connect().execute("DELETE FROM fragment")


# -- harvesting --------------------------------------------------------------


@pytest.fixture()
def job(tmp_path):
    """A job directory with CBH:<SMILES> systems carrying energies."""
    molsystem = pytest.importorskip("molsystem")
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    db = molsystem.SystemDB(filename=f"file:{job_dir / 'seamm.db'}")
    model = "MOPAC:SQM@PM7"
    for smiles, E in (("CCF", -262.0), ("C", -50.0), ("[CH3]", 150.0)):
        system = db.create_system(name=f"CBH:{smiles}", make_current=True)
        c = system.create_configuration(name="initial", make_current=True)
        c.from_smiles(smiles)
        p = c.properties
        name = f"total energy#Thermochemistry#{model}"
        p.add(name, "float", units="E_h", noerror=True)
        p.put(name, E / 2625.4996394799)
        if smiles == "C":
            for n, v in (
                (f"ZPE#Thermochemistry#{model}", 115.0),
                (f"H(298.15 K)#Thermochemistry#{model}", -40.0),
            ):
                p.add(n, "float", units="kJ/mol", noerror=True)
                p.put(n, v)
            n = f"settings fingerprint#Thermochemistry#{model}"
            p.add(n, "str", units=None, noerror=True)
            p.put(n, "mopac:deadbeef")
    # A species named with a non-canonical SMILES (as Read Structure names them)
    system = db.create_system(name="CBH:C(F)", make_current=True)
    c = system.create_configuration(name="initial", make_current=True)
    c.from_smiles("CF")
    c.properties.add(
        f"total energy#Thermochemistry#{model}", "float", units="E_h", noerror=True
    )
    c.properties.put(f"total energy#Thermochemistry#{model}", -230.0 / 2625.4996394799)
    # A system that is not a CBH species (named by SMILES, as Read Structure
    # does): skipped unless all systems are harvested
    system = db.create_system(name="O", make_current=True)
    system.create_configuration(name="initial", make_current=True).from_smiles("O")
    db.close()
    (job_dir / "job_data.json").write_text(json.dumps({"job id": 42}))
    return job_dir


def test_harvest(job, store):
    counts = harvest_job(job, store)
    assert counts["inserted"] == 4 and counts["skipped"] == []
    records = {r["species"]: r for r in store.records()}
    assert set(records) == {"CCF", "C", "[CH3]", "CF"}  # C(F) canonicalized
    assert records["CCF"]["E"] == pytest.approx(-262.0)  # E_h converted
    assert not records["CCF"]["verified"]  # no fingerprint stored
    assert records["CCF"]["model_chemistry"] == "MOPAC:SQM@PM7"
    assert records["CCF"]["job"].startswith("42:")
    assert records["[CH3]"]["multiplicity"] == 2
    c = records["C"]
    assert c["verified"] and c["energy_fingerprint"] == "mopac:deadbeef"
    assert c["ZPE"] == 115.0 and c["H"] == -40.0 and c["T"] == 298.15
    assert c["geometry_hash"] and c["source"].startswith("total energy#")
    # Harvesting again changes nothing
    counts = harvest_job(job, store)
    assert counts == {
        "inserted": 0,
        "duplicate": 4,
        "logged": 0,
        "superseded": 0,
        "skipped": [],
    }
    # and the step's lookup finds them
    assert store.get("CCF", model_chemistry="MOPAC:SQM@PM7")["E"] == pytest.approx(
        -262.0
    )
    assert store.get("C", energy_fingerprint="mopac:deadbeef")["ZPE"] == 115.0


def test_harvest_all_systems(job, store):
    counts = harvest_job(job, store, all_systems=True)
    assert counts["inserted"] == 4
    assert "O: no property matching" in counts["skipped"][0]


def test_harvest_missing_job(store, tmp_path):
    with pytest.raises(FileNotFoundError):
        harvest_job(tmp_path, store)


def test_cli_roundtrip(job, tmp_path, capsys):
    from seamm_thermochemistry.fragment_cli import main

    path = str(tmp_path / "cli.db")
    assert main(["--store", path, "harvest", str(job)]) == 0
    out = capsys.readouterr().out
    assert "4 inserted" in out
    assert (
        main(["--store", path, "lookup", "C", "--fingerprint", "mopac:deadbeef"]) == 0
    )
    assert "ZPE = 115.000" in capsys.readouterr().out
    assert main(["--store", path, "lookup", "CCCC"]) == 1
    assert main(["--store", path, "list"]) == 0
    assert "[unverified]" in capsys.readouterr().out
    assert main(["--store", path, "duplicates"]) == 0


def test_geometry_hash():
    h1 = F.geometry_hash(["O", "H"], [[0, 0, 0], [0.96, 0, 0]])
    h2 = F.geometry_hash(["O", "H"], [[0, 0, 0], [0.96004, 0, 0]])
    h3 = F.geometry_hash(["O", "H"], [[0, 0, 0], [0.97, 0, 0]])
    assert h1 == h2 != h3


def test_supersede_keeps_the_old_values(store):
    """A value known to be wrong is replaced only explicitly; the old values
    stay in the superseded table with the date and the reason."""
    old = _record(ZPE=10.0, H=-90.0, frequency_fingerprint="orca:freq1", job="1")
    assert store.put(old) == "inserted"
    better = _record(
        E=-103.0, ZPE=10.2, H=-92.0, frequency_fingerprint="orca:freq1", job="7"
    )
    # Without supersede a differing value is only logged
    assert store.put(better) == "logged"
    assert store.get("CC", energy_fingerprint="orca:abc123")["E"] == -100.0
    with pytest.raises(ValueError):
        store.supersede(better, "")
    assert store.supersede(better, "lowest conformer from B2PLYP geometries") == (
        "superseded"
    )
    assert len(store) == 1  # same key, same record
    got = store.get("CC", energy_fingerprint="orca:abc123")
    assert got["E"] == -103.0 and got["job"] == "7"
    history = store.superseded("CC")
    assert len(history) == 1
    h = history[0]
    assert h["E"] == -100.0 and h["job"] == "1" and h["current_E"] == -103.0
    assert h["reason"] == "lowest conformer from B2PLYP geometries"
    assert h["superseded_on"]
    # The same values again: nothing to replace
    assert store.supersede(better, "again") == "unchanged"
    assert len(store.superseded()) == 1
    # A new key is just inserted
    assert store.supersede(_record("CCC", E=-150.0), "new") == "inserted"
    assert store.superseded("CCC") == []


def test_supersede_keeps_values_the_new_record_lacks(store):
    """A corrected energy without a ZPE or H keeps the stored ones; a plain
    dict is normalized like FragmentRecord.make."""
    store.put(_record(ZPE=10.0, H=-90.0, frequency_fingerprint="orca:freq1"))
    new = _record(E=-105.0, ZPE=10.0, H=-90.0, frequency_fingerprint="orca:freq1")
    new["ZPE"] = None
    new["H"] = None
    assert store.supersede(new, "corrected energy only") == "superseded"
    got = store.get("CC", energy_fingerprint="orca:abc123")
    assert got["E"] == -105.0 and got["ZPE"] == 10.0 and got["H"] == -90.0
    # A plain dict with a None frequency fingerprint finds the same key
    plain = {
        "species": "CCC",
        "charge": 0,
        "multiplicity": 1,
        "energy_fingerprint": "orca:abc123",
        "frequency_fingerprint": None,
        "E": -1.0,
    }
    store.put(_record("CCC", E=-2.0))
    assert store.supersede(plain, "plain dict") == "superseded"
    assert len(store) == 2
    with pytest.raises(ValueError):
        store.supersede({"species": "C", "E": 1.0}, "no fingerprint")


def test_superseded_on_an_old_store(tmp_path):
    """A store made before the superseded table existed reads as nothing
    superseded, and gains the table when opened for writing."""
    path = tmp_path / "old.db"
    db = sqlite3.connect(path)
    db.executescript(F._SCHEMA.split("CREATE TABLE IF NOT EXISTS superseded")[0])
    db.commit()
    db.close()
    assert FragmentStore(path).superseded() == []
    with FragmentStore(path, create=True) as s:
        s.put(_record(E=-1.0))
        assert s.supersede(_record(E=-2.0), "test") == "superseded"
        assert len(s.superseded()) == 1


def test_harvest_supersede(job, store):
    harvest_job(job, store)
    # A corrected job: CCF's energy changed
    from molsystem import SystemDB

    db = SystemDB(filename=str(job / "seamm.db"))
    for system in db.systems:
        if system.name == "CBH:CCF":
            c = system.configuration
            name = next(k for k in c.properties.get("total energy*").keys())
            c.properties.put(name, -263.0 / 2625.4996394799)
    db.close()
    counts = harvest_job(job, store)
    assert counts["logged"] == 1 and counts["superseded"] == 0
    counts = harvest_job(job, store, supersede="corrected conformer")
    assert counts["superseded"] == 1 and counts["duplicate"] == 3
    assert store.get("CCF", model_chemistry="MOPAC:SQM@PM7")["E"] == pytest.approx(
        -263.0
    )
    assert store.superseded("CCF")[0]["E"] == pytest.approx(-262.0)


def test_cli_supersede(job, tmp_path, capsys):
    from seamm_thermochemistry.fragment_cli import main

    path = str(tmp_path / "cli.db")
    assert main(["--store", path, "harvest", str(job)]) == 0
    assert main(["--store", path, "superseded"]) == 0
    assert "Nothing superseded." in capsys.readouterr().out
    assert main(["--store", path, "harvest", str(job), "--supersede", "why"]) == 0
    assert "0 superseded" in capsys.readouterr().out
