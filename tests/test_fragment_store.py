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
    # A system that is not a CBH species
    system = db.create_system(name="water", make_current=True)
    system.create_configuration(name="initial", make_current=True).from_smiles("O")
    db.close()
    (job_dir / "job_data.json").write_text(json.dumps({"job id": 42}))
    return job_dir


def test_harvest(job, store):
    counts = harvest_job(job, store)
    assert counts["inserted"] == 3 and counts["skipped"] == []
    records = {r["species"]: r for r in store.records()}
    assert set(records) == {"CCF", "C", "[CH3]"}
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
    assert counts == {"inserted": 0, "duplicate": 3, "logged": 0, "skipped": []}
    # and the step's lookup finds them
    assert store.get("CCF", model_chemistry="MOPAC:SQM@PM7")["E"] == pytest.approx(
        -262.0
    )
    assert store.get("C", energy_fingerprint="mopac:deadbeef")["ZPE"] == 115.0


def test_harvest_missing_job(store, tmp_path):
    with pytest.raises(FileNotFoundError):
        harvest_job(tmp_path, store)


def test_cli_roundtrip(job, tmp_path, capsys):
    from seamm_thermochemistry.fragment_cli import main

    path = str(tmp_path / "cli.db")
    assert main(["--store", path, "harvest", str(job)]) == 0
    out = capsys.readouterr().out
    assert "3 inserted" in out
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
