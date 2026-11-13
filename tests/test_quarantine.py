"""Rows quarantined by hand (model chemistry "UNUSABLE: <reason>", unverified)
are never replaced automatically, only by an explicit supersede with a reason
(2026-11-13, after a trial harvest replaced the quarantined [CH3-] saddle); they
are not conformers, and the store logs what each harvest or write-back did."""

import pytest

pytest.importorskip("rdkit")
molsystem = pytest.importorskip("molsystem")

from seamm_thermochemistry import (  # noqa: E402
    EnergyRecord,
    FragmentRecord,
    FragmentStore,
    GeometryRecord,
    harvest_configuration,
)
from seamm_thermochemistry import fragment_cli  # noqa: E402

from tests.test_writeback import _embedded, _species  # noqa: E402

FP = "fp:r2scan"
REASON = "UNUSABLE: saddle point (687i), job 5684"


@pytest.fixture()
def store(tmp_path):
    with FragmentStore(tmp_path / "fragments.db", create=True) as s:
        yield s


def _quarantine(store, table, where="1 = 1"):
    db = store._connect(create=True)
    db.execute(
        f"UPDATE {table} SET verified = 0, model_chemistry = ? WHERE {where}",
        (REASON,),
    )
    db.commit()


def test_a_put_never_replaces_a_quarantined_fragment_record(store):
    record = FragmentRecord.make("[CH3-]", charge=-1, E=-1.0, energy_fingerprint="e")
    store.put(record)
    _quarantine(store, "fragment")
    new = FragmentRecord.make(
        "[CH3-]", charge=-1, E=-2.0, energy_fingerprint="e", model_chemistry="M@X"
    )
    assert store.put(new) == "quarantined"
    assert store.get("[CH3-]", charge=-1, energy_fingerprint="e") is None
    db = store._connect()
    row = db.execute("SELECT E, model_chemistry FROM fragment").fetchone()
    assert row["E"] == -1.0 and row["model_chemistry"] == REASON
    (h,) = store.history("fragment")
    assert h["kind"] == "quarantined" and "kept the quarantine: saddle" in h["reason"]
    assert store.superseded("[CH3-]") == []
    # A person releases it, with a reason
    assert store.supersede(new, "recomputed at the D3h minimum") == "superseded"
    got = store.get("[CH3-]", charge=-1, energy_fingerprint="e")
    assert got["E"] == -2.0 and got["verified"] and got["model_chemistry"] == "M@X"
    (old,) = store.superseded("[CH3-]")
    assert old["E"] == -1.0
    assert "recomputed" in old["reason"] and "released the quarantine: saddle" in (
        old["reason"]
    )
    assert store.put(new) == "duplicate"


def test_a_put_never_replaces_a_quarantined_geometry_or_energy_row(store):
    symbols, xyz = _embedded("C")
    g = GeometryRecord.make(
        "C", symbols, xyz, E=-1.0, ZPE=1.0, H=2.0, frequency_fingerprint=FP
    )
    store.put_geometry(g)
    e = EnergyRecord.make(g["geometry_hash"], -5.0, energy_fingerprint="e")
    store.put_energy(e)
    _quarantine(store, "geometry")
    _quarantine(store, "energy")
    assert store.geometries("C", frequency_fingerprint=FP) == []
    outcome, row = store.put_geometry(GeometryRecord(g, E=-1.5))
    assert outcome == "quarantined"
    assert store.geometries("C", frequency_fingerprint=FP) == []
    (h,) = store.history("geometry")
    assert h["kind"] == "quarantined" and "saddle" in h["reason"]
    new_energy = EnergyRecord.make(g["geometry_hash"], -6.0, energy_fingerprint="e")
    assert store.put_energy(new_energy)[0] == "quarantined"
    assert store.energy(g["geometry_hash"], energy_fingerprint="e") is None
    assert store.superseded_rows("C") == []
    # Released explicitly, with a reason
    outcome, _ = store.supersede_geometry(GeometryRecord(g, E=-1.5), "checked by hand")
    assert outcome == "superseded"
    assert store.supersede_energy(new_energy, "checked by hand")[0] == "superseded"
    (got,) = store.geometries("C", frequency_fingerprint=FP)
    assert got["E"] == -1.5 and got["verified"]
    assert store.energy(g["geometry_hash"], energy_fingerprint="e")["E"] == -6.0
    rows = store.superseded_rows("C")
    assert {d["table"] for d in rows} == {"geometry", "energy"}
    geometry = next(d for d in rows if d["table"] == "geometry")
    assert geometry["old"]["E"] == -1.0 and geometry["new"]["E"] == -1.5
    assert geometry["old"]["model_chemistry"].startswith("UNUSABLE:")
    assert "released the quarantine: saddle" in geometry["reason"]
    assert store.superseded_rows("CC") == []


def test_a_quarantined_row_cannot_supersede(store):
    record = FragmentRecord.make("C", E=-1.0, energy_fingerprint="e")
    store.put(record)
    bad = FragmentRecord.make(
        "C", E=-2.0, energy_fingerprint="e", model_chemistry=REASON
    )
    with pytest.raises(ValueError):
        store.supersede(bad, "no")


def test_the_harvest_keeps_the_ch3_minus_quarantine(store, tmp_path):
    """The 11-13 trial harvest of ChemAI 5684/5685: [CH3-], quarantined on 11-11
    as the D3h 687i saddle, came back from a job that had it as the saddle again
    and replaced the quarantine. A harvest now refuses the saddle, and a valid
    row with the key leaves the quarantine in place, reported."""
    db = molsystem.SystemDB(filename=f"file:{tmp_path / 'seamm.db'}")
    xyz = _embedded("[CH3-]")[1]

    def harvest(frequencies, **kwargs):
        return harvest_configuration(
            store,
            _species(db, "[CH3-]", xyz, -104000.0, frequencies=frequencies),
            "[CH3-]",
            energy="total energy#Thermochemistry#*",
            job="5684:/jobs/5684",
            **kwargs,
        )

    counts = harvest((300.0, 1500.0))  # stored first, as it was before 11-11
    assert counts["inserted"] == 1 and counts["geometry"] == {"inserted": 1}
    for table in ("fragment", "geometry", "energy"):
        _quarantine(store, table)

    # The saddle, by an ordinary harvest (not a write-back): refused
    counts = harvest((-687.0, 300.0, 1500.0))
    assert counts["inserted"] == counts["superseded"] == counts["quarantined"] == 0
    assert counts["geometry"] == {} and counts["energy"] == {}
    assert "imaginary mode of 687i" in counts["skipped"][0]
    # ... also with an explicit supersede: a release must pass every check
    counts = harvest((-687.0, 300.0, 1500.0), supersede="recomputed")
    assert counts["superseded"] == 0 and counts["geometry"] == {}

    # A valid row for the key: the quarantine stays, and the attempt is reported
    counts = harvest((310.0, 1500.0))
    assert counts["quarantined"] == 1 and counts["superseded"] == 0
    assert counts["geometry"] == {"quarantined": 1}
    assert counts["energy"] == {"quarantined": 1}
    assert any("quarantined" in s for s in counts["skipped"])
    rows = store._connect().execute("SELECT model_chemistry FROM fragment").fetchall()
    assert [r[0] for r in rows] == [REASON]
    assert store.geometries("[CH3-]", charge=-1, frequency_fingerprint=FP) == []
    assert store.superseded("[CH3-]") == [] and store.superseded_rows("[CH3-]") == []
    outcomes = {r["outcome"] for r in store.write_log(job="5684")}
    assert {"refused", "quarantined"} <= outcomes
    db.close()


def test_quarantined_rows_are_not_conformers(store, tmp_path):
    """The dissociated Cl2O2 conformer was the lowest: quarantined, it must not
    make a valid new conformer look 'not lower'."""
    db = molsystem.SystemDB(filename=f"file:{tmp_path / 'seamm.db'}")
    gauche = _embedded("CCCC", 65.0)[1]
    anti = _embedded("CCCC", 180.0)[1]
    eclipsed = _embedded("CCCC", 120.0)[1]

    def write(xyz, E):
        return harvest_configuration(
            store,
            _species(db, "CCCC", xyz, E),
            "CCCC",
            energy="total energy#Thermochemistry#*",
            write_back=True,
            job="1",
        )

    write(gauche, -500000.0)
    write(anti, -500010.0)  # lower and distinct: auto:
    _quarantine(store, "geometry", "conformer LIKE 'auto:%'")
    # Lower than the valid gauche one, higher than the quarantined anti one
    counts = write(eclipsed, -500002.0)
    assert counts["geometry"] == {"inserted": 1}
    db.close()


def test_the_write_log(store, tmp_path, capsys):
    db = molsystem.SystemDB(filename=f"file:{tmp_path / 'seamm.db'}")
    good = _species(db, "CO", _embedded("CO")[1], -300000.0)
    saddle = _species(db, "CC", _embedded("CC")[1], -400000.0, frequencies=(-90.0,))
    for configuration, smiles in ((good, "CO"), (saddle, "CC")):
        harvest_configuration(
            store,
            configuration,
            smiles,
            energy="total energy#Thermochemistry#*",
            write_back=True,
            job="5700:/jobs/5700",
        )
    db.close()
    rows = store.write_log(job="5700")
    outcomes = {(r["species"], r["table_name"], r["outcome"]) for r in rows}
    assert ("CO", "geometry", "inserted") in outcomes
    assert ("CO", "energy", "inserted") in outcomes
    assert ("CO", "fragment", "inserted") in outcomes
    (refused,) = [r for r in rows if r["outcome"] == "refused"]
    assert refused["species"] == "CC" and "90i" in refused["reason"]
    assert all(r["kind"] == "write-back" for r in rows)
    assert store.write_log(species="CO") and store.write_log(job="9999") == []
    code = fragment_cli.main(
        ["--store", str(store.path), "writeback-log", "--job", "5700"]
    )
    out = capsys.readouterr().out
    assert code == 0 and "refused: an imaginary mode of 90i" in out


def test_cli_superseded_shows_geometry_rows(store, capsys):
    symbols, xyz = _embedded("C")
    g = GeometryRecord.make(
        "C", symbols, xyz, E=-1.0, ZPE=1.0, H=2.0, frequency_fingerprint=FP
    )
    store.put_geometry(g)
    _quarantine(store, "geometry")
    store.supersede_geometry(GeometryRecord(g, E=-1.5), "checked by hand")
    assert fragment_cli.main(["--store", str(store.path), "superseded", "C"]) == 0
    out = capsys.readouterr().out
    assert "C geometry row 1" in out and "UNUSABLE: saddle" in out
    assert "released the quarantine" in out


def test_release_candidates_and_release(store, tmp_path, capsys):
    """A recomputed value for a quarantined species that passed the harvest's
    checks is listed as a release candidate; one reasoned command releases it."""
    db = molsystem.SystemDB(filename=f"file:{tmp_path / 'seamm.db'}")
    xyz = _embedded("[CH3-]")[1]

    def harvest(E):
        return harvest_configuration(
            store,
            _species(db, "[CH3-]", xyz, E),
            "[CH3-]",
            energy="total energy#Thermochemistry#*",
            job="5684:/jobs/5684",
            write_back=True,
        )

    harvest(-104000.0)
    for table in ("fragment", "geometry", "energy"):
        _quarantine(store, table)
    assert store.release_candidates() == []
    harvest(-104010.0)
    candidates = store.release_candidates("[CH3-]")
    assert {c["table"] for c in candidates} == {"fragment", "geometry", "energy"}
    for c in candidates:
        assert c["old"]["E"] == pytest.approx(-104000.0)
        assert c["new"]["E"] == pytest.approx(-104010.0)
        assert c["quarantine"].startswith("saddle") and c["job"] == "5684:/jobs/5684"
    assert store.release_candidates("C") == []

    path = str(store.path)
    store.close()
    assert fragment_cli.main(["--store", path, "release-candidates"]) == 0
    out = capsys.readouterr().out
    assert out.count("quarantined (saddle point (687i), job 5684)") == 3
    with pytest.raises(SystemExit):  # a reason is required
        fragment_cli.main(["--store", path, "release", str(candidates[0]["id"])])
    ids = [str(c["id"]) for c in candidates]
    code = fragment_cli.main(
        ["--store", path, "release", *ids, "--reason", "recomputed, checked"]
    )
    assert code == 0
    with FragmentStore(path) as s:
        assert s.release_candidates() == []
        got = s.get("[CH3-]", charge=-1, energy_fingerprint=FP)
        assert got["E"] == pytest.approx(-104010.0) and got["verified"]
        (g,) = s.geometries("[CH3-]", charge=-1, frequency_fingerprint=FP)
        assert g["E"] == pytest.approx(-104010.0)
        (old,) = s.superseded("[CH3-]")
        assert "recomputed, checked (release candidate" in old["reason"]
        assert "released the quarantine: saddle" in old["reason"]
        with pytest.raises(ValueError):
            s.release(candidates[0]["id"], "again")
    db.close()
