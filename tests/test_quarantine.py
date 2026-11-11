"""Rows quarantined by hand (model chemistry "UNUSABLE: <reason>", unverified)
give way to valid new rows with the same key, are not conformers, and the store
logs what each harvest or write-back did (P5b follow-up, 2026-11-11)."""

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


def test_a_valid_fragment_record_replaces_a_quarantined_one(store):
    record = FragmentRecord.make("[CH3-]", charge=-1, E=-1.0, energy_fingerprint="e")
    store.put(record)
    _quarantine(store, "fragment")
    new = FragmentRecord.make(
        "[CH3-]", charge=-1, E=-2.0, energy_fingerprint="e", model_chemistry="M@X"
    )
    assert store.put(new) == "superseded"
    got = store.get("[CH3-]", charge=-1, energy_fingerprint="e")
    assert got["E"] == -2.0 and got["verified"] and got["model_chemistry"] == "M@X"
    (old,) = store.superseded("[CH3-]")
    assert old["E"] == -1.0 and "replaced a quarantined record: saddle" in old["reason"]
    # A second, valid put is an ordinary duplicate
    assert store.put(new) == "duplicate"


def test_a_valid_geometry_or_energy_row_replaces_a_quarantined_one(store):
    symbols, xyz = _embedded("C")
    g = GeometryRecord.make(
        "C", symbols, xyz, E=-1.0, ZPE=1.0, H=2.0, frequency_fingerprint=FP
    )
    store.put_geometry(g)
    store.put_energy(
        EnergyRecord.make(g["geometry_hash"], -5.0, energy_fingerprint="e")
    )
    _quarantine(store, "geometry")
    _quarantine(store, "energy")
    assert store.geometries("C", frequency_fingerprint=FP) == []
    outcome, row = store.put_geometry(GeometryRecord(g, E=-1.5))
    assert outcome == "superseded"
    (got,) = store.geometries("C", frequency_fingerprint=FP)
    assert got["E"] == -1.5 and got["verified"]
    (h,) = store.history("geometry")
    assert h["kind"] == "superseded" and "quarantined row: saddle" in h["reason"]
    outcome, _ = store.put_energy(
        EnergyRecord.make(g["geometry_hash"], -6.0, energy_fingerprint="e")
    )
    assert outcome == "superseded"
    assert store.energy(g["geometry_hash"], energy_fingerprint="e")["E"] == -6.0


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
