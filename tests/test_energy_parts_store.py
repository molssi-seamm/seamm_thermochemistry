"""Energy parts (composite phase 2 §2) in the fragment store: harvest and
write-back keep a configuration's 'energy parts' property as the energy row's
byproducts, and a duplicate that brings parts adds them to a row without any."""

import pytest

pytest.importorskip("rdkit")
molsystem = pytest.importorskip("molsystem")

from seamm_thermochemistry.geometry_store import EnergyRecord  # noqa: E402

from tests.test_writeback import (  # noqa: E402,F401
    KJ,
    LEVEL,
    _embedded,
    _species,
    _write,
    db,
    store,
)

PARTS = {
    "parts version": 1,
    "units": "E_h",
    "parts": {"HF": -76.05, "CCSD correlation": -0.26, "(T)": -0.008},
    "total": ["HF", "CCSD correlation", "(T)"],
    "factors": {},
    "reference": "RHF",
    "frozen core orbitals": 1,
    "relativistic": "none",
    "diagnostics": {"T1": 0.007},
    "sources": {"HF": "property file"},
}


def _energy_rows(store):  # noqa: F811
    rows = store._connect().execute("SELECT * FROM energy").fetchall()
    return [dict(r) for r in rows]


def test_harvest_keeps_the_parts(db, store):  # noqa: F811
    symbols, xyz = _embedded("CCO")
    c = _species(db, "CCO", xyz, -406500.0)
    name = f"energy parts#Thermochemistry#{LEVEL}"
    c.properties.add(name, "json", noerror=True)
    c.properties.put(name, PARTS)
    _write(store, c, "CCO")
    (row,) = _energy_rows(store)
    import json

    assert json.loads(row["byproducts"]) == PARTS


def test_a_duplicate_with_parts_adds_them(store):  # noqa: F811
    plain = EnergyRecord.make("g:1", -1000.0, energy_fingerprint="fp:x")
    assert store.put_energy(plain)[0] == "inserted"
    with_parts = EnergyRecord.make(
        "g:1", -1000.0 + 1e-9, energy_fingerprint="fp:x", byproducts=PARTS
    )
    outcome, row_id = store.put_energy(with_parts)
    assert outcome == "parts added"
    row = store.energy("g:1", energy_fingerprint=["fp:x"])
    assert row["byproducts"] == PARTS
    # Once there, a further duplicate changes nothing
    assert store.put_energy(with_parts)[0] == "duplicate"
    history = store._connect().execute(
        "SELECT kind FROM row_history WHERE table_name = 'energy'"
    )
    assert [r[0] for r in history.fetchall()] == ["parts added"]
