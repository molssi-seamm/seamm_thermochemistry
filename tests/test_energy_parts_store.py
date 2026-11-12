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


def test_parts_are_not_added_when_the_energy_differs(store):  # noqa: F811
    """Not the same calculation (energy off by more than 1e-10 E_h): the row is
    left without parts and the attempt is logged."""
    plain = EnergyRecord.make("g:2", -1000.0, energy_fingerprint="fp:y")
    store.put_energy(plain)
    off = EnergyRecord.make(
        "g:2", -1000.0 + 1e-6, energy_fingerprint="fp:y", byproducts=PARTS
    )
    assert store.put_energy(off)[0] == "duplicate"
    row = store.energy("g:2", energy_fingerprint=["fp:y"])
    assert row["byproducts"] is None
    history = store._connect().execute(
        "SELECT kind, reason FROM row_history WHERE table_name = 'energy'"
    )
    ((kind, reason),) = history.fetchall()
    assert kind == "parts refused" and "1e-10" in reason


def test_parts_need_the_same_fingerprint(store):  # noqa: F811
    """Another settings fingerprint is another row, never this one's parts."""
    store.put_energy(EnergyRecord.make("g:3", -1000.0, energy_fingerprint="fp:a"))
    other = EnergyRecord.make(
        "g:3", -1000.0, energy_fingerprint="fp:b", byproducts=PARTS
    )
    assert store.put_energy(other)[0] == "inserted"
    assert store.energy("g:3", energy_fingerprint=["fp:a"])["byproducts"] is None
