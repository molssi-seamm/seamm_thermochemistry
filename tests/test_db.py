"""Tests for the ThermoDB schema and helper methods, against a tmp SQLite file."""

import pytest

from seamm_thermochemistry import ThermoDB


@pytest.fixture
def db(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        db.add_element(
            1,
            "H",
            multiplicity=2,
            term_symbol="2S1/2",
            standard_state="1/2 H2(g)",
            dfH0_0K=216.034,
            dfH0_298K=218.0,
            h298_minus_h0_atom=6.197,
            h298_minus_h0_std_state=8.468,
            s298_gas=114.72,
            reference="https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H",
        )
        db.add_element(6, "C", dfH0_0K=711.19, dfH0_298K=716.67)
        db.add_element(8, "O", dfH0_0K=246.79, dfH0_298K=249.18)
        yield db


def test_schema_creates_on_open(tmp_path):
    path = tmp_path / "fresh.db"
    assert not path.exists()
    with ThermoDB(path):
        pass
    assert path.exists()


def test_add_and_get_element(db):
    h = db.get_element(symbol="H")
    assert h["atomic_number"] == 1
    assert h["dfH0_0K"] == pytest.approx(216.034)
    assert db.get_element(atomic_number=1)["symbol"] == "H"


def test_get_element_unknown_returns_none(db):
    assert db.get_element(symbol="Xx") is None


def test_dfH0_0K_vs_298K(db):
    assert db.dfH0("H", at_0K=True) == pytest.approx(216.034)
    assert db.dfH0("H", at_0K=False) == pytest.approx(218.0)


def test_add_element_upsert_does_not_duplicate(db):
    db.add_element(1, "H", dfH0_0K=999.0)
    assert len(db.elements()) == 3
    assert db.get_element(symbol="H")["dfH0_0K"] == pytest.approx(999.0)


def test_add_atom_energy_requires_known_element(db):
    with pytest.raises(KeyError):
        db.add_atom_energy("Xx", "gaussian", "CBS-QB3", -100.0)


def test_add_and_get_atom_energy_roundtrip(db):
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0, units="kJ/mol")
    value = db.get_atom_energy("H", "gaussian", "CBS-QB3")
    assert value == pytest.approx(-1312.0)


def test_get_atom_energy_missing_returns_none(db):
    assert db.get_atom_energy("H", "gaussian", "CBS-QB3") is None


def test_add_atom_energy_upsert(db):
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1313.0)
    assert db.get_atom_energy("H", "gaussian", "CBS-QB3") == pytest.approx(-1313.0)
    assert len(db.list_methods("gaussian")) == 1


def test_atom_energy_unit_conversion(db):
    # 1 hartree = 2625.499639 kJ/mol
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -0.5, units="hartree")
    kJ = db.get_atom_energy("H", "gaussian", "CBS-QB3", units="kJ/mol")
    assert kJ == pytest.approx(-0.5 * 2625.499639)
    hartree = db.get_atom_energy("H", "gaussian", "CBS-QB3", units="hartree")
    assert hartree == pytest.approx(-0.5)


def test_atom_energy_correction_is_applied(db):
    db.add_atom_energy("H", "gaussian", "G4", -1300.0, correction=-2.5)
    assert db.get_atom_energy("H", "gaussian", "G4") == pytest.approx(-1302.5)


def test_ref_type_distinguishes_atom_vs_element_phase(db):
    db.add_atom_energy(
        "H", "vasp", "PBE", -55.0, ref_type="atom", settings="encut=700eV"
    )
    db.add_atom_energy(
        "H", "vasp", "PBE", -60.0, ref_type="element_phase", settings="encut=700eV"
    )
    atom = db.get_atom_energy(
        "H", "vasp", "PBE", ref_type="atom", settings="encut=700eV"
    )
    phase = db.get_atom_energy(
        "H", "vasp", "PBE", ref_type="element_phase", settings="encut=700eV"
    )
    assert atom == pytest.approx(-55.0)
    assert phase == pytest.approx(-60.0)


def test_get_reference_energies(db):
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
    db.add_atom_energy("C", "gaussian", "CBS-QB3", -98700.0)
    energies = db.get_reference_energies("gaussian", "CBS-QB3")
    assert energies == {"H": pytest.approx(-1312.0), "C": pytest.approx(-98700.0)}


def test_missing(db):
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
    assert db.missing("gaussian", "CBS-QB3", ["H", "C", "O"]) == ["C", "O"]


def test_coverage(db):
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
    db.add_atom_energy("O", "gaussian", "CBS-QB3", -98700.0)
    n, max_z = db.coverage("gaussian", "CBS-QB3")
    assert n == 2
    assert max_z == 8
    assert db.coverage("gaussian", "G4") == (0, None)


def test_dump_elements_csv_roundtrip(db, tmp_path):
    out = tmp_path / "elements.csv"
    db.dump_elements_csv(out)
    text = out.read_text()
    assert "atomic_number" in text.splitlines()[0]
    assert "H" in text


def test_dump_atom_energies_csv(db, tmp_path):
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
    out = tmp_path / "atom_energies.csv"
    db.dump_atom_energies_csv(out)
    text = out.read_text()
    assert "symbol" in text.splitlines()[0]
    assert "CBS-QB3" in text


def test_read_only_rejects_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        ThermoDB(tmp_path / "nope.db", read_only=True)


def test_batch_commits_once_at_the_end(db):
    with db.batch():
        db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
        # Uncommitted writes are still visible on the same connection...
        assert db.get_atom_energy("H", "gaussian", "CBS-QB3") == pytest.approx(-1312.0)
    # ...and still there (committed) after the block exits.
    assert db.get_atom_energy("H", "gaussian", "CBS-QB3") == pytest.approx(-1312.0)


def test_batch_restores_autocommit_after_exception(db):
    with pytest.raises(RuntimeError):
        with db.batch():
            db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
            raise RuntimeError("boom")
    assert db._autocommit is True
    # add_atom_energy after the failed batch commits immediately again.
    db.add_atom_energy("O", "gaussian", "CBS-QB3", -98700.0)
    assert db.get_atom_energy("O", "gaussian", "CBS-QB3") == pytest.approx(-98700.0)
