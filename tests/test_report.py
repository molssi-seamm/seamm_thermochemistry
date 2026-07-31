"""Tests for report.format_report: the shared detailed formation-energy
report (DB version, per-atom breakdown table, references)."""

import pytest

from seamm_thermochemistry import ThermoDB, format_report


@pytest.fixture
def db(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        db.add_element(
            1,
            "H",
            dfH0_0K=100.0,
            dfH0_298K=110.0,
            standard_state="1/2 H2(g)",
            reference="http://example.org/H",
            reference_note="TEST",
        )
        db.add_element(
            8,
            "O",
            dfH0_0K=200.0,
            dfH0_298K=210.0,
            standard_state="1/2 O2(g)",
            reference="http://example.org/O",
            reference_note="TEST",
        )
        db.set_s298_std_state("H", 50.0, reference="http://example.org/std_H")
        db.set_s298_std_state("O", 80.0, reference="http://example.org/std_O")
        db.add_atom_energy(
            "H", "test", "M", -10.0, ref_type="atom", source="h_atom.out"
        )
        db.add_atom_energy(
            "O", "test", "M", -70.0, ref_type="atom", source="o_atom.out"
        )
        yield db


H2O = {"H": 2, "O": 1}
SYSTEM_ENERGY = -100.0


def test_report_header_and_db_provenance_unpublished(db):
    text = format_report(H2O, db, "test", "M", name="water", level_label="test/M")
    assert "Thermochemistry of water with test/M" in text
    assert str(db.path) in text
    assert "unpublished working copy" in text
    assert "DOI" not in text


def test_report_db_provenance_prefers_doi(db):
    db.set_doi("10.5281/zenodo.21612188")
    text = format_report(H2O, db, "test", "M", name="water", level_label="test/M")
    assert "https://doi.org/10.5281/zenodo.21612188" in text
    assert "unpublished" not in text


def test_report_without_system_energy_has_no_sections(db):
    text = format_report(H2O, db, "test", "M")
    assert "Atomization Energy" not in text
    assert "Energy of Formation" not in text


def test_report_atomization_and_dfe0_sections(db):
    text = format_report(H2O, db, "test", "M", system_energy=SYSTEM_ENERGY)
    assert "Atomization Energy" in text
    assert "964.3093" not in text  # sanity: not leaking some other number
    assert "H(g)" in text and "O(g)" in text
    assert "h_atom.out" in text and "o_atom.out" in text
    assert "Energy of Formation (0 K" in text
    assert "http://example.org/H" in text
    # Numbers match the library function directly.
    from seamm_thermochemistry import atomization_energy, formation_energy

    atomization = atomization_energy(H2O, SYSTEM_ENERGY, db, "test", "M")
    dfe0 = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M")
    assert f"{atomization:.4f}" in text
    assert f"{dfe0:.4f}" in text


def test_report_dfht_dfgt_sections_when_given(db):
    text = format_report(
        H2O,
        db,
        "test",
        "M",
        system_energy=SYSTEM_ENERGY,
        system_enthalpy=-95.0,
        system_gibbs_energy=-96.0,
        temperature=300.0,
    )
    assert "Enthalpy of Formation" in text
    assert "Gibbs Energy of Formation" in text
    assert "S298_std_state" in text
    assert "http://example.org/std_H" in text

    from seamm_thermochemistry import formation_enthalpy, formation_gibbs_energy

    dfh = formation_enthalpy(H2O, -95.0, 300.0, db, "test", "M")
    dfg = formation_gibbs_energy(H2O, -96.0, 300.0, db, "test", "M")
    assert f"{dfh:.4f}" in text
    assert f"{dfg:.4f}" in text


def test_report_missing_atom_energy_reports_cleanly_no_exception():
    from seamm_thermochemistry import ThermoDB as TDB
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        with TDB(Path(tmp) / "empty.db") as empty_db:
            text = format_report(
                {"H": 1, "Cl": 1}, empty_db, "test", "M", system_energy=-5.0
            )
    assert "Cannot calculate the atomization energy" in text


def test_report_missing_gibbs_data_does_not_block_other_sections(db, tmp_path):
    # A DB with dfH0 but no s298_std_state for one element: DfE0/DfHT still
    # come through; only the Gibbs section reports the gap.
    with ThermoDB(tmp_path / "no_gibbs.db") as db2:
        db2.add_element(1, "H", dfH0_0K=100.0)
        db2.add_atom_energy("H", "test", "M", -10.0)
        text = format_report(
            {"H": 2},
            db2,
            "test",
            "M",
            system_energy=-5.0,
            system_enthalpy=-4.0,
            system_gibbs_energy=-4.5,
            temperature=300.0,
        )
    assert "Atomization Energy" in text
    assert "Energy of Formation (0 K" in text
    assert "Enthalpy of Formation" in text
    assert "Cannot calculate the Gibbs energy of formation" in text
