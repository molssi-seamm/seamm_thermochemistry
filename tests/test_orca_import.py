"""Tests for import_orca_atom_results: the S^2 vetting gate and the
duplicate-handling policy (unchanged if close, reported-not-overwritten if
not, unless force=True).
"""

import csv

import pytest

from seamm_thermochemistry import ThermoDB
from seamm_thermochemistry.importers import import_orca_atom_results

HEADER = [
    "Atomic Number",
    "Element",
    "Multiplicity",
    "E DFT@PBE0/def2-SV(P) (kJ/mol)",
    "S^2 DFT@PBE0/def2-SV(P)",
    "E DFT@PBE0/def2-SVP (kJ/mol)",
    "S^2 DFT@PBE0/def2-SVP",
]


def _write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(HEADER)
        writer.writerows(rows)
    return path


@pytest.fixture
def db(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        db.add_element(1, "H", multiplicity=2)
        db.add_element(2, "He", multiplicity=1)
        db.add_element(7, "N", multiplicity=4)
        db.add_element(65, "Tb", multiplicity=6)
        yield db


def test_clean_values_are_added(db, tmp_path):
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1308.22", "0.750001", "-1308.22", "0.750001"]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["added"]) == 2
    assert not summary["rejected"]
    assert not summary["spin_warnings"]
    assert db.get_atom_energy(
        "H", "orca", "PBE0", settings="def2-SV(P)"
    ) == pytest.approx(-1308.22)
    assert db.get_atom_energy(
        "H", "orca", "PBE0", settings="def2-SVP"
    ) == pytest.approx(-1308.22)


def test_blank_cell_is_skipped_not_an_error(db, tmp_path):
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1308.22", "0.750001", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["added"]) == 1
    assert db.get_atom_energy("H", "orca", "PBE0", settings="def2-SVP") is None


def test_singlet_with_blank_s2_is_fine(db, tmp_path):
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["2", "He", "1", "-7576.856", "", "-7576.856", ""]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["added"]) == 2
    assert not summary["rejected"]


def test_singlet_with_nonzero_s2_is_rejected(db, tmp_path):
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["2", "He", "1", "-7576.856", "0.5", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["rejected"]) == 1
    assert db.get_atom_energy("He", "orca", "PBE0", settings="def2-SV(P)") is None


def test_open_shell_missing_s2_is_rejected(db, tmp_path):
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1308.22", "", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["rejected"]) == 1
    assert "no S^2" in summary["rejected"][0]["reason"]


def test_moderately_off_s2_added_with_warning(db, tmp_path):
    # N: multiplicity 4, expected S^2 = 3.75. 10% off -> between the
    # default warn (2%) and reject (20%) thresholds.
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["7", "N", "4", "-142933.719", "4.125", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["added"]) == 1
    assert len(summary["spin_warnings"]) == 1
    assert db.get_atom_energy(
        "N", "orca", "PBE0", settings="def2-SV(P)"
    ) == pytest.approx(-142933.719)


def test_wildly_off_s2_is_rejected(db, tmp_path):
    # N: expected 3.75; 10.0 is far beyond the 20% reject threshold.
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["7", "N", "4", "-142933.719", "10.0", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["rejected"]) == 1
    assert db.get_atom_energy("N", "orca", "PBE0", settings="def2-SV(P)") is None


def test_tb_like_deviation_is_configurable(db, tmp_path):
    # Tb: multiplicity 6, expected S^2 = 8.75. The real ~20% deviation seen
    # for Tb in production sits right at the default reject boundary --
    # confirm both defaults (reject, just over 20%) and a looser threshold
    # (accept) work.
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["65", "Tb", "6", "-2160319.331", "10.6", "", ""]],  # ~21% off
    )
    summary = import_orca_atom_results(db, csv_path)
    assert len(summary["rejected"]) == 1

    summary2 = import_orca_atom_results(db, csv_path, reject_threshold=0.25)
    assert len(summary2["spin_warnings"]) == 1
    assert len(summary2["added"]) == 1


def test_duplicate_close_value_is_unchanged(db, tmp_path):
    db.add_atom_energy(
        "H", "orca", "PBE0", -1308.22, ref_type="atom", settings="def2-SV(P)"
    )
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1308.2201", "0.750001", "", ""]],  # tiny numeric noise
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["unchanged"]) == 1
    assert not summary["added"]
    assert not summary["conflicts"]
    assert db.get_atom_energy(
        "H", "orca", "PBE0", settings="def2-SV(P)"
    ) == pytest.approx(-1308.22)


def test_duplicate_conflict_is_not_overwritten_by_default(db, tmp_path):
    db.add_atom_energy(
        "H", "orca", "PBE0", -1308.22, ref_type="atom", settings="def2-SV(P)"
    )
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1400.00", "0.750001", "", ""]],  # genuinely different
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["conflicts"]) == 1
    assert summary["conflicts"][0]["delta"] == pytest.approx(-1400.00 - (-1308.22))
    # Original value is untouched.
    assert db.get_atom_energy(
        "H", "orca", "PBE0", settings="def2-SV(P)"
    ) == pytest.approx(-1308.22)


def test_duplicate_conflict_overwritten_with_force(db, tmp_path):
    db.add_atom_energy(
        "H", "orca", "PBE0", -1308.22, ref_type="atom", settings="def2-SV(P)"
    )
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1400.00", "0.750001", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path, force=True)

    assert len(summary["updated"]) == 1
    assert not summary["conflicts"]
    assert db.get_atom_energy(
        "H", "orca", "PBE0", settings="def2-SV(P)"
    ) == pytest.approx(-1400.00)


def test_dry_run_writes_nothing(db, tmp_path):
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1308.22", "0.750001", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path, dry_run=True)

    assert len(summary["added"]) == 1  # still classified/reported
    assert db.get_atom_energy("H", "orca", "PBE0", settings="def2-SV(P)") is None


def test_dry_run_does_not_overwrite_with_force(db, tmp_path):
    db.add_atom_energy(
        "H", "orca", "PBE0", -1308.22, ref_type="atom", settings="def2-SV(P)"
    )
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1400.00", "0.750001", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path, force=True, dry_run=True)

    assert len(summary["updated"]) == 1  # reported as what *would* happen
    assert db.get_atom_energy(
        "H", "orca", "PBE0", settings="def2-SV(P)"
    ) == pytest.approx(
        -1308.22
    )  # but nothing actually changed


def test_multiplicity_mismatch_is_reported_not_blocking(db, tmp_path):
    # DB says H is a doublet (2); CSV claims multiplicity 4 for H.
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "4", "-1308.22", "3.75", "", ""]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["multiplicity_mismatches"]) == 1
    assert summary["multiplicity_mismatches"][0] == {
        "symbol": "H",
        "csv_multiplicity": 4,
        "db_multiplicity": 2,
    }
    # Still processed using the CSV's own multiplicity for the S^2 gate.
    assert len(summary["added"]) == 1


def test_two_method_basis_pairs_in_one_row_both_parsed(db, tmp_path):
    csv_path = _write_csv(
        tmp_path / "results.csv",
        [["1", "H", "2", "-1308.22", "0.750001", "-1309.73", "0.750002"]],
    )
    summary = import_orca_atom_results(db, csv_path)

    assert len(summary["added"]) == 2
    methods = {(e["method"], e["basis"]) for e in summary["added"]}
    assert methods == {("PBE0", "def2-SV(P)"), ("PBE0", "def2-SVP")}
