"""Smoke tests for the legacy-file importers, against Paul's real master
files. Skipped (not failed) when those files or openpyxl aren't present, so
this suite stays green in an environment without Paul's ~/Downloads.
"""

from pathlib import Path

import pytest

from seamm_thermochemistry import ThermoDB, formation_energy

openpyxl = pytest.importorskip("openpyxl")
from seamm_thermochemistry import importers  # noqa: E402

REFERENCE_XLSX = Path("~/Downloads/Atom Reference Energies and States.xlsx").expanduser()
VASP_XLSX = Path("~/Downloads/VASP element_energies.xlsx").expanduser()

needs_reference_xlsx = pytest.mark.skipif(
    not REFERENCE_XLSX.exists(), reason=f"{REFERENCE_XLSX} not present"
)
needs_vasp_xlsx = pytest.mark.skipif(
    not VASP_XLSX.exists(), reason=f"{VASP_XLSX} not present"
)


@needs_reference_xlsx
def test_import_reference_xlsx(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        imported = importers.import_reference_xlsx(db, REFERENCE_XLSX)
        assert "H" in imported
        assert "O" in imported
        assert len(imported) >= 100

        h = db.get_element(symbol="H")
        assert h["dfH0_0K"] == pytest.approx(216.034)
        assert h["dfH0_298K"] == pytest.approx(218.0)
        assert h["standard_state"] == "1/2 H2(g)"

        # Elements past Kr have no experimental DfH0 today -- reflect that
        # honestly rather than hiding it.
        assert db.get_element(symbol="U")["dfH0_298K"] is None


@needs_vasp_xlsx
def test_import_vasp_workbook(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        summary = importers.import_vasp_workbook(db, VASP_XLSX)
        assert len(summary["elements"]) >= 90
        assert summary["atom_energy"] > 0
        assert summary["element_phase"] > 0

        # PBE@700 is one of the populated functionals.
        n, max_z = db.coverage("vasp", "PBE", ref_type="atom", settings="encut=700eV")
        assert n >= 90
        n, max_z = db.coverage(
            "vasp", "PBE", ref_type="element_phase", settings="encut=700eV"
        )
        assert n >= 90


@needs_vasp_xlsx
def test_vasp_atom_and_element_phase_are_independent(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        importers.import_vasp_workbook(db, VASP_XLSX)
        atom = db.get_atom_energy(
            "H", "vasp", "PBE-D3BJ", ref_type="atom", settings="encut=700eV"
        )
        phase = db.get_atom_energy(
            "H", "vasp", "PBE-D3BJ", ref_type="element_phase", settings="encut=700eV"
        )
        assert atom is not None
        assert phase is not None
        assert atom != phase


@needs_reference_xlsx
@needs_vasp_xlsx
def test_formation_energy_against_vasp_testing_sheet(tmp_path):
    """Cross-check against the numbers already worked out by hand in the
    VASP workbook's "Testing" sheet: PBE E(H2O) = -14.247 eV.

    This only exercises the `anchor=False` (materials/DfE0) path, since the
    reference workbook's experimental DfH0 doesn't extend usefully here --
    it is a schema/arithmetic smoke test, not a chemistry validation.
    """
    with ThermoDB(tmp_path / "test.db") as db:
        importers.import_reference_xlsx(db, REFERENCE_XLSX)
        importers.import_vasp_workbook(db, VASP_XLSX, import_elements=False)

        h2o_energy_eV = -14.247
        result = formation_energy(
            {"H": 2, "O": 1},
            h2o_energy_eV,
            db,
            "vasp",
            "PBE",
            ref_type="element_phase",
            settings="encut=700eV",
            anchor=False,
            units="eV",
        )
        # Water's formation energy relative to the elements should be
        # negative (stable) and of a physically sane magnitude (a couple of
        # eV, not tens of eV).
        assert -5.0 < result < 0.0
