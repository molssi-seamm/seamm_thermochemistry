#!/usr/bin/env python3
"""One-off: build data/thermochemistry.db from Paul's legacy master files.

Not part of the test suite (the files live outside the repo, in
~/Downloads) -- run by hand to (re)generate the prototype database and the
reviewable CSV snapshots.
"""

from pathlib import Path

from seamm_thermochemistry import ThermoDB, DEFAULT_DB_PATH
from seamm_thermochemistry import importers
from seamm_thermochemistry import derive_dfH0_0K

REFERENCE_XLSX = Path("~/Downloads/Atom Reference Energies and States.xlsx").expanduser()
VASP_XLSX = Path("~/Downloads/VASP element_energies.xlsx").expanduser()

GAUSSIAN_CSV = Path(
    "~/Work/SEAMM/gaussian_step/gaussian_step/data/atom_energies.csv"
).expanduser()
PSI4_CSV = Path("~/Work/SEAMM/psi4_step/psi4_step/data/atom_energies.csv").expanduser()

# A small slice, not the full ~5000-column composite-method grid -- enough
# to prove the wide-CSV importer end to end.
GAUSSIAN_METHODS = ["CBS-QB3", "G4"]
PSI4_METHODS = ["CBS-QB3", "G4"]


def main():
    if DEFAULT_DB_PATH.exists():
        DEFAULT_DB_PATH.unlink()

    with ThermoDB(DEFAULT_DB_PATH) as db:
        print(f"Building {DEFAULT_DB_PATH} ...")

        imported = importers.import_reference_xlsx(db, REFERENCE_XLSX)
        print(f"  experimental reference data: {len(imported)} elements")

        filled = derive_dfH0_0K(db)
        print(f"  derived DfH0(0K) for {len(filled)} elements: {filled}")

        summary = importers.import_vasp_workbook(db, VASP_XLSX, import_elements=False)
        print(
            f"  VASP: {len(summary['elements'])} elements, "
            f"{summary['atom_energy']} atom energies, "
            f"{summary['element_phase']} element-phase energies"
        )

        if GAUSSIAN_CSV.exists():
            summary = importers.import_wide_method_csv(
                db, GAUSSIAN_CSV, "gaussian", methods=GAUSSIAN_METHODS,
                import_elements=False,
            )
            print(f"  Gaussian ({GAUSSIAN_METHODS}): {summary['n_energies']} energies")

        if PSI4_CSV.exists():
            summary = importers.import_wide_method_csv(
                db, PSI4_CSV, "psi4", methods=PSI4_METHODS, import_elements=False,
            )
            print(f"  Psi4 ({PSI4_METHODS}): {summary['n_energies']} energies")

        print("\nMethods now in the database:")
        for code, method, ref_type, settings in db.list_methods():
            n, max_z = db.coverage(code, method, ref_type=ref_type, settings=settings)
            print(f"  {code:10s} {method:15s} {ref_type:14s} {settings:16s} "
                  f"{n:3d} elements, up to Z={max_z}")

        csv_dir = DEFAULT_DB_PATH.parent
        db.dump_elements_csv(csv_dir / "elements_snapshot.csv")
        db.dump_atom_energies_csv(csv_dir / "atom_energies_snapshot.csv")
        print(f"\nWrote review snapshots to {csv_dir}/")


if __name__ == "__main__":
    main()
