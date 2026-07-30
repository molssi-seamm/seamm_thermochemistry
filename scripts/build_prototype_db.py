#!/usr/bin/env python3
"""One-off: build data/thermochemistry.db from Paul's legacy master files.

Not part of the test suite (the files live outside the repo, in
~/Downloads) -- run by hand to (re)generate the prototype database and the
reviewable CSV snapshots.
"""

from collections import Counter
from pathlib import Path
import time

from seamm_thermochemistry import ThermoDB, DEFAULT_DB_PATH
from seamm_thermochemistry import importers
from seamm_thermochemistry import derive_dfH0_0K

REFERENCE_XLSX = Path(
    "~/Downloads/Atom Reference Energies and States.xlsx"
).expanduser()
VASP_XLSX = Path("~/Downloads/VASP element_energies.xlsx").expanduser()

GAUSSIAN_CSV = Path(
    "~/Work/SEAMM/gaussian_step/gaussian_step/data/atom_energies.csv"
).expanduser()
PSI4_CSV = Path("~/Work/SEAMM/psi4_step/psi4_step/data/atom_energies.csv").expanduser()

# Full composite-method/basis grid (~5000+ columns each) -- pass an explicit
# list here instead of None for a quick partial rebuild during development.
GAUSSIAN_METHODS = None
PSI4_METHODS = None


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
            t0 = time.perf_counter()
            summary = importers.import_wide_method_csv(
                db,
                GAUSSIAN_CSV,
                "gaussian",
                methods=GAUSSIAN_METHODS,
                import_elements=False,
            )
            dt = time.perf_counter() - t0
            print(
                f"  Gaussian: {len(summary['methods'])} methods, "
                f"{summary['n_energies']} energies ({dt:.1f}s)"
            )

        if PSI4_CSV.exists():
            t0 = time.perf_counter()
            summary = importers.import_wide_method_csv(
                db,
                PSI4_CSV,
                "psi4",
                methods=PSI4_METHODS,
                import_elements=False,
            )
            dt = time.perf_counter() - t0
            print(
                f"  Psi4: {len(summary['methods'])} methods, "
                f"{summary['n_energies']} energies ({dt:.1f}s)"
            )

        # One line per (code, ref_type) rather than per method -- with the
        # full composite-method grid there are thousands of methods.
        methods = db.list_methods()
        print(f"\n{len(methods)} (code, method, ref_type, settings) combinations:")
        by_code_type = Counter((code, ref_type) for code, _, ref_type, _ in methods)
        for (code, ref_type), n_methods in sorted(by_code_type.items()):
            print(f"  {code:10s} {ref_type:14s} {n_methods:5d} methods")

        csv_dir = DEFAULT_DB_PATH.parent
        db.dump_elements_csv(csv_dir / "elements_snapshot.csv")
        db.dump_atom_energies_csv(csv_dir / "atom_energies_snapshot.csv")
        print(f"\nWrote review snapshots to {csv_dir}/")


if __name__ == "__main__":
    main()
