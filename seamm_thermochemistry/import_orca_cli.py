#!/usr/bin/env python3
"""Import one or more ORCA atom-energy results CSVs into thermochemistry.db.

Installed as the ``seamm-thermochemistry-import-orca`` console script.

Usage
-----
    seamm-thermochemistry-import-orca results1.csv results2.csv ...
    seamm-thermochemistry-import-orca --dry-run results.csv
    seamm-thermochemistry-import-orca --force results.csv   # overwrite conflicts

See `seamm_thermochemistry.importers.import_orca_atom_results` for the CSV
shape expected and the vetting/duplicate-handling rules -- this module is
just the CLI + pretty-printing around that function.
"""

import argparse
from pathlib import Path

from .db import ThermoDB, DEFAULT_DB_PATH
from .importers import import_orca_atom_results


def _fmt(entry):
    bits = [entry["symbol"], f"{entry['method']}/{entry['basis']}"]
    if "energy" in entry:
        bits.append(f"{entry['energy']:.3f} kJ/mol")
    return " ".join(bits)


def report(csv_path, summary):
    print(f"\n{csv_path}")
    print(
        f"  added={len(summary['added'])} unchanged={len(summary['unchanged'])} "
        f"updated={len(summary['updated'])} conflicts={len(summary['conflicts'])} "
        f"rejected={len(summary['rejected'])} "
        f"spin_warnings={len(summary['spin_warnings'])}"
    )

    for entry in summary["rejected"]:
        if "method" in entry:
            print(f"  REJECTED  {_fmt(entry)} -- {entry['reason']}")
        else:
            print(f"  REJECTED  {entry['symbol']} -- {entry['reason']}")

    for entry in summary["spin_warnings"]:
        print(f"  SPIN?     {_fmt(entry)} -- {entry['reason']} (added anyway)")

    for entry in summary["conflicts"]:
        print(
            f"  CONFLICT  {_fmt(entry)} -- existing={entry['existing']:.3f} "
            f"new={entry['energy']:.3f} (delta={entry['delta']:+.3f}) "
            "-- NOT overwritten, rerun with --force to update"
        )

    for entry in summary["updated"]:
        print(
            f"  UPDATED   {_fmt(entry)} -- existing={entry['existing']:.3f} "
            f"new={entry['energy']:.3f} (delta={entry['delta']:+.3f})"
        )

    for entry in summary["multiplicity_mismatches"]:
        print(
            f"  MULT?     {entry['symbol']} -- CSV says {entry['csv_multiplicity']}, "
            f"database says {entry['db_multiplicity']}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_files", nargs="+", help="ORCA results CSV file(s)")
    parser.add_argument(
        "--db", default=None, help="Database path (default: DEFAULT_DB_PATH)"
    )
    parser.add_argument("--code", default="orca")
    parser.add_argument("--ref-type", default="atom")
    parser.add_argument(
        "--warn-threshold",
        type=float,
        default=0.02,
        help="Relative S^2 deviation above which to flag for review (default 2%%)",
    )
    parser.add_argument(
        "--reject-threshold",
        type=float,
        default=0.20,
        help="Relative S^2 deviation above which to skip entirely (default 20%%)",
    )
    parser.add_argument(
        "--energy-rel-tol",
        type=float,
        default=1e-6,
        help="Relative tolerance for treating a duplicate as unchanged",
    )
    parser.add_argument(
        "--energy-abs-tol",
        type=float,
        default=0.01,
        help="Absolute tolerance (kJ/mol) for treating a duplicate as unchanged",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing values that conflict, instead of just reporting them",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Classify and report everything without writing to the database",
    )
    args = parser.parse_args()

    db_path = Path(args.db).expanduser() if args.db else DEFAULT_DB_PATH

    totals = {
        "added": 0,
        "updated": 0,
        "unchanged": 0,
        "conflicts": 0,
        "rejected": 0,
        "spin_warnings": 0,
    }
    with ThermoDB(db_path) as db:
        for csv_file in args.csv_files:
            summary = import_orca_atom_results(
                db,
                csv_file,
                code=args.code,
                ref_type=args.ref_type,
                warn_threshold=args.warn_threshold,
                reject_threshold=args.reject_threshold,
                energy_rel_tol=args.energy_rel_tol,
                energy_abs_tol=args.energy_abs_tol,
                force=args.force,
                dry_run=args.dry_run,
            )
            report(csv_file, summary)
            for key in totals:
                totals[key] += len(summary[key])

    print("\n" + "=" * 60)
    label = "Would have" if args.dry_run else "Did"
    print(
        f"{label}: add {totals['added']}, update {totals['updated']}, "
        f"leave {totals['unchanged']} unchanged."
    )
    if totals["conflicts"]:
        print(
            f"{totals['conflicts']} conflicting value(s) NOT overwritten "
            "-- rerun with --force if that's wrong, or review them above."
        )
    if totals["rejected"]:
        print(f"{totals['rejected']} cell(s) rejected on the S^2 gate -- see above.")
    if totals["spin_warnings"]:
        print(
            f"{totals['spin_warnings']} cell(s) added despite an elevated S^2 "
            "-- worth a look."
        )


if __name__ == "__main__":
    main()
