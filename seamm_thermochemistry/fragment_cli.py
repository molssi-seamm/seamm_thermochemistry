# -*- coding: utf-8 -*-

"""``seamm-thermochemistry-fragments``: the fragment-energy store from the
command line::

    seamm-thermochemistry-fragments harvest JOB_DIR [JOB_DIR ...] [--store FILE]
    seamm-thermochemistry-fragments list [SMILES] [--store FILE]
    seamm-thermochemistry-fragments lookup SMILES [--charge Q] [--multiplicity M]
                                           [--model-chemistry MC] [--fingerprint FP]
    seamm-thermochemistry-fragments duplicates [--store FILE]
    seamm-thermochemistry-fragments harvest JOB_DIR ... --supersede "REASON"
    seamm-thermochemistry-fragments superseded [SMILES] [--store FILE]

``harvest`` is the single writer: run it on the host after a job has
finished, never on a compute node.
"""

import argparse
import sys

from . import fragment_store as F


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="The fragment-energy store of CBH species."
    )
    parser.add_argument(
        "--store",
        help="The store file (default: fragments.db next to the ATcT cache).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    harvest = sub.add_parser(
        "harvest",
        help="Promote the species energies of finished jobs into the store "
        "(insert-or-ignore, with provenance).",
    )
    harvest.add_argument("jobs", nargs="+", help="Job directories.")
    harvest.add_argument("--energy", default="total energy#*#*")
    harvest.add_argument("--zpe", default="ZPE#*#*")
    harvest.add_argument("--enthalpy", default="H(298.15 K)#*#*")
    harvest.add_argument("--fingerprint", default="settings fingerprint#*#*")
    harvest.add_argument("--prefix", default="CBH:", help="System-name prefix.")
    harvest.add_argument(
        "--all-systems",
        action="store_true",
        help="Harvest every system, its name (less the prefix) being the SMILES.",
    )
    harvest.add_argument(
        "--supersede",
        metavar="REASON",
        help="Replace stored values that differ (a corrected conformer, say), keeping "
        "the old ones in the superseded table with this reason. Without it a "
        "differing value is only logged.",
    )

    lst = sub.add_parser("list", help="List the records, or those of one species.")
    lst.add_argument("species", nargs="?", help="Canonical SMILES.")

    lookup = sub.add_parser("lookup", help="The best record for a species.")
    lookup.add_argument("species", help="Canonical SMILES.")
    lookup.add_argument("--charge", type=int, default=0)
    lookup.add_argument("--multiplicity", type=int, default=1)
    lookup.add_argument("--model-chemistry")
    lookup.add_argument("--fingerprint", help="The energy fingerprint.")

    sub.add_parser("duplicates", help="Show the duplicate log.")

    sup = sub.add_parser("superseded", help="Show the superseded values.")
    sup.add_argument("species", nargs="?", help="Canonical SMILES.")

    args = parser.parse_args(argv)

    try:
        if args.command == "harvest":
            store = F.FragmentStore(args.store, create=True)
            with store:
                for job in args.jobs:
                    counts = F.harvest_job(
                        job,
                        store,
                        energy=args.energy,
                        zpe=args.zpe,
                        enthalpy=args.enthalpy,
                        fingerprint=args.fingerprint,
                        name_prefix=args.prefix,
                        all_systems=args.all_systems,
                        supersede=args.supersede,
                    )
                    print(
                        f"{job}: {counts['inserted']} inserted, "
                        f"{counts['duplicate']} already there, "
                        f"{counts['logged']} differing (logged), "
                        f"{counts['superseded']} superseded"
                    )
                    for reason in counts["skipped"]:
                        print(f"    skipped {reason}")
            print(f"Store: {store.path}")
        elif args.command == "list":
            store = F.FragmentStore(args.store)
            records = store.records(args.species)
            if not records:
                print(f"No records in {store.path}")
            for r in records:
                _print(r)
        elif args.command == "lookup":
            store = F.FragmentStore(args.store)
            r = store.get(
                args.species,
                charge=args.charge,
                multiplicity=args.multiplicity,
                energy_fingerprint=args.fingerprint,
                model_chemistry=args.model_chemistry,
            )
            if r is None:
                print("No record.")
                return 1
            _print(r)
        elif args.command == "superseded":
            store = F.FragmentStore(args.store)
            rows = store.superseded(args.species)
            if not rows:
                print("Nothing superseded.")
            for d in rows:
                old = ", ".join(
                    f"{k} {d[k]:.3f}" for k in ("E", "ZPE", "H") if d.get(k) is not None
                )
                new = ", ".join(
                    f"{k} {d['current_' + k]:.3f}"
                    for k in ("E", "ZPE", "H")
                    if d.get("current_" + k) is not None
                )
                print(
                    f"{d['species']} on {d['superseded_on']}: {old} (from {d['job']}) "
                    f"replaced by {new}; reason: {d['reason']}"
                )
        elif args.command == "duplicates":
            store = F.FragmentStore(args.store)
            rows = store.duplicates()
            if not rows:
                print("No differing duplicates logged.")
            for d in rows:
                print(
                    f"{d['species']} {d['model_chemistry']}: stored E {d['stored_E']} "
                    f"vs {d['E']} (max difference {d['difference']:.3f} kJ/mol) "
                    f"from {d['job']} on {d['date']}"
                )
    except (F.sqlite3.Error, FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


def _print(r):
    flag = "" if r["verified"] else " [unverified]"
    values = ", ".join(
        f"{k} = {r[k]:.3f}" for k in ("E", "ZPE", "H") if r.get(k) is not None
    )
    print(
        f"{r['species']} (q={r['charge']}, 2S+1={r['multiplicity']}) "
        f"{r['model_chemistry']}{flag}: {values} kJ/mol; from {r['job']} "
        f"on {r['date']}"
    )


if __name__ == "__main__":
    sys.exit(main())
