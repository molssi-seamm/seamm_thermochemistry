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
    seamm-thermochemistry-fragments writeback-log [--job JOB] [--species SMILES]
    seamm-thermochemistry-fragments duplicates --layers [--official FILE]
    seamm-thermochemistry-fragments migrate [--store FILE]
    seamm-thermochemistry-fragments merge SOURCE.db ... [--store FILE]
    seamm-thermochemistry-fragments export --out FILE [--label TEXT]
    seamm-thermochemistry-fragments official build --version V --out FILE
                                           [--previous FILE] SOURCE.db ...

``--store`` is the local store (``fragments.db``), the one every command
writes; ``--official`` is the official store (``fragments-official.db``, read
only), which ``lookup`` and ``geometries`` read after the local one. Run
``harvest`` on the host after a job has finished, never on a compute node.
"""

import argparse
import sys

from . import fragment_store as F
from . import store_tools as T
from .layered_store import LayeredStore


def _print_counts(counts):
    for table in ("fragment", "geometry", "energy"):
        tally = counts.get(table) or {}
        if tally:
            print(
                f"    {table} rows: "
                + ", ".join(f"{n} {k}" for k, n in sorted(tally.items()))
            )
    if counts.get("history"):
        print(f"    {counts['history']} superseded values copied as history")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="The fragment-energy store of CBH species."
    )
    parser.add_argument(
        "--store",
        help="The local store file (default: fragments.db next to the ATcT cache).",
    )
    parser.add_argument(
        "--official",
        help="The official store file (default: fragments-official.db beside the "
        "local one).",
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
        "--species-map",
        metavar="CSV",
        help="Find the systems by name in this map (columns: job, system name, "
        "species; optionally configuration name, charge, multiplicity, conformer "
        "label), which gives each one's species, instead of the system-name prefix.",
    )
    harvest.add_argument(
        "--tables",
        choices=("all", "geometry"),
        default="all",
        help="'geometry' fills only the geometry and energy tables (e.g. several "
        "conformers of a species), not the fragment records.",
    )
    harvest.add_argument(
        "--supersede",
        metavar="REASON",
        help="Replace stored values that differ (a corrected conformer, say), keeping "
        "the old ones in the superseded table with this reason. Without it a "
        "differing value is only logged.",
    )
    harvest.add_argument(
        "--no-structure-check",
        action="store_true",
        help="Do not refuse structures whose bonds, perceived again from the "
        "coordinates, are not those of their species, or that have atoms too close "
        "together.",
    )
    harvest.add_argument(
        "--force",
        action="store_true",
        help="With --supersede, also replace values that differ by less than the "
        "duplicate tolerance (0.1 kJ/mol), e.g. to apply a correction to every "
        "record of a level of theory. Identical values are left alone.",
    )

    lst = sub.add_parser("list", help="List the records, or those of one species.")
    lst.add_argument("species", nargs="?", help="Canonical SMILES.")

    lookup = sub.add_parser("lookup", help="The best record for a species.")
    lookup.add_argument("species", help="Canonical SMILES.")
    lookup.add_argument("--charge", type=int, default=0)
    lookup.add_argument("--multiplicity", type=int, default=1)
    lookup.add_argument("--model-chemistry")
    lookup.add_argument("--fingerprint", help="The energy fingerprint.")

    geo = sub.add_parser(
        "geometries",
        help="The geometries of a species and the energies on each.",
    )
    geo.add_argument("species", help="Canonical SMILES.")
    geo.add_argument("--charge", type=int, default=0)
    geo.add_argument("--multiplicity", type=int, default=1)
    geo.add_argument("--model-chemistry", help="The frequency level.")
    geo.add_argument("--fingerprint", help="The frequency fingerprint.")
    geo.add_argument(
        "--include-imaginary",
        action="store_true",
        help="Also geometries with imaginary modes.",
    )

    dup = sub.add_parser("duplicates", help="Show the duplicate log.")
    dup.add_argument(
        "--layers",
        action="store_true",
        help="Instead, the keys whose values differ between the local and the "
        "official store.",
    )

    sub.add_parser(
        "migrate",
        help="Bring the local store up to date (the meta table; model chemistries "
        "stored bare with their program), after a dated backup.",
    )

    merge = sub.add_parser(
        "merge",
        help="Add the rows of other stores (another installation's copy, a "
        "contribution, a campaign store) to the local store.",
    )
    merge.add_argument("sources", nargs="+", help="Store files.")
    merge.add_argument("--label", help="A name for the sources, for the history.")

    export = sub.add_parser(
        "export",
        help="Write the local rows the official store lacks, or has with other "
        "values, to a contribution file for the next official version.",
    )
    export.add_argument("--out", required=True, help="The new contribution file.")
    export.add_argument("--label", help="Free text recorded in the file.")

    official = sub.add_parser("official", help="Official stores (maintainers).")
    official_sub = official.add_subparsers(dest="official_command", required=True)
    build = official_sub.add_parser(
        "build",
        help="Build a new official store from the previous one and contributions.",
    )
    build.add_argument("sources", nargs="*", help="Contribution or store files.")
    build.add_argument("--version", required=True, help="Its version (CalVer).")
    build.add_argument("--out", required=True, help="The new file.")
    build.add_argument("--previous", help="The previous official store.")
    build.add_argument("--source", help="Where it is published, e.g. a DOI.")
    build.add_argument(
        "--frequency-fingerprint",
        action="append",
        default=[],
        help="Keep only geometries with this frequency fingerprint (repeat for "
        "several), and what is on them.",
    )
    build.add_argument(
        "--energy-fingerprint",
        action="append",
        default=[],
        help="Keep only energies with this fingerprint (repeat for several).",
    )
    build.add_argument(
        "--exclude",
        metavar="CSV",
        help="Species states to leave out, with the reason: columns species, "
        "charge, multiplicity, reason.",
    )
    build.add_argument(
        "--no-atct",
        action="store_true",
        help="Do not check multiplicities against ATcT's spin states.",
    )

    wlog = sub.add_parser(
        "writeback-log",
        help="What harvests and jobs' write-back added to or refused from the store.",
    )
    wlog.add_argument("--job", help="Only this job (part of its label, e.g. 5689).")
    wlog.add_argument("--species", help="Only this species (canonical SMILES).")

    sup = sub.add_parser("superseded", help="Show the superseded values.")
    sup.add_argument("species", nargs="?", help="Canonical SMILES.")

    args = parser.parse_args(argv)
    if getattr(args, "force", False) and not args.supersede:
        parser.error("--force needs --supersede REASON")

    try:
        if args.command == "harvest":
            species_map = (
                F.load_species_map(args.species_map) if args.species_map else None
            )
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
                        force=args.force,
                        species_map=species_map,
                        tables=args.tables,
                        structure_check=not args.no_structure_check,
                    )
                    print(
                        f"{job}: {counts['inserted']} inserted, "
                        f"{counts['duplicate']} already there, "
                        f"{counts['logged']} differing (logged), "
                        f"{counts['superseded']} superseded"
                    )
                    for table in ("geometry", "energy"):
                        tally = counts.get(table) or {}
                        if tally:
                            print(
                                f"    {table} rows: "
                                + ", ".join(
                                    f"{n} {k}" for k, n in sorted(tally.items())
                                )
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
            store = LayeredStore(args.store, args.official)
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
            for note in store.notes:
                print(f"    note: {note}")
        elif args.command == "geometries":
            store = LayeredStore(args.store, args.official)
            rows = store.geometries(
                args.species,
                charge=args.charge,
                multiplicity=args.multiplicity,
                frequency_fingerprint=args.fingerprint,
                model_chemistry=args.model_chemistry,
                include_imaginary=args.include_imaginary,
            )
            if not rows:
                print("No geometry.")
                return 1
            for g in rows:
                flag = "" if g["verified"] else " [unverified]"
                imaginary = (
                    f", {g['n_imaginary']} imaginary" if g.get("n_imaginary") else ""
                )
                label = f" conformer '{g['conformer']}'" if g["conformer"] else ""
                print(
                    f"{g['species']}{label} at {_level(g)}{flag}: "
                    f"E = {g['E']:.3f}, ZPE = {g['ZPE']}, H = {g['H']} kJ/mol"
                    f"{imaginary}; {g['geometry_hash']} from {g['job']} "
                    f"[{store.layer_name(g['layer'])}]"
                )
                for e in store.energies(g["geometry_hash"]):
                    eflag = "" if e["verified"] else " [unverified]"
                    print(
                        f"    {_level(e)}{eflag}: E = {e['E']:.3f} kJ/mol "
                        f"from {e['job']} [{store.layer_name(e['layer'])}]"
                    )
            for note in store.notes:
                print(f"    note: {note}")
        elif args.command == "writeback-log":
            store = F.FragmentStore(args.store)
            rows = store.write_log(job=args.job, species=args.species)
            if not rows:
                print("Nothing logged.")
            for r in rows:
                what = r["table_name"] or "species"
                conformer = f" conformer '{r['conformer']}'" if r["conformer"] else ""
                row = f" row {r['row']}" if r["row"] is not None else ""
                reason = f": {r['reason']}" if r["reason"] else ""
                print(
                    f"{r['date']} {r['kind']} {r['job']}: {r['species']} "
                    f"(q={r['charge']}, 2S+1={r['multiplicity']}){conformer} "
                    f"{what}{row} {r['outcome']}{reason}"
                )
        elif args.command == "superseded":
            store = F.FragmentStore(args.store)
            rows = store.superseded(args.species)
            others = store.superseded_rows(args.species)
            if not rows and not others:
                print("Nothing superseded.")

            def values(d):
                return ", ".join(
                    f"{k} {d[k]:.3f}" for k in ("E", "ZPE", "H") if d.get(k) is not None
                )

            for d in others:
                label = f" conformer '{d['conformer']}'" if d["conformer"] else ""
                where = f"{d['table']} row {d['row']}"
                print(
                    f"{d['species']}{label} {where} on {d['date']}: "
                    f"{values(d['old'])} [{d['old']['model_chemistry']}] (from "
                    f"{d['old']['job']}) replaced by {values(d['new'])} "
                    f"[{d['new']['model_chemistry']}] (from {d['new']['job']}); "
                    f"reason: {d['reason']}"
                )
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
        elif args.command == "duplicates" and args.layers:
            rows = T.layer_conflicts(args.store, args.official)
            if not rows:
                print("The local and official stores agree on every shared key.")
            for d in rows:
                key = ", ".join(str(k) for k in d["key"])
                how = (
                    "the official value is used: the local one was superseded"
                    if d["superseded"]
                    else "the local value is used"
                )
                print(
                    f"{d['table']} ({key}): local {d['local']}, official "
                    f"{d['official']}, {d['difference']:.2f} kJ/mol apart; {how}"
                )
        elif args.command == "migrate":
            result = T.migrate_store(args.store)
            if result["backup"] is None:
                print("The store is up to date; nothing changed.")
            else:
                old, new = result["schema"]
                changed = ", ".join(
                    f"{n} {t} rows" for t, n in result["changed"].items()
                )
                print(
                    f"Migrated (schema {old} -> {new}; names normalized: "
                    f"{changed or 'none'}). Backup: {result['backup']}"
                )
        elif args.command == "merge":
            counts = T.merge_stores(
                F.fragment_store_path(args.store), args.sources, label=args.label
            )
            _print_counts(counts)
        elif args.command == "export":
            counts = T.export_local(
                args.out, args.store, args.official, label=args.label
            )
            print(
                f"Exported to {args.out}: {counts['fragment']} fragment, "
                f"{counts['geometry']} geometry and {counts['energy']} energy rows "
                f"({counts['differing']} differing from the official store)."
            )
        elif args.command == "official" and args.official_command == "build":
            atct = None
            if not args.no_atct:
                try:
                    from .atct import ATcT

                    atct = ATcT()
                except Exception as e:
                    print(f"Note: no ATcT for the multiplicity check ({e}).")
            path, digest, counts = T.build_official(
                args.out,
                args.version,
                args.sources,
                previous=args.previous,
                source=args.source,
                frequency_fingerprints=args.frequency_fingerprint,
                energy_fingerprints=args.energy_fingerprint,
                exclude=T.load_exclusions(args.exclude) if args.exclude else None,
                atct=atct,
            )
            _print_counts(counts)
            print(f"Built {path} (version {args.version}), sha256 {digest}")
            print(f"Report: {path.with_name(path.name + '.report.json')}")
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
    except (F.sqlite3.Error, FileNotFoundError, FileExistsError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


def _level(row):
    """The row's model chemistry as people write it: PROGRAM:model when the
    program is known (rows store the bare model and the program)."""
    model, program = F.split_model_chemistry(
        row.get("model_chemistry"), row.get("program")
    )
    return f"{program}:{model}" if program else f"{model}"


def _print(r):
    flag = "" if r["verified"] else " [unverified]"
    values = ", ".join(
        f"{k} = {r[k]:.3f}" for k in ("E", "ZPE", "H") if r.get(k) is not None
    )
    print(
        f"{r['species']} (q={r['charge']}, 2S+1={r['multiplicity']}) "
        f"{_level(r)}{flag}: {values} kJ/mol; from {r['job']} "
        f"on {r['date']}"
    )


if __name__ == "__main__":
    sys.exit(main())
