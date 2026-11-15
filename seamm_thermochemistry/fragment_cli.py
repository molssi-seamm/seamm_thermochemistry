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

``lookup`` and ``geometries`` with ``--fingerprint`` answer as a calculation
would: a verified row only by its fingerprint. Without it they browse: the rows
at ``--model-chemistry``, verified ones too (listed first), or the levels the
species is stored at, with their fingerprints; when nothing matches they say
why (another charge or multiplicity, another level, quarantined rows). ``list``
and ``lookup`` print each record's fingerprints.
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
        "--require-converged",
        action="store_true",
        help="Refuse species with no record of whether their optimization "
        "converged (results from before it was recorded are taken by default).",
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

    lookup = sub.add_parser(
        "lookup",
        help="The records of a species at a level, or the levels it is stored at.",
        description="With --fingerprint, the record a calculation with that "
        "energy fingerprint would use. Without, the species' records at "
        "--model-chemistry (verified ones first), or the levels to choose from "
        "when it is stored at several. 'list SMILES' shows every record.",
    )
    lookup.add_argument("species", help="Canonical SMILES.")
    lookup.add_argument("--charge", type=int, default=0)
    lookup.add_argument("--multiplicity", type=int, default=1)
    lookup.add_argument(
        "--model-chemistry",
        help="The level, e.g. 'ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP'.",
    )
    lookup.add_argument("--fingerprint", help="The energy fingerprint.")

    geo = sub.add_parser(
        "geometries",
        help="The geometries of a species and the energies on each.",
        description="With --fingerprint, the geometries a calculation with that "
        "frequency fingerprint would use. Without, the species' geometries at "
        "--model-chemistry (verified ones first), or the levels to choose from "
        "when it is stored at several.",
    )
    geo.add_argument("species", help="Canonical SMILES.")
    geo.add_argument("--charge", type=int, default=0)
    geo.add_argument("--multiplicity", type=int, default=1)
    geo.add_argument(
        "--model-chemistry", help="The frequency level, e.g. 'ORCA:CMP@R2SCAN-3C'."
    )
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
        "--keep-job-paths",
        action="store_true",
        help="Keep each job's directory in the provenance; by default only the job "
        "id is kept, since an official store is published.",
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

    cand = sub.add_parser(
        "release-candidates",
        help="New values kept out by a row quarantined by hand, which passed the "
        "harvest's checks: candidates for 'release'.",
    )
    cand.add_argument("species", nargs="?", help="Canonical SMILES.")

    rel = sub.add_parser(
        "release",
        help="Release quarantined rows: replace each with a release candidate's "
        "values, keeping the quarantined ones as superseded.",
    )
    rel.add_argument("ids", nargs="+", type=int, help="Release candidate ids.")
    rel.add_argument("--reason", required=True, help="Why (recorded with each).")

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
                        require_converged=args.require_converged,
                    )
                    print(
                        f"{job}: {counts['inserted']} inserted, "
                        f"{counts['duplicate']} already there, "
                        f"{counts['logged']} differing (logged), "
                        f"{counts['superseded']} superseded"
                        + (
                            f", {counts['quarantined']} kept quarantined"
                            if counts.get("quarantined")
                            else ""
                        )
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
                    if counts.get("no convergence record"):
                        print(
                            f"    {counts['no convergence record']} species with no "
                            "record of whether their optimization converged (taken)"
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
            if args.fingerprint is None:
                return _browse(store, args, "fragment")
            r = store.get(
                args.species,
                charge=args.charge,
                multiplicity=args.multiplicity,
                energy_fingerprint=args.fingerprint,
                model_chemistry=args.model_chemistry,
            )
            if r is None:
                print(
                    f"No record of {_state(args)} with the energy fingerprint "
                    f"{args.fingerprint}."
                )
                _browse(store, args, "fragment", levels_only=True)
                return 1
            _print(r, store)
            for note in store.notes:
                print(f"    note: {note}")
        elif args.command == "geometries":
            store = LayeredStore(args.store, args.official)
            if args.fingerprint is None:
                return _browse(store, args, "geometry")
            rows = store.geometries(
                args.species,
                charge=args.charge,
                multiplicity=args.multiplicity,
                frequency_fingerprint=args.fingerprint,
                model_chemistry=args.model_chemistry,
                include_imaginary=args.include_imaginary,
            )
            if not rows:
                print(
                    f"No geometry of {_state(args)} with the frequency fingerprint "
                    f"{args.fingerprint}."
                )
                _browse(store, args, "geometry", levels_only=True)
                return 1
            for g in rows:
                _print_geometry(g, store)
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
        elif args.command == "release-candidates":
            store = F.FragmentStore(args.store)
            rows = store.release_candidates(args.species)
            if not rows:
                print("No release candidates.")

            def values(d):
                return ", ".join(
                    f"{k} {d[k]:.3f}" for k in ("E", "ZPE", "H") if d.get(k) is not None
                )

            for d in rows:
                label = f" conformer '{d['conformer']}'" if d["conformer"] else ""
                print(
                    f"{d['id']}: {d['species']} (q={d['charge']}, "
                    f"2S+1={d['multiplicity']}){label} {d['table']} row {d['row']}, "
                    f"quarantined ({d['quarantine']}): {values(d['old'])}; new "
                    f"{values(d['new'])} [{d['new']['model_chemistry']}] from "
                    f"{d['job']} on {d['date']}"
                )
        elif args.command == "release":
            with F.FragmentStore(args.store, create=True) as store:
                for i in args.ids:
                    store.release(i, args.reason)
                    print(f"Released with candidate {i}.")
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
                job_ids_only=not args.keep_job_paths,
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


def _layer(row, store):
    """' [local]' or ' [official <version>]' for a row read through layers."""
    if store is None or "layer" not in row:
        return ""
    return f" [{store.layer_name(row['layer'])}]"


def _print(r, store=None):
    flag = "" if r["verified"] else " [unverified]"
    values = ", ".join(
        f"{k} = {r[k]:.3f}" for k in ("E", "ZPE", "H") if r.get(k) is not None
    )
    print(
        f"{r['species']} (q={r['charge']}, 2S+1={r['multiplicity']}) "
        f"{_level(r)}{flag}: {values} kJ/mol; from {r['job']} "
        f"on {r['date']}{_layer(r, store)}"
    )
    fingerprints = [
        f"{what} {r[key]}"
        for what, key in (
            ("energy", "energy_fingerprint"),
            ("frequency", "frequency_fingerprint"),
        )
        if r.get(key)
    ]
    if fingerprints:
        print(f"    fingerprints: {'; '.join(fingerprints)}")


def _print_geometry(g, store):
    flag = "" if g["verified"] else " [unverified]"
    imaginary = f", {g['n_imaginary']} imaginary" if g.get("n_imaginary") else ""
    label = f" conformer '{g['conformer']}'" if g["conformer"] else ""
    print(
        f"{g['species']}{label} (q={g['charge']}, 2S+1={g['multiplicity']}) at "
        f"{_level(g)}{flag}: E = {_kj(g['E'])}, ZPE = {_kj(g['ZPE'])}, "
        f"H = {_kj(g['H'])} kJ/mol{imaginary}; {g['geometry_hash']} from {g['job']}"
        f"{_layer(g, store)}"
    )
    if g.get("frequency_fingerprint"):
        print(f"    frequency fingerprint: {g['frequency_fingerprint']}")
    for e in store.energies(g["geometry_hash"]):
        eflag = "" if e["verified"] else " [unverified]"
        if F.is_unusable(e["model_chemistry"]):
            continue
        print(
            f"    {_level(e)}{eflag}: E = {e['E']:.3f} kJ/mol from {e['job']}"
            f"{_layer(e, store)}; energy fingerprint {e['energy_fingerprint']}"
        )


def _state(args):
    return f"{args.species} (q={args.charge}, 2S+1={args.multiplicity})"


def _level_key(row, table):
    """What tells one level from another: the model chemistry, whether the row
    is verified, and its fingerprints (energy and frequency for records,
    frequency for geometries)."""
    if table == "fragment":
        fingerprints = (
            ("energy", row.get("energy_fingerprint") or ""),
            ("frequency", row.get("frequency_fingerprint") or ""),
        )
    else:
        fingerprints = (("frequency", row.get("frequency_fingerprint") or ""),)
    return (_level(row), bool(row["verified"]), fingerprints)


def _kj(value):
    return "None" if value is None else f"{value:.3f}"


def _browse(store, args, table, levels_only=False):
    """``lookup`` or ``geometries`` without a fingerprint: every usable row of
    the species at the charge and multiplicity, at the level asked for, or the
    levels to choose from; and, when nothing matches, why.

    A calculation matches a verified row only by its fingerprint and an
    unverified one by its model chemistry (:meth:`FragmentStore.get`,
    :meth:`FragmentStore.geometries`); people browsing the store usually know
    the model chemistry alone, so here it matches verified rows too, listed
    first. With `levels_only`, only the levels are listed (after a fingerprint
    matched nothing). Returns the exit code.
    """
    what = "record" if table == "fragment" else "geometry"
    species = args.species
    if table == "fragment":
        rows = store.all_records(species)
    else:
        rows = store.all_geometries(species)
    if not rows:
        if not levels_only:
            print(f"No {what} of {species} in {_stores(store)}.")
        return 1
    state = [
        r
        for r in rows
        if r["charge"] == args.charge and r["multiplicity"] == args.multiplicity
    ]
    if not state:
        others = sorted({(r["charge"], r["multiplicity"]) for r in rows})
        print(
            f"No {what} of {_state(args)}; it is stored with "
            + ", ".join(f"--charge {q} --multiplicity {m}" for q, m in others)
            + "."
        )
        return 1
    usable = [r for r in state if not F.is_unusable(r["model_chemistry"])]
    quarantined = [r for r in state if F.is_unusable(r["model_chemistry"])]
    hidden = []
    if table == "geometry" and not args.include_imaginary:
        hidden = [r for r in usable if r.get("n_imaginary")]
        usable = [r for r in usable if not r.get("n_imaginary")]

    def notes():
        for r in quarantined:
            print(
                f"    note: {what} {r['id']}{_layer(r, store)} is quarantined: "
                f"{F.quarantine_reason(r['model_chemistry'])}"
            )
        if hidden:
            print(
                f"    note: {len(hidden)} geometr{'y' if len(hidden) == 1 else 'ies'}"
                " with imaginary modes not shown (--include-imaginary shows them)"
            )

    levels = {}
    for r in usable:
        levels.setdefault(_level_key(r, table), []).append(r)

    def list_levels(heading):
        print(heading)
        for (level, verified, fingerprints), group in levels.items():
            flag = "" if verified else " [unverified]"
            layers = []
            for r in group:
                name = _layer(r, store).strip()
                if name not in layers:
                    layers.append(name)
            n = len({(r["layer"], r.get("conformer", "")) for r in group})
            if table == "geometry" and n > len(layers):
                count = f", {n} conformers"
            elif n > len(layers):
                count = f", {n} records"
            else:
                count = ""
            print(f"    {level}{flag} {' '.join(layers)}{count}")
            for what, fingerprint in fingerprints:
                if fingerprint:
                    print(f"        {what} fingerprint {fingerprint}")

    if levels_only:
        if levels:
            list_levels(f"{_state(args)} is stored at:")
        notes()
        return 1

    if args.model_chemistry:
        matching = [
            r
            for r in usable
            if F.same_model_chemistry(
                args.model_chemistry, r["model_chemistry"], None, r.get("program")
            )
        ]
        if not matching:
            print(f"No {what} of {_state(args)} at {args.model_chemistry}.")
            if levels:
                list_levels("It is stored at:")
            notes()
            return 1
        matching.sort(key=lambda r: not r["verified"])  # verified first
    elif len(levels) > 1:
        list_levels(
            f"{_state(args)} is stored at {len(levels)} levels; choose one with "
            "--model-chemistry or --fingerprint:"
        )
        notes()
        return 1
    else:
        matching = usable
    if not matching:
        print(f"No usable {what} of {_state(args)}.")
        notes()
        return 1
    for r in matching:
        if table == "fragment":
            _print(r, store)
        else:
            _print_geometry(r, store)
    notes()
    for note in store.notes:
        print(f"    note: {note}")
    return 0


def _stores(store):
    paths = [str(store.local.path)]
    if store.official is not None:
        paths.append(str(store.official.path))
    return " or ".join(paths)


if __name__ == "__main__":
    sys.exit(main())
