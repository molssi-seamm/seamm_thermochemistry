# -*- coding: utf-8 -*-

"""``seamm-thermochemistry-references``: the supplementary reference library
from the command line::

    seamm-thermochemistry-references sources [--library FILE]
    seamm-thermochemistry-references list [--library FILE]
    seamm-thermochemistry-references show SMILES [--charge Q] [--multiplicity M]
                                      [--local FILE] [--official FILE]
    seamm-thermochemistry-references check [--library FILE]
    seamm-thermochemistry-references accept MEASUREMENT --by WHO --reason TEXT
                                      [--library FILE]
    seamm-thermochemistry-references accept-file ACCEPTS.csv --library FILE
    seamm-thermochemistry-references add-curated DIR --library FILE
    seamm-thermochemistry-references webbook fetch SPECIES.csv --cache DIR
    seamm-thermochemistry-references webbook import --cache DIR --library FILE
    seamm-thermochemistry-references offsets --library FILE [--atct-version V]

``add-curated`` (maintainers) reads hand-curated CSVs (printed JANAF, CODATA,
D0, W4-17, derived values) into a library file being built for publication;
it always names its file. ``webbook fetch`` saves NIST WebBook pages for the
species in a CSV (smiles, charge, multiplicity) into a cache, a few seconds
apart and logged; ``webbook import`` reads the cache into a library, offline.
``offsets`` measures each source against ATcT (rms over shared species).

``--library`` is one library file, by default the local layer
(``references.db``); ``show`` reads both layers, local and official
(``references-official.db``). ``accept`` releases an estimate or a held
species, recording who and why.
"""

import argparse
import sys

from . import references as R


def _fmt(value, digits=2):
    return "" if value is None else f"{value:.{digits}f}"


def _library(path, write=False):
    path = R.local_references_path(path)
    if not write and not path.exists():
        raise SystemExit(f"No reference library at {path}")
    return R.ReferenceLibrary(path, create=write)


def cmd_sources(args):
    library = _library(args.library)
    print(f"{'rank':>4}  {'name':<16} {'kind':<11} {'offset':>7} {'n':>4}  citation")
    for s in library.sources():
        print(
            f"{s['rank']:>4}  {s['name']:<16} {s['kind']:<11} "
            f"{_fmt(s['offset_kj']):>7} {s['offset_n'] or '':>4}  {s['citation'] or ''}"
        )


def cmd_list(args):
    library = _library(args.library)
    db = library._connect()
    for row in db.execute("SELECT * FROM species ORDER BY formula, smiles"):
        best = library.best(row["id"])
        value = (
            f"{best['dfh_298k']:9.2f} ± {_fmt(best['expanded_uncertainty']):>5} "
            f"{best['source_name']}"
            if best is not None and best["dfh_298k"] is not None
            else "  (none chosen)"
        )
        print(f"{row['formula']:<10} {row['smiles']:<30} {row['charge']:>2} {value}")


def cmd_show(args):
    layers = R.LayeredReferences(args.local, args.official)
    rows = layers.candidates(
        args.smiles, args.charge, args.multiplicity, include_held=True
    )
    if not rows:
        print(f"No reference values for {args.smiles} (charge {args.charge}).")
        return 1
    for r in rows:
        mark = "*" if r["preferred"] else " "
        print(
            f"{mark} {r['layer']:<8} {r['source_name']:<14} "
            f"ΔfH298 {_fmt(r['dfh_298k']):>9}  ΔfH0 {_fmt(r['dfh_0k']):>9}  "
            f"± {_fmt(r['uncertainty'])} ({r['unc_level']}) -> 95 % "
            f"{_fmt(r['expanded_uncertainty'])}, total {_fmt(r['total_uncertainty'])}"
        )
        state = "" if r["state"] == "chosen" else f"{r['state'].upper()}, "
        print(f"    {r['citation']}  [{state}{r['status']}, {r['rule']}]")
        for key in ("flag", "override_reason", "converted", "decision_reason"):
            if r.get(key):
                print(f"    {key.replace('_', ' ')}: {r[key]}")
    return 0


def cmd_check(args):
    library = _library(args.library)
    problems = library.check()
    for p in problems:
        print(f"{p['smiles']} (measurement {p['measurement']}): {p['problem']}")
    if not problems:
        print("Every derived value agrees with its derivation.")
    return 1 if problems else 0


def cmd_accept(args):
    library = _library(args.library, write=True)
    ref = library.accept(args.measurement, by=args.by, reason=args.reason)
    print(f"Accepted measurement {args.measurement} (reference {ref}).")
    return 0


def cmd_accept_file(args):
    library = R.ReferenceLibrary(args.library, create=True)
    refs = library.accept_file(args.file)
    print(f"Applied {len(refs)} accepts.")
    return 0


def cmd_add_curated(args):
    from .curated import import_curated

    library = R.ReferenceLibrary(args.library, create=True, layer=args.layer)
    report, problems = import_curated(library, args.directory)
    counts = {}
    for filename, smiles, outcome, detail in report:
        counts[outcome] = counts.get(outcome, 0) + 1
        if outcome != "added":
            print(f"{outcome:>8}  {filename}: {smiles}  {detail}")
    print(", ".join(f"{n} {k}" for k, n in sorted(counts.items())))
    for p in problems:
        print(
            f"PROBLEM  {p['smiles']} (measurement {p['measurement']}): {p['problem']}"
        )
    return 1 if problems else 0


def cmd_webbook(args):
    import csv

    from . import webbook

    if args.action == "fetch":
        if not args.species:
            raise SystemExit("webbook fetch needs a CSV of species")
        with open(args.species, newline="", encoding="utf-8") as fd:
            species = [
                (
                    r["smiles"],
                    int(r.get("charge") or 0),
                    int(r.get("multiplicity") or 1),
                )
                for r in csv.DictReader(fd)
            ]
        rows = webbook.fetch(species, args.cache, delay=args.delay)
        found = sum(1 for r in rows if r["file"])
        print(f"{len(rows)} new species looked up, {found} found.")
        return 0
    if not args.library:
        raise SystemExit("webbook import needs --library")
    library = R.ReferenceLibrary(args.library, create=True, layer=args.layer)
    counts = {}
    for smiles, outcome, detail in webbook.import_webbook(library, args.cache):
        counts[outcome] = counts.get(outcome, 0) + 1
        if outcome == "refused":
            print(f" refused  {smiles}  {detail}")
    print(", ".join(f"{n} {k}" for k, n in sorted(counts.items())))
    return 0


def cmd_offsets(args):
    from .atct import ATcT

    library = R.ReferenceLibrary(args.library, create=True)
    atct = ATcT(version=args.atct_version, base=args.atct_base)
    result = R.compute_offsets(library, atct)
    print(f"Offsets against ATcT ver. {atct.version} (rms of source - ATcT):")
    for name, (rms, n, diffs) in result.items():
        used = (
            "used"
            if n >= R.OFFSET_MIN_SHARED
            else f"not used (< {R.OFFSET_MIN_SHARED})"
        )
        print(f"  {name:<16} n = {n:>3}  rms = {_fmt(rms)} kJ/mol  {used}")
        if args.verbose:
            for smiles, d in diffs:
                print(f"      {smiles:<28} {d:+8.2f}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="seamm-thermochemistry-references",
        description="The supplementary reference library of enthalpies of formation.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("sources", help="List the sources, with their offsets.")
    p.add_argument("--library")
    p.set_defaults(func=cmd_sources)

    p = sub.add_parser("list", help="List the species and their chosen values.")
    p.add_argument("--library")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("show", help="Every value for one species, both layers.")
    p.add_argument("smiles")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--multiplicity", type=int)
    p.add_argument("--local")
    p.add_argument("--official")
    p.set_defaults(func=cmd_show)

    p = sub.add_parser("check", help="Check every derived value.")
    p.add_argument("--library")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("accept", help="Accept an estimate or a held value.")
    p.add_argument("measurement", type=int)
    p.add_argument("--by", required=True)
    p.add_argument("--reason", required=True)
    p.add_argument("--library")
    p.set_defaults(func=cmd_accept)

    p = sub.add_parser("accept-file", help="Apply accepts kept in a CSV (maintainers).")
    p.add_argument("file")
    p.add_argument("--library", required=True)
    p.set_defaults(func=cmd_accept_file)

    p = sub.add_parser("add-curated", help="Read curated CSVs (maintainers).")
    p.add_argument("directory")
    p.add_argument("--library", required=True)
    p.add_argument("--layer", default="official", choices=("official", "local"))
    p.set_defaults(func=cmd_add_curated)

    p = sub.add_parser("webbook", help="NIST WebBook pages (maintainers).")
    p.add_argument("action", choices=("fetch", "import"))
    p.add_argument("species", nargs="?")
    p.add_argument("--cache", required=True)
    p.add_argument("--library")
    p.add_argument("--layer", default="official", choices=("official", "local"))
    p.add_argument("--delay", type=float, default=3.0)
    p.set_defaults(func=cmd_webbook)

    p = sub.add_parser("offsets", help="Measure each source against ATcT.")
    p.add_argument("--library", required=True)
    p.add_argument("--atct-version")
    p.add_argument("--atct-base")
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(func=cmd_offsets)

    args = parser.parse_args(argv)
    try:
        return args.func(args) or 0
    except (R.ReferenceError, KeyError, FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
