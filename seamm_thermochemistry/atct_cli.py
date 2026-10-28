# -*- coding: utf-8 -*-

"""``seamm-thermochemistry-atct``: manage the local copy of the Active
Thermochemical Tables and look entries up from the command line::

    seamm-thermochemistry-atct download [--from FILE] [--base DIR]
    seamm-thermochemistry-atct list [--base DIR]
    seamm-thermochemistry-atct lookup SMILES [--charge Q] [--multiplicity M]
                                      [--conformer LABEL] [--version V]
    seamm-thermochemistry-atct index [--version V]   (rebuild the index)
"""

import argparse
import sys

from . import atct as A


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="The local copy of the Active Thermochemical Tables (ATcT)."
    )
    parser.add_argument(
        "--base",
        help="The directory holding the ATcT versions (default: the installation's "
        "Parameters/thermochemistry/ATcT).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    download = sub.add_parser(
        "download",
        help="Download the complete tables (/api/v1/all/) as a new version. An "
        "existing version is never replaced.",
    )
    download.add_argument(
        "--from",
        dest="source",
        help="A local copy of the /all/ download, e.g. saved from a browser, "
        "instead of fetching it from atct.anl.gov.",
    )
    download.add_argument(
        "--user-agent",
        default=A.USER_AGENT,
        help="The User-Agent header to send (default: %(default)s).",
    )

    sub.add_parser("list", help="List the cached versions.")

    lookup = sub.add_parser("lookup", help="Find the ATcT entry for a molecule.")
    lookup.add_argument("smiles", help="The SMILES, with any double-bond stereo.")
    lookup.add_argument("--charge", type=int, default=0)
    lookup.add_argument("--multiplicity", type=int)
    lookup.add_argument("--conformer", help="A conformer or isomer label to select.")
    lookup.add_argument("--version", help="The TN version (default: newest).")
    lookup.add_argument(
        "--all", action="store_true", help="Show every matching entry instead."
    )

    index = sub.add_parser("index", help="Rebuild the SMILES index of a version.")
    index.add_argument("--version", help="The TN version (default: newest).")

    args = parser.parse_args(argv)

    try:
        if args.command == "download":
            directory = A.download_all(
                args.base, user_agent=args.user_agent, source=args.source
            )
            print(f"ATcT version {directory.name} is now in {directory}")
            tables = A.ATcT(version=directory.name, base=args.base)
            index = tables.index
            print(
                f"Indexed {sum(len(v) for v in index['stereo'].values())} entries by "
                f"SMILES; {len(index['unparsed'])} could not be parsed by RDKit "
                "(they can still be used by ATcT ID)."
            )
        elif args.command == "list":
            base = A.atct_directory(args.base)
            versions = A.ATcT.available(base)
            if not versions:
                print(f"No ATcT data in {base}")
            for version in versions:
                tables = A.ATcT(version=version, base=base)
                meta = tables.metadata
                print(
                    f"ver. {version}: {len(tables.entries)} entries, downloaded "
                    f"{meta.get('downloaded', '?')} from {meta.get('source', '?')}"
                )
        elif args.command == "lookup":
            tables = A.ATcT(version=args.version, base=args.base)
            if args.all:
                for entry in tables.candidates(args.smiles, gas_only=False):
                    print(entry.label())
            else:
                entry = tables.lookup(
                    args.smiles,
                    charge=args.charge,
                    multiplicity=args.multiplicity,
                    conformer=args.conformer,
                )
                print(entry.label())
                print(f"ΔfH°(0 K) = {entry['dfH0_0K']} {entry['units']}")
                print(f"TN version {tables.version}; {tables.citation}")
        elif args.command == "index":
            tables = A.ATcT(version=args.version, base=args.base)
            path = tables.directory / "index.json"
            path.unlink(missing_ok=True)
            index = tables.index
            print(
                f"Rebuilt the index of ver. {tables.version}: "
                f"{len(index['stereo'])} canonical SMILES, "
                f"{len(index['unparsed'])} unparsed entries."
            )
    except (A.ATcTError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
