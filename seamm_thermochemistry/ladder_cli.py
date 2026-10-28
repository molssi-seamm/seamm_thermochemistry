# -*- coding: utf-8 -*-

"""``seamm-thermochemistry-ladder``: laddered reference values for CBH
reactions from the fragment store (CBH-ANL)::

    seamm-thermochemistry-ladder build NAME --energy LEVEL [--frequencies LEVEL]
        [--zpe-scale 0.9863] [--atct-version V] [--store FILE] SMILES [SMILES ...]
    seamm-thermochemistry-ladder show NAME [--store FILE]
    seamm-thermochemistry-ladder list [--store FILE]

LEVEL is a settings fingerprint or a model chemistry string, as stored by
``seamm-thermochemistry-fragments harvest``.
"""

import argparse
import sys

from . import ladder as L


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Laddered (CBH-ANL) reference values from the fragment store."
    )
    parser.add_argument("--store", help="The fragment store file.")
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build", help="Build and save a named ladder.")
    build.add_argument("name", help="The name of the reference set.")
    build.add_argument(
        "species", nargs="+", help="SMILES of the species to ladder (targets too)."
    )
    build.add_argument(
        "--energy",
        required=True,
        help="Fingerprint or model chemistry of the energies.",
    )
    build.add_argument(
        "--frequencies",
        help="Fingerprint or model chemistry of the ZPE and thermal corrections "
        "(default: the same records as the energies).",
    )
    build.add_argument("--zpe-scale", type=float, default=L.DEFAULT_ZPE_SCALE)
    build.add_argument("--atct-version", help="The TN version (default newest).")
    build.add_argument("--atct-base", help="The ATcT directory.")
    build.add_argument("--protect", default="default")
    build.add_argument(
        "--fetch-covariance",
        action="store_true",
        help="Fetch missing covariance blocks from atct.anl.gov.",
    )

    show = sub.add_parser("show", help="Show a ladder.")
    show.add_argument("name")

    sub.add_parser("list", help="List the ladders.")

    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            ladder = L.build_ladder(
                args.species,
                name=args.name,
                energy_level=args.energy,
                frequency_level=args.frequencies,
                zpe_scale=args.zpe_scale,
                atct_version=args.atct_version,
                atct_base=args.atct_base,
                store_path=args.store,
                protect=args.protect,
                fetch_covariance=args.fetch_covariance,
            )
            path = ladder.save(args.store)
            print(ladder.summary())
            for s, why in ladder.data["missing"].items():
                print(f"  not laddered: {s}: {why}")
            if ladder.data["unverified"]:
                print(
                    "  unverified (no settings fingerprint): "
                    + ", ".join(ladder.data["unverified"])
                )
            print(f"Saved as {path}")
        elif args.command == "show":
            ladder = L.Ladder.load(args.name, args.store)
            print(ladder.summary())
        elif args.command == "list":
            directory = L.ladder_path("x", args.store).parent
            names = sorted(p.stem for p in directory.glob("*.json"))
            if not names:
                print(f"No ladders in {directory}")
            for name in names:
                print(name)
    except (L.LadderError, FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
