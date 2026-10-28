#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The drift of CBH-n enthalpies of formation against ATcT's own values.

Elliott et al. (JPCA 127, 1512 (2023)) found that CBH-2 results referenced
to ATcT values for the fragments drift upward along the n-alkanes, by about
(n-3) times the propane/ethane reference mismatch, while "laddered" values
do not. This script reproduces that comparison for any set of targets:

    cbh_drift.py energies.csv [--rung 2] [--protect none] [--kcal]

``energies.csv`` has the columns ``smiles`` and ``dfH0_0K`` (kJ/mol, or
kcal/mol with --kcal): the enthalpies of formation at 0 K of every target
and fragment at ONE consistent level of theory (e.g. ANL0 or G4 atomization
values), from which the reaction enthalpies follow exactly. The ATcT values
come from the local cache (seamm-thermochemistry-atct). The reference
uncertainty (±ref) uses the cached covariance block when there is one,
otherwise quadrature, and says which; ±drift is the ATcT uncertainty of the
drift itself, references minus the target, which is correlated with them
(``seamm-thermochemistry-atct covariance`` caches the block).
"""

import argparse
import csv
import sys

from seamm_thermochemistry import ATcT, cbh, reaction_uncertainty
from seamm_thermochemistry.atct import canonical_smiles

KCAL = 4.184


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("energies", help="CSV with columns smiles, dfH0_0K")
    parser.add_argument("--rung", type=int, default=2)
    parser.add_argument("--protect", default="none")
    parser.add_argument("--kcal", action="store_true", help="values in kcal/mol")
    parser.add_argument("--version", help="ATcT version (default newest)")
    args = parser.parse_args(argv)

    values = {}
    with open(args.energies, newline="") as fd:
        for row in csv.DictReader(line for line in fd if not line.startswith("#")):
            can = canonical_smiles(row["smiles"])
            v = float(row["dfH0_0K"])
            values[can] = v * KCAL if args.kcal else v

    tables = ATcT(version=args.version)
    print(f"ATcT ver. {tables.version}; energies from {args.energies}\n")
    header = (
        f"{'target':14s} {'DrH(0K)':>9s} {'DfH CBH':>9s} {'ATcT':>9s} {'drift':>7s} "
        f"{'±ref':>6s} {'±drift':>6s} method"
    )
    print(header)
    for can, own_value in values.items():
        try:
            own = tables.lookup(can)
        except Exception:
            continue
        try:
            rxn = cbh.cbh_reactions(can, rungs=[args.rung], protect=args.protect)[
                args.rung
            ]
        except cbh.RungUnavailable as e:
            print(f"{can:14s} {e}")
            continue
        if rxn.trivial:
            continue
        missing = [s for s in rxn.coefficients if s not in values]
        if missing:
            print(f"{can:14s} no energy for {missing}")
            continue
        drH = sum(c * values[s] for s, c in rxn.coefficients.items())
        try:
            refs = {s: tables.lookup(s) for s in rxn.references}
        except Exception as e:
            print(f"{can:14s} ATcT: {e}")
            continue
        ref = sum(c * refs[s]["dfH0_0K"] for s, c in rxn.references.items())
        dfH = ref - drH
        ids = [refs[s]["id"] for s in rxn.references]
        cov = tables.covariance(ids, fetch=False)
        u, used = reaction_uncertainty(
            list(rxn.references.values()),
            [refs[s]["uncertainty_298K"] or 0.0 for s in rxn.references],
            cov,
        )
        # The drift = references - target: its uncertainty with the target
        # included, so that the correlation with the target is counted
        udrift, used2 = reaction_uncertainty(
            list(rxn.references.values()) + [-1.0],
            [refs[s]["uncertainty_298K"] or 0.0 for s in rxn.references]
            + [own["uncertainty_298K"] or 0.0],
            tables.covariance(ids + [own["id"]], fetch=False),
        )
        method = "covariance" if used else "quadrature"
        if used != used2:
            method += " (±drift " + ("covariance" if used2 else "quadrature") + ")"
        print(
            f"{can:14s} {drH:9.2f} {dfH:9.2f} {own['dfH0_0K']:9.2f} "
            f"{dfH - own['dfH0_0K']:+7.2f} {u:6.2f} {udrift:6.2f} {method}   {rxn}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
