#!/usr/bin/env python3
"""Publish the computed reference values (references-computed.db) on Zenodo.

The computed layer holds enthalpies of formation computed with SEAMM's
high-level (tier-A) recipes for CBH reference species that neither ATcT nor
the experimental reference library has; see ``seamm_thermochemistry.computed``.
It is published as its own record, apart from the experimental library (Paul,
2026-11-15), and versioned on its own as values are added.

It is built with ``seamm-thermochemistry-references add-computed VALUES.csv
--library references-computed.db`` and published here, in two steps:

    draft    Create the Zenodo deposit (``--new-deposit`` the first time,
             ``--new-version`` after that), stamp its prereserved DOI and the
             version into the file's meta, write the build report beside it,
             upload both and set the metadata. Nothing is public yet: review
             the draft on Zenodo.
    publish  Publish a reviewed draft (``--deposit ID``), registering its DOI.

The token is the [ZENODO] section of ~/.seamm.d/seammrc (seamm_util.Zenodo).
The file is stamped in place: run ``draft`` on a copy built for publication.
"""

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from seamm_thermochemistry.computed import LAYER, SOURCE
from seamm_thermochemistry.references import SCHEMA_VERSION, ReferenceLibrary
from seamm_util import Zenodo

#: The published record's concept id, stable across versions; None until the
#: first version is published
CONCEPT_ID = None
#: The experimental reference library's concept DOI
LIBRARY_CONCEPT_DOI = "10.5281/zenodo.24175073"

FILENAME = "references-computed.db"
REPORT = FILENAME + ".report.json"

TITLE = "SEAMM Thermochemistry Computed Reference Values (enthalpies of formation)"
KEYWORDS = [
    "SEAMM",
    "thermochemistry",
    "enthalpy of formation",
    "connectivity-based hierarchy",
    "CCSD(T)-F12",
    "computed reference data",
]
CREATORS = [
    {
        "name": "Saxe, Paul",
        "affiliation": "MolSSI, Virginia Tech",
        "orcid": "0000-0002-8641-9448",
    }
]
LICENSE = "cc-by-4.0"
RELATED = [
    {
        "identifier": LIBRARY_CONCEPT_DOI,
        "relation": "isSupplementTo",
        "resource_type": "dataset",
    },
    {
        "identifier": "https://github.com/molssi-seamm/seamm_thermochemistry",
        "relation": "isSupplementedBy",
        "resource_type": "software",
    },
    {
        "identifier": "https://github.com/molssi-seamm/thermochemistry_step",
        "relation": "isSupplementedBy",
        "resource_type": "software",
    },
]


def description(report):
    rows = "".join(
        f"<li><b>{v['name'] or v['smiles']}</b> ({v['smiles']}): "
        f"ΔfH°(298.15 K) = {v['dfH298_kJmol']:.2f} ± {v['uncertainty_kJmol']:.1f} "
        f"kJ/mol (95 %); {(v['route'] or '').removeprefix('none: ')}</li>"
        for v in report["values"]
    )
    recipes = ", ".join(sorted({v["recipe"] for v in report["values"]}))
    return (
        "<p>Gas-phase enthalpies of formation <b>computed</b> with the high-level "
        "(tier-A) recipes of SEAMM (Simulation Environment for Atomistic and "
        "Molecular Simulations), for reference species of the connectivity-based "
        "hierarchy (CBH) thermochemistry that neither the Active Thermochemical "
        "Tables (ATcT) nor the experimental SEAMM Thermochemistry Reference "
        f"Library (doi:{LIBRARY_CONCEPT_DOI}) has. They are kept apart from the "
        "experimental values.</p>"
        "<p>How SEAMM uses them: only when computed references are allowed (the "
        "default), and always after every experimental value, so a species with an "
        "experimental value never takes a computed one. Results that use one are "
        "labelled as using computed references. The values are references for "
        "other species' reactions, not benchmarks or results for the species "
        "themselves.</p>"
        "<p>Each value records its recipe, the SEAMM job and iteration that "
        "computed it, the geometry level and conformer, how its H298 - H0 was "
        "obtained, the route to the enthalpy of formation, and how its uncertainty "
        "was set (from the recipe, before the value was known). ΔfH°(298.15 K) is "
        "primary; ΔfH°(0 K) is derived from it with the species' own H298 - H0 "
        "and the CODATA enthalpies of the elements.</p>"
        f"<p>This version ({report['version']}): {len(report['values'])} species, "
        f"computed with {recipes}:</p><ul>{rows}</ul>"
        "<p>Files: references-computed.db (SQLite, the layout of the reference "
        "library: tables source, species, measurement, reference; each value's "
        "provenance is in measurement.notes as JSON) and "
        "references-computed.db.report.json (the build report).</p>"
    )


def build_report(path, version, doi):
    library = ReferenceLibrary(path)
    db = library._connect()
    values = []
    for r in db.execute(
        "SELECT sp.smiles, sp.name, sp.formula, sp.charge, sp.multiplicity, m.* "
        "FROM reference r JOIN species sp ON sp.id = r.species "
        "JOIN measurement m ON m.id = r.measurement "
        "WHERE r.state = 'chosen' ORDER BY sp.formula, sp.smiles"
    ):
        provenance = json.loads(r["notes"] or "{}")
        values.append(
            {
                "smiles": r["smiles"],
                "name": r["name"],
                "formula": r["formula"],
                "charge": r["charge"],
                "multiplicity": r["multiplicity"],
                "dfH298_kJmol": r["dfh_298k"],
                "dfH0_kJmol": r["dfh_0k"],
                "uncertainty_kJmol": r["uncertainty"],
                "uncertainty_level": r["unc_level"],
                "H298_H0_kJmol": r["h298_h0"],
                "recipe": r["method"],
                "job": provenance.get("job"),
                "iteration": provenance.get("iteration"),
                "route": provenance.get("atomization"),
                "computed_on": provenance.get("computed_on"),
            }
        )
    held = db.execute(
        "SELECT COUNT(*) FROM measurement WHERE status = 'held'"
    ).fetchone()[0]
    report = {
        "file": FILENAME,
        "version": version,
        "doi": doi,
        "schema": SCHEMA_VERSION,
        "layer": library.layer,
        "created": library.meta.get("created"),
        "source": SOURCE,
        "values": values,
        "replaced_values": held,
    }
    library.close()
    return report


def draft(args):
    path = Path(args.library).expanduser().resolve()
    if path.name != FILENAME:
        raise SystemExit(f"The file must be named {FILENAME}, not {path.name}")
    zenodo = Zenodo(use_sandbox=args.sandbox)
    if args.deposit is not None:
        record = zenodo.get_deposit_record(args.deposit)
        if record.submitted:
            raise SystemExit(f"Deposit {args.deposit} is already published")
        for name in (FILENAME, REPORT):
            if name in record.files():
                record.remove_file(name)
    elif args.new_version:
        if CONCEPT_ID is None:
            raise SystemExit("No CONCEPT_ID yet: use --new-deposit the first time.")
        latest = zenodo.get_latest_public_record(CONCEPT_ID)
        record = zenodo.add_version(latest["id"])
        for name in (FILENAME, REPORT):
            if name in record.files():
                record.remove_file(name)
    else:
        record = zenodo.create_record()
    doi = record.doi

    library = ReferenceLibrary(path, create=True)  # an existing file, written
    if library.layer != LAYER:
        raise SystemExit(f"{path} is a {library.layer} library, not the computed one")
    library.set_meta("version", args.version)
    library.set_meta("doi", doi)
    library.close()

    report = build_report(path, args.version, doi)
    data = path.read_bytes()
    report["sha256"] = hashlib.sha256(data).hexdigest()
    report["md5"] = hashlib.md5(data).hexdigest()
    report_path = path.with_name(REPORT)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")

    record.add_file(path, binary=True)
    record.add_file(report_path)
    record.metadata.update(
        {
            "title": TITLE,
            "upload_type": "dataset",
            "description": description(report),
            "creators": CREATORS,
            "keywords": KEYWORDS,
            "license": LICENSE,
            "version": args.version,
            "related_identifiers": RELATED,
        }
    )
    record.update_metadata()
    print(f"DOI (reserved): {doi}")
    print(f"sha256 {report['sha256']}  md5 {report['md5']}")
    print(f"Draft, deposit id {record['id']}: {record.data['links']['html']}")
    print("Review it on Zenodo, then: publish --deposit ID")


def publish(args):
    zenodo = Zenodo(use_sandbox=args.sandbox)
    record = zenodo.get_deposit_record(args.deposit)
    record.publish()
    print(f"Published: {record.doi}")
    print(f"Concept id (set CONCEPT_ID to it): {record['conceptrecid']}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--sandbox", action="store_true", help="sandbox.zenodo.org")
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("draft", help="create and fill a draft")
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--new-deposit", action="store_true")
    mode.add_argument("--new-version", action="store_true")
    mode.add_argument("--deposit", type=int, help="fill this unpublished draft")
    p.add_argument("--library", required=True, help=f"the {FILENAME} to publish")
    p.add_argument(
        "--version",
        default=f"{datetime.date.today():%Y.%-m.%-d}",
        help="the version (default today's date, YYYY.M.D)",
    )
    p.set_defaults(func=draft)
    p = commands.add_parser("publish", help="publish a reviewed draft")
    p.add_argument("--deposit", required=True, type=int)
    p.set_defaults(func=publish)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
