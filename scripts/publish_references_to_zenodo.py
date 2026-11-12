#!/usr/bin/env python3
"""Publish the official reference library (references-official.db) on Zenodo.

The library holds enthalpies of formation that ATcT lacks, for the reference
species of SEAMM's CBH thermochemistry; see
``seamm_thermochemistry.references``. It is built with the maintainer commands
of ``seamm-thermochemistry-references`` (add-curated, webbook import,
accept-file, offsets, check) and published here, in two steps:

    draft    Create the Zenodo deposit (``--new-deposit`` the first time,
             ``--new-version`` after that), stamp its prereserved DOI and the
             version into the file's meta, write the build report beside it,
             upload both and set the metadata. Nothing is public yet: review
             the draft on Zenodo.
    publish  Publish a reviewed draft (``--deposit ID``), registering its DOI.

Use ``--sandbox`` first (sandbox.zenodo.org, no real DOI). The tokens are the
[ZENODO] and [SANDBOX] sections of ~/.seamm.d/seammrc (seamm_util.Zenodo).

The file is stamped in place: run ``draft`` on a copy built for publication.
A sandbox run stamps a sandbox DOI; the real draft overwrites it.
"""

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from seamm_thermochemistry.references import SCHEMA_VERSION, ReferenceLibrary
from seamm_util import Zenodo

#: The published record's concept id, stable across versions: None until the
#: first version exists
CONCEPT_ID = None

FILENAME = "references-official.db"
REPORT = FILENAME + ".report.json"

TITLE = "SEAMM Thermochemistry Reference Library"
KEYWORDS = [
    "SEAMM",
    "thermochemistry",
    "enthalpy of formation",
    "connectivity-based hierarchy",
    "isodesmic reactions",
    "reference data",
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

#: How each source is credited, kept in the file's source table (licence)
CREDITS = {
    "JANAF-1998": "NIST Standard Reference Database 13 (NIST-JANAF Thermochemical "
    "Tables); values from the printed 4th edition, cited by table, with credit to "
    "NIST.",
    "NBS-82": "The NBS Tables of Chemical Thermodynamic Properties (J. Phys. Chem. "
    "Ref. Data 11, Suppl. 2, 1982); values cited by table and page, with credit to "
    "NIST.",
    "WebBook-SRD69": "NIST Chemistry WebBook, NIST Standard Reference Database 69; "
    "each value cites its original measurement, with credit to NIST.",
    "CODATA-1989": "CODATA Key Values for Thermodynamics (1989); published data, "
    "cited.",
    "D0": "Spectroscopic dissociation energies from the cited literature.",
    "W4-17": "W4-17 (Karton, Sylvetsky and Martin, J. Comput. Chem. 38, 2063, "
    "2017); published data, cited.",
    "derived": "Derived from the cited values of the other sources; each "
    "derivation is in the file.",
}


def description(report):
    sources = "".join(
        f"<li><b>{s['name']}</b>: {s['citation']} ({s['chosen']} values"
        + (
            f"; offset from ATcT {s['offset_kJmol']:.2f} kJ/mol rms over "
            f"{s['offset_n']} species"
            if s["offset_kJmol"] is not None
            else ""
        )
        + ")</li>"
        for s in report["sources"]
    )
    return (
        "<p>Enthalpies of formation of gas-phase species that the Active "
        "Thermochemical Tables (ATcT) lack, for the reference species of the "
        "connectivity-based hierarchy (CBH) thermochemistry in SEAMM (Simulation "
        "Environment for Atomistic and Molecular Simulations). SEAMM's "
        "Thermochemistry step takes each reference species from ATcT first and from "
        "this library only when ATcT lacks it; results that use it are labelled "
        "as mixed.</p>"
        "<p>The library keeps every value each source gives, with its original "
        "citation and its uncertainty as stated; the value chosen per species and "
        "source and the rule that chose it; any decision by a person (who, when "
        "and why); and the derivation of each derived value. Rules: an estimate "
        "is never chosen without a recorded decision; a later source replaces an "
        "earlier one only when three times as precise and in agreement; values "
        "whose measurements disagree beyond three times their combined uncertainty "
        "are held; ions are in ATcT's stationary-electron convention. Each "
        "source's rms offset from ATcT, measured over species both have, is "
        "recorded and added to its uncertainties.</p>"
        f"<p>This version ({report['version']}): {report['species_chosen']} "
        f"species with a chosen value and {report['values']} values in all, from "
        f"these sources, in rank order:</p><ul>{sources}</ul>"
        f"<p>Held pending a computed arbitration: "
        f"{', '.join(h['name'] or h['smiles'] for h in report['held'])}.</p>"
        "<p>The values are reproduced from, and cite, their original sources: "
        "please cite those as well as this record. The NIST-JANAF tables (SRD 13), "
        "the NBS tables and the NIST Chemistry WebBook (SRD 69) are credited to "
        "the National Institute of Standards and Technology.</p>"
        "<p>Files: references-official.db (SQLite; tables source, species, "
        "measurement, reference, derivation, calibration) and "
        "references-official.db.report.json (the build report).</p>"
    )


def build_report(path, version, doi):
    library = ReferenceLibrary(path)
    db = library._connect()
    sources = []
    for row in db.execute("SELECT * FROM source ORDER BY rank"):
        chosen = db.execute(
            "SELECT COUNT(*) FROM reference WHERE source = ? AND state = 'chosen'",
            (row["id"],),
        ).fetchone()[0]
        sources.append(
            {
                "name": row["name"],
                "kind": row["kind"],
                "rank": row["rank"],
                "citation": row["citation"],
                "credit": row["licence"],
                "chosen": chosen,
                "offset_kJmol": row["offset_kj"],
                "offset_n": row["offset_n"],
                "offset_atct": row["offset_atct"],
                "offset_full_kJmol": row["offset_full_kj"],
                "offset_excluded": row["offset_excluded"],
            }
        )
    held = [
        {
            "smiles": r["smiles"],
            "name": r["name"],
            "source": r["source"],
            "reason": r["decision_reason"] or r["rule"],
        }
        for r in db.execute(
            "SELECT sp.smiles, sp.name, s.name AS source, r.decision_reason, r.rule "
            "FROM reference r JOIN species sp ON sp.id = r.species "
            "JOIN source s ON s.id = r.source WHERE r.state = 'held' "
            "ORDER BY sp.smiles"
        )
    ]
    report = {
        "file": FILENAME,
        "version": version,
        "doi": doi,
        "schema": SCHEMA_VERSION,
        "layer": library.layer,
        "created": library.meta.get("created"),
        "species": db.execute("SELECT COUNT(*) FROM species").fetchone()[0],
        "species_chosen": db.execute(
            "SELECT COUNT(DISTINCT species) FROM reference WHERE state = 'chosen'"
        ).fetchone()[0],
        "values": db.execute("SELECT COUNT(*) FROM measurement").fetchone()[0],
        "sources": sources,
        "held": held,
        "check": library.check(),
    }
    library.close()
    return report


def draft(args):
    path = Path(args.library).expanduser().resolve()
    if path.name != FILENAME:
        raise SystemExit(f"The file must be named {FILENAME}, not {path.name}")
    zenodo = Zenodo(use_sandbox=args.sandbox)
    if args.deposit is not None:
        # An unpublished draft made earlier, e.g. by a run that stopped
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
    if library.layer != "official":
        raise SystemExit(f"{path} is a {library.layer} library, not the official one")
    problems = library.check()
    if problems:
        raise SystemExit(f"The library does not check: {problems}")
    library.set_meta("version", args.version)
    library.set_meta("doi", doi)
    db = library._writer()
    for name, credit in CREDITS.items():
        db.execute("UPDATE source SET licence = ? WHERE name = ?", (credit, name))
    db.commit()
    library.close()

    report = build_report(path, args.version, doi)
    data = path.read_bytes()
    report["sha256"] = hashlib.sha256(data).hexdigest()
    report["md5"] = hashlib.md5(data).hexdigest()
    report_path = path.with_name(REPORT)
    report_path.write_text(json.dumps(report, indent=2) + "\n")

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
    print(f"Concept id (the installer's): {record['conceptrecid']}")


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
