#!/usr/bin/env python3
"""Publish seamm_thermochemistry's reference database to Zenodo.

Creates a new Zenodo deposition, uploads data/thermochemistry.db, sets
metadata, and (only with --publish) publishes it -- registering a DOI.
Run with --sandbox first against sandbox.zenodo.org: no real DOI, fully
reversible, the right way to validate this before ever touching
production Zenodo. Requires a [ZENODO]/[SANDBOX] token in
~/.seamm.d/seammrc (see seamm_util.zenodo.Zenodo).

Metadata/creator/license conventions mirror the existing SEAMM_packages.json
Zenodo record (id 7789854) so this fits the project's established house
style rather than inventing a new one.
"""

import argparse

from seamm_thermochemistry import DEFAULT_DB_PATH
from seamm_util import Zenodo

TITLE = "SEAMM Thermochemistry Reference Database"
DESCRIPTION = (
    "<p>Shared atomic reference-energy database for SEAMM (Simulation "
    "Environment for Atomistic and Molecular Simulations). Powers "
    "physically meaningful energy/enthalpy-of-formation reporting in "
    "SEAMM's electronic-structure and periodic-DFT plug-ins (gaussian_step, "
    "psi4_step, vasp_step), replacing each plug-in's own arbitrary "
    "total-energy zero with a reference to the elements in their standard "
    "states.</p>"
    "<p>This version: experimental 0&nbsp;K/298&nbsp;K reference data for "
    "102 elements; computed atomic reference energies for the full "
    "Gaussian (5192 composite-method/basis combinations) and Psi4 (4674) "
    "grids, up to 36 elements each, and VASP (PBE, PBE-D3BJ, r2SCAN, "
    "r2SCAN-D3BJ @700&nbsp;eV, both isolated-atom and standard-state-phase "
    "conventions) across up to 94 elements. Coverage is expected to grow in "
    "future versions.</p>"
)
KEYWORDS = ["SEAMM", "thermochemistry", "reference energies", "formation energy"]
CREATORS = [
    {
        "name": "Saxe, Paul",
        "affiliation": "MolSSI, Virginia Tech",
        "orcid": "0000-0002-8641-9448",
    }
]
LICENSE = "cc-by-4.0"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sandbox", action="store_true", help="Use sandbox.zenodo.org, not production"
    )
    parser.add_argument(
        "--publish",
        action="store_true",
        help="Actually publish (registers the DOI in production). Without "
        "this, a draft is created/updated for review on Zenodo.",
    )
    args = parser.parse_args()

    zenodo = Zenodo(use_sandbox=args.sandbox)
    record = zenodo.create_record()

    record.title = TITLE
    record.description = DESCRIPTION
    record.creators = CREATORS
    record.upload_type = "dataset"
    record.keywords = KEYWORDS
    record.metadata["license"] = LICENSE

    record.add_file(DEFAULT_DB_PATH, binary=True)

    if args.publish:
        record.publish()
        print(f"Published: {record.doi}")
        print(
            f"Concept id (use with get_latest_public_record): {record['conceptrecid']}"
        )
    else:
        record.update_metadata()
        print(f"Draft created (not published): {record.data['links']['html']}")
        print("Review it on Zenodo, then re-run with --publish, or publish by hand.")


if __name__ == "__main__":
    main()
