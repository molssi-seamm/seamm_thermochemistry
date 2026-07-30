#!/usr/bin/env python3
"""Publish or update seamm_thermochemistry's reference database on Zenodo.

Two modes, and you must pick one explicitly (no silent default, since
picking the wrong one either fragments the dataset into an unrelated
second deposit or clobbers the wrong thing):

    --new-deposit   First-ever release. Creates a brand-new Zenodo record.
    --new-version   Every release after that. Creates a new *version* of
                    the existing record (CONCEPT_ID below) -- same concept
                    DOI, new version DOI, and get_latest_public_record()
                    keeps resolving to whichever is newest automatically.

Only uploads/updates data/thermochemistry.db, sets metadata, and (only
with --publish) publishes it -- registering the DOI. Without --publish, a
draft is created/updated for you to review on Zenodo first.

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

# The published record's concept id -- stable across all versions. Set once
# the first version exists; https://zenodo.org/records/21612188, concept
# DOI 10.5281/zenodo.21612187.
CONCEPT_ID = 21612187

DATABASE_FILENAME = "thermochemistry.db"

TITLE = "SEAMM Thermochemistry Reference Database"
DESCRIPTION = (
    "<p>Shared atomic reference-energy database for SEAMM (Simulation "
    "Environment for Atomistic and Molecular Simulations). Powers "
    "physically meaningful energy/enthalpy-of-formation reporting in "
    "SEAMM's electronic-structure and periodic-DFT plug-ins (gaussian_step, "
    "psi4_step, vasp_step, orca_step), replacing each plug-in's own "
    "arbitrary total-energy zero with a reference to the elements in their "
    "standard states.</p>"
    "<p>This version: experimental 0&nbsp;K/298&nbsp;K reference data for "
    "102 elements; computed atomic reference energies for the full "
    "Gaussian (5192 composite-method/basis combinations) and Psi4 (4674) "
    "grids, up to 36 elements each; VASP (PBE, PBE-D3BJ, r2SCAN, "
    "r2SCAN-D3BJ @700&nbsp;eV, both isolated-atom and standard-state-phase "
    "conventions) across up to 94 elements; and ORCA (REVDSD-PBEP86-D4/2021, "
    "PWLDA, VWN, VWN3, each across the 16-basis def2 family: SV(P)/SVP/"
    "TZVP/TZVP(-f)/QZVPP/SVPD/TZVPD/TZVPPD/QZVPD/QZVPPD and the ma-def2 "
    "augmented variants) across up to 85 elements, each value vetted "
    "against its expected spin before inclusion. Coverage is expected to "
    "keep growing in future versions.</p>"
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


def _set_metadata(record):
    record.title = TITLE
    record.description = DESCRIPTION
    record.creators = CREATORS
    record.upload_type = "dataset"
    record.keywords = KEYWORDS
    record.metadata["license"] = LICENSE


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--new-deposit",
        action="store_true",
        help="First-ever release: create a brand-new Zenodo record.",
    )
    mode.add_argument(
        "--new-version",
        action="store_true",
        help=f"Create a new version of the existing record (concept id {CONCEPT_ID}).",
    )
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

    if args.new_version:
        latest = zenodo.get_latest_public_record(CONCEPT_ID)
        record = zenodo.add_version(latest["id"])
        # The new draft starts as a copy of the previous version, files
        # included -- remove the old database before adding the refreshed
        # one, or Zenodo ends up with two files of the same name.
        if DATABASE_FILENAME in record.files():
            record.remove_file(DATABASE_FILENAME)
    else:
        record = zenodo.create_record()

    _set_metadata(record)
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
