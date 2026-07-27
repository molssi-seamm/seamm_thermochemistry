# -*- coding: utf-8 -*-

"""Installer for seamm_thermochemistry: fetches the shared atomic
reference-energy database from its Zenodo record.

This package has no external executable and no conda environment -- the
only "installation" step is downloading `thermochemistry.db`, which is not
shipped in the pip package (see the design doc at
~/Sites/reference-energy/2026-07-24_reference-energy/: the database is
regenerated from external master files and will keep growing). It lives
outside the package instead, in a configured SEAMM-root directory --
the same pattern used for DFTB+'s Slater-Koster parameter sets
(~/SEAMM/Parameters/slako) -- registered in seamm.ini's [thermochemistry]
section as `database-path`.

Because there is no executable/environment to manage, this installer does
NOT call `super().check()` / `super().install()` -- those are entirely
about the executables/conda-environment machinery `InstallerBase` provides
for a typical plug-in, which doesn't apply here. `check`/`install`/
`update`/`uninstall` are implemented directly instead.
"""

import hashlib
import logging
from pathlib import Path

import seamm_installer
from seamm_util import Zenodo

logger = logging.getLogger(__name__)


class Installer(seamm_installer.InstallerBase):
    """Install/update/remove the seamm_thermochemistry reference database."""

    # https://zenodo.org/records/21612188 (DOI 10.5281/zenodo.21612188),
    # concept DOI 10.5281/zenodo.21612187 -- this id resolves to whichever
    # version is newest via get_latest_public_record, so it never needs to
    # change when a new version of thermochemistry.db is published.
    zenodo_concept_id = 21612187

    database_filename = "thermochemistry.db"

    def __init__(self, logger=logger):
        super().__init__(logger=logger)
        logger.debug("Initializing the seamm_thermochemistry installer.")

        self.section = "thermochemistry"
        # No executable, no conda environment: this package is data-only.
        self.executables = []
        self.environment = None
        self.environment_file = None

    def _configured_database_path(self):
        """The configured database path, or the default install location."""
        data = self.configuration.get_values(self.section)
        if "database-path" in data and data["database-path"] != "":
            return Path(data["database-path"]).expanduser().resolve()
        return self.root / "Parameters" / "thermochemistry" / self.database_filename

    def check(self):
        """Check that the reference database is installed, offering to
        install it if not.

        Returns
        -------
        bool
            True if the database is present (installing it first if asked
            and needed), False if it's missing and installation was
            declined.
        """
        path = self._configured_database_path()

        if path.exists():
            return True

        if self.options.yes or self.ask_yes_no(
            "The thermochemistry reference database is not installed at "
            f"'{path}'.\nDownload it from Zenodo now?",
            default="yes",
        ):
            self.install()
            return True

        return False

    def install(self):
        """Download the reference database from Zenodo and register its
        location in seamm.ini's [thermochemistry] section."""
        if self.zenodo_concept_id is None:
            raise RuntimeError(
                "seamm_thermochemistry has no Zenodo record configured yet "
                "(Installer.zenodo_concept_id is None) -- nothing to download."
            )

        path = self._configured_database_path()
        path.parent.mkdir(parents=True, exist_ok=True)

        print("Getting the thermochemistry reference database from Zenodo...")
        zenodo = Zenodo()
        record = zenodo.get_latest_public_record(self.zenodo_concept_id)
        record.download_file(self.database_filename, path)

        self._verify_checksum(record, path)

        print(f"Done! Installed to {path}.")

        if not self.configuration.section_exists(self.section):
            self.configuration.add_section(self.section)
        self.configuration.set_value(self.section, "database-path", str(path))
        self.configuration.save()

    def update(self):
        """Re-download the reference database (a newer Zenodo version may
        be out -- get_latest_public_record always resolves to the current
        one, so this is just install() again)."""
        path = self._configured_database_path()
        if not path.exists():
            print(
                f"The thermochemistry database is not installed at {path}; "
                "installing it fresh instead of updating."
            )
        self.install()

    def uninstall(self):
        """Remove the installed reference database and clear the config."""
        data = self.configuration.get_values(self.section)
        if "database-path" not in data or data["database-path"] == "":
            print("The thermochemistry database is not installed; nothing to do.")
            return

        path = Path(data["database-path"]).expanduser().resolve()
        if path.exists():
            print(f"Deleting the thermochemistry database at {path}.")
            path.unlink()

        self.configuration.set_value(self.section, "database-path", "")
        self.configuration.save()
        print("Done!")

    def _verify_checksum(self, record, path):
        """Verify the downloaded file's md5 against Zenodo's manifest.

        Zenodo's public records API reports a "checksum" (format
        "md5:<hex>") per file -- confirmed against a live record. A silent
        corruption in a shared reference-energy dataset is exactly the
        kind of error that should be loud, not discovered downstream in
        someone's formation-energy numbers.
        """
        expected = None
        for entry in record["files"]:
            if entry.get("key") == self.database_filename:
                checksum = entry.get("checksum", "")
                if checksum.startswith("md5:"):
                    expected = checksum[len("md5:") :]
                break

        if expected is None:
            logger.warning(
                "No md5 checksum found in the Zenodo manifest; skipping "
                "integrity check."
            )
            return

        digest = hashlib.md5(path.read_bytes()).hexdigest()
        if digest != expected:
            path.unlink(missing_ok=True)
            raise RuntimeError(
                f"Checksum mismatch downloading {self.database_filename}: "
                f"expected {expected}, got {digest}. The download was "
                "removed; please try again."
            )
