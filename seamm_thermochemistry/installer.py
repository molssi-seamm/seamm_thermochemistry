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

        # --force for update and uninstall: replace or delete a database with
        # local changes
        for command in ("update", "uninstall"):
            if command in self.subparser:
                self.subparser[command].add_argument(
                    "--force",
                    action="store_true",
                    help="Replace (update) or delete (uninstall) the database even "
                    "if it has local changes; update keeps a backup.",
                )

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

    def _shared_report(self):
        """In an installation that shares the default installation's codes and
        data, never download: report the database it will use. Returns True if
        this installation shares (and so nothing else should be done)."""
        if not getattr(self, "shared_codes", False):
            return False
        path = self._configured_database_path()
        if path.exists():
            print(f"    This installation uses the shared reference database {path}.")
        else:
            print(
                "    This installation shares the default installation's data, which "
                f"has no reference database ({path}). Install it there first."
            )
        return True

    # ---- the published versions and the local copy ---------------------------

    @staticmethod
    def _md5(path):
        return hashlib.md5(Path(path).read_bytes()).hexdigest()

    def _expected_md5(self, record):
        """The md5 Zenodo lists for the database in this record, or None."""
        for entry in record["files"]:
            if entry.get("key") == self.database_filename:
                checksum = entry.get("checksum", "")
                if checksum.startswith("md5:"):
                    return checksum[len("md5:") :]
        return None

    def _published_md5s(self, record):
        """The md5 of the database in every published version, as a set.

        Asks Zenodo for all versions of the record; if that fails, just the latest.
        """
        result = set()
        latest = self._expected_md5(record)
        if latest:
            result.add(latest)
        try:
            import requests

            url = f"https://zenodo.org/api/records/{record['id']}/versions"
            response = requests.get(
                url, params={"allversions": "true", "size": 25}, timeout=30
            )
            if response.status_code == 200:
                for hit in response.json().get("hits", {}).get("hits", []):
                    for entry in hit.get("files", []):
                        checksum = entry.get("checksum", "")
                        if entry.get("key") == self.database_filename and (
                            checksum.startswith("md5:")
                        ):
                            result.add(checksum[len("md5:") :])
        except Exception as e:
            logger.debug(f"Could not list the published versions: {e}")
        return result

    def _local_state(self, path, record):
        """ "published" if the local file is a published version (or the copy
        this installer recorded), else "changed"."""
        local = self._md5(path)
        recorded = self.configuration.get_values(self.section).get("database-md5", "")
        if local == recorded or local in self._published_md5s(record):
            return "published"
        return "changed"

    def _latest_record(self):
        if self.zenodo_concept_id is None:
            raise RuntimeError(
                "seamm_thermochemistry has no Zenodo record configured yet "
                "(Installer.zenodo_concept_id is None) -- nothing to download."
            )
        return Zenodo().get_latest_public_record(self.zenodo_concept_id)

    def _download(self, record, path):
        """Download the latest database to `path`, safely.

        The file is downloaded and verified beside `path`, the existing database
        (if any) is backed up, and only then is the new file moved into place, so a
        failed download never touches the database.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name + ".download")
        partial.unlink(missing_ok=True)
        print("Getting the thermochemistry reference database from Zenodo...")
        record.download_file(self.database_filename, partial)
        self._verify_checksum(record, partial)
        if path.exists():
            from datetime import datetime

            stamp = datetime.now().strftime("%Y-%m-%d-%H%M%S")
            backup = path.with_name(f"{path.name}.bak-{stamp}")
            path.replace(backup)
            print(f"    The previous database was saved as {backup.name}.")
        partial.replace(path)
        self._record(path)
        print(f"Done! Installed to {path}.")

    def _register_path(self, path):
        """Register the database's location (only)."""
        if not self.configuration.section_exists(self.section):
            self.configuration.add_section(self.section)
        self.configuration.set_value(self.section, "database-path", str(path))
        self.configuration.save()

    def _record(self, path):
        """Register the database's location and the md5 of the copy just installed.

        Only for a file this installer downloaded: the md5 marks it as unmodified,
        so a later update may replace it.
        """
        self._register_path(path)
        self.configuration.set_value(self.section, "database-md5", self._md5(path))
        self.configuration.save()

    # ---- the installer's commands ----------------------------------------------

    def install(self):
        """Download the reference database from Zenodo and register its
        location in seamm.ini's [thermochemistry] section.

        An existing database is never overwritten here; see `update`.
        """
        if self._shared_report():
            return
        path = self._configured_database_path()
        if path.exists():
            print(f"The thermochemistry database is already installed at {path}.")
            self._register_path(path)  # make sure it is used; its md5 is not recorded
            self.update()
            return
        self._download(self._latest_record(), path)

    def update(self):
        """Bring the database up to the latest published version -- unless it has
        local changes.

        The database is replaced only if the local file is exactly one of the
        published versions (or the copy this installer installed). A database that
        has been edited, e.g. with curated reference energies, is kept as it is and
        reported; ``update --force`` replaces it, keeping a backup.
        """
        if self._shared_report():
            return
        path = self._configured_database_path()
        record = self._latest_record()
        if not path.exists():
            print(
                f"The thermochemistry database is not installed at {path}; "
                "installing it fresh instead of updating."
            )
            self._download(record, path)
            return
        latest = self._expected_md5(record)
        if latest is not None and self._md5(path) == latest:
            print(f"The thermochemistry database at {path} is up to date.")
            return
        force = bool(getattr(self.options, "force", False)) if self.options else False
        if self._local_state(path, record) == "changed" and not force:
            print(
                f"The thermochemistry database at {path} has local changes, so it "
                "was not replaced by the published version. To replace it anyway "
                "(keeping a backup), run 'seamm-thermochemistry-installer update "
                "--force'."
            )
            return
        self._download(record, path)

    def uninstall(self):
        """Remove the installed reference database and clear the config.

        A database with local changes is left in place (only the configuration is
        cleared) unless ``uninstall --force`` is given.
        """
        if self._shared_report():
            return
        data = self.configuration.get_values(self.section)
        if "database-path" not in data or data["database-path"] == "":
            print("The thermochemistry database is not installed; nothing to do.")
            return

        path = Path(data["database-path"]).expanduser().resolve()
        force = bool(getattr(self.options, "force", False)) if self.options else False
        if path.exists():
            keep = False
            if not force:
                try:
                    keep = self._local_state(path, self._latest_record()) == "changed"
                except Exception:
                    keep = True  # cannot tell: keep it
            if keep:
                print(
                    f"The thermochemistry database at {path} has local changes, so it "
                    "was left in place; 'uninstall --force' deletes it."
                )
            else:
                print(f"Deleting the thermochemistry database at {path}.")
                path.unlink()

        self.configuration.set_value(self.section, "database-path", "")
        self.configuration.set_value(self.section, "database-md5", "")
        self.configuration.save()
        print("Done!")

    def _verify_checksum(self, record, path):
        """Verify the downloaded file's md5 against Zenodo's manifest.

        Zenodo's public records API reports a "checksum" (format "md5:<hex>") per
        file. `path` is the freshly downloaded file, which is removed on a
        mismatch; the installed database is never touched here.
        """
        expected = self._expected_md5(record)
        if expected is None:
            logger.warning(
                "No md5 checksum found in the Zenodo manifest; skipping "
                "integrity check."
            )
            return

        digest = self._md5(path)
        if digest != expected:
            path.unlink(missing_ok=True)
            raise RuntimeError(
                f"Checksum mismatch downloading {self.database_filename}: "
                f"expected {expected}, got {digest}. The download was "
                "removed; please try again."
            )
