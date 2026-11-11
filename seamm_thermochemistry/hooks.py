# -*- coding: utf-8 -*-

"""The JobServer's ``org.molssi.seamm.jobserver.job_finished`` hook: harvest a
finished job's CBH species into the host's fragment store.

The Thermochemistry step writes the species it computes into the fragment store
itself when the store is on a local disk of the machine running the job
(write-back). When it cannot -- a job that ran on a cluster, or a store on a
network filesystem -- it leaves the marker ``harvest_fragments.json`` in the
job's directory. When the job comes back and the JobServer finalizes it, this
hook harvests such a job into the store on the JobServer's host, once: the
marker is then renamed ``harvest_fragments.done.json``, with the tallies.

The hook never changes the job or its status, and any error is logged, not
raised.
"""

import json
import logging
from pathlib import Path
import shutil
import sqlite3
import tempfile

logger = logging.getLogger(__name__)

#: The marker the step leaves in the job's directory
MARKER = "harvest_fragments.json"
DONE = "harvest_fragments.done.json"


def write_marker(job_directory, species, reason):
    """Add species to the job's marker (the step calls this when it cannot write
    the store). Never raises."""
    try:
        path = Path(job_directory) / MARKER
        data = {"reason": reason, "species": []}
        if path.exists():
            data = json.loads(path.read_text())
        data["reason"] = reason
        for s in species:
            if s not in data["species"]:
                data["species"].append(s)
        partial = path.with_suffix(".json.partial")
        partial.write_text(json.dumps(data, indent=1))
        partial.replace(path)
    except Exception as e:  # pragma: no cover
        logger.warning(f"Could not write the fragment-store marker: {e}")


def harvest_finished_job(job_id=None, wdir=None, state=None, logger=logger, **kwargs):
    """Harvest a finished job's species into this host's fragment store if the
    job asked for it (the marker), whatever its final state: the species of a
    job's finished Loop iterations are valid even if a later step failed.

    Returns the tallies, or None when there was nothing to do or it failed.
    """
    try:
        if wdir is None:
            return None
        wdir = Path(wdir)
        marker = wdir / MARKER
        if not marker.exists() or not (wdir / "seamm.db").exists():
            return None
        from .fragment_store import FragmentStore, harvest_job, store_is_writable

        writable, why = store_is_writable()
        if not writable:
            logger.info(
                f"Job {job_id}: not harvesting its fragments here: the fragment "
                f"store is {why}."
            )
            return None
        # Harvest a copy: the job's database may be in WAL mode, and the
        # dashboard may be reading it
        with tempfile.TemporaryDirectory(prefix="harvest-") as tmp:
            tmp = Path(tmp)
            source = sqlite3.connect(f"file:{wdir / 'seamm.db'}?mode=ro", uri=True)
            target = sqlite3.connect(tmp / "seamm.db")
            try:
                source.backup(target)
            finally:
                target.close()
                source.close()
            if (wdir / "job_data.json").exists():
                shutil.copy2(wdir / "job_data.json", tmp / "job_data.json")
            with FragmentStore(create=True) as store:
                counts = harvest_job(tmp, store)
        # The job is identified by its own directory, not the copy
        summary = {
            "job": job_id,
            "state": state,
            "store": str(FragmentStore().path),
            "inserted": counts["inserted"],
            "duplicate": counts["duplicate"],
            "logged": counts["logged"],
            "geometry": counts["geometry"],
            "energy": counts["energy"],
            "skipped": counts["skipped"],
        }
        (wdir / DONE).write_text(json.dumps(summary, indent=1))
        marker.unlink()
        logger.info(
            f"Job {job_id}: harvested its fragments into the store: geometry "
            f"{counts['geometry']}, energy {counts['energy']}"
            + (f", {len(counts['skipped'])} skipped" if counts["skipped"] else "")
        )
        return counts
    except Exception as e:
        logger.warning(f"Job {job_id}: harvesting its fragments failed: {e}")
        return None
