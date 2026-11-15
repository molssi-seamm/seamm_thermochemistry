"""Species whose optimization did not converge never reach the store (Paul,
2026-11-15): the Thermochemistry step records "optimization converged" on each
species, the harvest refuses "no", and the JobServer hook harvests only the
species the job listed as storable. Results from before the record existed are
taken, and counted."""

import json

import pytest

pytest.importorskip("rdkit")
molsystem = pytest.importorskip("molsystem")

from seamm_thermochemistry import FragmentStore, harvest_configuration  # noqa: E402
from seamm_thermochemistry import fragment_store as F  # noqa: E402
from seamm_thermochemistry import hooks  # noqa: E402

from tests.test_writeback import (  # noqa: E402,F401
    LEVEL,
    _embedded,
    _species,
    db,
    returned_job,
    store,
)


def _converged(configuration, value):
    name = f"{F.CONVERGED_PROPERTY}#Thermochemistry#{LEVEL}"
    configuration.properties.add(name, "str", noerror=True)
    configuration.properties.put(name, value)
    return configuration


def _harvest(target, configuration, smiles, **kwargs):
    return harvest_configuration(
        target,
        configuration,
        smiles,
        energy="total energy#Thermochemistry#*",
        job="9:/jobs/9",
        **kwargs,
    )


def test_an_unconverged_species_is_refused(db, store):  # noqa: F811
    c = _converged(_species(db, "CO", _embedded("CO")[1], -300000.0), "no")
    counts = _harvest(store, c, "CO")
    assert counts["geometry"] == {} and counts["energy"] == {}
    assert counts["inserted"] == 0
    assert "its optimization did not converge" in counts["skipped"][0]


@pytest.mark.parametrize("value", ["yes", "not optimized"])
def test_a_converged_species_is_stored(db, store, value):  # noqa: F811
    c = _converged(_species(db, "CO", _embedded("CO")[1], -300000.0), value)
    counts = _harvest(store, c, "CO")
    assert counts["geometry"] == {"inserted": 1}
    assert counts["no convergence record"] == 0


def test_old_results_are_taken_and_counted(db, store):  # noqa: F811
    c = _species(db, "CO", _embedded("CO")[1], -300000.0)
    counts = _harvest(store, c, "CO")
    assert counts["geometry"] == {"inserted": 1}
    assert counts["no convergence record"] == 1
    # unless a record is required
    counts = _harvest(
        FragmentStore(store.path.parent / "other.db", create=True),
        c,
        "CO",
        require_converged=True,
    )
    assert counts["geometry"] == {}
    assert "no record of whether its optimization converged" in counts["skipped"][0]


@pytest.fixture()
def remote_job(returned_job):  # noqa: F811
    """A returned job with CO (converged) and CF (not converged)."""
    job, path = returned_job
    job_db = molsystem.SystemDB(filename=f"file:{job / 'seamm.db'}")
    (co,) = [s for s in job_db.get_systems() if s.name == "CBH:CO"]
    _converged(co.configuration, "yes")
    _converged(_species(job_db, "CF", _embedded("CF")[1], -350000.0), "no")
    job_db.close()
    return job, path


def _stored(path):
    with FragmentStore(path) as fragments:
        rows = fragments._connect().execute("SELECT DISTINCT species FROM geometry")
        return {r[0] for r in rows}


def test_the_hook_harvests_only_the_listed_species(remote_job):
    """The step lists only the species it may store; the hook takes those and
    nothing else, even a converged one (CF would be refused anyway)."""
    job, path = remote_job
    hooks.write_marker(job, ["CO"], "the store is on gpfs", convergence_recorded=True)
    counts = hooks.harvest_finished_job(job_id=42, wdir=job, state="finished")
    assert _stored(path) == {"CO"}
    assert any("not in the job's list" in s for s in counts["skipped"])
    done = json.loads((job / hooks.DONE).read_text())
    assert done["no convergence record"] == 0


def test_an_unconverged_fragment_from_a_remote_job_is_refused(remote_job):
    """Even when listed (an older step listed every species), CF is refused."""
    job, path = remote_job
    hooks.write_marker(job, ["CO", "CF"], "the store is on gpfs")
    counts = hooks.harvest_finished_job(job_id=42, wdir=job, state="finished")
    assert _stored(path) == {"CO"}
    assert any(
        "CBH:CF: refused, its optimization did not converge" in s
        for s in counts["skipped"]
    )


def test_a_marker_without_a_list_harvests_all_with_the_checks(remote_job):
    job, path = remote_job
    (job / hooks.MARKER).write_text(json.dumps({"reason": "old step"}))
    hooks.harvest_finished_job(job_id=42, wdir=job, state="finished")
    assert _stored(path) == {"CO"}


def test_the_step_records_convergence_so_no_record_is_refused(
    returned_job,  # noqa: F811
):
    """CO from a step that records convergence, but without the record: refused."""
    job, path = returned_job
    hooks.write_marker(job, ["CO"], "x", convergence_recorded=True)
    counts = hooks.harvest_finished_job(job_id=42, wdir=job, state="finished")
    assert not path.exists() or _stored(path) == set()
    assert any("no record" in s for s in counts["skipped"])
