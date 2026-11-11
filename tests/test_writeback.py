"""Write-back to the fragment store (P5b, see
docs/developer_guide/campaigns/2026-11-10): the structure guard, imaginary
modes refused, new conformers only when lower and structurally distinct, where
the store may be written, and the JobServer hook that harvests returning jobs."""

import json

import numpy as np
import pytest

pytest.importorskip("rdkit")
molsystem = pytest.importorskip("molsystem")

from seamm_thermochemistry import (  # noqa: E402
    FragmentStore,
    check_structure,
    harvest_configuration,
)
from seamm_thermochemistry import fragment_store as F  # noqa: E402
from seamm_thermochemistry import hooks  # noqa: E402

LEVEL = "ORCA:CMP@R2SCAN-3C"
KJ = 2625.4996394799


def _embedded(smiles, dihedral=None, seed=7):
    """(symbols, coordinates) of an RDKit embedding, optionally with the first
    heavy-atom dihedral set (degrees)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol)
    conformer = mol.GetConformer()
    if dihedral is not None:
        rdMolTransforms.SetDihedralDeg(conformer, 0, 1, 2, 3, float(dihedral))
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    return symbols, conformer.GetPositions().tolist()


def _pf6(close=False):
    d = 1.606
    xyz = [[0, 0, 0]] + [
        [s * d if i == k else 0.0 for i in range(3)] for k in range(3) for s in (1, -1)
    ]
    if close:  # two fluorines 0.5 Å apart, as from the old SMILES embedding
        xyz[2] = [xyz[1][0] - 0.1, 0.49, 0.0]
    return ["P"] + ["F"] * 6, xyz


# -- the structure guard -------------------------------------------------------


def test_the_species_itself_passes():
    assert check_structure(*_embedded("CCO"), "CCO", 0, 1) is None
    assert check_structure(*_pf6(), "F[P-](F)(F)(F)(F)F", -1, 1) is None


def test_a_rearranged_structure_is_refused():
    # Dimethyl ether's coordinates filed as ethanol: the same atoms, other bonds
    problem = check_structure(*_embedded("COC"), "CCO", 0, 1)
    assert problem and "not those of CCO" in problem


def test_ligands_on_top_of_each_other_are_refused():
    # Two fluorines 0.5 Å apart read as bonded, so the bonds do not match; a
    # contact just outside bonding distance is caught by the contact check
    problem = check_structure(*_pf6(close=True), "F[P-](F)(F)(F)(F)F", -1, 1)
    assert problem and "not those of" in problem
    symbols, xyz = _pf6()
    # One F turned towards +y: F...F 1.8 Å, too long to read as a bond but
    # below the 2.0 Å floor for two fluorines on the same atom
    angle = np.radians(158.2)
    xyz[2] = [1.606 * np.cos(angle), 1.606 * np.sin(angle), 0.0]
    problem = check_structure(symbols, xyz, "F[P-](F)(F)(F)(F)F", -1, 1)
    assert problem and "too close" in problem


def test_charge_multiplicity_and_atoms_must_fit():
    symbols, xyz = _embedded("[CH3]")
    assert check_structure(symbols, xyz, "[CH3]", 0, 2) is None
    assert "multiplicity" in check_structure(symbols, xyz, "[CH3]", 0, 1)
    assert "charge" in check_structure(*_pf6(), "F[P-](F)(F)(F)(F)F", 0, 1)
    assert "atoms are not" in check_structure(*_embedded("CO"), "CCO", 0, 1)


# -- write-back -----------------------------------------------------------------


@pytest.fixture()
def db(tmp_path):
    db = molsystem.SystemDB(filename=f"file:{tmp_path / 'seamm.db'}")
    yield db
    db.close()


@pytest.fixture()
def store(tmp_path):
    with FragmentStore(tmp_path / "fragments.db", create=True) as s:
        yield s


def _species(db, smiles, coordinates, E, ZPE=200.0, frequencies=(300.0, 1500.0)):
    """A computed species as the Thermochemistry step leaves it."""
    system = db.create_system(name=f"CBH:{smiles}", make_current=True)
    c = system.create_configuration(name="computed", make_current=True)
    c.from_smiles(smiles)
    c.atoms.set_coordinates(coordinates, fractionals=False)
    p = c.properties
    s = f"#Thermochemistry#{LEVEL}"

    def put(name, kind, value, units=None):
        p.add(name, kind, units=units, noerror=True)
        p.put(name, value)

    put("total energy" + s, "float", E / KJ, "E_h")
    put("ZPE" + s, "float", ZPE, "kJ/mol")
    put("H(298.15 K)" + s, "float", E + ZPE + 15.0, "kJ/mol")
    put("settings fingerprint" + s, "str", "fp:r2scan")
    put("frequency fingerprint" + s, "str", "fp:r2scan")
    put("frequencies" + s, "json", list(frequencies))
    put("calculation" + s, "json", {"program": "ORCA", "version": "6.1.1"})
    return c


def _write(store, configuration, smiles):
    return harvest_configuration(
        store,
        configuration,
        smiles,
        energy="total energy#Thermochemistry#*",
        write_back=True,
        job="9:/jobs/9",
    )


def _conformers(store, smiles="CCCC"):
    rows = store.geometries(smiles, frequency_fingerprint="fp:r2scan")
    return sorted((r["conformer"], round(r["E"], 3)) for r in rows)


def test_write_back_rows(db, store):
    counts = _write(store, _species(db, "CO", _embedded("CO")[1], -300000.0), "CO")
    assert counts["geometry"] == {"inserted": 1} and counts["energy"] == {"inserted": 1}
    (g,) = store.geometries("CO", frequency_fingerprint="fp:r2scan")
    assert (g["model_chemistry"], g["program"]) == ("CMP@R2SCAN-3C", "ORCA")
    assert g["job"] == "9:/jobs/9"


def test_imaginary_modes_are_refused(db, store):
    c = _species(db, "CO", _embedded("CO")[1], -300000.0, frequencies=(-55.0, 300.0))
    counts = _write(store, c, "CO")
    assert counts["geometry"] == {} and counts["energy"] == {}
    assert "imaginary mode of 55i" in counts["skipped"][0]
    # Noise below 20i is not an imaginary mode
    c = _species(db, "CO", _embedded("CO")[1], -300000.0, frequencies=(-12.0, 300.0))
    assert _write(store, c, "CO")["geometry"] == {"inserted": 1}


def test_a_rearranged_species_is_not_written(db, store):
    c = _species(db, "CCO", _embedded("CCO")[1], -400000.0)
    c.atoms.set_coordinates(_embedded("COC")[1], fractionals=False)
    counts = _write(store, c, "CCO")
    assert counts["geometry"] == {} and "refused" in counts["skipped"][0]


def test_a_lower_distinct_conformer_is_added(db, store):
    anti = _embedded("CCCC", 180.0)[1]
    gauche = _embedded("CCCC", 65.0)[1]
    # The gauche conformer first, then the anti, 3 kJ/mol lower
    _write(store, _species(db, "CCCC", gauche, -500000.0), "CCCC")
    counts = _write(store, _species(db, "CCCC", anti, -500003.0), "CCCC")
    assert counts["geometry"] == {"inserted": 1}
    labels = _conformers(store)
    assert labels[0] == ("", -500000.0) and labels[1][0].startswith("auto:")
    # The step takes the lowest
    best = min(
        store.geometries("CCCC", frequency_fingerprint="fp:r2scan"),
        key=lambda r: r["E"] + r["ZPE"],
    )
    assert best["E"] == pytest.approx(-500003.0)


def test_the_same_minimum_again_is_only_logged(db, store):
    anti = _embedded("CCCC", 180.0)[1]
    _write(store, _species(db, "CCCC", anti, -500000.0), "CCCC")
    moved = (np.asarray(anti) + 0.001).tolist()  # another hash, the same minimum
    counts = _write(store, _species(db, "CCCC", moved, -500000.12), "CCCC")
    assert counts["geometry"] == {"logged": 1}
    assert len(_conformers(store)) == 1


def test_a_higher_conformer_is_only_logged(db, store):
    _write(store, _species(db, "CCCC", _embedded("CCCC", 180.0)[1], -500003.0), "CCCC")
    counts = _write(
        store, _species(db, "CCCC", _embedded("CCCC", 65.0)[1], -500000.0), "CCCC"
    )
    assert counts["geometry"] == {"logged": 1}
    assert len(_conformers(store)) == 1


def test_conformer_rmsd_is_symmetry_aware():
    symbols, xyz = _embedded("C")
    # Methane with two hydrogens swapped: the same structure
    swapped = [xyz[0], xyz[2], xyz[1], xyz[3], xyz[4]]
    assert F.conformer_rmsd(symbols, xyz, symbols, swapped) == pytest.approx(
        0.0, abs=1e-6
    )


# -- where the store may be written --------------------------------------------


def test_the_store_is_written_only_on_a_local_filesystem(monkeypatch, tmp_path):
    monkeypatch.setattr(F, "filesystem_type", lambda path: ("zfs", "/D3"))
    monkeypatch.setattr("seamm_thermochemistry.db._SEAMM_INI_PATH", tmp_path / "no")
    assert F.store_is_writable(tmp_path / "f.db")[0]
    monkeypatch.setattr(F, "filesystem_type", lambda path: ("gpfs", "/projects"))
    writable, why = F.store_is_writable(tmp_path / "f.db")
    assert not writable and "gpfs" in why
    monkeypatch.setattr(F, "filesystem_type", lambda path: (None, None))
    assert not F.store_is_writable(tmp_path / "f.db")[0]
    # seamm.ini overrides
    ini = tmp_path / "seamm.ini"
    ini.write_text("[thermochemistry]\nfragment-store-writable = yes\n")
    monkeypatch.setattr("seamm_thermochemistry.db._SEAMM_INI_PATH", ini)
    assert F.store_is_writable(tmp_path / "f.db")[0]


def test_filesystem_type_of_a_real_path(tmp_path):
    kind, point = F.filesystem_type(tmp_path)
    assert kind and point and str(tmp_path.resolve()).startswith(point.rstrip("/"))


# -- the JobServer hook ----------------------------------------------------------


@pytest.fixture()
def returned_job(tmp_path, monkeypatch):
    """A job directory as it comes back from a cluster, with the marker, and
    the store redirected to a temporary file."""
    job = tmp_path / "Job_000042"
    job.mkdir()
    db = molsystem.SystemDB(filename=f"file:{job / 'seamm.db'}")
    _species(db, "CO", _embedded("CO")[1], -300000.0)
    db.close()
    (job / "job_data.json").write_text(
        "!MolSSI job_data 1.0\n" + json.dumps({"job id": 42})
    )
    path = tmp_path / "store" / "fragments.db"
    real = F.fragment_store_path
    monkeypatch.setattr(F, "fragment_store_path", lambda p=None: real(p) if p else path)
    monkeypatch.setattr(F, "store_is_writable", lambda p=None: (True, "test"))
    return job, path


def test_the_hook_harvests_a_marked_job_once(returned_job):
    job, path = returned_job
    hooks.write_marker(job, ["CO"], "the store is on gpfs")
    counts = hooks.harvest_finished_job(job_id=42, wdir=job, state="finished")
    assert counts["geometry"] == {"inserted": 1}
    assert not (job / hooks.MARKER).exists()
    done = json.loads((job / hooks.DONE).read_text())
    assert done["job"] == 42 and done["geometry"] == {"inserted": 1}
    with FragmentStore(path) as store:
        (g,) = store.geometries("CO", frequency_fingerprint="fp:r2scan")
        assert g["job"].startswith("42:")
    # Once only: the marker is gone
    assert hooks.harvest_finished_job(job_id=42, wdir=job, state="finished") is None


def test_the_hook_ignores_unmarked_jobs_and_never_raises(returned_job, monkeypatch):
    job, _ = returned_job
    assert hooks.harvest_finished_job(job_id=42, wdir=job, state="finished") is None
    hooks.write_marker(job, ["CO"], "x")
    monkeypatch.setattr(
        F, "harvest_job", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom"))
    )
    assert hooks.harvest_finished_job(job_id=42, wdir=job, state="error") is None
    assert (job / hooks.MARKER).exists()  # kept for a later try


def test_the_hook_does_not_write_a_network_store(returned_job, monkeypatch):
    job, path = returned_job
    hooks.write_marker(job, ["CO"], "x")
    monkeypatch.setattr(F, "store_is_writable", lambda p=None: (False, "on gpfs"))
    assert hooks.harvest_finished_job(job_id=42, wdir=job) is None
    assert not path.exists() and (job / hooks.MARKER).exists()


def test_write_marker_accumulates(tmp_path):
    hooks.write_marker(tmp_path, ["C", "CC"], "a")
    hooks.write_marker(tmp_path, ["CC", "O"], "b")
    data = json.loads((tmp_path / hooks.MARKER).read_text())
    assert data == {"reason": "b", "species": ["C", "CC", "O"]}


@pytest.mark.parametrize("smiles", ["C=C(F)F", "FC(F)(F)F", "FC(F)F", "FN(F)F"])
def test_fluorines_on_one_atom_are_not_too_close(smiles):
    """Real geminal fluorines are 2.12-2.17 Å apart, below the 2.2 Å molsystem
    uses for fluorines that do not share an atom (1,1-difluoroethylene was
    refused by the first official build's dry run)."""
    assert check_structure(*_embedded(smiles), smiles, 0, 1) is None


def test_the_hook_never_writes_a_newer_store(returned_job):
    """A store from newer code: the hook leaves it alone and keeps the marker,
    so a JobServer running newer code can harvest the job later."""
    import sqlite3

    from seamm_thermochemistry.geometry_store import SCHEMA_VERSION

    job, path = returned_job
    with FragmentStore(path, create=True):
        pass
    db = sqlite3.connect(path)
    db.execute(f"PRAGMA user_version = {SCHEMA_VERSION + 1}")
    db.commit()
    db.close()
    before = path.read_bytes()
    hooks.write_marker(job, ["CO"], "x")
    assert hooks.harvest_finished_job(job_id=42, wdir=job, state="finished") is None
    assert (job / hooks.MARKER).exists()
    assert path.read_bytes() == before
