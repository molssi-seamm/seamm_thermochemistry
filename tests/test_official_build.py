"""The official build's sanity pass and selection (P5e, see
docs/developer_guide/campaigns/2026-11-10, §6): saddle points, structures that
are not their species, and wrong spin states are refused unless excluded with a
reason; only the rows reachable at the chosen levels are kept; and the report is
written beside the file."""

import json

import pytest

pytest.importorskip("rdkit")

from seamm_thermochemistry import (  # noqa: E402
    EnergyRecord,
    FragmentRecord,
    FragmentStore,
    GeometryRecord,
    build_official,
)
from seamm_thermochemistry import fragment_cli  # noqa: E402
from seamm_thermochemistry import store_tools as T  # noqa: E402

GEO = "fp:r2scan"
REVDSD = "fp:revdsd"


def _xyz(smiles, seed=3):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    if mol.GetNumAtoms() > 2:
        AllChem.MMFFOptimizeMolecule(mol)
    return [a.GetSymbol() for a in mol.GetAtoms()], mol.GetConformer().GetPositions()


def _put(
    store, smiles, *, charge=0, multiplicity=1, fp=GEO, frequencies=None, xyz=None
):
    symbols, coordinates = xyz or _xyz(smiles)
    g = GeometryRecord.make(
        smiles,
        symbols,
        coordinates,
        charge=charge,
        multiplicity=multiplicity,
        E=-1000.0,
        ZPE=50.0,
        H=-940.0,
        frequency_fingerprint=fp,
        frequencies=frequencies or [500.0, 1500.0],
    )
    store.put_geometry(g)
    store.put_energy(
        EnergyRecord.make(g["geometry_hash"], -1100.0, energy_fingerprint=REVDSD)
    )
    store.put_energy(
        EnergyRecord.make(g["geometry_hash"], -1200.0, energy_fingerprint="fp:other")
    )
    store.put(
        FragmentRecord.make(
            smiles,
            charge=charge,
            multiplicity=multiplicity,
            E=-1000.0,
            ZPE=50.0,
            H=-940.0,
            energy_fingerprint=fp,
            frequency_fingerprint=fp,
        )
    )
    return g


class FakeEntry(dict):
    def __init__(self, id, value, multiplicity=None):
        super().__init__(id=id, dfH0_298K=value)
        self.multiplicity = multiplicity
        self.is_plain = multiplicity is None
        self.is_gas = True


class FakeATcT:
    """ATcT with methylene: plain = triplet (391.6), singlet labelled (429.1)."""

    version = "9.999"

    def candidates(self, smiles, charge=None):
        if smiles == "[CH2]":
            return [
                FakeEntry("ch2*0", 391.615),
                FakeEntry("ch2*1", 391.615, 3),
                FakeEntry("ch2*2", 429.14, 1),
            ]
        return []


@pytest.fixture()
def contribution(tmp_path):
    path = tmp_path / "contribution.db"
    with FragmentStore(path, create=True, layer="contribution") as s:
        _put(s, "C")
        _put(s, "[CH3]", multiplicity=2)  # a doublet: kept
        _put(s, "[NH4+]", charge=1)  # an ion: kept
        _put(s, "CC", fp="fp:subproduct")  # another level: deselected
    return path


# -- the checks -------------------------------------------------------------------


def test_ground_states():
    assert T.ground_state_multiplicity("[BH]") == (
        1,
        "singlet ground state (cbh.SINGLET_GROUND_STATES)",
    )
    assert T.ground_state_multiplicity("Cl[Si]Cl")[0] == 1
    assert T.ground_state_multiplicity("[CH2]", 0, FakeATcT()) == (
        3,
        "ATcT 9.999 (ch2*1)",
    )
    assert T.ground_state_multiplicity("[CH3]", 0, FakeATcT()) == (None, None)


def test_sanity_check(tmp_path):
    with FragmentStore(tmp_path / "s.db", create=True) as s:
        _put(s, "C")
        _put(s, "[CH3-]", charge=-1, frequencies=[-687.0, 500.0])  # a saddle
        _put(s, "[BH]", multiplicity=3)  # the wrong spin state
        _put(s, "[CH2]", multiplicity=1)  # the singlet, not ground
        _put(s, "CCO", xyz=_xyz("COC"))  # rearranged
        problems = T.sanity_check(s, atct=FakeATcT())
    assert set(problems) == {
        ("[CH3-]", -1, 1),
        ("[BH]", 0, 3),
        ("[CH2]", 0, 1),
        ("CCO", 0, 1),
    }
    assert "687i" in problems[("[CH3-]", -1, 1)][0]
    assert "ground state's (1, singlet ground state" in problems[("[BH]", 0, 3)][0]
    assert "ATcT 9.999" in problems[("[CH2]", 0, 1)][0]
    assert "not those of CCO" in problems[("CCO", 0, 1)][0]


# -- the build --------------------------------------------------------------------


def test_selection_and_report(contribution, tmp_path):
    out = tmp_path / "official.db"
    path, digest, _ = build_official(
        out,
        "2026.11.12",
        [contribution],
        frequency_fingerprints=[GEO],
        energy_fingerprints=[REVDSD],
        source="10.5281/zenodo.0",
    )
    with FragmentStore(path, immutable=True) as s:
        db = s._connect()
        species = {r[0] for r in db.execute("SELECT species FROM geometry")}
        assert species == {"C", "[CH3]", "[NH4+]"}
        assert {r[0] for r in db.execute("SELECT energy_fingerprint FROM energy")} == {
            REVDSD
        }
        assert {r[0] for r in db.execute("SELECT species FROM fragment")} == species
        assert s.meta()["version"] == "2026.11.12" and s.layer == "official"
    report = json.loads((tmp_path / "official.db.report.json").read_text())
    assert report["rows"] == {"fragment": 3, "geometry": 3, "energy": 3}
    assert report["species"] == 3 and report["sha256"] == digest
    assert report["deselected"]["geometry"] == 1 and report["checked"]


def test_problems_stop_the_build_until_excluded(contribution, tmp_path):
    with FragmentStore(contribution, create=True) as s:
        _put(s, "[CH3-]", charge=-1, frequencies=[-687.0, 500.0])
        _put(s, "[BH]", multiplicity=3)
    out = tmp_path / "official.db"
    with pytest.raises(ValueError, match="2 species states fail") as e:
        build_official(out, "2026.11.12", [contribution], frequency_fingerprints=[GEO])
    assert "[CH3-]" in str(e.value) and "--exclude" in str(e.value)
    assert not out.exists()
    exclusions = tmp_path / "exclude.csv"
    exclusions.write_text(
        "species,charge,multiplicity,reason\n"
        "[CH3-],-1,1,saddle point in job 5684 (687i)\n"
        "[BH],0,3,triplet; BH is a singlet\n"
    )
    exclude = T.load_exclusions(exclusions)
    build_official(
        out, "2026.11.12", [contribution], frequency_fingerprints=[GEO], exclude=exclude
    )
    report = json.loads((tmp_path / "official.db.report.json").read_text())
    excluded = {(e["species"], e["reason"]): e["rows"] for e in report["excluded"]}
    assert all(e["conformer"] is None for e in report["excluded"])
    assert excluded[("[CH3-]", "saddle point in job 5684 (687i)")] > 0
    with FragmentStore(out, immutable=True) as s:
        assert s.geometries("[BH]", multiplicity=3, frequency_fingerprint=GEO) == []


def test_an_exclusion_needs_a_reason(tmp_path):
    path = tmp_path / "x.csv"
    path.write_text("species,charge,multiplicity,reason\nC,0,1,\n")
    with pytest.raises(ValueError, match="needs a reason"):
        T.load_exclusions(path)


def test_cli_build(contribution, tmp_path, capsys):
    out = tmp_path / "official.db"
    code = fragment_cli.main(
        [
            "official",
            "build",
            "--version",
            "2026.11.12",
            "--out",
            str(out),
            "--frequency-fingerprint",
            GEO,
            "--energy-fingerprint",
            REVDSD,
            "--no-atct",
            str(contribution),
        ]
    )
    text = capsys.readouterr().out
    assert code == 0 and "Report:" in text and out.exists()


def test_singlet_carbenes_and_borylene_pass_the_structure_check():
    """Radical electrons may pair: only the parity is fixed by the SMILES, so the
    singlet ground states of BH and SiCl2 are not refused (the 1694 set had them
    as triplets); the ground-state check is the official build's."""
    from seamm_thermochemistry import check_structure

    symbols, xyz = ["B", "H"], [[0, 0, 0], [0, 0, 1.23]]
    assert check_structure(symbols, xyz, "[BH]", 0, 1) is None
    assert check_structure(symbols, xyz, "[BH]", 0, 3) is None
    assert "multiplicity" in check_structure(symbols, xyz, "[BH]", 0, 2)


def test_one_conformer_can_be_excluded(tmp_path):
    """A rearranged conformer (here: dimethyl ether's coordinates filed as a
    conformer of ethanol) is left out, the species' good conformer kept."""
    path = tmp_path / "c.db"
    with FragmentStore(path, create=True, layer="contribution") as s:
        good = _put(s, "CCO")
        bad = GeometryRecord.make(
            "CCO",
            *_xyz("COC"),
            conformer="bad",
            E=-1010.0,
            ZPE=50.0,
            H=-950.0,
            frequency_fingerprint=GEO,
            frequencies=[500.0, 1500.0],
        )
        s.put_geometry(bad)
        s.put_energy(
            EnergyRecord.make(bad["geometry_hash"], -1110.0, energy_fingerprint=REVDSD)
        )
    out = tmp_path / "official.db"
    with pytest.raises(ValueError, match="conformer 'bad'"):
        build_official(out, "v", [path], frequency_fingerprints=[GEO])
    exclusions = tmp_path / "x.csv"
    exclusions.write_text(
        "species,charge,multiplicity,conformer,reason\n"
        "CCO,0,1,bad,rearranged: dimethyl ether\n"
    )
    build_official(
        out,
        "v",
        [path],
        frequency_fingerprints=[GEO],
        exclude=T.load_exclusions(exclusions),
    )
    with FragmentStore(out, immutable=True) as s:
        (g,) = s.geometries("CCO", frequency_fingerprint=GEO)
        assert g["geometry_hash"] == good["geometry_hash"]
        assert s.energies(bad["geometry_hash"]) == []


def test_quarantined_rows_are_left_out_and_not_renamed(tmp_path):
    """Rows marked by hand with the model chemistry "UNUSABLE: <reason>" (and
    unverified) are never found, migrate leaves their text alone, export skips
    them and the official build drops them, listed in its report."""
    import sqlite3

    path = tmp_path / "local.db"
    with FragmentStore(path, create=True) as s:
        _put(s, "C")
        _put(s, "CC")
    reason = "UNUSABLE: dissociated into two CH3 (C...C 2.4 Å) @ job 5689"
    db = sqlite3.connect(path)
    for table in ("geometry", "fragment"):
        db.execute(
            f"UPDATE {table} SET verified = 0, model_chemistry = ? "
            "WHERE species = 'CC'",
            (reason,),
        )
    db.commit()
    db.close()
    T.migrate_store(path)
    db = sqlite3.connect(path)
    assert db.execute(
        "SELECT model_chemistry, program FROM geometry WHERE species = 'CC'"
    ).fetchone() == (reason, None)
    db.close()
    with FragmentStore(path) as s:
        assert s.geometries("CC", frequency_fingerprint=GEO) == []
    contribution = tmp_path / "c.db"
    counts = T.export_local(contribution, path, tmp_path / "none.db")
    assert counts["geometry"] == 1  # C only
    out = tmp_path / "official.db"
    build_official(out, "v", [path], frequency_fingerprints=[GEO])
    report = json.loads((tmp_path / "official.db.report.json").read_text())
    tables = {d["table"] for d in report["quarantined"]}
    assert tables == {"geometry", "fragment"}
    assert report["quarantined"][0]["reason"].startswith("dissociated")
    with FragmentStore(out, immutable=True) as s:
        assert (
            s._connect()
            .execute("SELECT COUNT(*) FROM geometry WHERE species = 'CC'")
            .fetchone()[0]
            == 0
        )


def test_triplet_ground_states():
    """O2, S2, SO and OCCO are triplets a SMILES writes closed-shell; a singlet row
    of them is refused (5689's singlet OCCO was a 419i saddle)."""
    for smiles in ("O=O", "S=S", "O=S", "O=C=C=O"):
        multiplicity, source = T.ground_state_multiplicity(smiles)
        assert multiplicity == 3 and "TRIPLET_GROUND_STATES" in source
    assert T.ground_state_multiplicity("[GeH2]")[0] == 1  # one of the six added


def test_a_singlet_dioxygen_row_is_refused(tmp_path):
    with FragmentStore(tmp_path / "s.db", create=True) as s:
        _put(s, "O=O", xyz=(["O", "O"], [[0, 0, 0], [0, 0, 1.21]]))
        _put(s, "O=O", multiplicity=3, xyz=(["O", "O"], [[0, 0, 0], [0, 0, 1.21]]))
        problems = T.sanity_check(s)
    assert set(problems) == {("O=O", 0, 1)}
    assert "ground state's (3, triplet" in problems[("O=O", 0, 1)][0]
