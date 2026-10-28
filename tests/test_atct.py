"""Tests for the local ATcT cache, the index, the selection rules and the
covariance handling, on the synthetic fixture in tests/data."""

import json
import math
from pathlib import Path

import pytest

rdkit = pytest.importorskip("rdkit")

from seamm_thermochemistry import atct as A  # noqa: E402

FIXTURE = Path(__file__).parent / "data" / "atct_fixture.json"


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    base = tmp_path_factory.mktemp("atct")
    A.download_all(base, source=FIXTURE)
    return base


@pytest.fixture(scope="module")
def tables(base):
    return A.ATcT(base=base)


# -- normalization ----------------------------------------------------------


def test_normalize_all_style_keys():
    raw = {
        "ATcT_TN_Version": "1.220",
        "ATcT_ID": "74-82-8*0",
        "Name": "Methane",
        "Formula": "CH4  (g)",
        "Delta_Hf_0K": "-66.550",
        "Delta_Hf_298K": "-74.513",
        "Delta_Hf298K_uncertainty": "0.043",
        "unit": "kJ/mol",
        "SMILES": "C",
        "CASRN": "74-82-8",
        "charge": "0",
    }
    e = A.normalize_entry(raw)
    assert e["version"] == "1.220"
    assert e["id"] == "74-82-8*0"
    assert e["formula"] == "CH4"
    assert e["phase"] == "g" and e["descriptor"] == ""
    assert e["dfH0_0K"] == pytest.approx(-66.55)
    assert e["dfH0_298K"] == pytest.approx(-74.513)
    assert e["uncertainty_298K"] == pytest.approx(0.043)
    assert e["exact"] is False
    assert e["charge"] == 0


def test_normalize_missing_unit_defaults_to_kJ_per_mol():
    # The real /all/ file has "unit": null for the 20 exact element entries
    e = A.normalize_entry({"Formula": "H2  (g)", "unit": None})
    assert e["units"] == "kJ/mol"


def test_normalize_species_route_keys():
    raw = {
        "ATcT_TN_Version": "1.220",
        "ATcT_ID": "1333-74-0*0",
        "Name": "Hydrogen",
        "Formula": "H2  (g)",
        "∆fH_0K": "0",
        "∆fH_298K": "0",
        "∆fH_298K_uncertainty": "exact",
        "units": "kJ/mol",
        "SMILES": "[H][H]",
        "InChI": "InChI=1S/H2/h1H",
        "InChI_Key": "UFHFLCQGNIYNRP-UHFFFAOYSA-N",
        "XYZ": "2\n\nH 0 0 0\nH 0 0 0.74\n",
        "charge": 0,
    }
    e = A.normalize_entry(raw)
    assert e["dfH0_0K"] == 0.0 and e["dfH0_298K"] == 0.0
    assert e["exact"] is True and e["uncertainty_298K"] == 0.0
    assert e["inchikey"] == "UFHFLCQGNIYNRP-UHFFFAOYSA-N"
    assert e["units"] == "kJ/mol"


@pytest.mark.parametrize(
    "formula, phase, descriptor",
    [
        ("H2O  (cr,l)", "cr,l", ""),
        ("CH2  (g, triplet)", "g", "triplet"),
        ("C2H2Cl2  (g, trans)", "g", "trans"),
        ("H2  (g, para)", "g", "para"),
        ("C6H5CH3 (g, ortho)", "g", "ortho"),
        ("NaCl (cr)", "cr", ""),
        ("C2H5NO2  (g, conformer Ip)", "g", "conformer Ip"),
    ],
)
def test_state_parsing(formula, phase, descriptor):
    e = A.normalize_entry({"Formula": formula})
    assert (e["phase"], e["descriptor"]) == (phase, descriptor)


# -- SMILES -----------------------------------------------------------------


def test_canonical_smiles_keeps_double_bond_stereo():
    assert A.canonical_smiles("C/C=C\\C") != A.canonical_smiles("C/C=C/C")
    assert A.canonical_smiles("C/C=C/C", stereo=False) == A.canonical_smiles("CC=CC")


def test_canonical_smiles_drops_a_single_stereocentre_keeps_two():
    # One centre: enantiomers have the same DfH, so the centre is dropped
    assert A.canonical_smiles("C[C@H](N)C(=O)O") == A.canonical_smiles("CC(N)C(=O)O")
    # Two centres: the relative configuration matters, so it is kept
    rr = A.canonical_smiles("C[C@H](O)[C@@H](O)C")
    rs = A.canonical_smiles("C[C@H](O)[C@H](O)C")
    assert rr != rs
    assert "@" in rr


def test_canonical_smiles_flat_drops_stereo_defining_hydrogen():
    assert A.canonical_smiles("[H]/N=C\\O", stereo=False) == "N=CO"
    assert A.canonical_smiles("[H]/N=C\\O") != "N=CO"


def test_canonical_smiles_unparseable():
    assert A.canonical_smiles("O=[N]=O") is None


# -- the cache on disk --------------------------------------------------------


def test_download_from_local_source(base):
    directory = base / "9.999"
    assert (directory / "atct_all.json").exists()
    meta = json.loads((directory / "metadata.json").read_text())
    assert meta["version"] == "9.999"
    assert meta["n_entries"] == 27
    assert len(meta["sha256"]) == 64
    assert "ver. 9.999 of the Thermochemical Network" in meta["citation"]
    # the file is the download, unmodified in content
    assert json.loads((directory / "atct_all.json").read_text()) == json.loads(
        FIXTURE.read_text()
    )


def test_download_refuses_to_overwrite(base):
    with pytest.raises(A.ATcTError, match="already"):
        A.download_all(base, source=FIXTURE)


def test_download_refuses_mixed_versions(tmp_path):
    data = json.loads(FIXTURE.read_text())
    data[0]["ATcT_TN_Version"] = "1.0"
    src = tmp_path / "mixed.json"
    src.write_text(json.dumps(data))
    with pytest.raises(A.ATcTError, match="mixes"):
        A.download_all(tmp_path / "cache", source=src)


def test_available_versions_sorted_numerically(tmp_path):
    for v in ("1.228", "1.9", "1.122q", "1.220"):
        d = tmp_path / v
        d.mkdir()
        (d / "atct_all.json").write_text("[]")
    assert A.ATcT.available(tmp_path) == ["1.9", "1.122q", "1.220", "1.228"]


def test_latest_version_is_default(tmp_path):
    A.download_all(tmp_path, source=FIXTURE)
    other = tmp_path / "0.001"
    other.mkdir()
    (other / "atct_all.json").write_text("[]")
    (other / "metadata.json").write_text("{}")
    assert A.ATcT(base=tmp_path).version == "9.999"
    assert A.ATcT(version="0.001", base=tmp_path).version == "0.001"
    with pytest.raises(FileNotFoundError):
        A.ATcT(version="1.0", base=tmp_path)


def test_missing_cache_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="No ATcT data"):
        A.ATcT(base=tmp_path / "nothing")


# -- the index -----------------------------------------------------------------


def test_index_built_and_saved(tables):
    index = tables.index
    assert (tables.directory / "index.json").exists()
    assert index["rdkit"] == rdkit.__version__
    assert "10102-44-0*0" in index["unparsed"]  # O=[N]=O
    assert index["stereo"][A.canonical_smiles("C")] == ["74-82-8*0"]
    # cis and trans 2-butene are different stereo keys but one flat key
    flat = index["flat"][A.canonical_smiles("CC=CC", stereo=False)]
    assert sorted(flat) == ["590-18-1*0", "624-64-6*0"]


# -- selection rules --------------------------------------------------------


def test_lookup_plain_entry(tables):
    e = tables.lookup("C")
    assert e["id"] == "74-82-8*0"
    assert e["dfH0_298K"] == pytest.approx(-74.513)
    assert e.is_plain and e.is_gas


def test_lookup_gas_only(tables):
    # Water has (g) and (cr,l); only the gas entry is a candidate
    assert tables.lookup("O")["id"] == "7732-18-5*0"
    assert [e["id"] for e in tables.candidates("O", gas_only=False)] == [
        "7732-18-5*0",
        "7732-18-5*3",
    ]


def test_lookup_charge(tables):
    assert tables.lookup("[OH-]", charge=-1)["id"] == "14280-30-9*0"
    with pytest.raises(A.MissingEntry):
        tables.lookup("[OH-]", charge=0)


def test_lookup_prefers_star0_over_conformers(tables):
    e = tables.lookup("CC(=O)O")
    assert e["id"] == "64-19-7*0"
    assert tables.lookup("NCC(=O)O")["id"] == "56-40-6*0"


def test_lookup_conformer_label(tables):
    assert tables.lookup("CC(=O)O", conformer="anti")["id"] == "64-19-7*2"
    assert tables.lookup("CC(=O)O", conformer="syn")["id"] == "64-19-7*1"
    with pytest.raises(A.MissingEntry, match="labelled"):
        tables.lookup("CC(=O)O", conformer="gauche")


def test_lookup_nuclear_spin_isomers_ignored(tables):
    # H2 has (g), (g, ortho) and (g, para): the plain entry is used
    e = tables.lookup("[H][H]")
    assert e["id"] == "1333-74-0*0"
    assert e["exact"] and e["uncertainty_298K"] == 0.0


def test_lookup_isotopologues_resolved_by_mass(tables):
    # D2 and HD share the SMILES [H][H] with H2: the natural-abundance mass wins
    assert tables.lookup("[H][H]")["name"] == "Hydrogen"
    ids = sorted(e["id"] for e in tables.candidates("[H][H]"))
    assert "7782-39-0*0" in ids and "13983-20-5*0" in ids


def test_lookup_spin_state(tables):
    # Plain CH2 is the triplet ground state; (g, singlet) is labelled
    assert tables.lookup("[CH2]", multiplicity=1)["id"] == "2465-56-7*1"
    assert tables.lookup("[CH2]", multiplicity=3)["id"] == "2465-56-7*0"
    assert tables.lookup("[CH2]")["id"] == "2465-56-7*0"


def test_lookup_cis_trans_need_stereo(tables):
    assert tables.lookup("C/C=C\\C")["id"] == "590-18-1*0"
    assert tables.lookup("C/C=C/C")["id"] == "624-64-6*0"
    with pytest.raises(A.AmbiguousEntry) as info:
        tables.lookup("CC=CC")
    ids = sorted(e["id"] for e in info.value.candidates)
    assert ids == ["590-18-1*0", "624-64-6*0"]
    assert "590-18-1*0" in str(info.value)


def test_lookup_mixture_entry_when_no_stereo_given(tables):
    # 1,2-dichloroethene has a plain (g) entry besides cis and trans: a
    # stereo-free query gets the plain one, a stereo query the isomer
    assert tables.lookup("ClC=CCl")["id"] == "540-59-0*0"
    assert tables.lookup("Cl/C=C\\Cl")["id"] == "156-59-2*0"
    assert tables.lookup("Cl/C=C/Cl")["id"] == "156-60-5*0"


def test_lookup_stereo_query_falls_back_to_unspecified_entry(tables):
    # A CBH fragment can carry stereo the tables never specify (the E/Z of an
    # N-H imine): the entry written without stereo is used
    assert tables.lookup("[H]/N=C\\O")["id"] == "75-12-7*9"
    assert tables.lookup("N=CO")["id"] == "75-12-7*9"
    # but where the tables do have the stereo, it must match
    assert tables.lookup("Cl/C=C\\Cl")["id"] == "156-59-2*0"


def test_lookup_enantiomer_ignores_single_centre(tables):
    assert tables.lookup("C[C@@H](N)C(=O)O")["id"] == "50-99-7*0"
    assert tables.lookup("CC(N)C(=O)O")["id"] == "50-99-7*0"


def test_lookup_by_id_and_missing(tables):
    assert tables.lookup(atct_id="64-19-7*2")["name"] == "Acetic acid"
    with pytest.raises(A.MissingEntry):
        tables.lookup(atct_id="0-00-0*0")
    with pytest.raises(A.MissingEntry):
        tables.lookup("CCCCCCCC")
    with pytest.raises(A.ATcTError, match="cannot parse"):
        tables.lookup("not a smiles")


def test_entry_label_and_groups(tables):
    e = tables.lookup("CC(=O)O", conformer="anti")
    assert e.casrn_group == "64-19-7" and e.index == 2
    assert "64-19-7*2" in e.label() and "± 0.7" in e.label()


# -- covariance --------------------------------------------------------------

IDS = ["593-53-3*0", "74-84-0*0", "74-82-8*0", "353-36-6*0"]
UNC = {"593-53-3*0": 0.29, "74-84-0*0": 0.14, "74-82-8*0": 0.043, "353-36-6*0": 0.35}


def _block(ids):
    """A plausible covariance block in the API's layout, rows in `ids` order."""
    n = len(ids)
    matrix = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            if i == j:
                matrix[i][j] = UNC[a] ** 2
            else:
                matrix[i][j] = 0.3 * UNC[a] * UNC[b]
    return {
        "species": [
            {"id": k, "atctid": _id, "label": _id} for k, _id in enumerate(ids)
        ],
        "species_atctids": ids,
        "units": "kJ/mol",
        "matrix": matrix,
    }


def test_covariance_fetched_reordered_and_cached(tables, monkeypatch):
    calls = []

    def fake_get(url, **kwargs):
        calls.append(url)
        # The API returns the rows in its own order
        return _block(["74-82-8*0", "74-84-0*0", "353-36-6*0", "593-53-3*0"])

    monkeypatch.setattr(A, "_get", fake_get)
    sigma = tables.covariance(IDS)
    assert len(calls) == 1 and "covariance/species" in calls[0]
    # In the requested order
    for i, _id in enumerate(IDS):
        assert math.sqrt(sigma[i][i]) == pytest.approx(UNC[_id])
    assert sigma[0][3] == pytest.approx(0.3 * UNC["593-53-3*0"] * UNC["353-36-6*0"])
    # A subset is served from the cache without another request
    sub = tables.covariance(["353-36-6*0", "74-82-8*0"])
    assert len(calls) == 1
    assert math.sqrt(sub[0][0]) == pytest.approx(0.35)
    assert (tables.directory / "covariance.json").exists()
    # and no fetch when asked not to
    assert tables.covariance(["7732-18-5*0"], fetch=False) is None


def test_covariance_rejects_wrong_scale(tables):
    block = _block(IDS)
    block["matrix"][0][0] = (UNC[IDS[0]] / 2) ** 2  # 1-sigma instead of 95 %
    with pytest.raises(A.ATcTError, match="95 %"):
        tables.parse_covariance(block)


def test_covariance_failure_falls_back_to_none(tables, monkeypatch):
    def failing(url, **kwargs):
        raise A.ATcTError("HTTP 403")

    monkeypatch.setattr(A, "_get", failing)
    assert tables.covariance(["7732-18-5*0", "14280-30-9*0"]) is None


def test_reaction_uncertainty():
    # Fluoroethane CBH-1: CH3CH2F + CH4 -> CH3F + C2H6; reference part
    # CH3F + C2H6 - CH4 with the fixture's uncertainties.
    coef = [1, 1, -1]
    unc = [0.29, 0.14, 0.043]
    quad, used = A.reaction_uncertainty(coef, unc)
    assert not used
    assert quad == pytest.approx(math.sqrt(0.29**2 + 0.14**2 + 0.043**2))
    cov = [
        [u * v * (1.0 if i == j else 0.3) for j, v in enumerate(unc)]
        for i, u in enumerate(unc)
    ]
    full, used = A.reaction_uncertainty(coef, unc, cov)
    assert used
    expected = 0.0
    for i in range(3):
        for j in range(3):
            expected += coef[i] * cov[i][j] * coef[j]
    assert full == pytest.approx(math.sqrt(expected))
    assert full != pytest.approx(quad)
