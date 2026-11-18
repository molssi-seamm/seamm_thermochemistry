"""E/Z isomers: one rule for ATcT, the reference library and the computed values
(seamm_thermochemistry.stereo). A query with a defined E/Z gets that isomer's
value or none; never the generic value of a species ATcT records without E/Z.

The ATcT entries are ATcT 1.228's (ethanimine, difluorodiazene, formamide enol
and diazene), with a made-up imine that ATcT would list only generically."""

import json

import pytest

rdkit = pytest.importorskip("rdkit")

from seamm_thermochemistry import atct as A  # noqa: E402
from seamm_thermochemistry import references as R  # noqa: E402
from seamm_thermochemistry.stereo import (  # noqa: E402
    descriptor_isomer,
    double_bond_isomers,
)


def _entry(atct_id, name, formula, h0, h298, unc, smiles, casrn, mass):
    return {
        "ATcT_TN_Version": "1.228",
        "ATcT_ID": atct_id,
        "Name": name,
        "Formula": formula,
        "Delta_Hf_0K": h0,
        "Delta_Hf_298K": h298,
        "Delta_Hf298K_uncertainty": unc,
        "unit": "kJ/mol",
        "mass": mass,
        "SMILES": smiles,
        "CASRN": casrn,
        "charge": "0",
    }


ENTRIES = [
    # Ethanimine: isomers only in the descriptor, with a generic entry
    _entry(
        "56003-81-7*0",
        "(1E)-Ethanimine",
        "CH3CHNH  (g, trans)",
        "55.79",
        "39.58",
        "0.72",
        "CC=N",
        "20729-41-3",
        "43.0678",
    ),
    _entry(
        "56003-82-8*0",
        "(1Z)-Ethanimine",
        "CH3CHNH  (g, cis)",
        "58.9",
        "44.2",
        "1.4",
        "CC=N",
        "20729-41-3",
        "43.0678",
    ),
    _entry(
        "20729-41-3*0",
        "Ethanimine",
        "CH3CHNH  (g)",
        "55.79",
        "41.16",
        "0.72",
        "CC=N",
        "20729-41-3",
        "43.0678",
    ),
    # Difluorodiazene
    _entry(
        "13812-43-6*2",
        "(Z)-Difluorodiazene",
        "FNNF  (g, cis)",
        "78.4",
        "72.9",
        "2.1",
        "FN=NF",
        "10578-16-2",
        "66.01029",
    ),
    _entry(
        "13776-62-0*1",
        "(E)-Difluorodiazene",
        "FNNF  (g, trans)",
        "84.2",
        "78.9",
        "2.3",
        "FN=NF",
        "10578-16-2",
        "66.01029",
    ),
    _entry(
        "10578-16-2*0",
        "Difluorodiazene",
        "FNNF  (g)",
        "78.4",
        "73.4",
        "2.1",
        "FN=NF",
        "10578-16-2",
        "66.01029",
    ),
    # Formamide enol (formimidic acid): anti/syn of N-H to O, then of O-H
    _entry(
        "60100-09-6*0",
        "Formamide enol",
        "NHCHOH  (g)",
        "-130.6",
        "-141.6",
        "1.3",
        "N=CO",
        "60100-09-6",
        "45.04066",
    ),
    _entry(
        "60100-09-6*1",
        "Formamide enol",
        "NHCHOH  (g, anti-syn)",
        "-130.6",
        "-141.7",
        "1.3",
        "N=CO",
        "60100-09-6",
        "45.04066",
    ),
    _entry(
        "60100-09-6*2",
        "Formamide enol",
        "NHCHOH  (g, anti-anti)",
        "-116.6",
        "-127.6",
        "2.1",
        "N=CO",
        "60100-09-6",
        "45.04066",
    ),
    _entry(
        "60100-09-6*3",
        "Formamide enol",
        "NHCHOH  (g, syn-anti)",
        "-115.9",
        "-126.6",
        "2.1",
        "N=CO",
        "60100-09-6",
        "45.04066",
    ),
    _entry(
        "60100-09-6*4",
        "Formamide enol",
        "NHCHOH  (g, syn-syn)",
        "-109.8",
        "-120.6",
        "2.1",
        "N=CO",
        "60100-09-6",
        "45.04066",
    ),
    # Diazene: the E/Z is in ATcT's SMILES
    _entry(
        "15626-43-4*0",
        "(E)-Diazene",
        "HNNH  (g, trans)",
        "207.27",
        "200.13",
        "0.27",
        "[H]/N=N/[H]",
        "15626-43-4",
        "30.02936",
    ),
    _entry(
        "15626-42-3*0",
        "(Z)-Diazene",
        "HNNH  (g, cis)",
        "228.92",
        "221.78",
        "0.56",
        "[H]/N=N\\[H]",
        "15626-42-3",
        "30.02936",
    ),
    # A made-up imine listed only generically
    _entry(
        "99999-99-9*0",
        "N-Methylethanimine",
        "CH3CHNCH3  (g)",
        "40.0",
        "25.0",
        "2.0",
        "CC=NC",
        "99999-99-9",
        "57.09",
    ),
]


@pytest.fixture(scope="module")
def tables(tmp_path_factory):
    base = tmp_path_factory.mktemp("atct_ez")
    source = base / "all.json"
    source.write_text(json.dumps(ENTRIES))
    A.download_all(base, source=source)
    return A.ATcT(base=base)


# -- the shared helpers -----------------------------------------------------------


@pytest.mark.parametrize(
    "smiles, isomers",
    [
        ("[H]/N=C/C", ["E"]),
        ("[H]/N=C\\C", ["Z"]),
        ("F/N=N/F", ["E"]),
        ("[H]/N=C/O", ["E"]),
        ("[H]/N=C(\\C)O", ["E"]),  # acetimidic acid, E as thermo keys it
        ("[H]/N=C(/C)O", ["Z"]),
        ("CC=N", []),
        ("C/C=C/C=C/C", ["E", "E"]),
    ],
)
def test_double_bond_isomers(smiles, isomers):
    assert double_bond_isomers(smiles) == isomers


@pytest.mark.parametrize(
    "descriptor, atct_id, isomer",
    [
        ("trans", None, "E"),
        ("cis", None, "Z"),
        ("g, trans", None, "E"),
        ("anti-syn", "60100-09-6*1", "E"),
        ("syn-anti", "60100-09-6*3", "Z"),
        ("trans-cis equilib", None, None),
        ("", None, None),
        ("gauche", None, None),
        ("anti", None, None),  # a conformer, not an E/Z isomer
    ],
)
def test_descriptor_isomer(descriptor, atct_id, isomer):
    assert descriptor_isomer(descriptor, atct_id) == isomer


# -- ATcT ------------------------------------------------------------------------


@pytest.mark.parametrize(
    "smiles, atct_id, dfh298",
    [
        ("[H]/N=C/C", "56003-81-7*0", 39.58),  # (E)-ethanimine, not the generic 41.16
        ("[H]/N=C\\C", "56003-82-8*0", 44.2),
        ("F/N=N/F", "13776-62-0*1", 78.9),  # (E)-FNNF, not the generic 73.4
        ("F/N=N\\F", "13812-43-6*2", 72.9),
        ("[H]/N=C/O", "60100-09-6*1", -141.7),  # E: lowest of anti-syn, anti-anti
        ("[H]/N=C\\O", "60100-09-6*3", -126.6),  # Z: lowest of syn-anti, syn-syn
        ("[H]/N=N/[H]", "15626-43-4*0", 200.13),  # E/Z in ATcT's SMILES
        ("[H]/N=N\\[H]", "15626-42-3*0", 221.78),
    ],
)
def test_lookup_gives_the_isomer(tables, smiles, atct_id, dfh298):
    entry = tables.lookup(smiles)
    assert entry["id"] == atct_id
    assert entry["dfH0_298K"] == pytest.approx(dfh298)


@pytest.mark.parametrize(
    "smiles, atct_id",
    [("CC=N", "20729-41-3*0"), ("FN=NF", "10578-16-2*0"), ("N=CO", "60100-09-6*0")],
)
def test_without_ez_the_generic_entry_as_before(tables, smiles, atct_id):
    assert tables.lookup(smiles)["id"] == atct_id


def test_no_isomer_entry_is_refused_not_the_generic(tables):
    """ATcT has the imine only generically: an (E) query gets nothing."""
    assert tables.lookup("CC=NC")["id"] == "99999-99-9*0"
    with pytest.raises(A.AmbiguousEntry, match="no entry labelled as the E isomer"):
        tables.lookup("C/C=N/C")


def test_an_explicit_conformer_label_still_wins(tables):
    assert tables.lookup("N=CO", conformer="syn-syn")["id"] == "60100-09-6*4"


# -- the reference library and the computed values: rule 1 ------------------------


@pytest.fixture
def lib(tmp_path):
    library = R.ReferenceLibrary(tmp_path / "references.db", create=True)
    library.add_source("WebBook-SRD69", "measured", 50)
    yield library
    library.close()


def test_library_values_are_per_isomer(lib):
    """The library keys the exact SMILES: an isomer's value is used only for
    that isomer, and a value recorded without E/Z for neither."""
    e = lib.add_species("[H]/N=C(\\C)O")
    z = lib.add_species("[H]/N=C(/C)O")
    plain = lib.add_species("CC=NC")
    for sid, value in ((e, -189.5), (z, -177.4), (plain, 25.0)):
        lib.add_measurement(sid, "WebBook-SRD69", dfh_298k=value, citation="test")
    assert [r["id"] for r in lib.species("[H]/N=C(\\C)O")] == [e]
    assert [r["id"] for r in lib.species("[H]/N=C(/C)O")] == [z]
    assert lib.species("C/C=N/C") == []  # no value for an isomer
    assert [r["id"] for r in lib.species("CC=NC")] == [plain]
