"""Model-chemistry names stored bare with a program column, and lookups by a set
of expected fingerprints (P5a of the fragment store, see
docs/developer_guide/campaigns/2026-11-10)."""

import pytest

from seamm_thermochemistry import (
    EnergyRecord,
    FragmentRecord,
    FragmentStore,
    GeometryRecord,
    split_model_chemistry,
)
from seamm_thermochemistry.geometry_store import same_model_chemistry

WATER = (["O", "H", "H"], [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]])


@pytest.fixture()
def store(tmp_path):
    with FragmentStore(tmp_path / "fragments.db", create=True) as s:
        yield s


@pytest.mark.parametrize(
    "text, program, expected",
    [
        ("ORCA:CMP@R2SCAN-3C", None, ("CMP@R2SCAN-3C", "ORCA")),
        ("CMP@R2SCAN-3C", "ORCA", ("CMP@R2SCAN-3C", "ORCA")),
        ("CMP@R2SCAN-3C", None, ("CMP@R2SCAN-3C", None)),
        ("ORCA:DFT@B3LYP/bse:def2-SVP", None, ("DFT@B3LYP/bse:def2-SVP", "ORCA")),
        ("DFT@B3LYP/bse:def2-SVP", None, ("DFT@B3LYP/bse:def2-SVP", None)),
        (
            "QC@CCSD(T)-F12D_RI/cc-pVTZ-F12",
            "ORCA",
            ("QC@CCSD(T)-F12D_RI/cc-pVTZ-F12", "ORCA"),
        ),
        ("orca:CMP@R2SCAN-3C", "ORCA", ("CMP@R2SCAN-3C", "ORCA")),
        # Fingerprints and other strings without '@' are left alone
        ("sha256:6b46c8", None, ("sha256:6b46c8", None)),
        (None, "ORCA", (None, "ORCA")),
    ],
)
def test_split_model_chemistry(text, program, expected):
    assert split_model_chemistry(text, program) == expected


def test_contradicting_program_is_logged_and_the_explicit_one_kept(caplog):
    assert split_model_chemistry("Gaussian:DFT@B3LYP/6-31G*", "ORCA") == (
        "DFT@B3LYP/6-31G*",
        "ORCA",
    )
    assert "keeping ORCA" in caplog.text


def test_same_model_chemistry():
    assert same_model_chemistry("ORCA:CMP@R2SCAN-3C", "CMP@R2SCAN-3C")
    assert same_model_chemistry("ORCA:CMP@R2SCAN-3C", "CMP@R2SCAN-3C", None, "orca")
    assert not same_model_chemistry("ORCA:CMP@R2SCAN-3C", "Psi4:CMP@R2SCAN-3C")
    assert not same_model_chemistry("ORCA:CMP@R2SCAN-3C", "ORCA:DFT@R2SCAN-3C")


def test_records_store_the_bare_model_and_the_program():
    g = GeometryRecord.make(
        "O", *WATER, model_chemistry="ORCA:CMP@R2SCAN-3C", frequency_fingerprint="f"
    )
    assert (g["model_chemistry"], g["program"]) == ("CMP@R2SCAN-3C", "ORCA")
    e = EnergyRecord.make("sha256:x", 1.0, model_chemistry="ORCA:HF", program="ORCA")
    assert (e["model_chemistry"], e["program"]) == ("ORCA:HF", "ORCA")  # no '@'
    e = EnergyRecord.make("sha256:x", 1.0, model_chemistry="ORCA:HF@def2-SVP")
    assert (e["model_chemistry"], e["program"]) == ("HF@def2-SVP", "ORCA")
    r = FragmentRecord.make("O", E=1.0, model_chemistry="MOPAC:SQM@PM7")
    assert (r["model_chemistry"], r["program"]) == ("SQM@PM7", "MOPAC")
    # The unverified key keeps the prefixed spelling of earlier stores
    assert r["energy_fingerprint"] == "unverified:MOPAC:SQM@PM7@?"


def test_one_level_whatever_the_spelling():
    # Optimized at "CMP@..." (a sub-step) with frequencies at "ORCA:CMP@...":
    # the same level, so no pairing error
    g = GeometryRecord.make(
        "O",
        *WATER,
        model_chemistry="ORCA:CMP@R2SCAN-3C",
        optimization_model_chemistry="CMP@R2SCAN-3C",
        frequency_fingerprint="f",
    )
    assert g["model_chemistry"] == "CMP@R2SCAN-3C"


@pytest.mark.parametrize("asked", ["ORCA:CMP@R2SCAN-3C", "CMP@R2SCAN-3C"])
def test_unverified_lookups_accept_either_spelling(store, asked):
    store.put_geometry(
        GeometryRecord.make(
            "O", *WATER, model_chemistry="CMP@R2SCAN-3C", program="ORCA", E=1.0
        )
    )
    (g,) = store.geometries("O", model_chemistry=asked)
    assert g["model_chemistry"] == "CMP@R2SCAN-3C"
    store.put_energy(
        EnergyRecord.make(g["geometry_hash"], 2.0, model_chemistry="ORCA:HF@SVP")
    )
    for name in ("ORCA:HF@SVP", "HF@SVP"):
        assert store.energy(g["geometry_hash"], model_chemistry=name)["E"] == 2.0
    assert store.energy(g["geometry_hash"], model_chemistry="Psi4:HF@SVP") is None
    store.put(FragmentRecord.make("O", E=3.0, model_chemistry="ORCA:HF@SVP"))
    assert store.get("O", model_chemistry=asked.replace("CMP@R2SCAN-3C", "HF@SVP"))[
        "E"
    ] == (3.0)


def test_rows_written_before_the_names_were_normalized(store):
    """Rows with a prefixed name and no program (as in ChemAI's store before
    migration) are still found, by either spelling, and '_' in a name is not a
    wildcard."""
    db = store._connect(create=True)
    for name in ("ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP", "ORCA:DFT@X"):
        db.execute(
            "INSERT INTO energy (geometry_hash, energy_fingerprint, verified, "
            "model_chemistry, program, E) VALUES (?, ?, 0, ?, NULL, ?)",
            ("sha256:g", f"unverified:{name}@?", name, len(name)),
        )
    db.commit()
    for asked in (
        "DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP",
        "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP",
    ):
        row = store.energy("sha256:g", model_chemistry=asked)
        assert row["model_chemistry"].endswith("D4_2021/def2-TZVPP")
    assert store.energy("sha256:g", model_chemistry="DFT@REVDSD-PBEP86-D4x2021") is (
        None
    )


def test_lookup_by_a_set_of_fingerprints(store):
    """A boron species' energy has the NoCOSX fingerprint; an ordinary one's has
    the default: the caller passes the fingerprints it expects for each."""
    for species, fp in (("O", "fp:cosx"), ("B", "fp:nocosx")):
        symbols = ["O", "H", "H"] if species == "O" else ["B", "H", "H", "H"]
        coordinates = WATER[1] if species == "O" else [[0, 0, 0]] + WATER[1]
        g = GeometryRecord.make(
            species,
            symbols,
            coordinates,
            E=1.0,
            ZPE=0.5,
            H=1.6,
            frequency_fingerprint="fp:r2scan",
            optimization_fingerprint="fp:r2scan",
        )
        store.put_geometry(g)
        store.put_energy(
            EnergyRecord.make(g["geometry_hash"], 10.0, energy_fingerprint=fp)
        )
        store.put(
            FragmentRecord.make(species, E=10.0, ZPE=0.5, H=10.6, energy_fingerprint=fp)
        )
    expected = ["fp:nocosx", "fp:cosx"]
    for species in ("O", "B"):
        (g,) = store.geometries(species, frequency_fingerprint={"fp:r2scan", "x"})
        assert store.energy(g["geometry_hash"], energy_fingerprint=expected)
        assert store.energy(g["geometry_hash"], energy_fingerprint="fp:other") is None
        (got,) = store.species_at(
            species, frequency_fingerprint=["fp:r2scan"], energy_fingerprint=expected
        )
        assert got["E"] == 10.0
        assert store.get(species, energy_fingerprint=expected)["E"] == 10.0
        # The fragment row has a ZPE, so its frequency key is not empty
        assert (
            store.get(species, energy_fingerprint=expected, frequency_fingerprint="")
            is None
        )
    # The target's fingerprint alone misses the other kind of species
    (g,) = store.geometries("O", frequency_fingerprint="fp:r2scan")
    assert store.energy(g["geometry_hash"], energy_fingerprint="fp:nocosx") is None
    # The geometry level itself, by a set containing its fingerprint
    (same,) = store.species_at(
        "O", frequency_fingerprint="fp:r2scan", energy_fingerprint=["a", "fp:r2scan"]
    )
    assert same["E"] == 1.0 and same["energy"] is None
