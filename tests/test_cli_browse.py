"""``lookup`` and ``geometries`` from the command line without a fingerprint:
people browsing the store know the model chemistry, not the fingerprint, so the
commands list the levels, match verified rows by model chemistry and say why
nothing matched (thermo, 2026-11-12 and 11-15: "No record." even for ethane)."""

import sqlite3

import pytest

from seamm_thermochemistry import (
    EnergyRecord,
    FragmentRecord,
    FragmentStore,
    GeometryRecord,
)
from seamm_thermochemistry import fragment_cli

REVDSD = "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
R2SCAN = "ORCA:CMP@R2SCAN-3C"
E_FP = "sha256:revdsd"
F_FP = "sha256:r2scan"
WATER = (["O", "H", "H"], [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]])
BENT = (["O", "H", "H"], [[0, 0, 0.2], [0, 0.7, -0.5], [0, -0.7, -0.5]])


def _geometry(xyz=WATER, conformer="", E=-200000.0):
    return GeometryRecord.make(
        "O",
        *xyz,
        conformer=conformer,
        E=E,
        ZPE=55.0,
        H=E + 65.0,
        frequency_fingerprint=F_FP,
        model_chemistry=R2SCAN,
        job="5735",
    )


@pytest.fixture()
def stores(tmp_path):
    local = tmp_path / "fragments.db"
    official = tmp_path / "fragments-official.db"
    with FragmentStore(local, create=True) as s:
        # Verified, at revDSD: found by a calculation only by its fingerprint
        s.put(
            FragmentRecord.make(
                "O",
                E=-199000.0,
                ZPE=55.0,
                H=-198935.0,
                energy_fingerprint=E_FP,
                model_chemistry=REVDSD,
                job="5735",
            )
        )
        # Unverified, MOPAC: found by its model chemistry
        s.put(FragmentRecord.make("O", E=-241.8, model_chemistry="MOPAC:PM6", job="1"))
        # A radical, stored at multiplicity 2
        s.put(
            FragmentRecord.make(
                "[OH]",
                multiplicity=2,
                E=-100.0,
                energy_fingerprint=E_FP,
                model_chemistry=REVDSD,
            )
        )
        # A quarantined record
        s.put(
            FragmentRecord.make(
                "O", E=-1.0, energy_fingerprint="sha256:bad", model_chemistry=REVDSD
            )
        )
        g = _geometry()
        s.put_geometry(g)
        s.put_energy(
            EnergyRecord.make(
                g["geometry_hash"],
                -199000.0,
                energy_fingerprint=E_FP,
                model_chemistry=REVDSD,
            )
        )
        s.put_geometry(_geometry(BENT, conformer="saddle", E=-199990.0))
    db = sqlite3.connect(local)
    with db:
        db.execute(
            "UPDATE fragment SET verified = 0, model_chemistry = "
            "'UNUSABLE: test quarantine' WHERE energy_fingerprint = 'sha256:bad'"
        )
        db.execute(
            "UPDATE geometry SET n_imaginary = 1, max_imaginary = 300.0 "
            "WHERE conformer = 'saddle'"
        )
    db.close()
    with FragmentStore(official, create=True, layer="official") as s:
        s.set_meta(version="2026.11.11")
        s.put(
            FragmentRecord.make(
                "O",
                E=-199000.0,
                ZPE=55.0,
                H=-198935.0,
                energy_fingerprint=E_FP,
                model_chemistry=REVDSD,
                job="5602",
            )
        )
    return local, official


@pytest.fixture()
def run(stores, capsys):
    local, official = stores

    def run(*args):
        code = fragment_cli.main(
            ["--store", str(local), "--official", str(official)] + list(args)
        )
        return code, capsys.readouterr().out

    return run


def test_lookup_without_a_level_lists_the_levels(run):
    code, out = run("lookup", "O")
    assert code == 1
    assert "is stored at 2 levels" in out
    assert f"{REVDSD} [local] [official 2026.11.11]" in out
    assert f"energy fingerprint {E_FP}" in out
    assert "MOPAC:PM6 [unverified] [local]" in out
    assert "is quarantined: test quarantine" in out


def test_lookup_by_model_chemistry_finds_verified_records(run):
    code, out = run("lookup", "O", "--model-chemistry", REVDSD)
    assert code == 0
    assert "E = -199000.000" in out and "[unverified]" not in out
    assert f"fingerprints: energy {E_FP}" in out
    assert "[local]" in out and "[official 2026.11.11]" in out
    # The bare name and the unverified level
    assert run("lookup", "O", "--model-chemistry", REVDSD.split(":")[1])[0] == 0
    code, out = run("lookup", "O", "--model-chemistry", "MOPAC:PM6")
    assert code == 0 and "E = -241.800" in out and "[unverified]" in out


def test_lookup_says_why_nothing_matched(run):
    code, out = run("lookup", "O", "--model-chemistry", "MOPAC:PM7")
    assert code == 1
    assert "No record of O (q=0, 2S+1=1) at MOPAC:PM7." in out
    assert "It is stored at:" in out and REVDSD in out
    code, out = run("lookup", "[OH]")
    assert code == 1 and "stored with --charge 0 --multiplicity 2" in out
    assert run("lookup", "[OH]", "--multiplicity", "2")[0] == 0  # one level
    code, out = run("lookup", "CCCC")
    assert code == 1 and "No record of CCCC in" in out


def test_lookup_by_fingerprint_is_unchanged(run):
    code, out = run("lookup", "O", "--fingerprint", E_FP)
    assert code == 0 and "E = -199000.000" in out and "[local]" in out
    code, out = run("lookup", "O", "--fingerprint", "sha256:other")
    assert code == 1
    assert "No record of O (q=0, 2S+1=1) with the energy fingerprint" in out
    assert "is stored at:" in out


def test_geometries_without_a_fingerprint(run):
    code, out = run("geometries", "O")  # one level: shown
    assert code == 0
    assert f"O (q=0, 2S+1=1) at {R2SCAN}: E = -200000.000, ZPE = 55.000" in out
    assert f"frequency fingerprint: {F_FP}" in out
    assert f"energy fingerprint {E_FP}" in out
    assert "1 geometry with imaginary modes not shown" in out
    code, out = run("geometries", "O", "--model-chemistry", R2SCAN)
    assert code == 0 and "E = -200000.000" in out
    code, out = run("geometries", "O", "--include-imaginary")
    assert code == 0 and "conformer 'saddle'" in out and "1 imaginary" in out
    code, out = run("geometries", "O", "--fingerprint", "sha256:other")
    assert code == 1 and "No geometry of O" in out and F_FP in out


def test_list_prints_fingerprints(stores, capsys):
    local, _ = stores
    assert fragment_cli.main(["--store", str(local), "list", "O"]) == 0
    assert f"fingerprints: energy {E_FP}" in capsys.readouterr().out
