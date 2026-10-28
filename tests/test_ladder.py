"""The CBH-ANL ladder on a synthetic fragment store and the synthetic ATcT
fixture: self-consistency, the recursion, the ATcT coefficients and
uncertainties, the measured method term, and what happens when a rung is
missing from the store."""

from pathlib import Path

import pytest

rdkit = pytest.importorskip("rdkit")

from seamm_thermochemistry import atct as A  # noqa: E402
from seamm_thermochemistry import ladder as L  # noqa: E402
from seamm_thermochemistry.fragment_store import FragmentRecord  # noqa: E402
from seamm_thermochemistry.fragment_store import FragmentStore  # noqa: E402

FIXTURE = Path(__file__).parent / "data" / "atct_fixture.json"
LEVEL = "ORCA:QC@CCSD(T)-F12/bse:cc-pVTZ-F12"
FREQ = "ORCA:DFT@B2PLYP-D3/bse:cc-pVTZ"


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    base = tmp_path_factory.mktemp("atct")
    A.download_all(base, source=FIXTURE)
    return base


@pytest.fixture()
def tables(base):
    return A.ATcT(base=base)


def _fill(store, h0, *, level=LEVEL, freq=None, zpe_scale=1.0, thermal=None):
    """Records whose E + scale*ZPE equals the given H0 per species."""
    for smiles, value in h0.items():
        zpe = 10.0
        E = value - zpe_scale * zpe
        if freq is None:
            store.put(
                FragmentRecord.make(
                    smiles,
                    multiplicity=1 + smiles.count("[CH"),
                    E=E,
                    ZPE=zpe,
                    H=E + zpe + (thermal or {}).get(smiles, 5.0),
                    model_chemistry=level,
                )
            )
        else:
            store.put(FragmentRecord.make(smiles, E=E, model_chemistry=level))
            store.put(
                FragmentRecord.make(
                    smiles,
                    E=E + 100.0,  # a different level's energy
                    ZPE=zpe,
                    H=E + 100.0 + zpe + (thermal or {}).get(smiles, 5.0),
                    model_chemistry=freq,
                )
            )


def test_self_consistent_ladder_reproduces_atct(tmp_path, base, tables):
    """High-level values equal to ATcT's own give laddered values equal to
    ATcT, a zero method term, and unit coefficients that reduce to the
    anchors with the reaction stoichiometry."""
    h0 = {s: tables.lookup(s)["dfH0_0K"] for s in ("[H][H]", "C", "CC", "CCC", "CCCC")}
    with FragmentStore(tmp_path / "f.db", create=True) as store:
        _fill(store, h0)
    ladder = L.build_ladder(
        ["CCCC"],
        name="test",
        energy_level=LEVEL,
        atct_base=base,
        store_path=tmp_path / "f.db",
    )
    assert set(ladder.species) == {"[H][H]", "C", "CC", "CCC", "CCCC"}
    assert ladder["C"]["rung"] == 0 and ladder["[H][H]"]["rung"] == 0
    assert ladder["CC"]["rung"] == 1 and ladder["CCC"]["rung"] == 2
    assert ladder["CCCC"]["rung"] == 3
    for s in ("CC", "CCC", "CCCC"):
        assert ladder[s]["H0"] == pytest.approx(h0[s], abs=1e-9)
        assert ladder[s]["atct_difference_0K"] == pytest.approx(0.0, abs=1e-9)
    assert ladder.method_uncertainty["n"] == 3
    assert ladder.method_uncertainty["rms_0K"] == pytest.approx(0.0, abs=1e-9)
    # Ethane via CBH-0: CC + H2 -> 2 CH4, so L(CC) = 2 ATcT(CH4) - ATcT(H2) + DrH
    c = ladder["CC"]["atct_coefficients"]
    assert c == {"74-82-8*0": pytest.approx(2.0), "1333-74-0*0": pytest.approx(-1.0)}
    # Propane via CBH-1: CCC + CH4 -> 2 CC: 2 L(CC) - ATcT(CH4) = 3 CH4 - 2 H2
    c = ladder["CCC"]["atct_coefficients"]
    assert c == {"74-82-8*0": pytest.approx(3.0), "1333-74-0*0": pytest.approx(-2.0)}
    # The ATcT uncertainty: H2 is exact, so 3 x methane's 0.043 in quadrature
    assert ladder["CCC"]["uncertainty"] == pytest.approx(3 * 0.043)
    assert ladder["CC"]["reaction"].startswith("CBH-0:")
    assert ladder["CCCC"]["reaction"].startswith("CBH-2:")
    # Saved and loaded again
    path = ladder.save(tmp_path / "f.db")
    assert path == tmp_path / "ladders" / "test.json"
    again = L.Ladder.load("test", tmp_path / "f.db")
    assert again["CCC"]["H0"] == pytest.approx(h0["CCC"])
    e = again.entry("CCC")
    assert e["source"] == "laddered" and e["dfH0_0K"] == pytest.approx(h0["CCC"])
    assert again.entry("C")["source"] == "ATcT"
    assert again.entry("CCCCC") is None
    assert "CBH-2 CCC" in again.summary()


def test_ladder_measures_the_method_term(tmp_path, base, tables):
    """A level whose reaction enthalpies are off by a constant per CBH-0 step
    gives laddered values drifting from ATcT, and the drift is measured."""
    h0 = {s: tables.lookup(s)["dfH0_0K"] for s in ("[H][H]", "C", "CC", "CCC")}
    h0["CC"] += 1.0  # ethane's CBH-0 reaction is 1 kJ/mol off at this level
    with FragmentStore(tmp_path / "f.db", create=True) as store:
        _fill(store, h0)
    ladder = L.build_ladder(
        ["CCC"],
        name="off",
        energy_level=LEVEL,
        atct_base=base,
        store_path=tmp_path / "f.db",
    )
    assert ladder["CC"]["atct_difference_0K"] == pytest.approx(1.0)
    # Propane's CBH-1 reaction is exact at this level, but it inherits the
    # ethane error twice (2 CC on the right): L(CCC) = ATcT(CCC) + 2 - 2 x 1?
    # No: L(CCC) = H(CCC) - 2 H(CC) + H(C) + 2 L(CC) - ATcT(C); with H = ATcT
    # except H(CC) = ATcT + 1 and L(CC) = ATcT(CC) + 1 the two cancel.
    assert ladder["CCC"]["atct_difference_0K"] == pytest.approx(0.0, abs=1e-9)
    mu = ladder.method_uncertainty
    assert mu["n"] == 2 and mu["max_0K"] == pytest.approx(1.0)
    assert mu["rms_0K"] == pytest.approx((1.0 / 2) ** 0.5)


def test_frequency_level_and_zpe_scale(tmp_path, base, tables):
    """ZPE and thermal corrections from another level's records, scaled."""
    h0 = {s: tables.lookup(s)["dfH0_0K"] for s in ("[H][H]", "C", "CC")}
    with FragmentStore(tmp_path / "f.db", create=True) as store:
        _fill(store, h0, freq=FREQ)
    ladder = L.build_ladder(
        ["CC"],
        name="z",
        energy_level=LEVEL,
        frequency_level=FREQ,
        zpe_scale=0.9,
        atct_base=base,
        store_path=tmp_path / "f.db",
    )
    # Every species' H0 is lowered by 0.1 x ZPE(10) = 1, and the CBH-0
    # reaction CC + H2 -> 2 CH4 conserves the number of molecules, so the
    # shift cancels in the laddered value
    assert ladder["CC"]["H0"] == pytest.approx(h0["CC"], abs=1e-9)
    assert ladder["CC"]["H298"] is not None
    assert ladder.data["frequency_level"] == FREQ and ladder.data["zpe_scale"] == 0.9


def test_missing_rung_is_reported(tmp_path, base, tables):
    h0 = {s: tables.lookup(s)["dfH0_0K"] for s in ("[H][H]", "C", "CCC")}
    with FragmentStore(tmp_path / "f.db", create=True) as store:
        _fill(store, h0)  # no ethane at this level
    ladder = L.build_ladder(
        ["CCC"],
        name="m",
        energy_level=LEVEL,
        atct_base=base,
        store_path=tmp_path / "f.db",
    )
    assert "CCC" not in ladder and "CC" in ladder.data["missing"]
    assert "needs CC" in ladder.data["missing"]["CCC"]
    assert ladder.data["unverified"] == []


def test_cli(tmp_path, base, tables, capsys):
    from seamm_thermochemistry.ladder_cli import main

    h0 = {s: tables.lookup(s)["dfH0_0K"] for s in ("[H][H]", "C", "CC")}
    with FragmentStore(tmp_path / "f.db", create=True) as store:
        _fill(store, h0)
    argv = ["--store", str(tmp_path / "f.db")]
    assert (
        main(
            argv + ["build", "demo", "CC", "--energy", LEVEL, "--atct-base", str(base)]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert "Ladder 'demo'" in out and "unverified" in out
    assert main(argv + ["show", "demo"]) == 0
    assert main(argv + ["list"]) == 0
    assert "demo" in capsys.readouterr().out
    assert main(argv + ["show", "nope"]) == 1
