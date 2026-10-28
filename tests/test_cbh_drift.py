"""Golden case: the drift of CBH-2 enthalpies of formation along the n-alkanes
when the fragments take ATcT values (Elliott et al., J. Phys. Chem. A 127, 1512
(2023), Table 7). Needs the real ATcT ver. 1.220 cache, which is not
distributed with the package, so the test is skipped without it.

The 0 K enthalpies at one level of theory come from Elliott's Table 7: ethane
(ANL1') and propane (ANL0+) reference values, n-butane to n-heptane CBH-2 at
CCSD(T)-F12b/cc-pVTZ-F12//B2PLYP-D3/cc-pVTZ. Any CBH-2 reaction among them
returns that method's reaction enthalpy. Referenced to ATcT 1.220 propane and
ethane, the results drift upward against ATcT's own alkane values by +0.38,
+1.18, +1.61 and +1.72 kJ/mol for C4 to C7 (science, 2026-10-28; quadrature
over the 298 K 95 % uncertainties). The drift exceeds the ATcT reference +
target uncertainties but sits inside CBH-ANL's own claimed 2σ, so it does not
say which side is wrong. To be rerun with the covariance block once cached.
"""

import pytest

rdkit = pytest.importorskip("rdkit")

from seamm_thermochemistry import ATcT, cbh, reaction_uncertainty  # noqa: E402
from seamm_thermochemistry.atct import canonical_smiles  # noqa: E402

# kJ/mol at 0 K, Elliott 2023 Table 7 (kcal/mol x 4.184)
ENERGIES = {
    "CC": -68.977,
    "CCC": -83.475,
    "CCCC": -98.784,
    "CCCCC": -114.098,
    "CCCCCC": -129.578,
    "CCCCCCC": -145.143,
}

EXPECTED_DRIFT = {4: 0.38, 5: 1.18, 6: 1.61, 7: 1.72}
EXPECTED_DRH = {4: 0.812, 5: 1.628, 6: 2.611, 7: 3.678}


@pytest.fixture(scope="module")
def tables():
    try:
        t = ATcT(version="1.220")
    except FileNotFoundError:
        pytest.skip("ATcT ver. 1.220 is not cached on this machine")
    return t


def test_alkane_cbh2_drift(tables):
    values = {canonical_smiles(s): v for s, v in ENERGIES.items()}
    for n in (4, 5, 6, 7):
        target = canonical_smiles("C" * n)
        rxn = cbh.cbh_reactions(target, rungs=[2], protect="none")[2]
        assert rxn.balanced and not rxn.trivial
        drH = sum(c * values[s] for s, c in rxn.coefficients.items())
        assert drH == pytest.approx(EXPECTED_DRH[n], abs=0.005)
        refs = {s: tables.lookup(s) for s in rxn.references}
        ref = sum(c * refs[s]["dfH0_0K"] for s, c in rxn.references.items())
        dfH = ref - drH
        own = tables.lookup(target)["dfH0_0K"]
        assert dfH - own == pytest.approx(EXPECTED_DRIFT[n], abs=0.05)
        u, used = reaction_uncertainty(
            list(rxn.references.values()),
            [refs[s]["uncertainty_298K"] for s in rxn.references],
            tables.covariance([refs[s]["id"] for s in rxn.references], fetch=False),
        )
        # The drift is outside the quadrature reference uncertainty for C5+
        if n >= 5:
            assert dfH - own > u
