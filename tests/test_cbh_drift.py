"""Golden case: the drift of CBH-2 enthalpies of formation along the n-alkanes
when the fragments take ATcT values (Elliott et al., J. Phys. Chem. A 127, 1512
(2023), Table 7). Needs the real ATcT ver. 1.220 cache, which is not
distributed with the package, so the test is skipped without it.

The 0 K enthalpies at one level of theory come from Elliott's Table 7: ethane
(ANL1') and propane (ANL0+) reference values, n-butane to n-heptane CBH-2 at
CCSD(T)-F12b/cc-pVTZ-F12//B2PLYP-D3/cc-pVTZ. Any CBH-2 reaction among them
returns that method's reaction enthalpy. Referenced to ATcT 1.220 propane and
ethane, the results drift upward against ATcT's own alkane values by +0.38,
+1.18, +1.61 and +1.72 kJ/mol for C4 to C7 (science, 2026-10-28). With the
ATcT covariance block for ethane to octane (cached 2026-10-28 from a browser
download) the 95 % uncertainty of the reference combination is 0.23, 0.32,
0.43 and 0.53 kJ/mol, smaller than the quadrature 0.32, 0.51, 0.70 and 0.89
because the alkane values are positively correlated, and the uncertainty of the
drift itself (references minus the target) is 0.18, 0.37, 0.47 and 0.63. The
drift exceeds those limits for all four, but sits inside CBH-ANL's own claimed
2σ, so it does not say which side is wrong. Without the block the quadrature
numbers are checked instead.
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
# 95 % uncertainty of the reference combination and of the drift, kJ/mol
EXPECTED_U_REF = {
    "covariance": {4: 0.23, 5: 0.32, 6: 0.43, 7: 0.53},
    "quadrature": {4: 0.32, 5: 0.51, 6: 0.70, 7: 0.89},
}
EXPECTED_U_DRIFT = {
    "covariance": {4: 0.18, 5: 0.37, 6: 0.47, 7: 0.63},
    "quadrature": {4: 0.37, 5: 0.59, 6: 0.77, 7: 0.99},
}


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
        ids = [refs[s]["id"] for s in rxn.references]
        coefficients = list(rxn.references.values())
        uncertainties = [refs[s]["uncertainty_298K"] for s in rxn.references]
        u, used = reaction_uncertainty(
            coefficients, uncertainties, tables.covariance(ids, fetch=False)
        )
        method = "covariance" if used else "quadrature"
        assert u == pytest.approx(EXPECTED_U_REF[method][n], abs=0.005)
        udrift, used = reaction_uncertainty(
            coefficients + [-1.0],
            uncertainties + [tables.lookup(target)["uncertainty_298K"]],
            tables.covariance(ids + [tables.lookup(target)["id"]], fetch=False),
        )
        method = "covariance" if used else "quadrature"
        assert udrift == pytest.approx(EXPECTED_U_DRIFT[method][n], abs=0.005)
        # The drift is outside the 95 % limits: of the references alone for
        # C5+ either way, and of the drift itself for C4 too with the covariance
        if n >= 5:
            assert dfH - own > u
        if used:
            assert dfH - own > udrift
