"""Golden tests for the CBH reaction generator against the reactions published
by Ramabhadran & Raghavachari (2011) and Raghavachari et al. (2023), checked
by hand in the design review (cbh_papers.md, section 2.5)."""

import re

import pytest

rdkit = pytest.importorskip("rdkit")

from seamm_thermochemistry import cbh  # noqa: E402
from seamm_thermochemistry.atct import canonical_smiles  # noqa: E402


def parse(text):
    """'A + 2 B -> C + D' as {canonical SMILES: coefficient}, with the target
    (first reactant) at -1, reactants negative, products positive."""
    left, right = text.split("->")
    result = {}
    for side, sign in ((left, -1), (right, 1)):
        for term in side.split(" + "):
            term = term.strip()
            if not term:
                continue
            match = re.match(r"^(\d+)\s+(.*)$", term)
            if match:
                n, smiles = int(match.group(1)), match.group(2)
            else:
                n, smiles = 1, term
            result[canonical_smiles(smiles)] = (
                result.get(canonical_smiles(smiles), 0) + sign * n
            )
    return result


def check(smiles, rung, expected, protect="none"):
    reactions = cbh.cbh_reactions(smiles, rungs=[rung], protect=protect)
    reaction = reactions[rung]
    assert reaction.balanced, (str(reaction), reaction.balance_errors)
    assert not reaction.trivial, str(reaction)
    assert reaction.coefficients == parse(expected), str(reaction)
    return reaction


# -- Ramabhadran & Raghavachari 2011 ---------------------------------------


def test_cyclohexanone():
    p = "O=C1CCCCC1"
    check(p, 0, f"{p} + 8 [H][H] -> 6 C + O")
    check(p, 1, f"{p} + 7 C -> 6 CC + C=O")
    check(p, 2, f"{p} + 6 CC -> 5 CCC + CC(C)=O")
    check(p, 3, f"{p} + 5 CCC + CC(C)=O -> 4 CCCC + 2 CCC(C)=O")


def test_cyclopentadiene():
    p = "C1=CCC=C1"
    check(p, 0, f"{p} + 7 [H][H] -> 5 C")
    check(p, 1, f"{p} + 5 C -> 3 CC + 2 C=C")
    check(p, 2, f"{p} + 3 CC + 2 C=C -> CCC + 4 C=CC")
    # The ring double bonds give cis-2-butene and s-cis (no stereo) butadiene
    check(p, 3, f"{p} + CCC + 4 C=CC -> 2 C=CCC + C=CC=C + 2 C/C=C\\C")


def test_oxazoline():
    p = "C1COC=N1"
    check(p, 0, f"{p} + 6 [H][H] -> 3 C + N + O")
    check(p, 1, f"{p} + 3 C + N + O -> CC + 2 CO + CN + C=N")
    # The ring C=N is cis, so the formimidic acid fragment keeps H and O cis
    # (the paper writes it without stereo)
    check(p, 2, f"{p} + CC + 2 CO + CN + C=N -> CCN + CCO + [H]/N=C\\O + CN=C + COC")


def test_heptanethiol():
    p = "CCCCCCCS"
    check(p, 0, f"{p} + 7 [H][H] -> 7 C + S")
    check(p, 1, f"{p} + 6 C -> 6 CC + CS")
    check(p, 2, f"{p} + 5 CC -> 5 CCC + CCS")
    check(p, 3, f"{p} + 4 CCC -> 4 CCCC + CCCS")


def test_propyl_pentenoate():
    p = "C=CCCC(=O)OCCC"
    check(p, 0, f"{p} + 11 [H][H] -> 8 C + 2 O")
    check(p, 1, f"{p} + O + 7 C -> 5 CC + C=C + C=O + 2 CO")
    check(p, 2, f"{p} + 4 CC + 2 CO -> C=CC + CC(=O)O + 3 CCC + COC + CCO")
    check(
        p,
        3,
        f"{p} + 2 CCC + COC + CCO + CC(=O)O -> "
        "C=CCC + CCCC + CCC(=O)O + CC(=O)OC + CCOC + CCCO",
    )


# -- Raghavachari et al. 2023 ----------------------------------------------


def test_methionine():
    p = "CSCCC(N)C(=O)O"
    check(p, 0, f"{p} + 9 [H][H] -> 5 C + N + 2 O + S")
    check(p, 1, f"{p} + 6 C + S -> 3 CC + CN + 2 CS + C=O + CO")
    check(p, 2, f"{p} + 3 CC + CS -> CSC + CCS + CCC + CC(C)N + CC(=O)O")
    check(p, 3, f"{p} + CCC + CCS + CC(C)N -> CCSC + CCCS + CCC(C)N + CC(N)C(=O)O")
    check(p, 4, f"{p} + CCCS + CCC(C)N -> CCCSC + CC(N)CCS + CCC(N)C(=O)O")


def test_methyl_3_butenoate():
    p = "C=CCC(=O)OC"
    check(p, 0, f"{p} + 8 [H][H] -> 5 C + 2 O")
    check(p, 2, f"{p} + 2 CC + CO -> C=CC + CCC + CC(=O)O + COC")


def test_dimethylpentene_isodesmic():
    check("CCC(C)=C(C)C", 1, "CCC(C)=C(C)C + 5 C -> C=C + 5 CC")


def test_methyl_butanoate():
    p = "COC(=O)CCC"
    check(p, 2, f"{p} + 2 CC + CO -> COC + CC(=O)O + 2 CCC")
    check(p, 3, f"{p} + CC(=O)O + CCC -> CC(=O)OC + CCC(=O)O + CCCC")


def test_methyl_butanoate_radical():
    # Debnath 2019 / 2023 Figure 10: the radical centre stays in every
    # fragment containing it
    p = "COC(=O)C[CH2]"
    check(p, 2, f"{p} + CC + CO -> COC + CC(=O)O + CC[CH2]")
    check(p, 3, f"{p} + CC(=O)O -> COC(C)=O + OC(=O)C[CH2]")


# -- the design note's examples ----------------------------------------------


def test_fluoroethane_and_trivial_rung():
    p = "CCF"
    reactions = cbh.cbh_reactions(p, rungs=[0, 1, 2, 3], protect="none")
    assert reactions[0].coefficients == parse(f"{p} + 2 [H][H] -> 2 C + F")
    assert reactions[1].coefficients == parse(f"{p} + C -> CF + CC")
    assert reactions[2].trivial
    assert reactions[3].trivial
    assert reactions[1].references == {canonical_smiles("C"): -1, "CF": 1, "CC": 1}
    assert "CBH-1" in str(reactions[1]) and "trivial" in str(reactions[2])


def test_chlorotoluene_protected_ring():
    p = "Cc1ccc(Cl)cc1"
    reactions = cbh.cbh_reactions(p, rungs=[0, 1, 2], protect="default")
    assert reactions[0].coefficients == parse(f"{p} + 2 [H][H] -> c1ccccc1 + C + Cl")
    assert reactions[1].coefficients == parse(
        f"{p} + c1ccccc1 -> Clc1ccccc1 + Cc1ccccc1"
    )
    assert reactions[2].trivial
    assert reactions[1].balanced
    assert any("aromatic" in s for s in reactions[1].protected)


def test_chlorotoluene_kekule():
    p = "Cc1ccc(Cl)cc1"
    check(p, 1, f"{p} + 8 C -> 3 C=C + 4 CC + CCl")
    check(p, 2, f"{p} + 3 C=C + 3 CC -> 4 C=CC + C=C(C)Cl + C=C(C)C")


def test_naphthalene_fused_ring_is_one_group():
    reactions = cbh.cbh_reactions("c1ccc2ccccc2c1", rungs=[0, 1], protect="default")
    # One super-node holding the whole molecule: every rung is trivial
    assert reactions[0].trivial and reactions[1].trivial
    assert len(reactions[0].protected) == 1
    assert "10 atoms" in reactions[0].protected[0]


def test_small_rings():
    # cyclopropane: CBH-1 fine; CBH-2 needs the ring protected
    assert cbh.cbh_reactions("C1CC1", rungs=[1], protect="none")[
        1
    ].coefficients == parse("C1CC1 + 3 C -> 3 CC")
    with pytest.raises(cbh.RungUnavailable, match="at least 4 members"):
        cbh.cbh_reactions("C1CC1", rungs=[2], protect="none")
    # protected by default: the ring is one node, every rung is trivial
    assert cbh.cbh_reactions("C1CC1", rungs=[1, 2])[1].trivial
    # cyclobutane: CBH-2 fine, CBH-3 unavailable
    assert cbh.cbh_reactions("C1CCC1", rungs=[2], protect="none")[
        2
    ].coefficients == parse("C1CCC1 + 4 CC -> 4 CCC")
    with pytest.raises(cbh.RungUnavailable):
        cbh.cbh_reactions("C1CCC1", rungs=[3], protect="none")


def test_nitro_group_is_protected():
    p = "C[N+](=O)[O-]"
    reactions = cbh.cbh_reactions(p, rungs=[0, 1], protect="default")
    # CH3-NO2 contracted: CBH-0 gives CH4 + HNO2 (the capped nitro group)
    assert reactions[0].coefficients == parse(f"{p} + [H][H] -> C + O=[NH+][O-]")
    assert reactions[1].trivial
    assert reactions[0].balanced
    # unprotected, the nitro N alone would be capped to a meaningless NH4+
    bad = cbh.cbh_reactions(p, rungs=[0], protect="none")[0]
    assert canonical_smiles("[NH4+]") in bad.coefficients


def test_stereo_kept_in_fragments():
    # trans-2-pentene: the CBH-3 bond fragment of C2=C3 keeps the trans
    reactions = cbh.cbh_reactions("C/C=C/CC", rungs=[1, 2, 3], protect="none")
    r3 = reactions[3]
    assert not r3.trivial
    assert canonical_smiles("C/C=C/CC") in r3.coefficients  # the target
    assert canonical_smiles("C/C=C/C") in r3.products
    assert canonical_smiles("C/C=C\\C") not in r3.products


def test_charge_and_radical_balance():
    reactions = cbh.cbh_reactions("C[NH3+]", rungs=[0, 1], protect="none")
    assert reactions[0].coefficients == parse("C[NH3+] + [H][H] -> C + [NH4+]")
    assert reactions[0].balanced
    reactions = cbh.cbh_reactions("[CH2]C", rungs=[0, 1], protect="none")
    assert reactions[0].coefficients == parse("[CH2]C + [H][H] -> C + [CH3]")
    assert reactions[1].trivial


def test_unique_fragments_across_targets():
    a = cbh.cbh_reactions("CCF", rungs=[1], protect="none")[1]
    b = cbh.cbh_reactions("CCCl", rungs=[1], protect="none")[1]
    c = cbh.cbh_reactions("CCF", rungs=[2], protect="none")[2]
    fragments = cbh.unique_fragments([a, b, c])
    assert set(fragments) == {canonical_smiles(s) for s in ("C", "CC", "CF", "CCl")}


def test_reaction_to_dict_and_bad_input():
    r = cbh.cbh_reactions("CCF", rungs=[1], protect="none")[1]
    d = r.to_dict()
    assert d["rung"] == 1 and d["balanced"] and not d["trivial"]
    with pytest.raises(cbh.CBHError):
        cbh.cbh_reactions("not smiles", rungs=[1])
    with pytest.raises(cbh.CBHError, match="SMARTS"):
        cbh.cbh_reactions("CCO", rungs=[1], protect=["[[["])


def _dihedral(p0, p1, p2, p3):
    import numpy as np

    p0, p1, p2, p3 = (np.asarray(p) for p in (p0, p1, p2, p3))
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return abs(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def test_lowest_conformer_is_all_anti():
    """The lowest MMFF conformer of the n-alkanes is all-anti (a single
    embedding of pentane or heptane is gauche)."""
    for n in (5, 6, 7, 8):
        symbols, xyz, energies, best = cbh.lowest_conformer("C" * n, n=20)
        assert len(energies) > 1 and energies[best] == min(energies)
        carbons = [i for i, s in enumerate(symbols) if s == "C"]
        for i in range(n - 3):
            d = _dihedral(*[xyz[carbons[i + k]] for k in range(4)])
            assert d > 165.0, f"C{n}: dihedral {i} is {d:.1f}"
    # Dimethyl carbonate: the cis-cis conformer (both methyls syn to C=O)
    symbols, xyz, energies, best = cbh.lowest_conformer("COC(=O)OC", n=20)
    heavy = [i for i, s in enumerate(symbols) if s != "H"]
    assert len(heavy) == 6
    # Radical and ion still embed (UFF fallback when MMFF has no parameters)
    symbols, xyz, energies, best = cbh.lowest_conformer("[CH2]CC", n=5)
    assert len(symbols) == 10
