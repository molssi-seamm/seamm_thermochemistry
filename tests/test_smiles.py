"""Reading SMILES that RDKit's valence check refuses (hypervalent atoms)."""

import pytest

pytest.importorskip("rdkit")

from seamm_thermochemistry.smiles import is_hypervalent, mol_from_smiles  # noqa: E402


def test_plain_smiles_are_unchanged():
    from rdkit import Chem

    for smiles in ("CCO", "c1ccccc1", "[CH3]", "O=[SH2]=O"):
        assert Chem.MolToSmiles(mol_from_smiles(smiles)) == Chem.MolToSmiles(
            Chem.MolFromSmiles(smiles)
        )


def test_hypervalent_halogens_have_no_hydrogens():
    for smiles, n in (("FCl(F)F", 4), ("FBr(F)(F)(F)F", 6)):
        mol = mol_from_smiles(smiles)
        assert mol.GetNumAtoms() == n
        assert all(a.GetTotalNumHs() == 0 for a in mol.GetAtoms())


def test_is_hypervalent():
    mol = mol_from_smiles("CS(C)(=O)=O")
    assert [a.GetSymbol() for a in mol.GetAtoms() if is_hypervalent(a)] == ["S"]
    for smiles in ("CSC", "CCl", "C[S+](C)C"):
        assert not any(is_hypervalent(a) for a in mol_from_smiles(smiles).GetAtoms())


def test_garbage_is_none():
    assert mol_from_smiles("C1CC(") is None


def test_only_hypervalent_elements_may_exceed():
    """ATcT writes NO2 as O=[N]=O: an over-valent nitrogen is not read as
    closed-shell O=N=O."""
    assert mol_from_smiles("O=[N]=O") is None
