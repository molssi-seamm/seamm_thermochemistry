# -*- coding: utf-8 -*-

"""SMILES that RDKit refuses but chemistry needs: hypervalent atoms.

RDKit's valence check allows a halogen one bond, so ClF3 (``FCl(F)F``) or BrF5
cannot be read with ``Chem.MolFromSmiles``, bracketed (``F[Cl](F)F``) or not.
:func:`mol_from_smiles` reads them with every sanitization step except the
valence check, so such an atom keeps the bonds it is written with and gets no
implicit hydrogens beyond them; RDKit writes them back unbracketed, which this
function and Open Babel both read correctly.
"""

__all__ = ["mol_from_smiles", "is_hypervalent", "NORMAL_VALENCE"]

#: The lowest (normal) valence of the main-group elements whose hypervalent
#: compounds are common: sulfur, selenium and tellurium (2), phosphorus and
#: arsenic (3), and the halogens chlorine, bromine and iodine (1).
NORMAL_VALENCE = {16: 2, 34: 2, 52: 2, 15: 3, 33: 3, 17: 1, 35: 1, 53: 1}


def _chem():
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    return Chem


def mol_from_smiles(smiles):
    """An RDKit molecule from a SMILES, or None if it cannot be read.

    A SMILES that RDKit reads as it is gives the same molecule as
    ``Chem.MolFromSmiles``. One that fails only RDKit's valence check on an
    element that is hypervalent in real compounds (S, Se, Te, P, As, Cl, Br, I:
    ClF3, BrF5) is read without that check; aromaticity, conjugation and
    hybridization are still perceived. Any other over-valent atom (O=[N]=O)
    still makes the SMILES unreadable.
    """
    Chem = _chem()
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        mol.UpdatePropertyCache(strict=False)
        flags = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
        failed = Chem.SanitizeMol(mol, flags, catchErrors=True)
    except Exception:
        return None
    if failed != Chem.SanitizeFlags.SANITIZE_NONE:
        return None
    # Only the elements that are hypervalent in real compounds may exceed
    # RDKit's valences: an over-valent N, C or O (ATcT writes NO2 as O=[N]=O)
    # is not what it seems, so such a SMILES is still unreadable
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in NORMAL_VALENCE:
            continue
        try:
            atom.UpdatePropertyCache(strict=True)
        except Exception:
            return None
    return mol


def is_hypervalent(atom):
    """Whether an RDKit atom is a neutral or anionic S, Se, Te, P, As, Cl, Br or I
    with more bonds than its normal valence (SF4, SF6, sulfonyl, phosphoryl,
    ClF3, PF6-, ...)."""
    normal = NORMAL_VALENCE.get(atom.GetAtomicNum())
    if normal is None or atom.GetFormalCharge() > 0:
        return False
    return atom.GetTotalValence() > normal
