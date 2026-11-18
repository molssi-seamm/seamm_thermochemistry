# -*- coding: utf-8 -*-

"""Double-bond (E/Z) isomers: one rule for every source of reference values.

The CBH code keys species by canonical SMILES, which carry the E/Z of each
stereogenic double bond. A value for one isomer is not a value for the other,
and the value of a mixture, or of a species recorded without saying which isomer
it is, is not the value of a particular isomer: for formimidic acid the isomers
differ by up to 21 kJ/mol. So, for every source (ATcT, the reference library and
the computed values):

1. A value recorded under the same SMILES, E/Z included, is that isomer's.
2. A query with a defined E/Z may also take a value recorded without E/Z when
   the source labels it as that isomer -- ATcT's descriptors "trans" and "E",
   "cis" and "Z", and for formimidic acid "anti-..." and "syn-..." (see
   :data:`ISOMER_WORDS`). If the source has several labelled values for the
   isomer (its conformers), the lowest is taken: by the reference convention,
   the reference value of conformers that share one SMILES (here one SMILES with
   its E/Z) is that of the species' ground conformer.
3. Otherwise a query with a defined E/Z gets no value from that source, and the
   next source is tried; the generic value is never used in its place.

Only E/Z isomers are covered. Conformer labels (syn, anti, gauche, ...) of
species without E/Z, and queries without E/Z, are handled as before. The
reference library and the computed values record a species under its exact
SMILES and have no isomer labels, so for them rule 1 is the whole rule.
"""

import re

#: The words of a descriptor that name the E/Z isomer, in general
ISOMER_WORDS = {"trans": "E", "e": "E", "cis": "Z", "z": "Z"}

#: Descriptors of particular species (by CAS group of the ATcT ID) whose words
#: are not cis/trans. Formimidic acid HN=CH-OH: the first word is the N-H relative
#: to the O across C=N (anti = E), the second the O-H conformer.
SPECIES_ISOMER_WORDS = {"60100-09-6": (0, {"anti": "E", "syn": "Z"})}


def _rdkit():
    from rdkit import Chem

    return Chem


def double_bond_isomers(smiles):
    """The E/Z (CIP) labels of the double bonds whose stereo the SMILES defines.

    Returns
    -------
    list of str
        "E" or "Z" for each such double bond; empty if there is none or RDKit
        cannot read the SMILES.
    """
    from .smiles import mol_from_smiles

    Chem = _rdkit()
    from rdkit.Chem import rdCIPLabeler

    mol = mol_from_smiles(smiles)
    if mol is None:
        return []
    stereo = (
        Chem.BondStereo.STEREOE,
        Chem.BondStereo.STEREOZ,
        Chem.BondStereo.STEREOCIS,
        Chem.BondStereo.STEREOTRANS,
    )
    bonds = [
        b
        for b in mol.GetBonds()
        if b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() in stereo
    ]
    if not bonds:
        return []
    rdCIPLabeler.AssignCIPLabels(mol)
    labels = []
    for bond in bonds:
        if bond.HasProp("_CIPCode"):
            labels.append(bond.GetProp("_CIPCode"))
    return labels


def descriptor_isomer(descriptor, atct_id=None):
    """The E/Z isomer that a descriptor (e.g. ATcT's "trans") names, or None.

    Parameters
    ----------
    descriptor : str
        The descriptor, e.g. "trans", "cis", "anti-syn".
    atct_id : str, optional
        The ATcT ID, for species whose descriptors use other words
        (:data:`SPECIES_ISOMER_WORDS`).

    Returns
    -------
    str or None
        "E", "Z", or None if the descriptor names no single isomer (e.g. an
        equilibrium mixture, or no descriptor).
    """
    if not descriptor or "equilib" in descriptor.lower():
        return None
    words = [w for w in re.split(r"[\s,;()\-]+", descriptor.lower()) if w]
    if atct_id is not None:
        group = str(atct_id).split("*")[0]
        if group in SPECIES_ISOMER_WORDS:
            position, mapping = SPECIES_ISOMER_WORDS[group]
            if len(words) > position:
                return mapping.get(words[position])
            return None
    found = {ISOMER_WORDS[w] for w in words if w in ISOMER_WORDS}
    return found.pop() if len(found) == 1 else None
