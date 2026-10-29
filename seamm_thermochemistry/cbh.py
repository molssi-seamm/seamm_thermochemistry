# -*- coding: utf-8 -*-

"""Connectivity-based hierarchy (CBH) reactions for a molecule.

The CBH of Ramabhadran and Raghavachari (J. Chem. Theory Comput. 7, 2094
(2011); Raghavachari et al., J. Chem. Theory Comput. 19, 3763 (2023)) builds,
from the hydrogen-suppressed molecular graph alone, one error-cancelling
reaction per rung:

=======  ==================================================  ================
Rung     Products (capped with H to the parent's valence)    Reactants
=======  ==================================================  ================
CBH-0    the hydride of every heavy atom                     H2, to balance H
CBH-1    every heavy-atom bond                               N0(atom) x (deg-1)
CBH-2    every heavy atom with its heavy neighbours          N0(bond) per bond
CBH-3    every bond with the neighbours of both atoms        N1(atom) x (deg-1)
CBH-4    every atom with two shells of neighbours            N1(bond) per bond
=======  ==================================================  ================

Here N_k(atom) is the atom with all heavy atoms within k bonds, N_k(bond) is
the bond's two atoms with all heavy atoms within k bonds of either, and deg
is the number of heavy neighbours. Odd rungs (2k+1) have N_k(bond) products
and N_k(atom) reactants with coefficient (deg - 1); even rungs (2k >= 2) have
N_k(atom) products and N_{k-1}(bond) reactants with coefficient 1. Identical
fragments are merged and species on both sides cancel. The target itself has
coefficient -1. This reproduces every published reaction we checked (see
tests/test_cbh.py).

Every fragment is the induced subgraph: it keeps the bonds, bond orders,
formal charges and unpaired electrons its atoms have in the parent, and every
cut bond is replaced by hydrogens that restore the atom's valence in the
parent. Double-bond stereochemistry and the relative configuration of
stereocentres are kept, so cis and trans fragments are different species.

Groups that should not be broken (aromatic ring systems, three- and
four-membered rings, charge-separated or hypervalent groups such as nitro and
sulfonyl, and whole inorganic ions) are contracted into single nodes of the
graph first ("coarse-grained" or "ring-preserved" CBH). The hierarchy on the
contracted graph is shallower: for 4-chlorotoluene with the ring protected,
CBH-1 is Cl-C6H4-CH3 + C6H6 -> C6H5Cl + C6H5CH3 and CBH-2 is already trivial.

Limitations: the simple overlap formula is valid only when every unprotected
ring of size r satisfies rung <= r - 2; such a rung is reported as
unavailable. Aromatic rings that are not protected are Kekulized, which makes
CBH-3 and higher depend on the Kekule structure chosen. A delocalized charge
or spin must be localized on one atom by the input SMILES.
"""

from collections import Counter, deque
import logging

from .atct import canonical_smiles

__all__ = [
    "CBHError",
    "RungUnavailable",
    "Reaction",
    "cbh_reactions",
    "unique_fragments",
    "DEFAULT_PROTECT",
    "PROTECT_SMARTS",
]

logger = logging.getLogger(__name__)

# SMARTS for groups whose hydrogen-capped fragments are meaningless: they are
# contracted into one node. Overlapping matches are merged.
PROTECT_SMARTS = {
    "nitro": "[N+](=O)[O-]",
    "nitro (hypervalent)": "N(=O)=O",
    "N-oxide": "[N+;!$(N=O)][O-]",
    "azide": "N=[N+]=[N-]",
    "diazo": "[C,N]=[N+]=[N-]",
    "diazonium": "[N+]#N",
    "sulfoxide": "[#16X3](=O)([#6,#7])[#6,#7]",
    "sulfonyl": "[#16X4](=O)(=O)",
    "sulfinyl (charge separated)": "[S+]([O-])",
    "phosphoryl": "[#15X4](=O)",
    "phosphonium": "[P+]",
    "PF6-": "[P-](F)(F)(F)(F)(F)F",
    "BF4-": "[B-](F)(F)(F)F",
    "perchlorate": "[Cl](=O)(=O)(=O)[O-]",
}

DEFAULT_PROTECT = ("aromatic", "small rings", "groups")


class CBHError(Exception):
    """A problem generating a CBH reaction."""


class RungUnavailable(CBHError):
    """The rung cannot be built for this molecule (small unprotected rings or
    parallel edges in the contracted graph)."""


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit import RDLogger
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "RDKit is needed for CBH reactions: install seamm-thermochemistry[cbh]"
        ) from e
    RDLogger.DisableLog("rdApp.*")
    return Chem


class Reaction:
    """One CBH reaction.

    Attributes
    ----------
    rung : int
    coefficients : {str: int}
        Stoichiometric coefficients by canonical SMILES (products positive,
        reactants negative, the target -1). Empty if the rung is trivial.
    fragments : {str: str}
        The fragment SMILES as generated (with stereo) by canonical SMILES.
    trivial : bool
        The net reaction is empty: the rung reproduces the target.
    balanced : bool
    balance_errors : [str]
    protected : [str]
        Descriptions of the groups contracted into single nodes.
    """

    def __init__(self, rung, target):
        self.rung = rung
        self.target = target
        self.coefficients = {}
        self.fragments = {}
        self.trivial = False
        self.balanced = True
        self.balance_errors = []
        self.protected = []

    @property
    def products(self):
        return {s: c for s, c in self.coefficients.items() if c > 0}

    @property
    def reactants(self):
        return {s: -c for s, c in self.coefficients.items() if c < 0}

    @property
    def references(self):
        """The species other than the target, with their coefficients."""
        return {s: c for s, c in self.coefficients.items() if s != self.target}

    def __str__(self):
        if self.trivial:
            return f"CBH-{self.rung}: trivial (reproduces the target)"

        def side(terms):
            out = []
            for smiles, c in terms.items():
                out.append(f"{c} {smiles}" if c != 1 else smiles)
            return " + ".join(out)

        return f"CBH-{self.rung}: {side(self.reactants)} -> {side(self.products)}"

    def to_dict(self):
        return {
            "rung": self.rung,
            "target": self.target,
            "coefficients": dict(self.coefficients),
            "fragments": dict(self.fragments),
            "trivial": self.trivial,
            "balanced": self.balanced,
            "balance errors": list(self.balance_errors),
            "protected": list(self.protected),
        }


# ---------------------------------------------------------------------------
# The contracted graph
# ---------------------------------------------------------------------------


class _Graph:
    """The hydrogen-suppressed graph with protected groups contracted.

    nodes : [frozenset of atom indices]
    adjacency : {node index: {node index: n_bonds}}
    """

    def __init__(self, mol, protect):
        Chem = _rdkit()
        self.mol = mol
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        groups, self.protected = _protected_groups(mol, protect)

        # Each heavy atom belongs to at most one group
        owner = {}
        nodes = []
        for group in groups:
            nodes.append(frozenset(group))
            for idx in group:
                owner[idx] = len(nodes) - 1
        for idx in heavy:
            if idx not in owner:
                nodes.append(frozenset([idx]))
                owner[idx] = len(nodes) - 1
        self.nodes = nodes
        self.owner = owner

        self.adjacency = {i: {} for i in range(len(nodes))}
        self.parallel = False
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a not in owner or b not in owner:
                continue
            i, j = owner[a], owner[b]
            if i == j:
                continue
            self.adjacency[i][j] = self.adjacency[i].get(j, 0) + 1
            self.adjacency[j][i] = self.adjacency[j].get(i, 0) + 1
            if self.adjacency[i][j] > 1:
                self.parallel = True
        self.edges = sorted(
            {tuple(sorted((i, j))) for i in self.adjacency for j in self.adjacency[i]}
        )

        # Ring sizes in the contracted graph, from the molecule's rings
        ring_info = mol.GetRingInfo()
        self.ring_sizes = []
        for ring in ring_info.AtomRings():
            nodes_in_ring = {owner[idx] for idx in ring if idx in owner}
            if len(nodes_in_ring) > 1:
                self.ring_sizes.append(len(nodes_in_ring))
        _ = Chem  # silence unused

    def degree(self, i):
        return len(self.adjacency[i])

    def neighbourhood(self, start, k):
        """The nodes within k steps of the nodes in `start`."""
        seen = set(start)
        frontier = deque((s, 0) for s in start)
        while frontier:
            node, d = frontier.popleft()
            if d == k:
                continue
            for nb in self.adjacency[node]:
                if nb not in seen:
                    seen.add(nb)
                    frontier.append((nb, d + 1))
        return seen

    def atoms(self, node_set):
        atoms = set()
        for i in node_set:
            atoms |= self.nodes[i]
        return atoms


def _protected_groups(mol, protect):
    """The atom groups to contract, merged where they overlap, and their
    descriptions."""
    Chem = _rdkit()
    if protect is None or protect == "none" or protect == ():
        return [], []
    if isinstance(protect, str):
        protect = (protect,)

    groups = []
    labels = []
    ring_info = mol.GetRingInfo()
    options = set()
    smarts = []
    for item in protect:
        if item in ("aromatic", "small rings", "groups", "all rings"):
            options.add(item)
        elif item == "default":
            options.update(DEFAULT_PROTECT)
        else:
            smarts.append(item)

    if "aromatic" in options or "all rings" in options:
        for ring in ring_info.BondRings():
            bonds = [mol.GetBondWithIdx(b) for b in ring]
            if "all rings" in options or all(b.GetIsAromatic() for b in bonds):
                atoms = set()
                for b in bonds:
                    atoms.add(b.GetBeginAtomIdx())
                    atoms.add(b.GetEndAtomIdx())
                groups.append(atoms)
                labels.append(f"{len(atoms)}-membered aromatic ring")
    if "small rings" in options:
        for ring in ring_info.AtomRings():
            if len(ring) <= 4:
                groups.append(set(ring))
                labels.append(f"{len(ring)}-membered ring")
    if "groups" in options:
        for name, pattern in PROTECT_SMARTS.items():
            query = Chem.MolFromSmarts(pattern)
            for match in mol.GetSubstructMatches(query):
                atoms = {
                    idx for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1
                }
                groups.append(atoms)
                labels.append(name)
    for pattern in smarts:
        query = Chem.MolFromSmarts(pattern)
        if query is None:
            raise CBHError(f"Invalid SMARTS to protect: '{pattern}'")
        for match in mol.GetSubstructMatches(query):
            atoms = {idx for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1}
            groups.append(atoms)
            labels.append(f"SMARTS {pattern}")

    # Merge overlapping groups (fused aromatic rings, nitro on a ring, ...)
    merged = []
    merged_labels = []
    for atoms, label in zip(groups, labels):
        hits = [i for i, g in enumerate(merged) if g & atoms]
        if not hits:
            merged.append(set(atoms))
            merged_labels.append([label])
        else:
            first = hits[0]
            for i in reversed(hits[1:]):
                merged[first] |= merged.pop(i)
                merged_labels[first].extend(merged_labels.pop(i))
            merged[first] |= atoms
            merged_labels[first].append(label)
    descriptions = []
    for atoms, names in zip(merged, merged_labels):
        counts = Counter(names)
        text = ", ".join(f"{n} x {k}" if n > 1 else k for k, n in counts.items())
        descriptions.append(f"{text} ({len(atoms)} atoms)")
    return merged, descriptions


# ---------------------------------------------------------------------------
# Fragments
# ---------------------------------------------------------------------------


def _fragment_smiles(mol_h, keep):
    """The SMILES of the fragment induced by the heavy atoms in `keep`, every
    cut bond replaced by hydrogens restoring the parent valence.

    `mol_h` is the Kekulized parent with explicit hydrogens. A cut neighbour
    is turned into a hydrogen in place, so the neighbour order of each atom,
    and with it the stereo tags, is unchanged; a cut multiple bond gets extra
    hydrogens.
    """
    Chem = _rdkit()
    rw = Chem.RWMol(mol_h)
    keep = set(keep)
    anchor = {}  # capping atom -> the kept atom it caps
    extra_bonds = []  # (kept atom, new H) for cut atoms claimed twice
    for bond in list(rw.GetBonds()):
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        for a, other in ((a1, a2), (a2, a1)):
            if (
                a.GetIdx() in keep
                and other.GetIdx() not in keep
                and (other.GetAtomicNum() > 1 or other.GetIdx() in anchor)
            ):
                order = int(round(bond.GetBondTypeAsDouble()))
                if other.GetIdx() in anchor:
                    # Already a cap for another kept atom (a ring closure,
                    # and so already turned into a hydrogen): this kept atom
                    # gets fresh hydrogens instead
                    extra_bonds.append((a.GetIdx(), other.GetIdx(), order))
                    continue
                anchor[other.GetIdx()] = a.GetIdx()
                other.SetAtomicNum(1)
                other.SetFormalCharge(0)
                other.SetNumRadicalElectrons(0)
                other.SetIsAromatic(False)
                other.SetNoImplicit(True)
                other.SetNumExplicitHs(0)
                other.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
                for _ in range(order - 1):
                    h = rw.AddAtom(Chem.Atom(1))
                    rw.AddBond(a.GetIdx(), h, Chem.BondType.SINGLE)
    for a_idx, cap_idx, order in extra_bonds:
        rw.RemoveBond(a_idx, cap_idx)
        for _ in range(order):
            h = rw.AddAtom(Chem.Atom(1))
            rw.AddBond(a_idx, h, Chem.BondType.SINGLE)
    # A capping atom keeps only its bond to its anchor
    for cap, a_idx in anchor.items():
        for nb in list(rw.GetAtomWithIdx(cap).GetNeighbors()):
            if nb.GetIdx() != a_idx:
                rw.RemoveBond(cap, nb.GetIdx())
    # A kept atom that was a stereocentre in the parent but now carries two
    # or more hydrogens is no longer one
    for idx in keep:
        atom = rw.GetAtomWithIdx(idx)
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            n_h = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1)
            if n_h >= 2:
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    keep_all = set(keep) | set(anchor)
    for idx in keep:
        for nb in rw.GetAtomWithIdx(idx).GetNeighbors():
            if nb.GetAtomicNum() == 1:
                keep_all.add(nb.GetIdx())
    for idx in sorted(
        (a.GetIdx() for a in rw.GetAtoms() if a.GetIdx() not in keep_all), reverse=True
    ):
        rw.RemoveAtom(idx)
    frag = rw.GetMol()
    Chem.SanitizeMol(frag)
    frag = Chem.RemoveHs(frag)
    # Keep the stereo flags carried over from the parent (cleanIt would drop
    # the cis/trans set for ring double bonds)
    Chem.AssignStereochemistry(frag, cleanIt=False, force=True)
    frag.SetProp("_StereochemDone", "1")
    return Chem.MolToSmiles(frag)


def _composition(smiles):
    """(element counts including H, charge, radical electrons, heavy-bond
    counts) of a SMILES."""
    Chem = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise CBHError(f"Cannot parse the fragment SMILES '{smiles}'")
    Chem.Kekulize(mol, clearAromaticFlags=True)
    elements = Counter()
    charge = 0
    radicals = 0
    for atom in mol.GetAtoms():
        elements[atom.GetSymbol()] += 1
        elements["H"] += atom.GetTotalNumHs()
        charge += atom.GetFormalCharge()
        radicals += atom.GetNumRadicalElectrons()
    bonds = Counter()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() > 1 and b.GetAtomicNum() > 1:
            key = tuple(sorted((a.GetSymbol(), b.GetSymbol()))) + (
                bond.GetBondTypeAsDouble(),
            )
            bonds[key] += 1
    return elements, charge, radicals, bonds


# ---------------------------------------------------------------------------
# The reactions
# ---------------------------------------------------------------------------


def _prepare(smiles):
    """The parent as a Kekulized RDKit mol with explicit H, plus the version
    with aromatic flags for ring perception."""
    Chem = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise CBHError(f"RDKit cannot parse the SMILES '{smiles}'")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    kek = Chem.Mol(mol)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    _ring_double_bonds_are_cis(kek)
    mol_h = Chem.AddHs(kek)
    return mol, mol_h


def _ring_double_bonds_are_cis(mol):
    """A double bond in a ring of up to seven atoms is cis with respect to the
    ring. SMILES carries no stereo for it, so set it, so that a fragment that
    cuts the ring open (2-butene from cyclopentadiene) keeps the geometry."""
    Chem = _rdkit()
    ring_info = mol.GetRingInfo()
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            continue
        sizes = [len(ring) for ring in ring_info.BondRings() if bond.GetIdx() in ring]
        if not sizes or min(sizes) > 7:
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        # The ring neighbours on each end (ring bonds other than this one)
        ring = min((r for r in ring_info.BondRings() if bond.GetIdx() in r), key=len)
        ring_atoms = set()
        for idx in ring:
            rb = mol.GetBondWithIdx(idx)
            ring_atoms.add(rb.GetBeginAtomIdx())
            ring_atoms.add(rb.GetEndAtomIdx())
        na = [
            n.GetIdx()
            for n in a.GetNeighbors()
            if n.GetIdx() != b.GetIdx() and n.GetIdx() in ring_atoms
        ]
        nb = [
            n.GetIdx()
            for n in b.GetNeighbors()
            if n.GetIdx() != a.GetIdx() and n.GetIdx() in ring_atoms
        ]
        if len(na) != 1 or len(nb) != 1:
            continue
        bond.SetStereoAtoms(na[0], nb[0])
        bond.SetStereo(Chem.BondStereo.STEREOCIS)


def cbh_reactions(smiles, rungs=(0, 1, 2, 3, 4), protect="default"):
    """The CBH reactions of a molecule.

    Parameters
    ----------
    smiles : str
        The target, with its double-bond stereochemistry, formal charges and
        radical centres. Aromatic rings may be written aromatic.
    rungs : iterable of int
        The rungs wanted.
    protect : "default" | "none" | sequence
        Which groups to contract into single nodes: "aromatic" (aromatic ring
        systems, fused rings merged), "small rings" (3- and 4-membered),
        "groups" (the charge-separated and hypervalent groups in
        `PROTECT_SMARTS`), "all rings", or SMARTS patterns. "default" is
        aromatic + small rings + groups.

    Returns
    -------
    {int: Reaction}
        A Reaction for each rung. A rung that cannot be built on this graph
        raises `RungUnavailable` unless it is above a trivial rung, in which
        case it is trivial too.
    """
    Chem = _rdkit()
    mol, mol_h = _prepare(smiles)
    target = canonical_smiles(smiles)
    graph = _Graph(mol, DEFAULT_PROTECT if protect == "default" else protect)

    if not graph.nodes:
        raise CBHError("The molecule has no heavy atoms")

    cache = {}

    def fragment(node_set):
        key = frozenset(node_set)
        if key not in cache:
            atoms = graph.atoms(key)
            generated = _fragment_smiles(mol_h, atoms)
            cache[key] = (canonical_smiles(generated), generated)
        return cache[key]

    parent_elements, parent_charge, parent_radicals, parent_bonds = _composition(smiles)

    results = {}
    trivial_from = None
    for rung in sorted(rungs):
        reaction = Reaction(rung, target)
        reaction.protected = list(graph.protected)
        results[rung] = reaction
        if trivial_from is not None and rung > trivial_from:
            reaction.trivial = True
            continue

        if rung > 0 and graph.parallel:
            raise RungUnavailable(
                f"CBH-{rung}: a group is bonded twice to the same protected group "
                "(parallel edges in the contracted graph); protect it too or use "
                "protect='none'"
            )
        too_small = [r for r in graph.ring_sizes if rung > r - 2]
        if rung > 1 and too_small:
            raise RungUnavailable(
                f"CBH-{rung} needs every unprotected ring to have at least "
                f"{rung + 2} members; this molecule has rings of size "
                f"{sorted(set(too_small))}. Protect them (small rings are "
                "protected by default) or use a lower rung."
            )

        terms = Counter()
        if rung == 0:
            for i in range(len(graph.nodes)):
                terms[fragment({i})] += 1
        elif rung % 2 == 1:
            k = (rung - 1) // 2
            for i, j in graph.edges:
                terms[fragment(graph.neighbourhood({i, j}, k))] += 1
            for i in range(len(graph.nodes)):
                coefficient = graph.degree(i) - 1
                if coefficient > 0:
                    terms[fragment(graph.neighbourhood({i}, k))] -= coefficient
        else:
            k = rung // 2
            for i in range(len(graph.nodes)):
                terms[fragment(graph.neighbourhood({i}, k))] += 1
            for i, j in graph.edges:
                terms[fragment(graph.neighbourhood({i, j}, k - 1))] -= 1

        # Collect by canonical SMILES; remember the generated form
        coefficients = Counter()
        for (can, generated), c in terms.items():
            coefficients[can] += c
            reaction.fragments.setdefault(can, generated)
        coefficients[target] -= 1
        coefficients = {s: c for s, c in coefficients.items() if c != 0}

        if all(c == 0 for s, c in coefficients.items() if s != target):
            # Everything cancelled: the rung reproduces the target (or the
            # hierarchy has saturated on a contracted graph with no edges)
            reaction.trivial = True
            trivial_from = rung
            continue
        if target not in coefficients:
            # The target cancelled against its own fragment but other species
            # remain: an inconsistent rung (should not happen)
            reaction.coefficients = coefficients
            reaction.balanced = False
            reaction.balance_errors.append(
                "the target cancels but other species remain"
            )
            continue

        if rung == 0:
            # H2 balances the hydrogens
            n_h = -parent_elements["H"]
            for s, c in coefficients.items():
                if s != target:
                    n_h += c * _composition(s)[0]["H"]
            if n_h % 2 != 0:
                raise CBHError(f"CBH-0 of {smiles}: an odd hydrogen imbalance of {n_h}")
            if n_h != 0:
                h2 = canonical_smiles("[H][H]")
                coefficients[h2] = coefficients.get(h2, 0) - n_h // 2
                reaction.fragments.setdefault(h2, "[H][H]")

        reaction.coefficients = dict(sorted(coefficients.items(), key=lambda x: -x[1]))
        _check_balance(
            reaction, parent_elements, parent_charge, parent_radicals, parent_bonds
        )
    _ = Chem
    return results


def _check_balance(reaction, elements, charge, radicals, bonds):
    """Elements, charge, unpaired electrons, and for rung >= 1 the heavy-atom
    bond types, must balance."""
    net_elements = Counter()
    net_charge = 0
    net_radicals = 0
    net_bonds = Counter()
    for smiles, c in reaction.coefficients.items():
        e, q, r, b = _composition(smiles)
        for k, v in e.items():
            net_elements[k] += c * v
        net_charge += c * q
        net_radicals += c * r
        for k, v in b.items():
            net_bonds[k] += c * v
    errors = []
    bad = {k: v for k, v in net_elements.items() if v != 0}
    if bad:
        errors.append(f"elements do not balance: {bad}")
    if net_charge != 0:
        errors.append(f"charge does not balance: {net_charge}")
    if net_radicals != 0:
        errors.append(f"unpaired electrons do not balance: {net_radicals}")
    if reaction.rung >= 1:
        bad = {k: v for k, v in net_bonds.items() if v != 0}
        if bad:
            errors.append(f"heavy-atom bond types do not balance: {bad}")
    reaction.balance_errors = errors
    reaction.balanced = not errors


def unique_fragments(reactions):
    """The union of the reference species over several reactions, with the
    generated (stereo) SMILES, so that a batch computes each fragment once.

    Parameters
    ----------
    reactions : iterable of Reaction

    Returns
    -------
    {str: str}
        canonical SMILES -> generated SMILES
    """
    result = {}
    for reaction in reactions:
        if reaction.trivial:
            continue
        for can in reaction.references:
            result.setdefault(can, reaction.fragments.get(can, can))
    return result


def lowest_conformer(smiles, n=20, seed=0xF00D, return_mol=False):
    """The lowest-energy conformer of `n` ETKDG embeddings optimized with
    MMFF94: the starting structure for a species whose lowest rotamer is
    wanted (all-anti n-alkanes, anti ethanol, cis-cis dimethyl carbonate), as
    CBH-ANL starts every species from its lowest rotamer. A single embedding
    of a long chain is often gauche.

    Returns (symbols, coordinates in Å, mmff energies of all embeddings in
    kcal/mol, lowest index), or with `return_mol` the RDKit molecule (with
    explicit hydrogens and the chosen conformer only) first.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        raise CBHError(f"RDKit cannot parse '{smiles}'")
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=max(1, int(n)), params=params))
    if not ids:
        params.useRandomCoords = True
        ids = list(
            AllChem.EmbedMultipleConfs(mol, numConfs=max(1, int(n)), params=params)
        )
    if not ids:
        raise CBHError(f"RDKit could not embed '{smiles}'")
    energies = {}
    if AllChem.MMFFHasAllMoleculeParams(mol):
        props = AllChem.MMFFGetMoleculeProperties(mol)
        for cid in ids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            ff.Minimize(maxIts=2000)
            energies[cid] = ff.CalcEnergy()
    else:  # no MMFF parameters (some radicals, ions): UFF
        for cid in ids:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            ff.Minimize(maxIts=2000)
            energies[cid] = ff.CalcEnergy()
    best = min(energies, key=energies.get)
    conf = mol.GetConformer(best)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    coordinates = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    if return_mol:
        keep = Chem.Mol(mol)
        keep.RemoveAllConformers()
        keep.AddConformer(Chem.Conformer(conf), assignId=True)
        return keep, symbols, coordinates, [energies[c] for c in ids], ids.index(best)
    return symbols, coordinates, [energies[c] for c in ids], ids.index(best)
