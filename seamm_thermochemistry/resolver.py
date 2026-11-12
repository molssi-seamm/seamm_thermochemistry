# -*- coding: utf-8 -*-

"""Resolving the reference species of CBH reactions: ATcT first, then the
supplementary reference library.

The order (the reference-library design, §2, and the 2026-11-12 decisions):

1. ATcT, the version installed. If every reference is in ATcT, nothing else
   is used and the result is "ATcT-only".
2. The user's own references (the local layer's ``user`` source), when the
   job enables them.
3. The library's sources by rank: printed JANAF, derived values, CODATA, D0,
   W4-17 (0 K only, used only without a 298.15 K value), the WebBook
   measurements.
4. Computed references, when enabled.

**One source per element.** For each element that ATcT does not have (its
gas-phase atom is not in ATcT), the references containing it are taken from
the first source, in rank order, that has all of them, so that the source's
shared error cancels in the reaction. A later source still replaces that
source's value for a species where the replacement rule says so (science,
2026-11-12): the element then has several sources, the result is labelled
"mixed: several sources" with the sources named, and a one-line note says so.
It never blocks. When no single source has them all, each reference takes its
own best value, with the same label.

The values come back in the layout of an ATcT entry (``dfH0_0K``,
``dfH0_298K``, ``uncertainty_298K`` at 95 %, ``id``), with
``atct_coefficients`` empty (they have no covariance with ATcT) and the
library's details (``source``, ``layer``, ``citation``, ...).
"""

import logging
import re

from .atct import AmbiguousEntry, MissingEntry
from .cbh import CBHError
from .references import LayeredReferences

logger = logging.getLogger(__name__)

#: The settings of "Supplementary references"
SUPPLEMENTARY = ("none", "compilations", "compilations and computed")

#: The result labels, in increasing order of importance
LABELS = (
    "ATcT-only",
    "mixed",
    "mixed: several sources",
    "computed references",
    "user references",
)

#: H298 - H0 of the elements' standard states, kJ/mol per formula unit of the
#: standard state, from CODATA Key Values (Cox, Wagman & Medvedev 1989)
CODATA_H298_H0 = {
    "H": (8.468, 2),  # H2(g)
    "C": (1.050, 1),  # graphite
    "N": (8.670, 2),  # N2(g)
    "O": (8.680, 2),  # O2(g)
    "F": (8.825, 2),  # F2(g)
    "Cl": (9.181, 2),  # Cl2(g)
    "Br": (24.52, 2),  # Br2(l)
    "S": (4.412, 1),  # rhombic
    "P": (5.360, 1),  # white
    "Si": (3.217, 1),
    "B": (1.222, 1),
    "Al": (4.540, 1),
    "Li": (4.632, 1),
    "Na": (6.460, 1),
    "K": (7.088, 1),
    "Be": (1.950, 1),
    "Mg": (4.998, 1),
    "Ca": (5.736, 1),
    "Ge": (4.636, 1),
}


def element_h298_h0():
    """{symbol: H298 - H0 of the element's standard state, per atom}, kJ/mol
    (CODATA; CODATA has no Ga, As or Se)."""
    return {k: v / n for k, (v, n) in CODATA_H298_H0.items()}


def formula_counts(formula):
    counts = {}
    for symbol, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        counts[symbol] = counts.get(symbol, 0) + (int(n) if n else 1)
    return counts


def dfh_0k_from_298k(dfh_298k, h298_h0_species, formula):
    """ΔfH(0 K) from ΔfH(298.15 K), the species' own H298 - H0 and the
    elements' (standard states, per atom)::

        ΔfH(0) = ΔfH(298) - [H298 - H0](species) + Σ n_E [H298 - H0](E)
    """
    table = element_h298_h0()
    total = 0.0
    for symbol, n in formula_counts(formula).items():
        if symbol not in table:
            raise KeyError(f"No H298 - H0 for the element {symbol}")
        total += n * table[symbol]
    return dfh_298k - h298_h0_species + total


class ReferenceResolver:
    """Resolve CBH reference species from ATcT and the reference library.

    Parameters
    ----------
    atct : ATcT
        The ATcT tables.
    supplementary : str
        "none" (ATcT only, today's behaviour), "compilations", or
        "compilations and computed" (the default, Paul's decision 4).
    use_user : bool
        Use the user's own references (the local layer's ``user`` kind).
    layers : LayeredReferences, optional
        The libraries; by default the installed ones.
    multiplicity : callable, optional
        multiplicity(smiles) for a reference species (the CBH convention).
    """

    def __init__(
        self,
        atct,
        *,
        supplementary="compilations and computed",
        use_user=False,
        layers=None,
        multiplicity=None,
    ):
        if supplementary not in SUPPLEMENTARY:
            raise ValueError(
                f"Supplementary references must be one of {SUPPLEMENTARY}, not "
                f"{supplementary!r}"
            )
        self.atct = atct
        self.supplementary = supplementary
        self.use_user = use_user
        self.multiplicity = multiplicity
        if supplementary == "none" and not use_user:
            self.layers = None
        else:
            self.layers = layers if layers is not None else LayeredReferences()
        self._atct_elements = {}
        self.notes = []
        if self.layers is not None:
            self.layers.set_in_atct(self._in_atct)

    # -- helpers --------------------------------------------------------------

    def _charge_mult(self, smiles):
        from .cbh import _composition

        _, charge, _, _ = _composition(smiles)
        mult = self.multiplicity(smiles) if self.multiplicity else None
        return charge, mult

    def atct_lookup(self, smiles):
        charge, mult = self._charge_mult(smiles)
        return self.atct.lookup(smiles, charge=charge, multiplicity=mult)

    def _in_atct(self, smiles, charge=0, multiplicity=None):
        try:
            self.atct.lookup(smiles, charge=charge, multiplicity=multiplicity)
            return True
        except (MissingEntry, AmbiguousEntry):
            return False

    def atct_has_element(self, symbol):
        """Whether ATcT has the element's gas-phase atom."""
        if symbol not in self._atct_elements:
            try:
                self.atct.lookup(f"[{symbol}]", charge=0)
                self._atct_elements[symbol] = True
            except (MissingEntry, AmbiguousEntry):
                self._atct_elements[symbol] = False
            except Exception:  # pragma: no cover
                self._atct_elements[symbol] = False
        return self._atct_elements[symbol]

    def _allowed(self, ref):
        kind = ref["source_kind"]
        if kind == "user":
            return self.use_user
        if self.supplementary == "none":
            return False
        if kind == "computed":
            return self.supplementary == "compilations and computed"
        return True

    def candidates(self, smiles):
        """The library's usable values for a species, best rank first."""
        if self.layers is None:
            return []
        charge, mult = self._charge_mult(smiles)
        return [
            r for r in self.layers.candidates(smiles, charge, mult) if self._allowed(r)
        ]

    @staticmethod
    def _best(candidates):
        """The preferred value among candidates (overrides applied within a
        layer), as LayeredReferences.lookup does."""
        allowed = {r["reference_id"]: r for r in candidates}
        for r in candidates:
            if r["preferred"]:
                return r
        # The preferred value may be filtered out (e.g. a computed one):
        # the best remaining by rank, at 298.15 K first
        at_298 = [r for r in candidates if r["dfh_298k"] is not None]
        pool = at_298 or list(allowed.values())
        return pool[0] if pool else None

    @staticmethod
    def _follow_override(candidates, start):
        """The value that replaces `start` through the chain of overrides."""
        current = start
        while True:
            later = [
                r
                for r in candidates
                if r["overrides"] == current["reference_id"]
                and r["layer"] == current["layer"]
            ]
            if not later:
                return current
            current = later[0]

    def _entry(self, smiles, ref):
        u = ref["total_uncertainty"]
        return {
            "id": f"{ref['source_name']}:{ref['layer']}:{ref['id']}",
            "smiles": smiles,
            "dfH0_0K": ref["dfh_0k"],
            "dfH0_298K": ref["dfh_298k"],
            "uncertainty_298K": u,
            "exact": False,
            "atct_coefficients": {},
            "source": ref["source_name"],
            "source_kind": ref["source_kind"],
            "layer": ref["layer"],
            "citation": ref["citation"],
            "rule": ref["rule"],
            "flag": ref["flag"],
            "override_reason": ref["override_reason"],
            "decision_reason": ref["decision_reason"],
            "unc_level": ref["unc_level"],
            "offset": ref["offset_kj"],
            "offset_from": ref.get("offset_from"),
            "formula": ref.get("formula"),
        }

    # -- resolving ----------------------------------------------------------------

    def resolve(self, species):
        """Resolve a set of reference species (canonical SMILES).

        Returns (entries, missing, label, notes): entries {smiles: entry},
        ATcT entries carrying ``atct_coefficients`` {id: 1} and ``source``
        "ATcT"; missing {smiles: reason}; the result label (LABELS); and
        one-line notes about mixed sources, flags and conversions.
        """
        entries = {}
        missing = {}
        notes = []
        library_species = []
        for smiles in species:
            try:
                entry = dict(self.atct_lookup(smiles))
                entry["atct_coefficients"] = {entry["id"]: 1.0}
                entry["source"] = "ATcT"
                entries[smiles] = entry
            except (MissingEntry, AmbiguousEntry) as e:
                library_species.append(smiles)
                missing[smiles] = str(e)
            except CBHError as e:
                missing[smiles] = str(e)

        # The library: one source per element ATcT lacks
        cands = {s: self.candidates(s) for s in library_species}
        assigned = {}
        uncovered = {}
        for smiles in library_species:
            if not cands[smiles]:
                continue
            for symbol in self._elements(smiles):
                if not self.atct_has_element(symbol):
                    uncovered.setdefault(symbol, []).append(smiles)
        several = {}
        for symbol in sorted(uncovered, key=lambda e: len(uncovered[e])):
            group = [s for s in uncovered[symbol] if s not in assigned]
            if not group:
                continue
            sources = None
            for s in group:
                names = {r["source_name"] for r in cands[s]}
                sources = names if sources is None else sources & names
            ranked = sorted(
                {
                    (r["source_rank"], r["source_name"])
                    for s in group
                    for r in cands[s]
                    if r["source_name"] in (sources or ())
                }
            )
            if ranked:
                source = ranked[0][1]
                for s in group:
                    own = [r for r in cands[s] if r["source_name"] == source]
                    # local layer first within the source
                    start = own[0]
                    chosen = self._follow_override(cands[s], start)
                    if chosen["dfh_298k"] is None and any(
                        r["dfh_298k"] is not None for r in cands[s]
                    ):
                        # a 0 K-only value is used only without a 298 K one
                        chosen = self._best(cands[s])
                    assigned[s] = chosen
            else:
                for s in group:
                    assigned[s] = self._best(cands[s])
        for smiles in library_species:
            if smiles not in assigned and cands[smiles]:
                assigned[smiles] = self._best(cands[smiles])

        for smiles, ref in assigned.items():
            if ref is None:
                continue
            entries[smiles] = self._entry(smiles, ref)
            missing.pop(smiles, None)
            if ref.get("flag"):
                notes.append(f"{smiles} ({ref['source_name']}): {ref['flag']}")

        # Several sources for one element?
        for symbol, group in uncovered.items():
            names = sorted({entries[s]["source"] for s in group if s in entries})
            if len(names) > 1:
                several[symbol] = names
                notes.append(
                    f"The references with {symbol} come from several sources: "
                    + ", ".join(names)
                    + " (one source per element is the aim; the values are "
                    "used as they are)."
                )

        label = self.label(entries, several)
        self.notes.extend(notes)
        return entries, missing, label, notes

    def _elements(self, smiles):
        from .references import elements

        return elements(smiles)

    @staticmethod
    def label(entries, several=None):
        """The result label for the references used."""
        kinds = {
            e.get("source_kind") for e in entries.values() if e["source"] != "ATcT"
        }
        if "user" in kinds:
            return "user references"
        if "computed" in kinds:
            return "computed references"
        if several:
            return (
                "mixed: several sources ("
                + "; ".join(f"{k}: {', '.join(v)}" for k, v in sorted(several.items()))
                + ")"
            )
        if kinds:
            return "mixed"
        return "ATcT-only"

    def versions(self):
        """The ATcT version and the libraries' versions and paths, for the
        result's ``references`` property."""
        result = {"ATcT": self.atct.version}
        if self.layers is not None:
            result["libraries"] = self.layers.versions
        result["supplementary"] = self.supplementary
        result["use user references"] = self.use_user
        return result
