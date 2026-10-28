# -*- coding: utf-8 -*-

"""Laddered reference values for CBH reactions (CBH-ANL, Elliott et al.
JPCA 127, 1512 (2023), eqs. 8-11).

Only the CBH-0 species (H2, CH4, H2O, NH3, HF, ..., and the radical CBH-0
species such as CH3 and OH) take their values from ATcT. Every other
reference species gets its value from its *own* CBH reaction one rung below
its size, computed at a high level of theory, anchored to the laddered values
of that reaction's species::

    L(M) = H(M) - Σ_i c_i H(i) + Σ_i c_i L(i)

where H are the high-level enthalpies (E + scaled ZPE at 0 K, or with the
thermal correction at 298.15 K), c_i the coefficients of M's CBH-(r-1)
reaction (r = M's own rung, the lowest rung at which M is its own fragment)
and L(i) the laddered values of the reaction's other species, recursively
down to the ATcT anchors. Because every laddered value is a linear
combination of ATcT anchors plus high-level reaction enthalpies, each one
carries its coefficients on the ATcT species, so a reaction's ATcT
uncertainty is cᵀΣc with the ATcT covariance exactly, and a *method* term
for the high-level reaction enthalpies is measured by comparing the laddered
values with ATcT for every ladder species that ATcT has.

A ladder is built from the fragment-energy store (the high-level energies and
the frequency-level ZPE and thermal corrections, harvested from finished
jobs) and saved as a named reference set, a JSON file in ``ladders/`` next to
the store, with the provenance of every value.
"""

import datetime
import json
import logging
import math
from pathlib import Path

from .atct import ATcT, AmbiguousEntry, MissingEntry, reaction_uncertainty
from .fragment_store import FragmentStore, fragment_store_path

logger = logging.getLogger(__name__)

__all__ = ["Ladder", "build_ladder", "ladder_path", "LadderError"]

T_REF = 298.15
DEFAULT_ZPE_SCALE = 0.9863  # CBH-ANL, B2PLYP-D3/cc-pVTZ harmonic ZPE


class LadderError(RuntimeError):
    pass


def ladder_path(name, store_path=None):
    """The file of the named ladder: ``ladders/<name>.json`` next to the
    fragment store."""
    if name.endswith(".json") and Path(name).expanduser().exists():
        return Path(name).expanduser()
    safe = "".join(c if c.isalnum() or c in "-_.+@" else "_" for c in name)
    return fragment_store_path(store_path).parent / "ladders" / f"{safe}.json"


def _lookup(tables, smiles):
    """ATcT entry of a species with the charge and multiplicity its SMILES
    implies, or None."""
    from . import cbh

    _, q, radicals, _ = cbh._composition(smiles)
    try:
        return tables.lookup(smiles, charge=q, multiplicity=radicals + 1)
    except (MissingEntry, AmbiguousEntry):
        return None


def _own_rung(smiles, protect, max_rung=6):
    """The rung at which the species is its own fragment (its CBH-n reaction
    is trivial), and the reaction of the highest available rung below it.

    Returns (rung, reaction) with reaction None for a rung-0 species.
    """
    from . import cbh

    below = None
    for n in range(0, max_rung + 1):
        try:
            reaction = cbh.cbh_reactions(smiles, rungs=[n], protect=protect)[n]
        except cbh.RungUnavailable:
            continue
        except cbh.CBHError as e:
            if "no heavy atoms" in str(e):  # H2: a CBH-0 species
                return 0, None
            raise
        if reaction.trivial:
            return n, below
        below = reaction
    raise LadderError(f"{smiles} is not its own fragment up to CBH-{max_rung}")


class Ladder:
    """A named set of laddered reference values.

    ``ladder[species]`` is a dict with the keys ``rung``, ``H0`` and ``H298``
    (kJ/mol, the laddered ΔfH° at 0 K and 298.15 K; one may be None),
    ``atct_coefficients`` ({ATcT ID: coefficient}), ``uncertainty`` (the 95 %
    ATcT part), ``used_covariance``, ``reaction``, ``provenance`` and, for
    species that ATcT has, ``atct_difference_0K`` / ``atct_difference_298K``.
    """

    def __init__(self, data):
        self.data = data
        self.name = data["name"]
        self.species = data["species"]
        self.atct_version = data["atct_version"]
        self.method_uncertainty = data.get("method_uncertainty", {})

    @classmethod
    def load(cls, name, store_path=None):
        path = ladder_path(name, store_path)
        if not path.exists():
            raise FileNotFoundError(f"No ladder '{name}' ({path})")
        with open(path, encoding="utf-8") as fd:
            return cls(json.load(fd))

    def save(self, store_path=None):
        path = ladder_path(self.name, store_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fd:
            json.dump(self.data, fd, indent=1, ensure_ascii=False)
        return path

    def __contains__(self, species):
        return species in self.species

    def __getitem__(self, species):
        return self.species[species]

    def __len__(self):
        return len(self.species)

    def entry(self, species):
        """The laddered value in the layout of an ATcT entry (``dfH0_0K``,
        ``dfH0_298K``, ``uncertainty_298K``, ``id``), with the extra keys
        ``atct_coefficients``, ``method_uncertainty`` and ``source``; or None.
        """
        if species not in self.species:
            return None
        v = self.species[species]
        return {
            "id": f"ladder:{self.name}:{species}",
            "smiles": species,
            "dfH0_0K": v["H0"],
            "dfH0_298K": v["H298"],
            "uncertainty_298K": v["uncertainty"],
            "exact": False,
            "units": "kJ/mol",
            "atct_coefficients": v["atct_coefficients"],
            "used_covariance": v["used_covariance"],
            "method_uncertainty": self.method_uncertainty,
            "source": "laddered" if v["rung"] > 0 else "ATcT",
            "rung": v["rung"],
        }

    def summary(self):
        lines = [
            f"Ladder '{self.name}': {len(self.species)} species on ATcT ver. "
            f"{self.atct_version}; energies {self.data['energy_level']}, "
            f"frequencies {self.data['frequency_level']} (ZPE x "
            f"{self.data['zpe_scale']}); built {self.data['created']}."
        ]
        mu = self.method_uncertainty
        if mu.get("n"):
            lines.append(
                f"Method term from {mu['n']} ladder species in ATcT: RMS "
                f"{mu['rms_0K']:.2f} (0 K), max |diff| {mu['max_0K']:.2f} kJ/mol"
                + (
                    f"; 298.15 K RMS {mu['rms_298K']:.2f}"
                    if mu.get("rms_298K") is not None
                    else ""
                )
            )
        for s, v in sorted(self.species.items(), key=lambda x: (x[1]["rung"], x[0])):
            h0 = f"{v['H0']:9.2f}" if v["H0"] is not None else "        -"
            h298 = f"{v['H298']:9.2f}" if v["H298"] is not None else "        -"
            diff = (
                f"  (ATcT {v['atct_difference_0K']:+.2f})"
                if v.get("atct_difference_0K") is not None
                else ""
            )
            lines.append(
                f"  CBH-{v['rung']} {s:16s} {h0} {h298} ± {v['uncertainty']:.2f}"
                f"{diff}  {v['reaction'] or 'ATcT anchor'}"
            )
        return "\n".join(lines)


def build_ladder(
    species,
    *,
    name,
    energy_level,
    frequency_level=None,
    zpe_scale=DEFAULT_ZPE_SCALE,
    atct_version=None,
    atct_base=None,
    store_path=None,
    protect="default",
    fetch_covariance=False,
    anchors_from_atct=True,
):
    """Build a ladder for `species` (SMILES; their lower-rung fragments are
    added as needed) from the fragment store.

    Parameters
    ----------
    energy_level : str
        The high-level energies: a settings fingerprint, or a model chemistry
        string (then unverified records are accepted, and reported).
    frequency_level : str, optional
        The level of the ZPE and thermal corrections (default: the energy
        level's own record).
    zpe_scale : float
        Scale factor for the harmonic ZPE (0.9863 for B2PLYP-D3/cc-pVTZ).
    """
    from . import cbh
    from .atct import canonical_smiles

    tables = ATcT(version=atct_version, base=atct_base)
    store = FragmentStore(store_path)
    if len(store) == 0:
        raise LadderError(f"The fragment store {store.path} is empty or missing")

    def record(smiles, level):
        _, q, radicals, _ = cbh._composition(smiles)
        return store.get(
            smiles,
            charge=q,
            multiplicity=radicals + 1,
            energy_fingerprint=level,
            model_chemistry=level,
        )

    def enthalpies(smiles):
        """(H0, H298, provenance) at the high level, or an error message."""
        e = record(smiles, energy_level)
        if e is None or e.get("E") is None:
            return None, None, f"no energy at '{energy_level}' in the store"
        f = e if frequency_level is None else record(smiles, frequency_level)
        if f is None:
            return None, None, f"no frequencies at '{frequency_level}' in the store"
        prov = {
            "energy": {k: e.get(k) for k in ("job", "configuration", "verified")},
        }
        if f is not e:
            prov["frequencies"] = {
                k: f.get(k) for k in ("job", "configuration", "verified")
            }
        H0 = H298 = None
        if f.get("ZPE") is not None:
            H0 = e["E"] + zpe_scale * f["ZPE"]
            if f.get("H") is not None and f.get("E") is not None:
                thermal = f["H"] - f["E"] - f["ZPE"]
                H298 = H0 + thermal
        elif f.get("H") is not None and f.get("E") is not None:
            H298 = e["E"] + (f["H"] - f["E"])
        if H0 is None and H298 is None:
            return None, None, "the store record has neither ZPE nor H(298.15 K)"
        return H0, H298, prov

    values = {}  # smiles -> dict
    missing = {}
    building = set()

    def ladder(smiles):
        if smiles in values:
            return values[smiles]
        if smiles in missing:
            return None
        if smiles in building:
            raise LadderError(f"Circular ladder at {smiles}")
        building.add(smiles)
        try:
            rung, reaction = _own_rung(smiles, protect)
            entry = _lookup(tables, smiles)
            if rung == 0 or reaction is None:
                if entry is None or not anchors_from_atct:
                    missing[smiles] = "a CBH-0 species that ATcT does not have"
                    return None
                v = {
                    "rung": 0,
                    "H0": entry["dfH0_0K"],
                    "H298": entry["dfH0_298K"],
                    "atct_coefficients": {entry["id"]: 1.0},
                    "atct_id": entry["id"],
                    "uncertainty": entry["uncertainty_298K"] or 0.0,
                    "used_covariance": True,
                    "reaction": None,
                    "provenance": {"ATcT": entry["id"]},
                }
                values[smiles] = v
                return v
            # The species' own high-level enthalpies and those of the reaction
            H0, H298, prov = enthalpies(smiles)
            if H0 is None and H298 is None:
                missing[smiles] = prov
                return None
            refs = {}
            for ref in reaction.references:
                lv = ladder(ref)
                if lv is None:
                    missing[smiles] = f"its CBH-{rung - 1} reaction needs {ref}: " + (
                        missing.get(ref, "")
                    )
                    return None
                refs[ref] = lv
            ref_H = {}
            for ref in reaction.references:
                h0, h298, p = enthalpies(ref)
                if h0 is None and h298 is None:
                    missing[smiles] = f"the reaction species {ref} has {p}"
                    return None
                ref_H[ref] = (h0, h298)
            # L(M) = H(M) - Σ c H(i) + Σ c L(i); coefficients on ATcT species
            coefficients = {}
            L0 = H0
            L298 = H298
            for ref, c in reaction.references.items():
                h0, h298 = ref_H[ref]
                if L0 is not None:
                    L0 = (
                        L0 - c * h0 + c * refs[ref]["H0"]
                        if h0 is not None and refs[ref]["H0"] is not None
                        else None
                    )
                if L298 is not None:
                    L298 = (
                        L298 - c * h298 + c * refs[ref]["H298"]
                        if h298 is not None and refs[ref]["H298"] is not None
                        else None
                    )
                for _id, cc in refs[ref]["atct_coefficients"].items():
                    coefficients[_id] = coefficients.get(_id, 0.0) + c * cc
            coefficients = {k: v for k, v in coefficients.items() if abs(v) > 1e-12}
            ids = list(coefficients)
            unc = [tables.entries[_id]["uncertainty_298K"] or 0.0 for _id in ids]
            u, used = reaction_uncertainty(
                list(coefficients.values()),
                unc,
                tables.covariance(ids, fetch=fetch_covariance) if ids else None,
            )
            v = {
                "rung": rung,
                "H0": L0,
                "H298": L298,
                "atct_coefficients": coefficients,
                "atct_id": entry["id"] if entry else None,
                "uncertainty": u,
                "used_covariance": used,
                "reaction": str(reaction),
                "provenance": prov,
                "high_level": {"H0": H0, "H298": H298},
            }
            if entry is not None:
                if L0 is not None and entry["dfH0_0K"] is not None:
                    v["atct_difference_0K"] = L0 - entry["dfH0_0K"]
                if L298 is not None:
                    v["atct_difference_298K"] = L298 - entry["dfH0_298K"]
            values[smiles] = v
            return v
        finally:
            building.discard(smiles)

    for s in species:
        ladder(canonical_smiles(s))

    # The method term: laddered vs ATcT over the ladder species ATcT has
    d0 = [v["atct_difference_0K"] for v in values.values() if "atct_difference_0K" in v]
    d298 = [
        v["atct_difference_298K"]
        for v in values.values()
        if "atct_difference_298K" in v
    ]
    method = {"n": len(d0)}
    if d0:
        method["rms_0K"] = math.sqrt(sum(d * d for d in d0) / len(d0))
        method["max_0K"] = max(abs(d) for d in d0)
        method["mean_0K"] = sum(d0) / len(d0)
    if d298:
        method["rms_298K"] = math.sqrt(sum(d * d for d in d298) / len(d298))
        method["max_298K"] = max(abs(d) for d in d298)

    unverified = sorted(
        s
        for s, v in values.items()
        if v["rung"] > 0 and not v["provenance"].get("energy", {}).get("verified", True)
    )
    data = {
        "name": name,
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "atct_version": tables.version,
        "atct_citation": tables.citation,
        "energy_level": energy_level,
        "frequency_level": frequency_level or energy_level,
        "zpe_scale": zpe_scale,
        "protect": protect,
        "store": str(store.path),
        "species": values,
        "missing": missing,
        "unverified": unverified,
        "method_uncertainty": method,
    }
    return Ladder(data)
