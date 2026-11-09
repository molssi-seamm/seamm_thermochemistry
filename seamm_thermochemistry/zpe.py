# -*- coding: utf-8 -*-

"""Scale factors for harmonic zero-point energies, with their sources.

A harmonic ZPE is too large; scaled by a factor fitted for the level of theory
it comes close to the anharmonic one. The factor belongs to the geometry and
frequency level (the method and basis of the Hessian), so it is looked up from
that level's model chemistry, e.g. ``ORCA:CMP@R2SCAN-3C`` or
``ORCA:DFT@PBE0-D4/def2-TZVP``. Each entry carries its source, and a flag when
it is used outside what it was fitted for (another basis, or larger molecules
than the fit set); the thermochemistry step prints both and stores them with its
results.

The store keeps ZPEs unscaled; the factor is applied when the values are used,
so a corrected factor needs no recomputation.
"""

from dataclasses import dataclass
import re

__all__ = ["ZPEFactor", "ZPE_FACTORS", "zpe_scale_factor", "split_level"]


@dataclass(frozen=True)
class ZPEFactor:
    """A ZPE scale factor for a level of theory: the factor, its source, a
    caution about its use (e.g. a basis set other than the fitted one) or "",
    and the method and basis it is for (upper case; basis "" for a composite
    method)."""

    factor: float
    source: str
    flag: str = ""
    method: str = ""
    basis: str = ""


#: The factors, by (method, basis), both upper case with the dispersion
#: correction left out of the method (it hardly changes the frequencies).
ZPE_FACTORS = {
    ("R2SCAN-3C", ""): ZPEFactor(
        0.9899,
        "fitted to the anharmonic ZPVEs of 22 W4-2006 molecules (SEAMM, 2026-11)",
        "fitted on molecules of up to 6 atoms; an extrapolation for larger ones",
        "R2SCAN-3C",
        "",
    ),
    ("PBE0", "DEF2-TZVP"): ZPEFactor(
        0.975,
        "Truhlar's database of frequency scale factors, version 5, PBE0/MG3S",
        "Truhlar PBE0/MG3S value, basis-mismatched (used for def2-TZVP)",
        "PBE0",
        "DEF2-TZVP",
    ),
    ("B2PLYP", "CC-PVTZ"): ZPEFactor(
        0.9863,
        "the CBH-ANL protocol's factor for B2PLYP-D3/cc-pVTZ harmonic ZPEs",
        "",
        "B2PLYP",
        "CC-PVTZ",
    ),
}

_DISPERSION = re.compile(r"-(D2|D3|D3BJ|D3ZERO|D3\(BJ\)|D3\(0\)|D4)$", re.I)


def split_level(model_chemistry):
    """(method, basis) of a model chemistry string, upper case, with the
    program and type prefix, a ``bse:`` source and the dispersion suffix of the
    method removed; basis "" for a composite method.

    >>> split_level("ORCA:DFT@PBE0-D4/bse:def2-TZVP")
    ('PBE0', 'DEF2-TZVP')
    >>> split_level("ORCA:CMP@R2SCAN-3C")
    ('R2SCAN-3C', '')
    """
    text = str(model_chemistry or "").strip()
    if "@" in text:
        text = text.split("@", 1)[1]
    method, _, basis = text.partition("/")
    if basis.lower().startswith("bse:"):
        basis = basis[4:]
    method = method.strip().upper()
    if not method.endswith("-3C"):
        method = _DISPERSION.sub("", method)
    return method, basis.strip().upper()


def zpe_scale_factor(model_chemistry):
    """The ZPE scale factor for a geometry/frequency level, or None if the
    table has none for it."""
    return ZPE_FACTORS.get(split_level(model_chemistry))
