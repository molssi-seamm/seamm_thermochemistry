# -*- coding: utf-8 -*-

"""Corrections of harmonic frequencies, one frequency at a time, with their sources.

Some methods get their frequencies wrong in a way a single scale factor cannot fix.
PM7, for one, puts the X-H stretches near a quarter of their real spread, so its
frequencies above 2500 cm^-1 need a steep linear map (Rozanska et al. 2014, eq 1).
Such a correction is a piecewise-linear function of the frequency, fitted for a
method (and basis), and is applied to each real frequency before the zero-point
energy and the thermal functions (H(T) - H(0), S, Cp) are computed.

A correction replaces the ZPE scale factor for its method: the two are never
applied together. Each entry carries its citation, the elements its fit set
contained, and the frequency window where the map is discontinuous or was not
validated, so the thermochemistry step can warn about frequencies there. Its
``tag`` (id, version and the range used) goes into the frequency fingerprint of the
corrected thermal data, so corrected and uncorrected data are never mixed.

The functions are looked up like the ZPE scale factors, by the (method, basis) of
the frequency level's model chemistry (:func:`seamm_thermochemistry.zpe.split_level`).
"""

from dataclasses import dataclass, field

from .zpe import split_level

__all__ = [
    "FrequencyCorrection",
    "FREQUENCY_CORRECTIONS",
    "RANGES",
    "frequency_correction",
]

#: How much of a correction to apply: all of it, or only the branches that are
#: needed (a branch marked ``optional`` is then left out)
RANGES = ("full", "high-frequency only")


@dataclass(frozen=True)
class Branch:
    """``nu' = slope * nu + intercept`` for ``low <= nu < high`` (cm^-1).
    ``optional``: the source says the branch may be left out."""

    low: float
    high: float
    slope: float = 1.0
    intercept: float = 0.0
    optional: bool = False

    def __call__(self, nu):
        return self.slope * nu + self.intercept

    @property
    def identity(self):
        return self.slope == 1.0 and self.intercept == 0.0


@dataclass(frozen=True)
class FrequencyCorrection:
    """A piecewise-linear correction of the harmonic frequencies of one method.

    Its fields:

    id, version : str, int
        Its name and version; both go into :attr:`tag`.
    method, basis : str
        The level it is for, upper case, as :func:`split_level` gives it.
    branches : tuple of Branch
        Covering 0 to infinity in order; frequencies on a boundary use the
        higher branch.
    citation : str
        The source, for printing.
    bibtex : str
        The source as a BibTeX entry, for the references.
    fitted_elements : tuple of str
        The elements of the fit set. Stretches of hydrogen bonded to other
        elements were not validated.
    unvalidated : tuple of (float, float)
        The windows (cm^-1) where the map is discontinuous or was not
        validated, as raw frequencies.
    valid_max : float or None
        The largest raw frequency (cm^-1) the fit covers. Higher frequencies are
        left as they are (a steep branch would blow them up), and the step warns.
    note : str
        Anything else to say when it is used.
    """

    id: str
    version: int
    method: str
    basis: str
    branches: tuple
    citation: str
    bibtex: str = ""
    fitted_elements: tuple = ()
    unvalidated: tuple = ()
    valid_max: float = None
    note: str = ""
    aliases: tuple = field(default=())

    def tag(self, range_="full"):
        """The text that identifies the corrected data, for the fingerprints."""
        _check_range(range_)
        return f"fc:{self.id}:v{self.version}:{range_}"

    def branch(self, nu):
        """The branch for a frequency (cm^-1)."""
        for branch in self.branches:
            if branch.low <= nu < branch.high:
                return branch
        return self.branches[-1]

    def __call__(self, frequencies, range_="full"):
        """The corrected frequencies (cm^-1), in the same order.

        ``range_`` "high-frequency only" leaves the optional branches out.
        Frequencies that are not positive (imaginary modes given as negative
        numbers, zero modes), and those above :attr:`valid_max`, are returned
        unchanged.
        """
        _check_range(range_)
        out = []
        for nu in frequencies:
            nu = float(nu)
            if nu <= 0.0 or (self.valid_max is not None and nu > self.valid_max):
                out.append(nu)
                continue
            branch = self.branch(nu)
            if branch.optional and range_ != "full":
                out.append(nu)
            else:
                out.append(branch(nu))
        return out

    def not_applicable(self, frequencies):
        """The raw frequencies (cm^-1) above the fit's range, left uncorrected."""
        if self.valid_max is None:
            return []
        return [float(nu) for nu in frequencies if float(nu) > self.valid_max]

    def in_unvalidated_windows(self, frequencies):
        """The positive frequencies (cm^-1) in an unvalidated window."""
        return [
            float(nu)
            for nu in frequencies
            if any(low <= float(nu) <= high for low, high in self.unvalidated)
        ]

    def describe(self, range_="full"):
        """The function, as text."""
        _check_range(range_)
        lines = []
        for b in self.branches:
            if b.high == float("inf"):
                where = f"{b.low:g} cm^-1 and above"
            elif b.low == 0.0:
                where = f"below {b.high:g} cm^-1"
            else:
                where = f"{b.low:g} to {b.high:g} cm^-1"
            if b.identity or (b.optional and range_ != "full"):
                what = "unchanged"
                if b.optional and range_ != "full" and not b.identity:
                    what += " (the optional branch is not applied)"
            else:
                sign = "-" if b.intercept < 0 else "+"
                what = f"{b.slope:g} nu {sign} {abs(b.intercept):g}"
            lines.append(f"{where}: {what}")
        if self.valid_max is not None:
            lines.append(
                f"above {self.valid_max:g} cm^-1 (outside the fit): unchanged, with a "
                "warning"
            )
        return lines


def _check_range(range_):
    if range_ not in RANGES:
        raise ValueError(
            f"The range of a frequency correction must be one of {', '.join(RANGES)}, "
            f"not {range_!r}"
        )


_ROZANSKA_2014 = """@article{Rozanska2014,
    author = {Rozanska, Xavier and Stewart, James J. P. and Ungerer, Philippe and
              Leblanc, Benoit and Freeman, Chris and Saxe, Paul and Wimmer, Erich},
    title = {High-Throughput Calculations of Molecular Properties in the MedeA
             Environment: Accuracy of PM7 in Predicting Vibrational Frequencies,
             Ideal Gas Entropies, Heat Capacities, and Gibbs Free Energies of
             Organic Molecules},
    journal = {Journal of Chemical \\& Engineering Data},
    volume = {59},
    number = {10},
    pages = {3136--3143},
    year = {2014},
    doi = {10.1021/je500201y}
}"""

#: The corrections, by (method, basis), upper case
FREQUENCY_CORRECTIONS = {
    ("PM7", ""): FrequencyCorrection(
        id="PM7-Rozanska2014",
        version=1,
        method="PM7",
        basis="",
        branches=(
            Branch(0.0, 125.0),
            Branch(125.0, 2500.0, 1.006, -27.7, optional=True),
            Branch(2500.0, float("inf"), 3.99, -7849.2),
        ),
        citation=(
            "Rozanska, Stewart, Ungerer, Leblanc, Freeman, Saxe and Wimmer, "
            "J. Chem. Eng. Data 59, 3136 (2014), eq 1"
        ),
        bibtex=_ROZANSKA_2014,
        fitted_elements=("C", "H", "N", "O"),
        # Raw frequencies around both switches, where a small change moves a
        # mode to the other branch: 124 -> 124 but 126 -> 99.1, and 2499 ->
        # 2486.3 but 2501 -> 2129.8. (PM7's C-H/N-H stretches lie at 2620-2850,
        # above this window.)
        unvalidated=((98.0, 125.0), (2400.0, 2600.0)),
        # The fit's experimental frequencies end at about 3700 cm^-1 (the paper's
        # Figure 2; raw PM7 there about 2850), i.e. a raw (3700 + 7849.2) / 3.99 =
        # 2894.5. The slope-3.99 branch must not be extrapolated beyond it: PM7's
        # H2 stretch at 5577 cm^-1 would become 14402 (experiment ~4401).
        valid_max=2900.0,
        note=(
            "Fitted to the experimental frequencies of 52 organic molecules (C, H, "
            "N, O; experimental frequencies up to about 3700 cm^-1, raw PM7 up to "
            "about 2850). The map is discontinuous at 125 cm^-1 (124 -> 124, "
            "126 -> 99.1) and at 2500 cm^-1 (2499 -> 2486.3, 2501 -> 2129.8). The "
            "source calls the 125-2500 cm^-1 branch optional (0.6 %). For PM7 only: "
            "other MOPAC "
            "methods (PM6, PM6-ORG, ...) need their own fit."
        ),
    ),
}


def frequency_correction(model_chemistry):
    """The frequency correction for a frequency level, or None if there is none.

    Only the exact (method, basis) of an entry matches: a correction fitted for
    PM7 is not used for PM6, PM6-ORG or any other method.
    """
    return FREQUENCY_CORRECTIONS.get(split_level(model_chemistry))
