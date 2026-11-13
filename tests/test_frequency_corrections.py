# -*- coding: utf-8 -*-

"""Per-frequency corrections (Rozanska et al. 2014, eq 1, for PM7)."""

import pytest

from seamm_thermochemistry import frequency_correction
from seamm_thermochemistry.frequency_corrections import RANGES


@pytest.fixture
def pm7():
    return frequency_correction("MOPAC:SQM@PM7")


def test_lookup_is_pm7_only(pm7):
    assert pm7 is not None and pm7.id == "PM7-Rozanska2014"
    for level in (
        "MOPAC:SQM@PM6",
        "MOPAC:SQM@PM6-ORG",
        "MOPAC:SQM@AM1",
        "ORCA:CMP@R2SCAN-3C",
    ):
        assert frequency_correction(level) is None


def test_the_published_equation(pm7):
    """eq 1: unchanged below 125; 1.006 nu - 27.7 to 2500; 3.99 nu - 7849.2 above."""
    raw = [50.0, 124.0, 126.0, 1000.0, 2499.0, 2501.0, 2800.0]
    expected = [50.0, 124.0, 99.056, 978.3, 2486.294, 2129.79, 3322.8]
    assert pm7(raw) == pytest.approx(expected)


def test_no_extrapolation_above_the_fit(pm7):
    """The fit ends at a raw ~2900 cm^-1 (3700 experimental); PM7's H2 at 5577
    would otherwise become 14402."""
    assert pm7.valid_max == 2900.0
    assert pm7([2899.0, 2901.0, 5577.0]) == pytest.approx([3717.81, 2901.0, 5577.0])
    assert pm7.not_applicable([2800.0, 2901.0, 5577.0]) == [2901.0, 5577.0]


def test_high_frequency_only_leaves_the_optional_branch_out(pm7):
    assert pm7([126.0, 1000.0, 2501.0], "high-frequency only") == pytest.approx(
        [126.0, 1000.0, 2129.79]
    )


def test_imaginary_and_zero_modes_are_unchanged(pm7):
    assert pm7([-200.0, 0.0]) == [-200.0, 0.0]


def test_the_discontinuities_are_flagged(pm7):
    """Raw frequencies near 125 and 2500 cm^-1, where the map jumps."""
    raw = [90.0, 100.0, 125.0, 126.0, 1500.0, 2399.0, 2400.0, 2600.0, 2601.0, 2750.0]
    assert pm7.in_unvalidated_windows(raw) == [100.0, 125.0, 2400.0, 2600.0]
    assert pm7.fitted_elements == ("C", "H", "N", "O")


def test_tag_and_ranges(pm7):
    assert pm7.tag() == "fc:PM7-Rozanska2014:v1:full"
    assert (
        pm7.tag("high-frequency only") == "fc:PM7-Rozanska2014:v1:high-frequency only"
    )
    assert RANGES == ("full", "high-frequency only")
    with pytest.raises(ValueError, match="must be one of"):
        pm7([1000.0], "partial")


def test_citation(pm7):
    assert "10.1021/je500201y" in pm7.bibtex
    assert "J. Chem. Eng. Data 59, 3136 (2014)" in pm7.citation
