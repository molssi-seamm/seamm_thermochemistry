"""ZPE scale factors by level of theory, with their sources and flags."""

import pytest

from seamm_thermochemistry import zpe_scale_factor
from seamm_thermochemistry.zpe import split_level


@pytest.mark.parametrize(
    "level, expected",
    [
        ("ORCA:CMP@R2SCAN-3C", ("R2SCAN-3C", "")),
        ("ORCA:DFT@PBE0-D4/def2-TZVP", ("PBE0", "DEF2-TZVP")),
        ("ORCA:DFT@PBE0/bse:def2-TZVP", ("PBE0", "DEF2-TZVP")),
        ("Gaussian:DFT@B2PLYP-D3/cc-pVTZ", ("B2PLYP", "CC-PVTZ")),
        ("ORCA:DFT@WB97X-V/def2-TZVP", ("WB97X-V", "DEF2-TZVP")),
        ("", ("", "")),
    ],
)
def test_split_level(level, expected):
    assert split_level(level) == expected


def test_factors():
    r2 = zpe_scale_factor("ORCA:CMP@R2SCAN-3C")
    assert r2.factor == 0.9899 and "W4-2006" in r2.source
    assert "extrapolation" in r2.flag
    pbe0 = zpe_scale_factor("ORCA:DFT@PBE0-D4/def2-TZVP")
    assert pbe0.factor == 0.975 and "basis-mismatched" in pbe0.flag
    assert zpe_scale_factor("ORCA:DFT@B2PLYP-D3/cc-pVTZ").factor == 0.9863
    # Not in the table: no factor (the caller warns and uses 1)
    assert zpe_scale_factor("ORCA:DFT@PBE0-D4/def2-SVP") is None
    assert zpe_scale_factor("MOPAC:SQM@PM7") is None
