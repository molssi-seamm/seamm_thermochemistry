"""Tests for derive_dfH0_0K, checked against the two elements the master
sheet gives directly at 0K (H, He) and against a diatomic-vs-monatomic
standard-state case (O vs C)."""

import pytest

from seamm_thermochemistry import ThermoDB, derive_dfH0_0K


@pytest.fixture
def db(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        # H: diatomic standard state (1/2 H2), master's own 0K value known.
        db.add_element(
            1,
            "H",
            standard_state="1/2 H2(g)",
            dfH0_298K=218.0,
            h298_minus_h0_atom=6.197,
            h298_minus_h0_std_state=8.468,
        )
        # He: monatomic standard state, master's own 0K value known (0).
        db.add_element(
            2,
            "He",
            standard_state="He(g)",
            dfH0_298K=0.0,
            h298_minus_h0_atom=6.197,
            h298_minus_h0_std_state=6.197,
        )
        # C: monatomic standard state (graphite); no master 0K value to
        # check against, but exercises the "no diatomic divisor" path for
        # an element actually used in organic-molecule formation energies.
        db.add_element(
            6,
            "C",
            standard_state="C(s,gr)",
            dfH0_298K=716.67,
            h298_minus_h0_atom=6.536,
            h298_minus_h0_std_state=1.05,
        )
        # Missing data -- must be skipped, not raise.
        db.add_element(3, "Li", standard_state="Li(s)")
        yield db


def test_derives_H_matching_master_value(db):
    filled = derive_dfH0_0K(db)
    assert "H" in filled
    # Master sheet's own tabulated value: 216.034
    assert db.get_element(symbol="H")["dfH0_0K"] == pytest.approx(216.034, abs=0.01)


def test_derives_He_matching_master_value(db):
    filled = derive_dfH0_0K(db)
    assert "He" in filled
    assert db.get_element(symbol="He")["dfH0_0K"] == pytest.approx(0.0, abs=0.01)


def test_derives_monatomic_standard_state_C(db):
    derive_dfH0_0K(db)
    # 716.67 - 6.536 + 1.05 (no /2: graphite is the monatomic formula unit)
    expected = 716.67 - 6.536 + 1.05
    assert db.get_element(symbol="C")["dfH0_0K"] == pytest.approx(expected)


def test_skips_element_with_missing_data(db):
    filled = derive_dfH0_0K(db)
    assert "Li" not in filled
    assert db.get_element(symbol="Li")["dfH0_0K"] is None


def test_does_not_overwrite_existing_by_default(tmp_path):
    with ThermoDB(tmp_path / "test2.db") as db:
        db.add_element(
            1,
            "H",
            dfH0_0K=999.0,  # a deliberately "wrong" existing value
            dfH0_298K=218.0,
            h298_minus_h0_atom=6.197,
            h298_minus_h0_std_state=8.468,
        )
        filled = derive_dfH0_0K(db)
        assert "H" not in filled
        assert db.get_element(symbol="H")["dfH0_0K"] == pytest.approx(999.0)


def test_overwrite_true_recomputes(tmp_path):
    with ThermoDB(tmp_path / "test3.db") as db:
        db.add_element(
            1,
            "H",
            dfH0_0K=999.0,
            dfH0_298K=218.0,
            h298_minus_h0_atom=6.197,
            h298_minus_h0_std_state=8.468,
        )
        filled = derive_dfH0_0K(db, overwrite=True)
        assert "H" in filled
        assert db.get_element(symbol="H")["dfH0_0K"] == pytest.approx(216.034, abs=0.01)
