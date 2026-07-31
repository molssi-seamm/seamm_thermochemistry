"""Tests pinning down the sign conventions in formation.py against
hand-worked numbers, and against the existing gaussian_step/psi4_step
(anchor=True) and vasp_step (anchor=False) conventions they must match.
"""

import pytest

from seamm_thermochemistry import (
    ThermoDB,
    atomization_energy,
    formation_energy,
    formation_enthalpy,
    formation_gibbs_energy,
    MissingReferenceData,
)


@pytest.fixture
def db(tmp_path):
    with ThermoDB(tmp_path / "test.db") as db:
        # Round numbers chosen so the arithmetic is easy to check by hand.
        db.add_element(1, "H", dfH0_0K=100.0, dfH0_298K=110.0)
        db.add_element(8, "O", dfH0_0K=200.0, dfH0_298K=210.0)
        db.set_s298_std_state("H", 50.0, reference="test")
        db.set_s298_std_state("O", 80.0, reference="test")

        # Gas-atom (molecular-code) convention.
        db.add_atom_energy("H", "test", "M", -10.0, ref_type="atom")
        db.add_atom_energy("O", "test", "M", -70.0, ref_type="atom")

        # Standard-state-phase (VASP) convention, same code/method label.
        db.add_atom_energy("H", "test", "M", -5.0, ref_type="element_phase")
        db.add_atom_energy("O", "test", "M", -60.0, ref_type="element_phase")
        yield db


H2O = {"H": 2, "O": 1}
SYSTEM_ENERGY = -100.0  # kJ/mol, arbitrary ZPE-free electronic energy


def test_atomization_energy(db):
    # sum(n * E_ref) - E(system) = (2*-10 + -70) - (-100) = 10.0
    assert atomization_energy(H2O, SYSTEM_ENERGY, db, "test", "M") == pytest.approx(
        10.0
    )


def test_atomization_energy_positive_for_bound_system(db):
    # A system less stable than its separated atoms would combine to a
    # different sign, but for a normal bound molecule (energy below the
    # separated atoms) atomization energy is positive.
    assert atomization_energy(H2O, SYSTEM_ENERGY, db, "test", "M") > 0


def test_formation_energy_anchor_0K(db):
    # sum(n * DfH0_0K) - atomization = (2*100 + 200) - 10 = 390.0
    result = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M", anchor_at_0K=True)
    assert result == pytest.approx(390.0)


def test_formation_energy_anchor_298K(db):
    # sum(n * DfH0_298K) - atomization = (2*110 + 210) - 10 = 420.0
    result = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M", anchor_at_0K=False)
    assert result == pytest.approx(420.0)


def test_formation_energy_matches_dfH0_minus_atomization(db):
    # By construction: formation(anchor=True) == sum(anchor) - atomization,
    # the same shape as gaussian_step's `DfH_at - data["H atomization"]`.
    atomization = atomization_energy(H2O, SYSTEM_ENERGY, db, "test", "M")
    anchor_sum = 2 * db.dfH0("H") + 1 * db.dfH0("O")
    result = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M")
    assert result == pytest.approx(anchor_sum - atomization)


def test_formation_energy_no_anchor_matches_vasp_DfE0_convention(db):
    # VASP: DfE0 = E(system) - sum(n_X * mu_X) = -atomization_energy(...)
    atomization = atomization_energy(
        H2O, SYSTEM_ENERGY, db, "test", "M", ref_type="element_phase"
    )
    result = formation_energy(
        H2O, SYSTEM_ENERGY, db, "test", "M", ref_type="element_phase", anchor=False
    )
    assert result == pytest.approx(-atomization)
    # And directly, by the VASP formula:
    mu_sum = 2 * (-5.0) + 1 * (-60.0)
    assert result == pytest.approx(SYSTEM_ENERGY - mu_sum)


def test_formation_energy_no_anchor_needs_no_experimental_data(tmp_path):
    # A DB with atom energies but zero experimental dfH0 data must still
    # work for anchor=False (that's the whole point of the materials
    # convention: no experimental anchor required).
    with ThermoDB(tmp_path / "no_anchor.db") as db2:
        db2.add_element(1, "H")  # dfH0_0K / dfH0_298K left as NULL
        db2.add_element(8, "O")
        db2.add_atom_energy("H", "test", "M", -5.0, ref_type="element_phase")
        db2.add_atom_energy("O", "test", "M", -60.0, ref_type="element_phase")
        result = formation_energy(
            H2O, SYSTEM_ENERGY, db2, "test", "M", ref_type="element_phase", anchor=False
        )
        assert result == pytest.approx(SYSTEM_ENERGY - (2 * -5.0 + -60.0))


def test_missing_atom_energy_raises(db):
    with pytest.raises(MissingReferenceData):
        formation_energy({"H": 1, "Cl": 1}, SYSTEM_ENERGY, db, "test", "M")


def test_missing_experimental_anchor_raises(db, tmp_path):
    with ThermoDB(tmp_path / "no_anchor.db") as db2:
        db2.add_element(1, "H")  # no dfH0 values at all
        db2.add_atom_energy("H", "test", "M", -10.0)
        with pytest.raises(MissingReferenceData):
            formation_energy({"H": 2}, -5.0, db2, "test", "M")


def test_formation_enthalpy_at_300K(db):
    # anchor=(2*100+200)=400; ref_sum=(2*-10+-70)=-90; n_atoms=3;
    # thermal/atom = 2.5*R*300 = 6.235846963614930 kJ/mol.
    # DfH(300K) = 400 - 3*6.235846963614930 - (-90) + (-95) = 376.29245911...
    result = formation_enthalpy(H2O, -95.0, 300.0, db, "test", "M")
    assert result == pytest.approx(376.29245910915523)


def test_formation_enthalpy_reduces_to_formation_energy_at_0K(db):
    # At T=0 the ideal-gas thermal term vanishes, so DfH(0) == DfE0 with
    # the same (here ZPE-free) system energy.
    dfh0 = formation_enthalpy(H2O, SYSTEM_ENERGY, 0.0, db, "test", "M")
    dfe0 = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M", anchor_at_0K=True)
    assert dfh0 == pytest.approx(dfe0)


def test_formation_enthalpy_trivial_single_atom_is_T_independent(db):
    # Decisive regression check for the sign of the (5/2)RT term: for the
    # degenerate "molecule" that IS just one isolated reference atom (no
    # atomization step at all), DfH(T) must reduce to the plain 0 K anchor
    # for *any* T -- H_system's own (5/2)RT growth with T is exactly
    # cancelled by the same atom's thermal term in the Hess cycle. A wrong
    # sign here would leave a spurious 2*(5/2)RT term that grows with T.
    R = 8.31446261815324e-3
    for T in (0.0, 300.0, 500.0):
        H_system = -10.0 + 2.5 * R * T  # H atom's own E_ref + (5/2)RT
        result = formation_enthalpy({"H": 1}, H_system, T, db, "test", "M")
        assert result == pytest.approx(100.0)  # == db.dfH0("H"), for every T


def test_formation_enthalpy_missing_anchor_raises(db, tmp_path):
    with ThermoDB(tmp_path / "no_anchor.db") as db2:
        db2.add_element(1, "H")  # no dfH0 at all
        db2.add_atom_energy("H", "test", "M", -10.0)
        with pytest.raises(MissingReferenceData):
            formation_enthalpy({"H": 2}, -5.0, 300.0, db2, "test", "M")


def test_formation_gibbs_energy_at_300K(db):
    # std_state_term = 300*(2*50+80)/1000 = 300*180/1000 = 54.0
    # DfG(300K) = 400 + 54.0 - 3*6.235846963614930 - (-90) + (-96)
    #           = 429.2924591091553
    result = formation_gibbs_energy(H2O, -96.0, 300.0, db, "test", "M")
    assert result == pytest.approx(429.2924591091553)


def test_formation_gibbs_energy_trivial_single_atom(db):
    # Same decisive check as the enthalpy case, extended to G: for the
    # degenerate "molecule" that IS just one isolated reference atom,
    # DfG(T) = DfH0(0K) - T*(S_eff - S298_std_state)/1000 for WHATEVER
    # entropy S_eff went into building the atom's own G(T) -- i.e. the
    # atom's own entropy must cancel out of the result exactly, leaving
    # only the standard-state-element entropy behind. S_eff=90.0 here is
    # an arbitrary placeholder (not read from the DB, not a real physical
    # value) -- the point is the result must not depend on its value
    # except through G_system itself.
    R = 8.31446261815324e-3
    s_eff = 90.0
    for T in (0.0, 300.0, 500.0):
        G_system = -10.0 + 2.5 * R * T - T * s_eff / 1000.0
        result = formation_gibbs_energy({"H": 1}, G_system, T, db, "test", "M")
        expected = 100.0 - T * (s_eff - 50.0) / 1000.0
        assert result == pytest.approx(expected)


def test_formation_gibbs_energy_reduces_to_formation_energy_at_0K(db):
    # At T=0 both the ideal-gas and standard-state-entropy terms vanish,
    # so DfG(0) == DfE0 with the same (here ZPE-inclusive) system energy --
    # the derivation's own sanity check.
    dfg0 = formation_gibbs_energy(H2O, SYSTEM_ENERGY, 0.0, db, "test", "M")
    dfe0 = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M", anchor_at_0K=True)
    assert dfg0 == pytest.approx(dfe0)


def test_formation_gibbs_energy_missing_std_state_entropy_raises(db, tmp_path):
    with ThermoDB(tmp_path / "no_std_state.db") as db2:
        # dfH0 present, but s298_std_state was never set.
        db2.add_element(1, "H", dfH0_0K=100.0)
        db2.add_atom_energy("H", "test", "M", -10.0)
        with pytest.raises(MissingReferenceData):
            formation_gibbs_energy({"H": 2}, -5.0, 300.0, db2, "test", "M")


def test_formation_gibbs_energy_missing_anchor_raises(db, tmp_path):
    with ThermoDB(tmp_path / "no_anchor.db") as db2:
        # s298_std_state present, but dfH0 was never set.
        db2.add_element(1, "H")
        db2.set_s298_std_state("H", 50.0)
        db2.add_atom_energy("H", "test", "M", -10.0)
        with pytest.raises(MissingReferenceData):
            formation_gibbs_energy({"H": 2}, -5.0, 300.0, db2, "test", "M")


def test_unit_consistency_eV_vs_kJmol(db):
    # 1 hartree = 2625.499639 kJ/mol; use eV throughout and cross check.
    db.add_atom_energy("H", "test", "N", -10.0 / 96.485332, ref_type="atom", units="eV")
    db.add_atom_energy("O", "test", "N", -70.0 / 96.485332, ref_type="atom", units="eV")
    result_eV = formation_energy(
        H2O, SYSTEM_ENERGY / 96.485332, db, "test", "N", units="eV"
    )
    result_kJ = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M", units="kJ/mol")
    assert result_eV * 96.485332 == pytest.approx(result_kJ, rel=1e-6)
