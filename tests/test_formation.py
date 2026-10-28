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
        # h298_minus_h0_std_state is per formula unit of the standard state:
        # H2 and O2 here, so 4.0 and 6.0 per atom.
        db.add_element(
            1,
            "H",
            standard_state="1/2 H2(g)",
            dfH0_0K=100.0,
            dfH0_298K=110.0,
            h298_minus_h0_std_state=8.0,
        )
        db.add_element(
            8,
            "O",
            standard_state="1/2 O2(g)",
            dfH0_0K=200.0,
            dfH0_298K=210.0,
            h298_minus_h0_std_state=12.0,
        )
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


def test_formation_enthalpy_at_298K(db):
    # anchor = sum n*(DfH0(0K) - h_std) = 2*(100-4) + (200-6) = 386;
    # ref_sum = (2*-10 + -70) = -90.
    # DfH(298.15) = 386 - (-90) + (-95) = 381.0
    result = formation_enthalpy(H2O, -95.0, 298.15, db, "test", "M")
    assert result == pytest.approx(381.0)


def test_formation_enthalpy_at_300K_scales_std_state_term(db):
    # Away from 298.15 K the standard-state H(T)-H(0) is scaled linearly:
    # h_std(300) = (2*4 + 6) * 300/298.15 = 14.086869...
    # DfH(300) = 400 - 14.086869 + 90 - 95 = 380.913131...
    scale = 300.0 / 298.15
    result = formation_enthalpy(H2O, -95.0, 300.0, db, "test", "M")
    assert result == pytest.approx(400.0 - 14.0 * scale + 90.0 - 95.0)


def test_formation_enthalpy_reduces_to_formation_energy_at_0K(db):
    # At T=0 the standard-state thermal term vanishes, so DfH(0) == DfE0 with
    # the same (here ZPE-free) system energy.
    dfh0 = formation_enthalpy(H2O, SYSTEM_ENERGY, 0.0, db, "test", "M")
    dfe0 = formation_energy(H2O, SYSTEM_ENERGY, db, "test", "M", anchor_at_0K=True)
    assert dfh0 == pytest.approx(dfe0)


def test_formation_enthalpy_single_atom_gives_DfH_of_the_atom_at_T(db):
    # The degenerate "molecule" that IS one gas atom: DfH(T) must be the
    # atom's own DfH(X, T) = DfH(X, 0K) + [H_atom(T)-H_atom(0)] -
    # [H_std(X, T)-H_std(X, 0)], which DOES depend on T (H: 216.03 at 0 K,
    # 218.00 at 298 K). The atom's own thermal enthalpy enters only through
    # H_system; the standard state's is subtracted by the function.
    R = 8.31446261815324e-3
    for T in (0.0, 298.15, 500.0):
        H_system = -10.0 + 2.5 * R * T  # H atom's own E_ref + (5/2)RT
        result = formation_enthalpy({"H": 1}, H_system, T, db, "test", "M")
        expected = 100.0 + 2.5 * R * T - 4.0 * T / 298.15
        assert result == pytest.approx(expected)


def test_formation_enthalpy_water_JANAF_298K_minus_0K(tmp_path):
    # Real JANAF data for H and O, as in the published database: the
    # difference DfH(298.15) - DfH(0) must equal the molecule's own
    # H(298)-H(0) minus the elements', and reproduce experiment:
    # water DfH(0) = -238.92, DfH(298.15) = -241.83 kJ/mol (JANAF), i.e.
    # -2.91; with H(298)-H(0) = 9.904 for H2O(g), 8.468 for H2 and 8.683
    # for O2: 9.904 - 8.468 - 8.683/2 = -2.906. The old (5/2)RT-per-atom
    # formula gave 9.904 - 3*6.197 = -8.69, 5.8 kJ/mol too negative.
    with ThermoDB(tmp_path / "janaf.db") as db2:
        db2.add_element(
            1,
            "H",
            standard_state="1/2 H2(g)",
            dfH0_0K=216.034,
            dfH0_298K=217.998,
            h298_minus_h0_atom=6.197,
            h298_minus_h0_std_state=8.468,
        )
        db2.add_element(
            8,
            "O",
            standard_state="1/2 O2(g)",
            dfH0_0K=246.795,
            dfH0_298K=249.18,
            h298_minus_h0_atom=6.725,
            h298_minus_h0_std_state=8.683,
        )
        db2.add_atom_energy("H", "test", "M", -10.0, ref_type="atom")
        db2.add_atom_energy("O", "test", "M", -70.0, ref_type="atom")

        E0 = -100.0  # electronic + ZPE, arbitrary
        dfh_0 = formation_enthalpy(H2O, E0, 0.0, db2, "test", "M")
        dfh_298 = formation_enthalpy(H2O, E0 + 9.904, 298.15, db2, "test", "M")
        assert dfh_298 - dfh_0 == pytest.approx(9.904 - 8.468 - 8.683 / 2)
        assert dfh_298 - dfh_0 == pytest.approx(-2.91, abs=0.01)


def test_std_state_thermal_enthalpy_halves_diatomics_only(tmp_path):
    from seamm_thermochemistry.formation import _std_state_thermal_enthalpy

    with ThermoDB(tmp_path / "std.db") as db2:
        db2.add_element(
            1, "H", standard_state="1/2 H2(g)", h298_minus_h0_std_state=8.468
        )
        db2.add_element(6, "C", standard_state="C(s,gr)", h298_minus_h0_std_state=1.05)
        db2.add_element(2, "He", standard_state="He(g)")
        assert _std_state_thermal_enthalpy(db2, "H", 298.15) == pytest.approx(4.234)
        assert _std_state_thermal_enthalpy(db2, "C", 298.15) == pytest.approx(1.05)
        assert _std_state_thermal_enthalpy(db2, "C", 0.0) == pytest.approx(0.0)
        assert _std_state_thermal_enthalpy(db2, "He", 298.15) is None


def test_formation_enthalpy_missing_anchor_raises(db, tmp_path):
    with ThermoDB(tmp_path / "no_anchor.db") as db2:
        db2.add_element(1, "H")  # no dfH0 at all
        db2.add_atom_energy("H", "test", "M", -10.0)
        with pytest.raises(MissingReferenceData):
            formation_enthalpy({"H": 2}, -5.0, 300.0, db2, "test", "M")


def test_formation_gibbs_energy_at_298K(db):
    # std_state_term = 298.15*(2*50+80)/1000 = 53.667
    # DfG(298.15) = 386 + 53.667 - (-90) + (-96) = 433.667
    result = formation_gibbs_energy(H2O, -96.0, 298.15, db, "test", "M")
    assert result == pytest.approx(386.0 + 298.15 * 180.0 / 1000.0 + 90.0 - 96.0)


def test_formation_gibbs_energy_single_atom(db):
    # The degenerate "molecule" that IS one gas atom: DfG(T) must be
    # DfH(X, 0K) + [H_atom(T)-H_atom(0)] - [H_std(X,T)-H_std(X,0)]
    # - T*(S_eff - S298_std_state)/1000 for WHATEVER entropy S_eff went
    # into the atom's own G(T): the atom's entropy enters only through
    # G_system, and only the standard-state element's is added back.
    R = 8.31446261815324e-3
    s_eff = 90.0
    for T in (0.0, 298.15, 500.0):
        G_system = -10.0 + 2.5 * R * T - T * s_eff / 1000.0
        result = formation_gibbs_energy({"H": 1}, G_system, T, db, "test", "M")
        expected = 100.0 + 2.5 * R * T - 4.0 * T / 298.15 - T * (s_eff - 50.0) / 1000.0
        assert result == pytest.approx(expected)


def test_formation_gibbs_energy_equals_DfH_minus_T_DfS(db):
    # With G_system = H_system - T*S_system, DfG = DfH - T*(S_sys - sum n S_std)
    T = 298.15
    H_system = -95.0
    S_system = 188.8  # J/mol/K
    G_system = H_system - T * S_system / 1000.0
    dfh = formation_enthalpy(H2O, H_system, T, db, "test", "M")
    dfg = formation_gibbs_energy(H2O, G_system, T, db, "test", "M")
    dfs = (S_system - (2 * 50.0 + 80.0)) / 1000.0
    assert dfg == pytest.approx(dfh - T * dfs)


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
