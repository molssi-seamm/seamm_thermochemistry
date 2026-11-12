"""Resolving CBH references from ATcT and the reference library."""

import pytest

pytest.importorskip("rdkit")

from seamm_thermochemistry import references as R  # noqa: E402
from seamm_thermochemistry import resolver as RS  # noqa: E402
from seamm_thermochemistry.atct import MissingEntry  # noqa: E402

J = "JANAF 4th ed."


class FakeATcT:
    """ATcT with C, H, O, F, Cl species and atoms, but no P, Na or K."""

    version = "test"
    values = {
        "C": -74.5,
        "[H][H]": 0.0,
        "O": -241.8,
        "Cl": -92.2,
        "F": -272.7,
        "ClCl": 0.0,
        "[C]": 716.9,
        "[H]": 218.0,
        "[O]": 249.2,
        "[F]": 79.4,
        "[Cl]": 121.3,
    }

    def lookup(self, smiles, charge=0, multiplicity=None, **kw):
        from seamm_thermochemistry.atct import canonical_smiles

        can = canonical_smiles(smiles)
        for key, value in self.values.items():
            if canonical_smiles(key) == can:
                return {
                    "id": f"atct-{key}",
                    "dfH0_0K": value,
                    "dfH0_298K": value,
                    "uncertainty_298K": 0.1,
                    "exact": False,
                }
        raise MissingEntry(smiles)


@pytest.fixture
def layers(tmp_path):
    lib = R.ReferenceLibrary(tmp_path / "official.db", create=True, layer="official")
    lib.add_source("JANAF-1998", "evaluated", 10)
    lib.add_source("D0", "measured", 30)
    lib.add_source("W4-17", "literature", 40)
    lib.add_source("WebBook-SRD69", "measured", 50)
    lib.add_source("computed", "computed", 60)

    def add(smiles, source, value, u=1.0, level="2 sigma", status="measured", v0=None):
        sid = lib.add_species(smiles)
        lib.add_measurement(
            sid,
            source,
            dfh_298k=value,
            dfh_0k=v0,
            uncertainty=u,
            unc_level=level,
            status=status,
            citation=f"{source} {smiles}",
        )
        lib.choose(sid, source)
        return sid

    # P: JANAF has PH3 and PCl3; the WebBook has PH3 and CH3PH2 (JANAF lacks it)
    add("P", "JANAF-1998", 5.4, 3.4)
    add("ClP(Cl)Cl", "JANAF-1998", -288.7, 10.8)
    add("P", "WebBook-SRD69", 6.0, 2.0)
    add("CP", "WebBook-SRD69", -20.0, 4.0)
    # Na: JANAF 124.26 +- 38.4 overridden by D0 138.59 +- 1.39
    add("[NaH]", "JANAF-1998", 124.26, 38.4)
    add("[NaH]", "D0", 138.59, 1.39)
    add("[Na]F", "JANAF-1998", -290.5, 2.0)
    # Al: W4-17 at 0 K only
    sid = lib.add_species("[AlH3]")
    lib.add_measurement(
        sid, "W4-17", dfh_0k=135.12, uncertainty=4.47, unc_level="mixed",
        status="computed", citation="W4-17 AlH3",
    )  # fmt: skip
    lib.choose(sid, "W4-17")
    # a computed reference
    add("[K]Cl", "computed", -214.0, 2.0, status="computed")
    lib.apply_overrides()
    lib.close()
    yield R.LayeredReferences(tmp_path / "local.db", tmp_path / "official.db")


def test_supplementary_none_is_atct_only(layers):
    r = RS.ReferenceResolver(FakeATcT(), supplementary="none", layers=layers)
    entries, missing, label, notes = r.resolve(["C", "P"])
    assert entries["C"]["source"] == "ATcT" and "P" in missing
    assert label == "ATcT-only"
    assert r.layers is None  # nothing is read


def test_atct_wins_and_label(layers):
    r = RS.ReferenceResolver(FakeATcT(), layers=layers)
    entries, missing, label, _ = r.resolve(["C", "O"])
    assert label == "ATcT-only" and not missing
    assert entries["C"]["atct_coefficients"] == {"atct-C": 1.0}


def test_one_source_per_element(layers):
    r = RS.ReferenceResolver(FakeATcT(), layers=layers)
    # PH3 alone: JANAF (rank 10)
    entries, _, label, _ = r.resolve(["C", "P"])
    assert entries["P"]["source"] == "JANAF-1998" and label == "mixed"
    assert entries["P"]["atct_coefficients"] == {}
    assert entries["P"]["uncertainty_298K"] == pytest.approx(3.4)
    # PH3 + CH3PH2: only the WebBook has both, so both come from it
    entries, _, label, notes = r.resolve(["P", "CP"])
    assert entries["P"]["source"] == entries["CP"]["source"] == "WebBook-SRD69"
    assert label == "mixed"
    # PH3 + PCl3 + CH3PH2: no single source; each its best; several sources
    entries, _, label, notes = r.resolve(["P", "ClP(Cl)Cl", "CP"])
    assert entries["P"]["source"] == "JANAF-1998"
    assert entries["CP"]["source"] == "WebBook-SRD69"
    assert label.startswith("mixed: several sources") and "P:" in label
    assert any("several sources" in n for n in notes)


def test_override_inside_one_source_per_element(layers):
    r = RS.ReferenceResolver(FakeATcT(), layers=layers)
    entries, _, label, notes = r.resolve(["[NaH]", "[Na]F"])
    # JANAF has both, but D0 replaces JANAF for NaH: two sources for Na
    assert entries["[NaH]"]["source"] == "D0"
    assert entries["[Na]F"]["source"] == "JANAF-1998"
    assert "Na: D0, JANAF-1998" in label
    assert "below a third" in entries["[NaH]"]["override_reason"]


def test_0k_only_and_computed(layers):
    r = RS.ReferenceResolver(FakeATcT(), layers=layers)
    entries, missing, label, _ = r.resolve(["[AlH3]", "[K]Cl"])
    assert entries["[AlH3]"]["dfH0_298K"] is None
    assert entries["[AlH3]"]["dfH0_0K"] == 135.12
    assert label == "computed references"
    r = RS.ReferenceResolver(FakeATcT(), supplementary="compilations", layers=layers)
    entries, missing, label, _ = r.resolve(["[K]Cl"])
    assert "[K]Cl" in missing and label == "ATcT-only"


def test_user_references(tmp_path, layers):
    local = R.ReferenceLibrary(tmp_path / "local.db", create=True)
    local.add_source("user", "user", 5)
    sid = local.add_species("P")
    local.add_measurement(
        sid, "user", dfh_298k=5.0, uncertainty=0.5, status="user", citation="mine"
    )
    local.choose(sid, "user")
    local.close()
    both = R.LayeredReferences(tmp_path / "local.db", tmp_path / "official.db")
    r = RS.ReferenceResolver(FakeATcT(), layers=both)
    assert r.resolve(["P"])[0]["P"]["source"] == "JANAF-1998"  # not enabled
    r = RS.ReferenceResolver(FakeATcT(), layers=both, use_user=True)
    entries, _, label, _ = r.resolve(["P"])
    assert entries["P"]["source"] == "user" and label == "user references"


def test_dfh_0k_from_298k():
    # H2O: 298 K -241.826, H298-H0 9.904; H 8.468/2 per atom, O 8.683/2
    value = RS.dfh_0k_from_298k(-241.826, 9.904, "H2O")
    table = RS.element_h298_h0()
    assert value == pytest.approx(-241.826 - 9.904 + 2 * table["H"] + table["O"])
    assert value == pytest.approx(-238.92, abs=0.05)
