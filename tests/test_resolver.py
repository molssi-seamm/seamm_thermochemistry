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


@pytest.mark.parametrize(
    "supplementary, found, label",
    [
        ("none", False, "ATcT-only"),
        ("compilations", False, "ATcT-only"),  # strictly experimental
        ("compilations and computed", True, "computed references"),
    ],
)
def test_w4_17_is_computed(layers, supplementary, found, label):
    """Paul, 2026-11-15: "W4-17 are computed values ... we should handle them
    as such" -- a literature source, but its values' status is computed."""
    r = RS.ReferenceResolver(FakeATcT(), supplementary=supplementary, layers=layers)
    entries, missing, got, _ = r.resolve(["C", "[AlH3]"])
    assert ("[AlH3]" in entries) == found and ("[AlH3]" in missing) != found
    assert got == label
    if found:
        assert entries["[AlH3]"]["source"] == "W4-17"
        assert entries["[AlH3]"]["computed"]


def test_an_experiment_beats_w4_17_whatever_the_rank(layers, tmp_path):
    """W4-17 ranks ahead of the WebBook, but a computed value comes after
    every experimental one."""
    lib = R.ReferenceLibrary(tmp_path / "official.db", create=True)
    sid = lib.add_species("[AlH3]")
    lib.add_measurement(
        sid, "WebBook-SRD69", dfh_298k=130.0, uncertainty=8.0, citation="exp"
    )
    lib.choose(sid, "WebBook-SRD69")
    lib.close()
    both = R.LayeredReferences(tmp_path / "local.db", tmp_path / "official.db")
    rows = both.candidates("[AlH3]")
    assert [r["source_name"] for r in rows] == ["WebBook-SRD69", "W4-17"]
    entries, _, label, _ = RS.ReferenceResolver(FakeATcT(), layers=both).resolve(
        ["[AlH3]"]
    )
    assert entries["[AlH3]"]["source"] == "WebBook-SRD69" and label == "mixed"


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


def test_unreadable_smiles_is_missing(layers):
    r = RS.ReferenceResolver(FakeATcT(), layers=layers)
    entries, missing, label, _ = r.resolve(["C", "?"])
    assert "?" in missing and entries["C"]["source"] == "ATcT"


def test_a_held_experiment_lets_the_computed_value_through(layers, tmp_path):
    """Science, 2026-11-15: a computed value beats an experiment for one
    species only when the experimental value is held, with the reason and
    who decided -- never by changing ranks."""
    lib = R.ReferenceLibrary(tmp_path / "official.db", create=True)
    sid = lib.add_species("[AlH3]")
    lib.add_measurement(
        sid, "WebBook-SRD69", dfh_298k=130.0, uncertainty=8.0, citation="exp"
    )
    lib.choose(sid, "WebBook-SRD69")
    lib.hold("[AlH3]", 0, "WebBook-SRD69", by="Paul", reason="W4-17 preferred")
    lib.close()
    both = R.LayeredReferences(tmp_path / "local.db", tmp_path / "official.db")
    r = RS.ReferenceResolver(FakeATcT(), layers=both)
    entries, _, label, _ = r.resolve(["[AlH3]"])
    assert entries["[AlH3]"]["source"] == "W4-17"
    assert label == "computed references"
    r = RS.ReferenceResolver(FakeATcT(), supplementary="compilations", layers=both)
    assert "[AlH3]" in r.resolve(["[AlH3]"])[1]  # held, and W4-17 is computed


def test_dfh_0k_from_298k_ga_as_se():
    """Ga, As and Se, which CODATA lacks (2026-11-16): AsH3 with NBS-82's own
    numbers, 66.44 at 298 K and 74.06 at 0 K, H298 - H0 10.22 (made-up check
    of the arithmetic: ΔfH0 = ΔfH298 - H(species) + Σ H(elements))."""
    table = RS.element_h298_h0()
    assert table["Ga"] == pytest.approx(5.561)
    assert table["As"] == pytest.approx(5.130)
    assert table["Se"] == pytest.approx(5.519)
    value = RS.dfh_0k_from_298k(66.44, 10.0, "AsH3")
    assert value == pytest.approx(66.44 - 10.0 + 5.130 + 3 * table["H"])
    RS.dfh_0k_from_298k(1.0, 1.0, "ClGaH2")  # no KeyError
    RS.dfh_0k_from_298k(1.0, 1.0, "H2Se")


@pytest.mark.parametrize("supplementary", ["compilations", "compilations and computed"])
def test_nea_ranks_ahead_of_nbs82(tmp_path, supplementary):
    """Paul, 2026-11-18: the OECD-NEA thermochemical database joins the
    hierarchy just behind CODATA, so its H2Se and SeF6 replace NBS-82's by
    rank, without holds, under both settings. Made-up values."""
    from seamm_thermochemistry import curated

    header = (
        "smiles,formula,charge,multiplicity,name,phase,dfH0_kJmol,dfH298_kJmol,"
        "unc_kJmol,unc_level,status,citation,notes\n"
    )
    (tmp_path / "nbs82.csv").write_text(
        header
        + '[SeH2],H2Se,0,1,H2Se,g,,30.0,,unstated,evaluated,"NBS-82 test",\n'
        + "F[Se](F)(F)(F)(F)F,F6Se,0,1,SeF6,g,,-1117.0,,unstated,evaluated,"
        + '"NBS-82 test",\n'
    )
    (tmp_path / "nea.csv").write_text(
        header
        + '[SeH2],H2Se,0,1,H2Se,g,,29.0,2.0,95 %,evaluated,"NEA test",\n'
        + "F[Se](F)(F)(F)(F)F,F6Se,0,1,SeF6,g,,-1118.0,0.4,95 %,evaluated,"
        + '"NEA test",\n'
    )
    library = R.ReferenceLibrary(
        tmp_path / "official.db", create=True, layer="official"
    )
    report, problems = curated.import_curated(library, tmp_path)
    assert problems == [] and all(r[2] == "added" for r in report)
    assert library.source("NEA-TDB-2005")["rank"] == 7
    assert curated.CURATED_SOURCES["codata.csv"][2] < 7
    assert 7 < library.source("NBS-82")["rank"]
    library.close()
    layers = R.LayeredReferences(tmp_path / "local.db", tmp_path / "official.db")
    r = RS.ReferenceResolver(FakeATcT(), supplementary=supplementary, layers=layers)
    entries, missing, label, _ = r.resolve(["[SeH2]", "F[Se](F)(F)(F)(F)F"])
    assert not missing and label == "mixed"
    assert entries["[SeH2]"]["source"] == "NEA-TDB-2005"
    assert entries["[SeH2]"]["dfH0_298K"] == 29.0
    assert entries["F[Se](F)(F)(F)(F)F"]["source"] == "NEA-TDB-2005"
    assert entries["F[Se](F)(F)(F)(F)F"]["uncertainty_298K"] == pytest.approx(0.4)
    layers.close()


def test_curated_ranks_follow_the_atomic_heat_hierarchy():
    """Paul, 2026-11-19: the library's sources in the order of the atomic
    heats, CODATA / NEA-TDB > printed JANAF > derived > Gurvich (TSIV) > NBS-82,
    then D0, the WebBook and the literature, and computed values last."""
    from seamm_thermochemistry import curated
    from seamm_thermochemistry import webbook

    rank = {v[0]: v[2] for v in curated.CURATED_SOURCES.values()}
    order = [
        "CODATA-1989",
        "NEA-TDB-2005",
        "JANAF-1998",
        "derived",
        "Gurvich-TSIV",
        "NBS-82",
        "D0",
    ]
    ranks = [rank[name] for name in order]
    assert ranks == sorted(ranks) and len(set(ranks)) == len(ranks)
    assert rank["D0"] < webbook.SOURCE_RANK < rank["Literature"] < R.COMPUTED_RANK
    assert rank["W4-17"] == R.COMPUTED_RANK
    assert curated.CURATED_SOURCES["gurvich.csv"][1] == "evaluated"


def _p_library(path, ranks):
    """A made-up library of P species (ATcT lacks P), as in v6: P4 in JANAF
    and CODATA, which agree with CODATA far more precise; PH3, PCl3 and POCl3
    in JANAF only; PCl5 derived. Plus a species in JANAF, Gurvich and NBS-82,
    and one in NEA-TDB and NBS-82."""
    lib = R.ReferenceLibrary(path, create=True, layer="official")
    for name, (kind, rank) in ranks.items():
        lib.add_source(name, kind, rank)

    def add(smiles, source, value, u, status="measured"):
        (sid,) = [s["id"] for s in lib.species(smiles)] or [lib.add_species(smiles)]
        lib.add_measurement(
            sid, source, dfh_298k=value, uncertainty=u, unc_level="2 sigma",
            status=status, citation=f"{source} {smiles}",
        )  # fmt: skip
        lib.choose(sid, source)

    add("P12P3P1P23", "JANAF-1998", 58.9, 4.2)
    add("P12P3P1P23", "CODATA-1989", 58.9, 0.3)
    add("P", "JANAF-1998", 5.4, 3.4)
    add("ClP(Cl)Cl", "JANAF-1998", -288.7, 10.8)
    add("O=P(Cl)(Cl)Cl", "JANAF-1998", -558.5, 2.9)
    add("ClP(Cl)(Cl)(Cl)Cl", "derived", -360.0, 5.0, status="derived")
    add("[Na]Cl", "JANAF-1998", -181.4, 2.1)
    add("[Na]Cl", "Gurvich-TSIV", -179.0, 4.0)
    add("[Na]Cl", "NBS-82", -176.0, 8.0)
    add("[Na]F", "NBS-82", -290.0, 8.0)
    add("[Na]F", "Gurvich-TSIV", -291.0, 4.0)
    add("[SeH2]", "NEA-TDB-2005", 29.0, 2.0)
    add("[SeH2]", "NBS-82", 29.7, 8.0)
    lib.apply_overrides()
    lib.close()


def test_reranking_changes_no_value(tmp_path):
    """The check behind Paul's 2026-11-19 decision. In library v6 the one
    species with both a JANAF and a CODATA value is P4, and the two agree
    (58.9 kJ/mol), so ranking CODATA and NEA-TDB above JANAF changes no value.
    Each species' best value is the same. But with one source per element, P4
    beside a JANAF-only P species now takes JANAF's value with JANAF's
    uncertainty, so the result is a single source ("mixed"). Before, CODATA
    overrode JANAF, which gave "mixed: several sources" and CODATA's smaller
    uncertainty. The library is rebuilt with set_rank, as for library v7."""
    from seamm_thermochemistry import curated

    old = {
        "JANAF-1998": ("evaluated", 10),
        "derived": ("derived", 11),
        "Gurvich-TSIV": ("evaluated", 15),
        "CODATA-1989": ("evaluated", 20),
        "NEA-TDB-2005": ("evaluated", 21),
        "NBS-82": ("evaluated", 25),
    }
    _p_library(tmp_path / "old.db", old)
    _p_library(tmp_path / "new.db", old)
    lib = R.ReferenceLibrary(tmp_path / "new.db", create=True, layer="official")
    for name, _, rank, *_ in curated.CURATED_SOURCES.values():
        if name in old:
            lib.set_rank(name, rank)
    lib.apply_overrides()
    assert lib.source("CODATA-1989")["rank"] < lib.source("JANAF-1998")["rank"]
    assert "rank 20 -> 6" in lib.source("CODATA-1989")["notes"]
    lib.close()

    groups = [
        ["P12P3P1P23"],
        ["P12P3P1P23", "P"],
        ["P12P3P1P23", "P", "ClP(Cl)Cl", "O=P(Cl)(Cl)Cl"],
        ["P", "ClP(Cl)(Cl)(Cl)Cl"],
        ["P12P3P1P23", "ClP(Cl)(Cl)(Cl)Cl"],
        ["[Na]Cl"],
        ["[Na]Cl", "[Na]F"],
        ["[SeH2]"],
    ]
    results = {}
    for name in ("old", "new"):
        layers = R.LayeredReferences(
            tmp_path / "local.db", tmp_path / f"{name}.db", tmp_path / "none.db"
        )
        r = RS.ReferenceResolver(
            FakeATcT(), supplementary="compilations", layers=layers
        )
        results[name] = []
        for group in groups:
            entries, missing, label, notes = r.resolve(group)
            results[name].append(
                (
                    {s: (e["source"], e["dfH0_298K"]) for s, e in entries.items()},
                    sorted(missing),
                    label,
                    notes,
                )
            )
        layers.close()
    # No value changes, in any group
    for (a, *_), (b, *_) in zip(results["old"], results["new"]):
        assert {s: v for s, (_, v) in a.items()} == {s: v for s, (_, v) in b.items()}
    # P4 alone: CODATA either way. Beside JANAF-only P species: JANAF now, a
    # single source, where CODATA's override made it several before
    for k, (old_r, new_r) in enumerate(zip(results["old"], results["new"])):
        if k in (1, 2):
            assert old_r[0]["P12P3P1P23"][0] == "CODATA-1989"
            assert old_r[2].startswith("mixed: several sources")
            assert new_r[0]["P12P3P1P23"][0] == "JANAF-1998"
            assert new_r[2] == "mixed"
        else:
            assert new_r == old_r, groups[k]
    assert results["new"][0][0]["P12P3P1P23"] == ("CODATA-1989", 58.9)
