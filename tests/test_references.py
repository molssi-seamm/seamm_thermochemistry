"""The supplementary reference library: schema, rules and layers.

Values are the reference-data pilot's (2026-11-12), so the tests also record
what the rules decide for real cases.
"""

import sqlite3

import pytest

pytest.importorskip("rdkit")

from seamm_thermochemistry import references as R  # noqa: E402
from seamm_thermochemistry import references_cli  # noqa: E402

JANAF = "NIST-JANAF Thermochemical Tables, 4th ed. (Chase 1998)"


@pytest.fixture
def lib(tmp_path):
    library = R.ReferenceLibrary(tmp_path / "references.db", create=True)
    library.add_source("JANAF-1998", "evaluated", 10, citation=JANAF)
    library.add_source("CODATA-1989", "evaluated", 20)
    library.add_source("D0", "measured", 30)
    library.add_source("W4-17", "literature", 40)
    library.add_source("WebBook-SRD69", "measured", 50)
    library.add_source("derived", "derived", 15)
    yield library
    library.close()


def test_schema_version_and_layer(lib):
    db = sqlite3.connect(lib.path)
    assert db.execute("PRAGMA user_version").fetchone()[0] == R.SCHEMA_VERSION
    assert lib.layer == "local"


def test_newer_schema_is_never_written(tmp_path):
    path = tmp_path / "r.db"
    R.ReferenceLibrary(path, create=True).close()
    db = sqlite3.connect(path)
    db.execute(f"PRAGMA user_version = {R.SCHEMA_VERSION + 1}")
    db.commit()
    db.close()
    with pytest.raises(R.NewerReferenceSchemaError):
        R.ReferenceLibrary(path, create=True)
    # but it can still be read
    assert R.ReferenceLibrary(path).sources() == []


def test_burcat_is_never_a_source(lib):
    with pytest.raises(R.ReferenceError, match="Burcat"):
        lib.add_source("Burcat-2024", "evaluated", 90)
    sid = lib.add_species("P", formula="PH3")
    with pytest.raises(R.ReferenceError, match="Burcat"):
        lib.add_measurement(
            sid, "JANAF-1998", dfh_298k=5.47, citation="Burcat & Ruscic, ANL-05/20"
        )


def test_no_p_species_from_online_janaf(lib):
    lib.add_source("JANAF-online", "evaluated", 11)
    ph3 = lib.add_species("P", formula="H3P")
    with pytest.raises(R.ReferenceError, match="printed edition"):
        lib.add_measurement(ph3, "JANAF-online", dfh_298k=22.89, citation=JANAF)
    # the printed edition is fine, and other elements from the online tables
    lib.add_measurement(
        ph3, "JANAF-1998", dfh_298k=5.47, uncertainty=1.7, citation=JANAF
    )
    hf = lib.add_species("F", formula="HF")
    lib.add_measurement(hf, "JANAF-online", dfh_298k=-272.55, citation=JANAF)


def test_identity_checks(lib):
    with pytest.raises(R.ReferenceError, match="does not match"):
        lib.add_species("ClP(Cl)(Cl)(Cl)Cl", formula="Cl3P")
    pcl5 = lib.add_species("ClP(Cl)(Cl)(Cl)Cl", formula="Cl5P", cas="10026-13-8")
    assert lib.add_species("P(Cl)(Cl)(Cl)(Cl)Cl") == pcl5  # same canonical key
    with pytest.raises(R.ReferenceError, match="CAS"):
        lib.add_species("ClP(Cl)Cl", cas="10026-13-8")


def test_ions_need_and_get_the_stationary_electron_convention(lib):
    hp = lib.add_species("[H+]", charge=1, multiplicity=1)
    with pytest.raises(R.ReferenceError, match="convention"):
        lib.add_measurement(hp, "JANAF-1998", dfh_298k=1536.246, citation=JANAF)
    m = lib.add_measurement(
        hp,
        "JANAF-1998",
        dfh_298k=1536.246,
        dfh_0k=1528.085,
        uncertainty=0.04,
        ion_convention="thermal-electron",
        citation=JANAF,
    )
    row = lib._measurement(m)
    # ATcT's H+ at 298.15 K is 1530.05 kJ/mol (stationary electron)
    assert row["dfh_298k"] == pytest.approx(1530.049, abs=0.002)
    assert row["dfh_0k"] == 1528.085
    assert "thermal-electron" in row["converted"]


def test_expanded_uncertainty_policy():
    assert R.expanded_uncertainty(1.0, "2 sigma") == 1.0
    assert R.expanded_uncertainty(1.0, "1 sigma") == 2.0
    assert R.expanded_uncertainty(1.0, "unstated") == 2.0
    assert R.expanded_uncertainty(None, "none") == R.UNSTATED_FLOOR == 4.0
    assert R.expanded_uncertainty(None, "unstated", "estimate") is None


def test_later_source_wins_nah_and_kh():
    # NaH: JANAF 124.26 ± 19.2 vs D0 138.59 ± 1.39 -> D0 (pilot)
    assert R.later_source_wins((124.26, 19.2), (138.59, 1.39))
    # KH: JANAF 123.01 ± 14.6 vs D0 132.54 ± 0.8 -> D0
    assert R.later_source_wins((123.01, 14.6), (132.54, 0.8))
    # precise but disagreeing: no
    assert not R.later_source_wins((124.26, 6.0), (138.59, 1.0))
    # agreeing but not precise enough: no
    assert not R.later_source_wins((124.26, 3.0), (125.0, 1.5))


def test_override_marks_the_winner_and_keeps_both(lib):
    nah = lib.add_species("[NaH]", formula="HNa")
    lib.add_measurement(
        nah,
        "JANAF-1998",
        dfh_298k=124.26,
        uncertainty=19.2,
        unc_level="2 sigma",
        citation=JANAF,
    )
    lib.add_measurement(
        nah,
        "D0",
        dfh_298k=138.59,
        dfh_0k=140.55,
        uncertainty=1.39,
        unc_level="2 sigma",
        citation="Stwalley, Zemke & Yang, JPCRD 20, 153 (1991)",
    )
    lib.choose(nah, "JANAF-1998")
    lib.choose(nah, "D0")
    assert lib.apply_overrides() == 1
    best = lib.best(nah)
    assert best["source_name"] == "D0"
    assert "below a third" in best["override_reason"]
    assert len(lib.references(nah)) == 2


def _webbook(lib, sid, value, unc, year, status="measured", cite=None):
    return lib.add_measurement(
        sid,
        "WebBook-SRD69",
        dfh_298k=value,
        uncertainty=unc,
        unc_level="unstated",
        status=status,
        year=year,
        citation=cite or f"measurement {year}",
    )


def test_choice_rule_review_first(lib):
    sid = lib.add_species("c1ccc2cc3ccccc3cc2c1", formula="C14H10")
    _webbook(lib, sid, 230.9, 3.3, 1970)
    review = _webbook(
        lib, sid, 225.7, 3.4, 2008, status="review", cite="Roux et al. 2008"
    )
    lib.choose(sid, "WebBook-SRD69")
    best = lib.best(sid)
    assert best["id"] == review and best["rule"] == "review"


def test_choice_rule_most_recent_with_uncertainty_and_spread_flag(lib):
    # a spread of 1.6x the combined uncertainties: flagged, still chosen
    sid = lib.add_species("c1c[nH]cn1", formula="C3H4N2")
    _webbook(lib, sid, 136.0, 1.9, 1987, cite="Other 1987")
    _webbook(lib, sid, 132.9, 0.6, 1987, cite="Jiménez et al. 1987")
    _webbook(lib, sid, 140.0, None, 1990, cite="no uncertainty")
    lib.choose(sid, "WebBook-SRD69")
    best = lib.best(sid)
    assert best["dfh_298k"] == 132.9  # same year: the smaller uncertainty
    assert best["rule"] == "most recent with an uncertainty"
    assert best["state"] == "chosen"
    assert "Other 1987" in best["flag"]
    assert best["spread"] == pytest.approx(7.1)


def test_large_spread_without_review_is_held_until_accepted(lib):
    sid = lib.add_species("c1c[nH]cn1", formula="C3H4N2")
    _webbook(lib, sid, 139.3, 1.9, 1987, cite="Guthrie & Pike 1987")
    m = _webbook(lib, sid, 132.9, 0.6, 1987, cite="Jiménez et al. 1987")
    lib.choose(sid, "WebBook-SRD69")
    assert lib.best(sid) is None
    held = lib.references(sid, include_held=True)
    assert held[0]["state"] == "held"
    lib.accept(m, by="science", reason="Jiménez is the more careful combustion")
    best = lib.best(sid)
    assert best["dfh_298k"] == 132.9 and best["accepted_by"] == "science"


def test_estimates_are_never_chosen_without_an_accept(lib):
    sid = lib.add_species("[BeH2]", formula="BeH2")
    m = lib.add_measurement(
        sid,
        "JANAF-1998",
        dfh_298k=125.5,
        unc_level="estimate",
        status="estimate",
        citation=JANAF,
    )
    assert lib.choose(sid, "JANAF-1998") is None
    assert lib.best(sid) is None
    lib.accept(m, by="paul", reason="test")
    assert lib.best(sid)["id"] == m
    assert lib.best(sid)["expanded_uncertainty"] is None


def _pcl5(lib, value):
    pcl3 = lib.add_species("ClP(Cl)Cl", formula="Cl3P")
    lib.add_measurement(
        pcl3, "JANAF-1998", dfh_298k=-288.70, uncertainty=5.4, citation=JANAF
    )
    lib.choose(pcl3, "JANAF-1998")
    pcl5 = lib.add_species("ClP(Cl)(Cl)(Cl)Cl", formula="Cl5P")
    m = lib.add_measurement(
        pcl5,
        "derived",
        dfh_298k=value,
        uncertainty=5.5,
        status="derived",
        citation="JANAF Cl-177 adopted reaction enthalpy with JANAF PCl3(g)",
    )
    lib.choose(pcl5, "derived")
    lib.add_derivation(
        m,
        [
            {"target": True, "coefficient": 1},
            {"smiles": "ClP(Cl)Cl", "source": "JANAF-1998", "coefficient": -1},
            {"element": True, "smiles": "ClCl", "coefficient": -1},
        ],
        -87.38,
        drh_uncertainty=1.0,
        citation="JANAF Cl-177: PCl3(g) + Cl2(g) = PCl5(g)",
    )
    return pcl5, m


def test_derivation_check_pcl5(lib):
    _pcl5(lib, -376.08)
    assert lib.check() == []


def test_derivation_check_catches_janaf_tabulated_pcl5(lib):
    _pcl5(lib, -360.18)
    problems = lib.check()
    assert len(problems) == 1
    assert problems[0]["recomputed"] == pytest.approx(-376.08)


def test_derived_value_without_derivation_is_a_problem(lib):
    sid = lib.add_species("O=P(Cl)(Cl)Cl", formula="Cl3OP")
    lib.add_measurement(
        sid, "derived", dfh_298k=-558.5, status="derived", citation=JANAF
    )
    assert lib.check()[0]["problem"] == "no derivation"


def test_offset_needs_five_shared_species(lib):
    lib.set_offset("CODATA-1989", 0.9, 2, "1.228")
    assert lib.source("CODATA-1989")["offset_kj"] is None
    lib.set_offset("WebBook-SRD69", 1.2, 12, "1.228")
    sid = lib.add_species("c1ccsc1", formula="C4H4S")
    _webbook(lib, sid, 114.9, 1.0, 1990)
    lib.choose(sid, "WebBook-SRD69")
    best = lib.best(sid)
    assert best["expanded_uncertainty"] == 2.0
    assert best["total_uncertainty"] == pytest.approx((2.0**2 + 1.2**2) ** 0.5)


def test_layers(tmp_path):
    official = R.ReferenceLibrary(
        tmp_path / "references-official.db", create=True, layer="official"
    )
    official.add_source("JANAF-1998", "evaluated", 10)
    sid = official.add_species("[LiH]", formula="HLi")
    official.add_measurement(
        sid, "JANAF-1998", dfh_298k=140.62, uncertainty=0.04, citation=JANAF
    )
    official.choose(sid, "JANAF-1998")
    official.close()

    local = R.ReferenceLibrary(tmp_path / "references.db", create=True)
    local.add_source("user", "user", 5)
    sid = local.add_species("[LiH]")
    local.add_measurement(
        sid, "user", dfh_298k=139.3, uncertainty=1.0, status="user", citation="mine"
    )
    local.choose(sid, "user")
    local.close()

    layers = R.LayeredReferences(
        tmp_path / "references.db", tmp_path / "references-official.db"
    )
    rows = layers.candidates("[LiH]")
    assert [r["layer"] for r in rows] == ["local", "official"]
    assert layers.lookup("[LiH]")["dfh_298k"] == 139.3
    assert layers.lookup("[NaH]") is None
    layers.close()

    # an official file is never used as the local layer
    with pytest.raises(ValueError):
        R.LayeredReferences(
            tmp_path / "references-official.db", tmp_path / "references-official.db"
        )


def test_cli(lib, capsys):
    _pcl5(lib, -376.08)
    lib.close()
    path = str(lib.path)
    assert references_cli.main(["check", "--library", path]) == 0
    assert references_cli.main(["sources", "--library", path]) == 0
    assert references_cli.main(["list", "--library", path]) == 0
    assert (
        references_cli.main(
            ["show", "ClP(Cl)(Cl)(Cl)Cl", "--local", path, "--official", path + ".none"]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert "Every derived value agrees" in out and "-376.08" in out
