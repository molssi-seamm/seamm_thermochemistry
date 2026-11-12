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
    assert best["dfh_298k"] == 132.9 and best["decided_by"] == "science"


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


def test_a_0k_only_value_is_used_only_without_a_298k_one(lib):
    sid = lib.add_species("c1ccsc1", formula="C4H4S")
    lib.add_measurement(
        sid,
        "W4-17",
        dfh_0k=126.28,
        uncertainty=2.7,
        unc_level="mixed",
        status="computed",
        citation="Karton et al. 2017",
    )
    lib.choose(sid, "W4-17")
    assert lib.best(sid)["source_name"] == "W4-17"
    _webbook(lib, sid, 114.9, 1.0, 1990)
    lib.choose(sid, "WebBook-SRD69")
    assert lib.best(sid)["source_name"] == "WebBook-SRD69"


def test_curated_reader(tmp_path):
    from pathlib import Path

    from seamm_thermochemistry import curated

    library = R.ReferenceLibrary(tmp_path / "r.db", create=True)
    report, problems = curated.import_curated(
        library, Path(__file__).parent / "data" / "curated"
    )
    outcome = {(f, s): o for f, s, o, _ in report}
    assert outcome[("janaf_printed.csv", "C#[O+]")] == "added"
    assert outcome[("janaf_printed.csv", "ClP(Cl)(Cl)(Cl)Cl")] == "held"
    assert outcome[("janaf_printed.csv", "FS(F)(F)(F)(F)F")] == "held"
    assert problems == []  # PCl5's derivation holds

    def best(smiles, charge=0):
        (row,) = library.species(smiles, charge)
        return library.best(row["id"])

    # JANAF's thermal-electron HCO+ in ATcT's stationary-electron convention
    assert best("C#[O+]", 1)["dfh_298k"] == pytest.approx(833.0 - R.ELECTRON_H298)
    assert best("C#[O+]", 1)["dfh_0k"] == 827.2
    assert best("[NaH]")["source_name"] == "D0"  # the override
    assert best("[LiH]")["source_name"] == "JANAF-1998"
    assert best("ClP(Cl)(Cl)(Cl)Cl")["dfh_298k"] == -376.08
    assert best("[BeH2]") is None  # an estimate
    assert best("FS(F)(F)(F)(F)F") is None  # held in every source
    assert best("[AlH3]")["dfh_298k"] is None  # W4-17, 0 K only
    assert best("[AlH3]")["unc_level"] == "mixed"
    p = library.species("[P]")[0]
    assert p["multiplicity"] == 4  # blank in CODATA: the ground state
    library.close()


def test_cli_add_curated(tmp_path, capsys):
    from pathlib import Path

    path = str(tmp_path / "official.db")
    data = str(Path(__file__).parent / "data" / "curated")
    assert references_cli.main(["add-curated", data, "--library", path]) == 0
    out = capsys.readouterr().out
    assert "held" in out and "FS(F)(F)(F)(F)F" in out
    assert R.ReferenceLibrary(path).layer == "official"


DFH = "&#916;<sub>f</sub>H&#176;<sub>gas</sub>"


def _page():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent / "data" / "webbook"))
    from make_page import page

    return page(
        "Testazole",
        "C<sub>3</sub>H<sub>4</sub>N<sub>2</sub>",
        "99999-99-9",
        [
            (
                "exp",
                DFH,
                "120.0",
                "kJ/mol",
                "N/A",
                "ref-1",
                "Derived, 1991",
                "from solid",
            ),
            (
                "exp",
                DFH,
                "131.0 &#177; 1.9",
                "kJ/mol",
                "Cm",
                "ref-2",
                "Older, 1987",
                "",
            ),
            (
                "exp",
                DFH,
                "130.0 &#177; 0.6",
                "kJ/mol",
                "Ccb",
                "ref-3",
                "Newer, 1987",
                "",
            ),
            ("cal", DFH, "128.", "kJ/mol", "Review", "ref-4", "Evaluator, 1998", ""),
            (
                "exp",
                "S&#176;<sub>gas</sub>",
                "260.",
                "J/mol*K",
                "N/A",
                "ref-1",
                "x",
                "",
            ),
        ],
        [
            ("ref-1", "Derived, 1991 Derived, D., J. Test 1, 1 (1991)."),
            ("ref-2", "Older, 1987 Older, O., J. Test 2, 2 (1987)."),
            ("ref-3", "Newer, 1987 Newer, N., J. Test 3, 3 (1987)."),
            ("ref-4", "Evaluator, 1998 Evaluator, E., Tables (1998)."),
        ],
    )


def test_webbook_parse_page():
    from seamm_thermochemistry import webbook

    page = webbook.parse_page(_page())
    assert page["found"] and page["cas"] == "99999-99-9"
    assert page["formula"] == "C3H4N2"
    assert len(page["rows"]) == 4  # the entropy row is not an enthalpy
    first = page["rows"][1]
    assert (first["value"], first["uncertainty"], first["method"]) == (131.0, 1.9, "Cm")
    assert first["year"] == 1987 and first["citation"].startswith("Older, 1987 Older")
    assert webbook.parse_value("128.") == (128.0, None)
    assert webbook.parse_value("N/A") == (None, None)


def test_webbook_import(tmp_path):
    import csv

    from seamm_thermochemistry import webbook

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "a.html").write_text(_page())
    with open(cache / webbook.INDEX, "w", newline="") as fd:
        w = csv.DictWriter(fd, webbook.INDEX_FIELDS)
        w.writeheader()
        w.writerow(dict(file="a.html", smiles="c1c[nH]cn1", charge=0, multiplicity=1))
        w.writerow(dict(file="", smiles="C1CC1", charge=0, multiplicity=1))
    library = R.ReferenceLibrary(tmp_path / "r.db", create=True)
    report = webbook.import_webbook(library, cache)
    assert ("C1CC1", "no page", "") in report
    (sp,) = library.species("c1c[nH]cn1")
    assert sp["cas"] == "99999-99-9"
    statuses = sorted(m["status"] for m in library.measurements(sp["id"]))
    assert statuses == ["held", "measured", "measured", "measured"]  # cal held
    best = library.best(sp["id"])
    # same year: the smaller uncertainty; the N/A value is a last resort
    assert best["dfh_298k"] == 130.0
    assert best["citation"].startswith("Newer")
    assert best["expanded_uncertainty"] == 1.2  # unstated level: doubled
    assert best["flag"] is None  # 1.0 apart, within hypot(0.6, 1.9)
    library.close()


def test_webbook_fetch_retries_without_stereo(tmp_path):
    import io

    from seamm_thermochemistry import webbook

    urls = []

    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    def opener(request, timeout=None):
        urls.append(request.full_url)
        if "%2Ft" in request.full_url:  # the stereo InChI: not found
            return Response(b"<html><title>Name Not Found</title></html>")
        return Response(_page().encode())

    rows = webbook.fetch([("C[C@H]1C[C@@H]1C", 0, 1)], tmp_path, delay=0, opener=opener)
    assert len(urls) == 2
    assert rows[0]["inchi_retry"] and "/t" not in rows[0]["inchi_retry"]
    assert (tmp_path / "fetch_log.csv").read_text().count("webbook.nist.gov") == 2
    # cached: nothing is fetched again
    assert webbook.fetch([("C[C@H]1C[C@@H]1C", 0, 1)], tmp_path, opener=opener) == []
    assert len(urls) == 2


def test_offsets_need_five_shared_species(lib):
    from seamm_thermochemistry.atct import MissingEntry

    atct_values = {}
    for i, smiles in enumerate(["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"]):
        sid = lib.add_species(smiles)
        lib.add_measurement(
            sid, "WebBook-SRD69", dfh_298k=-10.0 * i, uncertainty=1.0, citation="x"
        )
        lib.choose(sid, "WebBook-SRD69")
        atct_values[lib._species_row(sid)["smiles"]] = -10.0 * i - (1 if i % 2 else -1)
    sid = lib.add_species("O")
    lib.add_measurement(sid, "CODATA-1989", dfh_298k=-241.0, citation="y")
    lib.choose(sid, "CODATA-1989")
    atct_values["O"] = -241.8

    class FakeATcT:
        version = "test"

        def lookup(self, smiles, charge=0, multiplicity=None):
            if smiles not in atct_values:
                raise MissingEntry(smiles)
            return {"dfH0_298K": atct_values[smiles]}

    result = R.compute_offsets(lib, FakeATcT())
    assert result["WebBook-SRD69"][:2] == (1.0, 6)
    assert lib.source("WebBook-SRD69")["offset_kj"] == 1.0
    assert result["CODATA-1989"][1] == 1
    assert lib.source("CODATA-1989")["offset_kj"] is None  # too few
    assert lib.source("CODATA-1989")["offset_atct"] == "test"


def test_accept_file_and_0k_rule(lib, tmp_path):
    sid = lib.add_species("Cc1ccon1", formula="C4H5NO")
    _webbook(lib, sid, 35.6, 0.67, 1978, cite="Hamilton et al. 1978")
    _webbook(lib, sid, 20.0, 3.0, 1940, cite="Tappi 1940")
    lib.choose(sid, "WebBook-SRD69")
    assert lib.best(sid) is None  # held: 5x the combined uncertainty
    accepts = tmp_path / "accepts.csv"
    accepts.write_text(
        "smiles,charge,source,dfh_298k,by,reason\n"
        'Cc1ccon1,0,WebBook-SRD69,35.6,science,"modern calorimetry preferred"\n'
    )
    lib.accept_file(accepts)
    assert lib.best(sid)["dfh_298k"] == 35.6
    assert lib.best(sid)["decision_reason"] == "modern calorimetry preferred"
    bad = tmp_path / "bad.csv"
    bad.write_text(
        "smiles,charge,source,dfh_298k,by,reason\nCc1ccon1,0,WebBook-SRD69,1.0,x,y\n"
    )
    with pytest.raises(R.ReferenceError, match="not one"):
        lib.accept_file(bad)

    aid = lib.add_species("[AlH3]", formula="AlH3")
    lib.add_measurement(
        aid, "W4-17", dfh_0k=135.12, uncertainty=4.47, unc_level="mixed",
        status="computed", citation="Karton et al. 2017",
    )  # fmt: skip
    lib.choose(aid, "W4-17")
    assert "0 K only" in lib.best(aid)["rule"]


def test_hold_by_decision(lib, tmp_path):
    sid = lib.add_species("c1cnoc1", formula="C3H3NO")
    _webbook(lib, sid, 82.0, 0.6, 1996, cite="Steele et al. 1996")
    _webbook(lib, sid, 78.58, 0.54, 1978, cite="McCormick & Hamilton 1978")
    lib.choose(sid, "WebBook-SRD69")
    assert lib.best(sid) is None  # 4.2x: held by the rule already
    decisions = tmp_path / "decisions.csv"
    decisions.write_text(
        "action,smiles,charge,source,dfh_298k,by,reason\n"
        "hold,c1cnoc1,0,WebBook-SRD69,,science,"
        '"4.2x combined; awaiting tier-A arbitration"\n'
    )
    lib.accept_file(decisions)
    (held,) = lib.references(sid, include_held=True)
    assert held["state"] == "held" and held["decided_by"] == "science"
    assert "arbitration" in held["decision_reason"]
    with pytest.raises(R.ReferenceError, match="not one"):
        lib.hold("C", 0, "WebBook-SRD69", by="x", reason="y")
