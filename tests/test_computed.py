"""The computed layer of the reference library (Paul, 2026-11-15): values
computed with SEAMM's tier-A recipes, in their own file, used only after every
experimental value. The values here are made up."""

import csv
import json

import pytest

pytest.importorskip("rdkit")

from seamm_thermochemistry import computed as C  # noqa: E402
from seamm_thermochemistry import references as R  # noqa: E402
from seamm_thermochemistry import references_cli  # noqa: E402
from seamm_thermochemistry import resolver as RS  # noqa: E402
from seamm_thermochemistry.atct import MissingEntry  # noqa: E402

FIELDS = (
    "smiles",
    "formula",
    "charge",
    "multiplicity",
    "name",
    "phase",
    "dfH0_kJmol",
    "dfH298_kJmol",
    "unc_kJmol",
    "unc_level",
    "basis",
    "status",
    "recipe",
    "job",
    "iteration",
    "model_chemistry",
    "geometry",
    "conformer",
    "H298_H0_kJmol",
    "thermal",
    "atomization",
    "atom_references",
    "computed_on",
    "citation",
    "notes",
)


def _row(smiles, value, job="5738", iteration="1", **kw):
    row = dict(
        smiles=smiles,
        charge=0,
        multiplicity=1,
        phase="g",
        dfH0_kJmol=value + 10.0,
        dfH298_kJmol=value,
        unc_kJmol=2.5,
        unc_level="95 %",
        basis="recipe estimate",
        status="computed (tier-A)",
        recipe="TEST:CMP@F12-recipe",
        job=job,
        iteration=iteration,
        geometry="r2SCAN-3c",
        conformer="lowest of 3",
        H298_H0_kJmol=15.0,
        thermal="harmonic, ZPE x 0.99",
        atomization="atomization, ATcT atoms",
        atom_references="ATcT test",
        computed_on="2026-11-15",
    )
    row.update(kw)
    return row


def _csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fd:
        w = csv.DictWriter(fd, FIELDS)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    return path


@pytest.fixture
def computed(tmp_path):
    library = R.ReferenceLibrary(
        tmp_path / "references-computed.db", create=True, layer=C.LAYER
    )
    yield library
    library.close()


def test_import_and_provenance(computed, tmp_path):
    report = C.import_computed(
        computed, _csv(tmp_path / "v.csv", [_row("CC1CCC1", 2.0)])
    )
    assert report == [("CC1CCC1", "added", "SEAMM job 5738, iteration 1")]
    (sp,) = computed.species("CC1CCC1")
    best = computed.best(sp["id"])
    assert best["dfh_298k"] == 2.0 and best["dfh_0k"] == 12.0
    assert best["source_name"] == C.SOURCE and best["source_kind"] == "computed"
    assert best["status"] == "computed" and best["expanded_uncertainty"] == 2.5
    assert best["method"] == "TEST:CMP@F12-recipe" and best["h298_h0"] == 15.0
    notes = json.loads(best["notes"])
    assert notes["conformer"] == "lowest of 3" and notes["job"] == "5738"


def test_reimport_is_idempotent_and_a_later_job_replaces(computed, tmp_path):
    path = _csv(tmp_path / "v.csv", [_row("CC1CCC1", 2.0)])
    C.import_computed(computed, path)
    assert C.import_computed(computed, path)[0][1] == "present"
    path = _csv(tmp_path / "w.csv", [_row("CC1CCC1", 3.0, job="5800")])
    assert C.import_computed(computed, path)[0][1] == "replaced"
    (sp,) = computed.species("CC1CCC1")
    assert computed.best(sp["id"])["dfh_298k"] == 3.0
    old = [m for m in computed.measurements(sp["id"]) if m["dfh_298k"] == 2.0]
    assert old[0]["status"] == "held"
    assert json.loads(old[0]["notes"])["replaced by"] == "SEAMM job 5800, iteration 1"


@pytest.mark.parametrize(
    "change, message",
    [
        ({"job": ""}, "no job"),
        ({"status": "computed"}, "status must be"),
        ({"unc_level": "1 sigma"}, "95 %"),
        ({"job": "/D3/someone/job5738"}, "SEAMM job id"),
        ({"notes": "see /Users/someone/run"}, "local path"),
        ({"geometry": "~/runs/geom.xyz"}, "local path"),
        ({"computed_on": "15 Nov"}, "YYYY-MM-DD"),
        ({"dfH0_kJmol": "", "dfH298_kJmol": ""}, "no enthalpy"),
    ],
)
def test_rows_are_checked(computed, tmp_path, change, message):
    report = C.import_computed(
        computed, _csv(tmp_path / "v.csv", [_row("CC1CCC1", 2.0, **change)])
    )
    assert report[0][1] == "refused" and message in report[0][2]
    assert computed.species("CC1CCC1") == []


def test_only_the_computed_layer_takes_computed_values(tmp_path):
    lib = R.ReferenceLibrary(tmp_path / "references.db", create=True)
    with pytest.raises(R.ReferenceError, match="computed"):
        C.import_computed(lib, _csv(tmp_path / "v.csv", [_row("CC1CCC1", 2.0)]))
    lib.close()


def _official(tmp_path):
    lib = R.ReferenceLibrary(
        tmp_path / "references-official.db", create=True, layer="official"
    )
    lib.add_source("WebBook-SRD69", "measured", 50)
    sid = lib.add_species("CC(F)(F)Cl")
    lib.add_measurement(
        sid, "WebBook-SRD69", dfh_298k=-536.2, uncertainty=2.0, citation="exp"
    )
    lib.choose(sid, "WebBook-SRD69")
    lib.close()


def _layers(tmp_path, rank=None):
    _official(tmp_path)
    lib = R.ReferenceLibrary(
        tmp_path / "references-computed.db", create=True, layer=C.LAYER
    )
    if rank is not None:  # a computed source ranked ahead of the experiments
        lib.add_source(C.SOURCE, "computed", rank)
    C.import_computed(
        lib,
        _csv(
            tmp_path / "v.csv",
            [_row("CC(F)(F)Cl", -530.0), _row("CC1CCC1", 2.0, iteration="2")],
        ),
    )
    lib.close()
    return R.LayeredReferences(
        tmp_path / "references.db", tmp_path / "references-official.db"
    )


@pytest.mark.parametrize("rank", [None, 1])
def test_experimental_values_always_win(tmp_path, rank):
    layers = _layers(tmp_path, rank)
    # found beside the official library, without being named
    assert [name for name, _ in layers.layers] == ["official", "computed"]
    rows = layers.candidates("CC(F)(F)Cl")
    assert [r["layer"] for r in rows] == ["official", "computed"]
    assert layers.lookup("CC(F)(F)Cl")["dfh_298k"] == -536.2
    only = layers.lookup("CC1CCC1")
    assert only["layer"] == "computed" and only["dfh_298k"] == 2.0
    layers.close()


def test_the_computed_layer_gives_only_computed_values(tmp_path):
    _official(tmp_path)
    lib = R.ReferenceLibrary(
        tmp_path / "references-computed.db", create=True, layer=C.LAYER
    )
    lib.add_source("WebBook-SRD69", "measured", 50)
    sid = lib.add_species("CC1CCC1")
    lib.add_measurement(sid, "WebBook-SRD69", dfh_298k=1.0, citation="stray")
    lib.choose(sid, "WebBook-SRD69")
    lib.close()
    layers = R.LayeredReferences(None, tmp_path / "references-official.db")
    assert layers.lookup("CC1CCC1") is None
    layers.close()


def test_a_wrong_file_is_not_read_as_the_computed_layer(tmp_path):
    _official(tmp_path)
    R.ReferenceLibrary(tmp_path / "references-computed.db", create=True).close()
    with pytest.raises(ValueError, match="not the computed one"):
        R.LayeredReferences(None, tmp_path / "references-official.db")


class FakeATcT:
    version = "test"
    values = {"C": -74.5, "CC": -84.0, "[C]": 716.9, "[H]": 218.0, "[F]": 79.4}

    def lookup(self, smiles, charge=0, multiplicity=None, **kw):
        from seamm_thermochemistry.atct import canonical_smiles

        for key, value in self.values.items():
            if canonical_smiles(key) == canonical_smiles(smiles):
                return {
                    "id": f"atct-{key}",
                    "dfH0_0K": value,
                    "dfH0_298K": value,
                    "uncertainty_298K": 0.1,
                    "exact": False,
                }
        raise MissingEntry(smiles)


def test_resolver(tmp_path):
    layers = _layers(tmp_path)
    r = RS.ReferenceResolver(FakeATcT(), supplementary="compilations", layers=layers)
    entries, missing, label, _ = r.resolve(["C", "CC(F)(F)Cl", "CC1CCC1"])
    assert entries["CC(F)(F)Cl"]["source"] == "WebBook-SRD69"
    assert "CC1CCC1" in missing and label == "mixed"

    r = RS.ReferenceResolver(FakeATcT(), layers=layers)  # and computed
    entries, missing, label, _ = r.resolve(["C", "CC(F)(F)Cl", "CC1CCC1"])
    assert entries["CC(F)(F)Cl"]["source"] == "WebBook-SRD69"  # experiment wins
    assert entries["CC1CCC1"]["source"] == C.SOURCE
    assert entries["CC1CCC1"]["source_kind"] == "computed"
    assert entries["CC1CCC1"]["layer"] == "computed"
    assert not missing and label == "computed references"
    layers.close()


def test_cli(tmp_path, capsys):
    _official(tmp_path)
    path = _csv(tmp_path / "v.csv", [_row("CC1CCC1", 2.0), _row("C1CC1", 1.0, job="")])
    db = tmp_path / "references-computed.db"
    assert references_cli.main(["add-computed", str(path), "--library", str(db)]) == 1
    out = capsys.readouterr().out
    assert "refused  C1CC1  no job" in out and "1 added, 1 refused" in out
    code = references_cli.main(
        [
            "show",
            "CC1CCC1",
            "--local",
            str(tmp_path / "none.db"),
            "--official",
            str(tmp_path / "references-official.db"),
        ]
    )
    out = capsys.readouterr().out
    assert code == 0 and "computed" in out and "SEAMM job 5738, iteration 1" in out
    assert "recipe TEST:CMP@F12-recipe" in out
