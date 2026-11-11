"""The fragment store in two layers (P5c, see
docs/developer_guide/campaigns/2026-11-10, §4): the meta table, lookups local
first by rules R1-R5, the official store read-only, and the tools -- migrate,
export, merge, the official build and the cross-layer conflicts."""

import sqlite3

import pytest

from seamm_thermochemistry import (
    EnergyRecord,
    FragmentRecord,
    FragmentStore,
    GeometryRecord,
    LayeredStore,
    build_official,
    export_local,
    layer_conflicts,
    merge_stores,
    migrate_store,
)
from seamm_thermochemistry import fragment_cli
from seamm_thermochemistry.geometry_store import SCHEMA_VERSION

WATER = (["O", "H", "H"], [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]])
MOVED = (["O", "H", "H"], [[0, 0, 0.12], [0, 0.76, -0.47], [0, -0.76, -0.47]])
FP = "fp:r2scan"


def _geometry(E=-200000.0, ZPE=55.0, xyz=WATER, conformer="", fp=FP, **kw):
    return GeometryRecord.make(
        "O",
        *xyz,
        conformer=conformer,
        E=E,
        ZPE=ZPE,
        H=E + 65.0,
        frequency_fingerprint=fp,
        model_chemistry="ORCA:CMP@R2SCAN-3C",
        **kw,
    )


def _energy(E, xyz=WATER, fp="fp:revdsd"):
    return EnergyRecord.make(
        _geometry(xyz=xyz)["geometry_hash"], E, energy_fingerprint=fp
    )


@pytest.fixture()
def paths(tmp_path):
    return tmp_path / "fragments.db", tmp_path / "fragments-official.db"


def _official(path, version="2026.12.1"):
    store = FragmentStore(path, create=True, layer="official")
    store.set_meta(version=version)
    return store


# -- meta and the layers -----------------------------------------------------------


def test_meta_and_layer(paths, tmp_path):
    local, official = paths
    with FragmentStore(local, create=True) as s:
        assert s.layer == "local" and s.schema_version() == SCHEMA_VERSION == 4
        assert "created" in s.meta()
    with _official(official) as s:
        assert s.layer == "official" and s.meta()["version"] == "2026.12.1"
    # A schema-2 store (no meta table) is local, and the writer adds the table
    old = tmp_path / "old.db"
    _schema2_store(old)
    assert FragmentStore(old).layer == "local" and FragmentStore(old).meta() == {}
    with FragmentStore(old, create=True) as s:
        assert s.meta()["layer"] == "local" and s.schema_version() == 4


def test_the_official_store_is_read_only(paths):
    local, official = paths
    _official(official).close()
    with FragmentStore(official, immutable=True) as s:
        with pytest.raises(PermissionError):
            s.put_geometry(_geometry())
    with LayeredStore(local, official, create_local=True) as layered:
        layered.put_geometry(_geometry())  # goes to the local store
    assert FragmentStore(local).n_geometries() == 1
    assert FragmentStore(official).n_geometries() == 0


def test_wrong_layers_are_refused(paths):
    local, official = paths
    _official(local).close()
    with pytest.raises(ValueError, match="official store, not a local one"):
        LayeredStore(local, official)
    local.unlink()
    FragmentStore(official, create=True).close()  # a local store in its place
    with LayeredStore(local, official) as layered:
        assert layered.official is None
        assert "not an official store" in layered.notes[0]


# -- R1-R5 -------------------------------------------------------------------------


def test_r1_local_first(paths):
    local, official = paths
    with FragmentStore(local, create=True) as s:
        s.put_geometry(_geometry(E=-200000.05))
    with _official(official) as s:
        s.put_geometry(_geometry())
    with LayeredStore(local, official) as layered:
        (g,) = layered.geometries("O", frequency_fingerprint=FP)
        assert g["layer"] == "local" and layered.notes == []


def test_r2_an_official_supersede_beats_the_old_local_value(paths):
    local, official = paths
    with FragmentStore(local, create=True) as s:
        s.put_geometry(_geometry(E=-200000.0))  # the old value, harvested earlier
    with _official(official) as s:
        s.put_geometry(_geometry(E=-200000.0))
        s.supersede_geometry(_geometry(E=-200003.0), reason="lower SCF solution")
    with LayeredStore(local, official) as layered:
        (g,) = layered.geometries("O", frequency_fingerprint=FP)
        assert g["layer"] == "official" and g["E"] == -200003.0
        assert "superseded" in layered.notes[0] and "lower SCF" in layered.notes[0]
        assert "official 2026.12.1" in layered.notes[0]


def test_r3_a_differing_local_value_is_used_and_noted(paths):
    local, official = paths
    with FragmentStore(local, create=True) as s:
        s.put_geometry(_geometry(E=-200001.0))
    with _official(official) as s:
        s.put_geometry(_geometry(E=-200000.0))
    with LayeredStore(local, official) as layered:
        (g,) = layered.geometries("O", frequency_fingerprint=FP)
        assert g["layer"] == "local" and g["E"] == -200001.0
        assert "differ" in layered.notes[0] and "1.00 kJ/mol" in layered.notes[0]


def test_r4_conformers_from_both_lowest_wins_ties_local(paths):
    local, official = paths
    with FragmentStore(local, create=True) as s:
        s.put_geometry(_geometry(E=-200000.0, conformer="a"))
    with _official(official) as s:
        s.put_geometry(_geometry(E=-200002.0, xyz=MOVED, conformer="b"))
    with LayeredStore(local, official) as layered:
        rows = layered.geometries("O", frequency_fingerprint=FP)
        assert {(r["conformer"], r["layer"]) for r in rows} == {
            ("a", "local"),
            ("b", "official"),
        }
        assert LayeredStore.lowest(rows)["layer"] == "official"
        # A tie within 0.1 kJ/mol goes to the local layer
        rows[0]["E"] = -200001.95 if rows[0]["layer"] == "local" else rows[0]["E"]
        rows[1]["E"] = -200001.95 if rows[1]["layer"] == "local" else rows[1]["E"]
        assert LayeredStore.lowest(rows)["layer"] == "local"


def test_verified_before_unverified_across_layers(paths):
    local, official = paths
    with FragmentStore(local, create=True) as s:
        s.put_geometry(
            GeometryRecord.make(
                "O", *WATER, E=-1.0, ZPE=1.0, H=2.0, model_chemistry="ORCA:CMP@X"
            )
        )
    with _official(official) as s:
        s.put_geometry(_geometry())
    with LayeredStore(local, official) as layered:
        rows = layered.geometries(
            "O", frequency_fingerprint=FP, model_chemistry="CMP@X"
        )
        assert [(r["layer"], r["verified"]) for r in rows] == [("official", True)]


def test_r5_a_local_energy_on_an_official_geometry(paths):
    """A user's revDSD single point on an official r2SCAN-3c geometry: the same
    geometry hash, so it is used, and the layers of each are given."""
    local, official = paths
    with _official(official) as s:
        s.put_geometry(_geometry())
    with FragmentStore(local, create=True) as s:
        s.put_energy(_energy(-200500.0))
    with LayeredStore(local, official) as layered:
        (got,) = layered.species_at(
            "O", frequency_fingerprint=FP, energy_fingerprint="fp:revdsd"
        )
        assert got["E"] == -200500.0 and got["layers"] == ("official", "local")


def test_r5_a_superseded_geometry_does_not_pair_with_the_old_energy(paths, tmp_path):
    """Science's condition (ii): a later official version supersedes the geometry
    with new coordinates, so its hash changes, and the local energy on the old
    one is no longer found."""
    local, official = paths
    with _official(official) as s:
        s.put_geometry(_geometry())
    with FragmentStore(local, create=True) as s:
        s.put_energy(_energy(-200500.0))
    newer = tmp_path / "v2.db"
    with _official(newer, "2027.1.1") as s:
        s.put_geometry(_geometry())
        s.supersede_geometry(_geometry(xyz=MOVED, E=-200002.0), reason="new minimum")
    with LayeredStore(local, newer) as layered:
        assert (
            layered.species_at(
                "O", frequency_fingerprint=FP, energy_fingerprint="fp:revdsd"
            )
            == []
        )


def test_fragment_table_and_energy_rules(paths):
    local, official = paths
    record = FragmentRecord.make("O", E=10.0, energy_fingerprint="fp:e")
    with FragmentStore(local, create=True) as s:
        s.put(record)
        s.put_energy(_energy(-1.0))
    with _official(official) as s:
        s.put(FragmentRecord.make("O", E=10.05, energy_fingerprint="fp:e"))
        s.put(FragmentRecord.make("C", E=20.0, energy_fingerprint="fp:e"))
        s.put_energy(_energy(-1.0))
    with LayeredStore(local, official) as layered:
        assert layered.get("O", energy_fingerprint="fp:e")["layer"] == "local"
        assert layered.get("C", energy_fingerprint="fp:e")["layer"] == "official"
        h = _geometry()["geometry_hash"]
        assert layered.energy(h, energy_fingerprint="fp:revdsd")["layer"] == "local"
        assert len(layered.energies(h)) == 2 and len(layered) == 3


# -- the tools -----------------------------------------------------------------------


def _populate(local, official):
    with FragmentStore(local, create=True) as s:
        s.put_geometry(_geometry())  # also in official, same values
        s.put_geometry(_geometry(E=-200009.0, xyz=MOVED, conformer="mine"))  # new
        s.put_energy(_energy(-200500.0))  # new
        s.put(FragmentRecord.make("O", E=10.0, energy_fingerprint="fp:e"))  # differs
    with _official(official) as s:
        s.put_geometry(_geometry())
        s.put(FragmentRecord.make("O", E=12.0, energy_fingerprint="fp:e"))


def test_export_only_what_official_lacks(paths, tmp_path):
    local, official = paths
    _populate(local, official)
    out = tmp_path / "contribution.db"
    counts = export_local(out, local, official, label="Paul's Mac")
    assert counts == {"fragment": 1, "geometry": 1, "energy": 1, "differing": 1}
    with FragmentStore(out) as s:
        assert s.layer == "contribution" and s.meta()["label"] == "Paul's Mac"
        assert s.meta()["official_version"] == "2026.12.1"
        (h,) = [h for h in s.history() if h["kind"] == "differs from official"]
        assert h["old"]["E"] == 12.0 and h["new"]["E"] == 10.0
    with pytest.raises(FileExistsError):
        export_local(out, local, official)


def test_merge_round_trip(paths, tmp_path):
    local, official = paths
    _populate(local, official)
    out = tmp_path / "contribution.db"
    export_local(out, local, official)
    target = tmp_path / "other-installation.db"
    counts = merge_stores(target, [out], label="from the Mac")
    assert counts["geometry"] == {"inserted": 1} and counts["energy"] == {"inserted": 1}
    with FragmentStore(target) as s:
        assert s.n_geometries() == 1 and len(s) == 1
        assert any(k.startswith("merged ") for k in s.meta())
    # Merging again: duplicates only
    counts = merge_stores(target, [out])
    assert counts["geometry"] == {"duplicate": 1}


def test_merge_carries_superseded_values(paths, tmp_path):
    local, _ = paths
    with FragmentStore(local, create=True) as s:
        s.put_geometry(_geometry(E=-200000.0))
        s.supersede_geometry(_geometry(E=-200003.0), reason="lower conformer")
    target = tmp_path / "t.db"
    counts = merge_stores(target, [local], label="ChemAI")
    assert counts["history"] == 1
    with FragmentStore(target) as s:
        (g,) = s.geometries("O", frequency_fingerprint=FP)
        (old,) = s.superseded_values("geometry", g["id"])
        assert old["E"] == -200000.0 and "lower conformer" in old["reason"]
        assert "ChemAI" in old["reason"]


def test_merge_refuses_an_official_target(paths):
    local, official = paths
    _populate(local, official)
    with pytest.raises(ValueError, match="official store"):
        merge_stores(official, [local])


def test_build_official(paths, tmp_path):
    local, official = paths
    _populate(local, official)
    contribution = tmp_path / "c.db"
    export_local(contribution, local, official)
    # The differing fragment value stops the build, and nothing is written
    out = tmp_path / "v2.db"
    with pytest.raises(ValueError, match="differ"):
        build_official(out, "2027.1.1", [contribution], previous=official)
    assert not out.exists()
    # Resolved by superseding in the contribution's source: build from the rest
    clean = tmp_path / "clean.db"
    with FragmentStore(clean, create=True, layer="contribution") as s:
        s.put_geometry(_geometry(E=-200009.0, xyz=MOVED, conformer="mine"))
    path, digest, counts = build_official(
        out, "2027.1.1", [clean], previous=official, source="10.5281/zenodo.1"
    )
    with FragmentStore(path, immutable=True) as s:
        meta = s.meta()
        assert meta["layer"] == "official" and meta["version"] == "2027.1.1"
        assert meta["source"] == "10.5281/zenodo.1" and s.n_geometries() == 2
    sidecar = (tmp_path / "v2.db.sha256").read_text()
    assert sidecar.split()[0] == digest and sidecar.split()[1] == "v2.db"


def test_layer_conflicts(paths):
    local, official = paths
    _populate(local, official)
    (c,) = layer_conflicts(local, official)
    assert c["table"] == "fragment" and c["difference"] == pytest.approx(2.0)
    assert c["local"]["E"] == 10.0 and not c["superseded"]


def _schema2_store(path):
    """A store as P1-P4 wrote it: prefixed model chemistry names, no program in
    the fragment table, no meta table (as ChemAI's fragments.db, 2026-11-10)."""
    with FragmentStore(path, create=True) as s:
        s.put(
            FragmentRecord.make(
                "O",
                E=1.0,
                energy_fingerprint="fp:a",
                model_chemistry="ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP",
            )
        )
        s.put_geometry(_geometry())
    db = sqlite3.connect(path)
    db.execute(
        "UPDATE fragment SET program = NULL, model_chemistry = "
        "'ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP'"
    )
    db.execute("UPDATE geometry SET model_chemistry = 'ORCA:CMP@R2SCAN-3C'")
    db.execute("DROP TABLE meta")
    db.execute("PRAGMA user_version = 2")
    db.commit()
    db.close()


def test_migrate(tmp_path):
    path = tmp_path / "fragments.db"
    _schema2_store(path)
    before = (
        sqlite3.connect(path)
        .execute("SELECT energy_fingerprint FROM fragment")
        .fetchall()
    )
    result = migrate_store(path)
    assert result["backup"].exists() and result["schema"] == (2, 4)
    assert result["changed"] == {"fragment": 1, "geometry": 1}
    db = sqlite3.connect(path)
    assert db.execute("SELECT model_chemistry, program FROM fragment").fetchall() == [
        ("DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP", "ORCA")
    ]
    assert db.execute("SELECT model_chemistry, program FROM geometry").fetchall() == [
        ("CMP@R2SCAN-3C", "ORCA")
    ]
    assert db.execute("SELECT energy_fingerprint FROM fragment").fetchall() == before
    db.close()
    assert FragmentStore(path).layer == "local"
    # Again: nothing to do, no backup
    again = migrate_store(path)
    assert again["backup"] is None and again["changed"] == {}


# -- the command line --------------------------------------------------------------


def test_cli(paths, tmp_path, capsys):
    local, official = paths
    _populate(local, official)
    base = ["--store", str(local), "--official", str(official)]

    def run(*args):
        assert fragment_cli.main(base + list(args)) == 0
        return capsys.readouterr().out

    out = run("geometries", "O", "--fingerprint", FP)
    assert "[local]" in out
    assert "local value is used" in run("duplicates", "--layers")
    assert "up to date" in run("migrate")
    contribution = tmp_path / "c.db"
    assert "1 differing" in run("export", "--out", str(contribution))
    target = tmp_path / "t.db"
    out = fragment_cli.main(["--store", str(target), "merge", str(contribution)])
    assert out == 0 and "inserted" in capsys.readouterr().out
    clean = tmp_path / "clean.db"
    FragmentStore(clean, create=True, layer="contribution").close()
    out = run(
        "official",
        "build",
        "--version",
        "2027.1.1",
        "--out",
        str(tmp_path / "v2.db"),
        "--previous",
        str(official),
        str(clean),
    )
    assert "version 2027.1.1" in out and "sha256" in out
    # Errors are reported, not raised
    assert fragment_cli.main(base + ["export", "--out", str(contribution)]) == 1
