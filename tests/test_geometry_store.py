"""The geometry and energy tables of the fragment store: keys by geometry hash,
the pairing and same-level guards, imaginary modes, duplicates and supersede,
combining an energy level with a frequency level, old stores, and harvest."""

import json
import sqlite3

import pytest

from seamm_thermochemistry import (
    EnergyRecord,
    FragmentStore,
    GeometryRecord,
    PairingError,
    harvest_job,
)
from seamm_thermochemistry import geometry_store as G

WATER = (["O", "H", "H"], [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]])
R2SCAN = "ORCA:DFT@R2SCAN-3C"
REVDSD = "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"


@pytest.fixture()
def store(tmp_path):
    with FragmentStore(tmp_path / "fragments.db", create=True) as s:
        yield s


def _geometry(conformer="", frequencies=(1600.0, 3650.0, 3750.0), **kw):
    kw.setdefault("frequency_fingerprint", "orca:r2scan")
    kw.setdefault("model_chemistry", R2SCAN)
    return GeometryRecord.make(
        "O",
        *WATER,
        conformer=conformer,
        E=-200000.0,
        ZPE=55.0,
        H=-199935.0,
        frequencies=frequencies,
        symmetry_number=2,
        moments=[0.6, 1.2, 1.8],
        **kw,
    )


def test_schema_version_and_tables(store):
    assert store.schema_version() == G.SCHEMA_VERSION
    tables = {
        r[0]
        for r in store._connect().execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }
    assert {"fragment", "geometry", "energy", "row_history"} <= tables


def test_old_store_is_migrated_by_the_writer(tmp_path):
    path = tmp_path / "fragments.db"
    db = sqlite3.connect(path)
    db.executescript(
        "CREATE TABLE fragment (id INTEGER PRIMARY KEY, species TEXT NOT NULL, "
        "charge INTEGER NOT NULL DEFAULT 0, multiplicity INTEGER NOT NULL DEFAULT 1, "
        "energy_fingerprint TEXT NOT NULL, frequency_fingerprint TEXT NOT NULL "
        "DEFAULT '', verified INTEGER NOT NULL DEFAULT 1, model_chemistry TEXT, "
        "frequency_model_chemistry TEXT, program TEXT, version TEXT, E REAL, ZPE "
        "REAL, H REAL, T REAL, geometry_hash TEXT, job TEXT, configuration TEXT, "
        "source TEXT, date TEXT, UNIQUE (species, charge, multiplicity, "
        "energy_fingerprint, frequency_fingerprint));"
    )
    db.close()
    # A reader of the old store sees no geometries, without an error
    with FragmentStore(path) as reader:
        assert reader.schema_version() == 0
        assert reader.geometries("O", model_chemistry=R2SCAN) == []
        assert reader.energy("sha256:x", model_chemistry=REVDSD) is None
    # The writer adds the tables
    with FragmentStore(path, create=True) as writer:
        assert writer.schema_version() == G.SCHEMA_VERSION
        assert writer.put_geometry(_geometry())[0] == "inserted"


def test_geometry_hash_is_the_qc_steps_fingerprint():
    from seamm_util.fingerprint import geometry_fingerprint

    g = _geometry()
    assert g["geometry_hash"] == geometry_fingerprint([8, 1, 1], WATER[1])
    assert g["geometry_hash"].startswith("sha256:")
    # Ghost-atom style symbols ("H:") map to the element
    assert G.hash_geometry(["O", "H:", "h"], WATER[1]) == g["geometry_hash"]


def test_geometry_roundtrip_and_conformers(store):
    assert store.put_geometry(_geometry())[0] == "inserted"
    assert store.put_geometry(_geometry(conformer="b"))[0] == "inserted"
    got = store.geometries("O", frequency_fingerprint="orca:r2scan")
    assert [g["conformer"] for g in got] == ["", "b"]
    g = got[0]
    assert g["verified"] and g["coordinates"] == [[0, 0, 0.117], *WATER[1][1:]]
    assert g["frequencies"] == [1600.0, 3650.0, 3750.0]
    assert g["n_imaginary"] == 0 and g["symmetry_number"] == 2
    assert g["moments"] == [0.6, 1.2, 1.8] and g["T"] == 298.15
    # Other charge, multiplicity or level: nothing
    assert store.geometries("O", charge=1, frequency_fingerprint="orca:r2scan") == []
    assert store.geometries("O", frequency_fingerprint="orca:pbe0") == []


def test_unverified_geometry_fallback(store):
    store.put_geometry(_geometry(frequency_fingerprint=None, version="6.1.1"))
    got = store.geometries("O", frequency_fingerprint="orca:r2scan")
    assert got == []  # a fingerprint asks for verified rows
    got = store.geometries(
        "O", frequency_fingerprint="orca:r2scan", model_chemistry=R2SCAN
    )
    assert len(got) == 1 and not got[0]["verified"]
    assert got[0]["frequency_fingerprint"] == f"unverified:{R2SCAN}@6.1.1"
    with pytest.raises(ValueError):
        GeometryRecord.make("O", *WATER)


def test_imaginary_modes(store):
    noise = _geometry(frequencies=(-15.0, 1600.0, 3650.0))
    assert noise["n_imaginary"] == 0 and noise["max_imaginary"] == 15.0
    saddle = _geometry(conformer="ts", frequencies=(-250.0, 1600.0, 3650.0))
    assert saddle["n_imaginary"] == 1
    store.put_geometry(noise)
    store.put_geometry(saddle)
    got = store.geometries("O", frequency_fingerprint="orca:r2scan")
    assert [g["conformer"] for g in got] == [""]
    got = store.geometries(
        "O", frequency_fingerprint="orca:r2scan", include_imaginary=True
    )
    assert [g["conformer"] for g in got] == ["", "ts"]


def test_one_level_for_geometry_and_frequencies():
    with pytest.raises(PairingError, match="optimized at"):
        _geometry(optimization_model_chemistry="ORCA:DFT@PBE0-D4/def2-TZVP")
    # The same level (the fingerprints may differ, e.g. a tighter SCF)
    g = _geometry(
        optimization_model_chemistry=R2SCAN, optimization_fingerprint="orca:other"
    )
    assert g["optimization_fingerprint"] == "orca:other"


def test_energy_rows_and_the_pairing_guard(store):
    _, gid = store.put_geometry(_geometry())
    h = _geometry()["geometry_hash"]
    e = EnergyRecord.make(
        h,
        -200500.0,
        energy_fingerprint="orca:revdsd",
        model_chemistry=REVDSD,
        calculation_geometry=h,
        byproducts={"MP2": -12.5},
    )
    assert store.put_energy(e) == ("inserted", 1)
    got = store.energy(h, energy_fingerprint="orca:revdsd")
    assert got["E"] == -200500.0 and got["geometry"] == gid
    assert got["byproducts"] == {"MP2": -12.5}
    with pytest.raises(PairingError, match="different geometry"):
        EnergyRecord.make(
            h, -1.0, energy_fingerprint="orca:revdsd", calculation_geometry="sha256:x"
        )
    # Unverified fallback by model chemistry
    store.put_energy(EnergyRecord.make(h, -200400.0, model_chemistry="ORCA:HF"))
    assert store.energy(h, model_chemistry="ORCA:HF")["E"] == -200400.0
    assert store.energy(h, energy_fingerprint="orca:nothing") is None
    assert len(store.energies(h)) == 2


def test_duplicates_are_logged_not_overwritten(store, caplog):
    h = _geometry()["geometry_hash"]
    store.put_energy(EnergyRecord.make(h, -100.0, energy_fingerprint="fp", job="1"))
    assert (
        store.put_energy(EnergyRecord.make(h, -100.05, energy_fingerprint="fp"))[0]
        == "duplicate"
    )
    assert (
        store.put_energy(
            EnergyRecord.make(h, -101.0, energy_fingerprint="fp", job="2")
        )[0]
        == "logged"
    )
    assert store.energy(h, energy_fingerprint="fp")["E"] == -100.0
    (entry,) = store.history("energy")
    assert entry["kind"] == "duplicate" and entry["difference"] == pytest.approx(1.0)
    assert entry["new"]["job"] == "2"
    assert "differing by 1.000" in caplog.text
    # The same for geometry rows
    store.put_geometry(_geometry())
    g = _geometry()
    g["ZPE"] = 57.0
    assert store.put_geometry(g)[0] == "logged"


def test_supersede(store):
    h = _geometry()["geometry_hash"]
    store.put_energy(EnergyRecord.make(h, -100.0, energy_fingerprint="fp"))
    with pytest.raises(ValueError):
        store.supersede_energy(
            EnergyRecord.make(h, -101.0, energy_fingerprint="fp"), ""
        )
    outcome, _ = store.supersede_energy(
        EnergyRecord.make(h, -101.0, energy_fingerprint="fp"), "wrong SCF solution"
    )
    assert outcome == "superseded"
    assert store.energy(h, energy_fingerprint="fp")["E"] == -101.0
    (entry,) = store.history("energy")
    assert entry["kind"] == "superseded" and entry["old"]["E"] == -100.0
    assert entry["reason"] == "wrong SCF solution"
    # Within the tolerance: unchanged unless forced
    small = EnergyRecord.make(h, -101.05, energy_fingerprint="fp")
    assert store.supersede_energy(small, "x")[0] == "unchanged"
    assert store.supersede_energy(small, "x", force=True)[0] == "superseded"
    # A geometry: a value the new row lacks keeps the stored one
    store.put_geometry(_geometry())
    g = _geometry()
    g["ZPE"], g["H"] = None, -199930.0
    assert store.supersede_geometry(g, "better thermal term")[0] == "superseded"
    (row,) = store.geometries("O", frequency_fingerprint="orca:r2scan")
    assert row["ZPE"] == 55.0 and row["H"] == -199930.0


def test_species_at_an_energy_level(store):
    store.put_geometry(_geometry(optimization_fingerprint="orca:r2scan"))
    h = _geometry()["geometry_hash"]
    store.put_energy(EnergyRecord.make(h, -200500.0, energy_fingerprint="orca:revdsd"))
    (got,) = store.species_at(
        "O", frequency_fingerprint="orca:r2scan", energy_fingerprint="orca:revdsd"
    )
    assert got["E"] == -200500.0 and got["ZPE"] == 55.0
    assert got["H - E"] == pytest.approx(65.0)
    # The frequency level itself, from the geometry row
    (same,) = store.species_at(
        "O", frequency_fingerprint="orca:r2scan", energy_fingerprint="orca:r2scan"
    )
    assert same["E"] == -200000.0 and same["energy"] is None
    # A level with no energy on this geometry: nothing
    assert (
        store.species_at(
            "O", frequency_fingerprint="orca:r2scan", energy_fingerprint="orca:f12"
        )
        == []
    )


def test_readonly_cannot_write(store, tmp_path):
    store.put_geometry(_geometry())
    with FragmentStore(tmp_path / "fragments.db") as reader:
        with pytest.raises(PermissionError):
            reader.put_geometry(_geometry(conformer="x"))
        with pytest.raises(PermissionError):
            reader.put_energy(
                EnergyRecord.make("sha256:x", 1.0, energy_fingerprint="f")
            )


# -- harvesting --------------------------------------------------------------


@pytest.fixture()
def job(tmp_path):
    """A job with a species at r2SCAN-3c (frequency level, with frequencies,
    symmetry number and moments) and a revDSD single point on the same
    geometry, plus one whose single point ran on another geometry."""
    molsystem = pytest.importorskip("molsystem")
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    db = molsystem.SystemDB(filename=f"file:{job_dir / 'seamm.db'}")
    for smiles in ("O", "C"):
        system = db.create_system(name=f"CBH:{smiles}", make_current=True)
        c = system.create_configuration(name="initial", make_current=True)
        c.from_smiles(smiles)
        p = c.properties
        ghash = G.hash_geometry(
            c.atoms.symbols, c.atoms.get_coordinates(fractionals=False)
        )

        def put(name, kind, value, units=None):
            p.add(name, kind, units=units, noerror=True)
            p.put(name, value)

        s = f"#ORCA#{R2SCAN}"
        put("total energy" + s, "float", -76.0, "E_h")
        put("ZPE" + s, "float", 55.0, "kJ/mol")
        put("H(298.15 K)" + s, "float", -76.0 * 2625.4996394799 + 65.0, "kJ/mol")
        put("settings fingerprint" + s, "str", "orca:r2scan")
        put("frequency fingerprint" + s, "str", "orca:r2scan")
        put("frequencies" + s, "json", [1600.0, 3650.0, 3750.0])
        put("symmetry number" + s, "int", 2)
        put("moments of inertia" + s, "json", [0.6, 1.2, 1.8])
        put("calculation" + s, "json", {"program": "ORCA", "version": "6.1.1"})
        s = f"#ORCA#{REVDSD}"
        put("total energy" + s, "float", -76.3, "E_h")
        put("settings fingerprint" + s, "str", "orca:revdsd")
        put(
            "calculation" + s,
            "json",
            {
                "program": "ORCA",
                "version": "6.1.1",
                "geometry fingerprint": ghash if smiles == "O" else "sha256:other",
            },
        )
    db.close()
    (job_dir / "job_data.json").write_text(json.dumps({"job id": 7}))
    return job_dir


def test_harvest_geometry_and_energy(job, store):
    counts = harvest_job(job, store, all_systems=False)
    assert counts["geometry"] == {"inserted": 2}
    assert counts["energy"] == {"inserted": 3}  # methane's revDSD refused
    (refused,) = counts["skipped"]
    assert refused.startswith("C (total energy#ORCA#") and "different geometry" in (
        refused
    )
    (g,) = store.geometries("O", frequency_fingerprint="orca:r2scan")
    assert g["frequencies"] == [1600.0, 3650.0, 3750.0]
    assert g["symmetry_number"] == 2 and g["moments"] == [0.6, 1.2, 1.8]
    assert g["program"] == "ORCA" and g["version"] == "6.1.1"
    assert g["optimization_fingerprint"] == "orca:r2scan"
    assert g["job"].startswith("7:")
    (got,) = store.species_at(
        "O", frequency_fingerprint="orca:r2scan", energy_fingerprint="orca:revdsd"
    )
    assert got["E"] == pytest.approx(-76.3 * 2625.4996394799)
    assert got["H - E"] == pytest.approx(65.0)
    # Harvesting again changes nothing
    again = harvest_job(job, store)
    assert again["geometry"] == {"duplicate": 2}
    assert again["energy"] == {"duplicate": 3}


def test_cli_geometries(job, tmp_path, capsys):
    from seamm_thermochemistry import fragment_cli

    path = str(tmp_path / "cli.db")
    assert fragment_cli.main(["--store", path, "harvest", str(job)]) == 0
    out = capsys.readouterr().out
    assert "geometry rows: 2 inserted" in out and "energy rows: 3 inserted" in out
    assert "different geometry" in out
    assert (
        fragment_cli.main(
            ["--store", path, "geometries", "O", "--fingerprint", "orca:r2scan"]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert f"O at {R2SCAN}" in out and REVDSD in out
    assert fragment_cli.main(["--store", path, "geometries", "CC"]) == 1


def test_harvest_with_a_species_map(tmp_path, store):
    """Systems named by id are found through a species map, which gives the
    species (with stereo), the configuration, the conformer label and the
    charge and multiplicity; 'geometry' fills only the geometry tables, so two
    conformers of one species go in side by side."""
    molsystem = pytest.importorskip("molsystem")
    from seamm_thermochemistry import load_species_map

    job_dir = tmp_path / "job"
    job_dir.mkdir()
    db = molsystem.SystemDB(filename=f"file:{job_dir / 'seamm.db'}")
    s = f"#ORCA#{R2SCAN}"
    for name, z in (("74-84-0_0", 0.0), ("74-84-0_1", 0.3), ("O_unlisted", 0.0)):
        system = db.create_system(name=name, make_current=True)
        first = system.create_configuration(name="initial", make_current=True)
        first.from_smiles("CC")
        final = system.create_configuration(
            name=f"{name} optimized", make_current=False
        )
        final.from_smiles("CC")
        xyz = final.atoms.get_coordinates(fractionals=False)
        xyz[0][2] += z
        final.atoms.set_coordinates(xyz, fractionals=False)
        for n, kind, value, units in (
            ("total energy" + s, "float", -79.0 - z, "E_h"),
            ("ZPE" + s, "float", 190.0, "kJ/mol"),
            (
                "H(298.15 K)" + s,
                "float",
                (-79.0 - z) * 2625.4996394799 + 200.0,
                "kJ/mol",
            ),
            ("settings fingerprint" + s, "str", "orca:r2scan", None),
            ("frequency fingerprint" + s, "str", "orca:r2scan", None),
        ):
            final.properties.add(n, kind, units=units, noerror=True)
            final.properties.put(n, value)
    db.close()
    (job_dir / "job_data.json").write_text(json.dumps({"job id": 5602}))
    path = tmp_path / "map.csv"
    path.write_text(
        "job,system name,configuration name,species,charge,multiplicity,"
        "conformer label\n"
        "5602,74-84-0_0,74-84-0_0 optimized,CC,0,1,74-84-0_0\n"
        "5602,74-84-0_1,74-84-0_1 optimized,CC,0,1,74-84-0_1\n"
        "5999,O_unlisted,x,O,0,1,\n"  # another job's: not used here
    )
    species_map = load_species_map(path)
    counts = harvest_job(job_dir, store, species_map=species_map, tables="geometry")
    assert counts["inserted"] == 0  # no fragment records
    assert counts["geometry"] == {"inserted": 2}
    assert counts["skipped"] == ["O_unlisted: not in the species map"]
    rows = store.geometries("CC", frequency_fingerprint="orca:r2scan")
    assert [r["conformer"] for r in rows] == ["74-84-0_0", "74-84-0_1"]
    # The named configuration (the optimized one), not the current one
    assert rows[1]["coordinates"][0][2] != rows[0]["coordinates"][0][2]
    assert len(store) == 0
    # A duplicate key in the map is an error
    path.write_text("job,system name,species\n5602,a,CC\n5602,a,CCC\n")
    with pytest.raises(ValueError, match="listed twice"):
        load_species_map(path)
