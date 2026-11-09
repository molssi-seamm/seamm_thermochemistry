# -*- coding: utf-8 -*-

"""The fragment-energy store: energies of CBH fragments (and targets) computed
earlier, so that a species is computed once per level of theory.

The store is one SQLite file, ``fragments.db`` next to the ATcT cache
(``Parameters/thermochemistry/``), and must live on a **local disk**: SQLite
locking is unreliable on GPFS/NFS, so the connection uses no memory mapping,
a busy timeout and the rollback journal, and nothing on a compute node ever
writes to it. Jobs store their species' energies as properties on their own
configurations (the systems named ``CBH:<SMILES>``); afterwards the single
writer, ``seamm-thermochemistry-fragments harvest <job>`` on the host,
promotes them with insert-or-ignore and full provenance.

A record is keyed by the species identity (canonical SMILES, charge and spin
multiplicity) and two opaque *fingerprints*: the one of the calculation that
gave the energy and the one of the calculation that gave the frequencies
(ZPE and thermal enthalpy), which can be at different levels. The QC step
owns the fingerprint (a hash of the generated input without the lines that do
not affect the energy, plus program and version); this module never
interprets it. The model chemistry string is kept for people to read. A
record without a fingerprint is keyed on the model chemistry and code
version instead and flagged *unverified*: it is used only when nothing
verified exists, and the caller is told.

E, ZPE and H(T) are separate fields in kJ/mol. A record whose key already
exists is ignored (first value wins); if its values differ from the stored
ones by more than 0.1 kJ/mol the attempt is logged in ``duplicate_log``, never
overwritten, which catches a different conformer or SCF solution. A value
known to be wrong (a wrong conformer, say) is replaced only explicitly, with
``supersede``: the old values move to the ``superseded`` table with the date
and the reason, and the record keeps its key and id, so nothing is lost and
the replacement is visible.

Reads are best effort: a missing, locked or corrupt store gives None and a
log message, and the caller computes the species.
"""

import datetime
import hashlib
import json
import logging
from pathlib import Path
import sqlite3

from .geometry_store import (
    GEOMETRY_SCHEMA,
    SCHEMA_VERSION,
    EnergyRecord,
    GeometryEnergyMixin,
    GeometryRecord,
    PairingError,
    hash_geometry,
)

logger = logging.getLogger(__name__)

__all__ = [
    "FragmentStore",
    "FragmentRecord",
    "GeometryRecord",
    "EnergyRecord",
    "PairingError",
    "fragment_store_path",
    "harvest_job",
    "load_species_map",
]

DUPLICATE_TOLERANCE = 0.1  # kJ/mol
UNVERIFIED = "unverified"

# Energy units to kJ/mol
_TO_KJ = {
    "kJ/mol": 1.0,
    "kj/mol": 1.0,
    "kcal/mol": 4.184,
    "E_h": 2625.4996394799,
    "hartree": 2625.4996394799,
    "Hartree": 2625.4996394799,
    "eV": 96.48533212331,
    None: 1.0,
    "": 1.0,
}

_SCHEMA = """
CREATE TABLE IF NOT EXISTS fragment (
    id INTEGER PRIMARY KEY,
    species TEXT NOT NULL,
    charge INTEGER NOT NULL DEFAULT 0,
    multiplicity INTEGER NOT NULL DEFAULT 1,
    energy_fingerprint TEXT NOT NULL,
    frequency_fingerprint TEXT NOT NULL DEFAULT '',
    verified INTEGER NOT NULL DEFAULT 1,
    model_chemistry TEXT,
    frequency_model_chemistry TEXT,
    program TEXT,
    version TEXT,
    E REAL,
    ZPE REAL,
    H REAL,
    T REAL,
    geometry_hash TEXT,
    job TEXT,
    configuration TEXT,
    source TEXT,
    date TEXT,
    UNIQUE (species, charge, multiplicity, energy_fingerprint, frequency_fingerprint)
);
CREATE INDEX IF NOT EXISTS fragment_species ON fragment (species, charge, multiplicity);
CREATE TABLE IF NOT EXISTS superseded (
    id INTEGER PRIMARY KEY,
    fragment INTEGER NOT NULL REFERENCES fragment (id),
    verified INTEGER,
    model_chemistry TEXT,
    frequency_model_chemistry TEXT,
    program TEXT,
    version TEXT,
    E REAL,
    ZPE REAL,
    H REAL,
    T REAL,
    geometry_hash TEXT,
    job TEXT,
    configuration TEXT,
    source TEXT,
    date TEXT,
    superseded_on TEXT,
    reason TEXT
);
CREATE TABLE IF NOT EXISTS duplicate_log (
    id INTEGER PRIMARY KEY,
    fragment INTEGER NOT NULL REFERENCES fragment (id),
    E REAL,
    ZPE REAL,
    H REAL,
    geometry_hash TEXT,
    job TEXT,
    configuration TEXT,
    source TEXT,
    date TEXT,
    difference REAL
);
"""


def fragment_store_path(path=None):
    """The store file: `path` if given, else ``fragments.db`` next to the ATcT
    cache (``<...>/Parameters/thermochemistry/fragments.db``)."""
    if path is not None:
        return Path(path).expanduser()
    from .atct import atct_directory

    return atct_directory().parent / "fragments.db"


def geometry_hash(symbols, coordinates, digits=4):
    """A hash of the geometry (symbols and coordinates rounded to `digits` Å),
    to tell apart records of the same species at different geometries."""
    text = ";".join(
        f"{s}:" + ",".join(f"{float(x):.{digits}f}" for x in xyz)
        for s, xyz in zip(symbols, coordinates)
    )
    return hashlib.sha256(text.encode()).hexdigest()[:16]


class FragmentRecord(dict):
    """One record of the store, a dict with the columns of ``fragment`` and
    ``verified`` as a bool."""

    @property
    def key(self):
        return (
            self["species"],
            self["charge"],
            self["multiplicity"],
            self["energy_fingerprint"],
            self["frequency_fingerprint"],
        )

    @staticmethod
    def make(
        species,
        *,
        charge=0,
        multiplicity=1,
        E=None,
        ZPE=None,
        H=None,
        T=298.15,
        energy_fingerprint=None,
        frequency_fingerprint=None,
        model_chemistry=None,
        frequency_model_chemistry=None,
        program=None,
        version=None,
        geometry_hash=None,
        job=None,
        configuration=None,
        source=None,
        date=None,
    ):
        """Build a record. Without an energy fingerprint the record is keyed on
        the model chemistry and code version and flagged unverified."""
        verified = bool(energy_fingerprint)
        if not energy_fingerprint:
            if not model_chemistry:
                raise ValueError(
                    "A record needs an energy fingerprint or a model chemistry."
                )
            energy_fingerprint = f"{UNVERIFIED}:{model_chemistry}@{version or '?'}"
        if frequency_fingerprint is None:
            frequency_fingerprint = ""
        if ZPE is None and H is None:
            frequency_fingerprint = ""
        elif not frequency_fingerprint:
            level = frequency_model_chemistry or model_chemistry
            frequency_fingerprint = f"{UNVERIFIED}:{level}@{version or '?'}"
        return FragmentRecord(
            species=species,
            charge=int(charge),
            multiplicity=int(multiplicity),
            energy_fingerprint=energy_fingerprint,
            frequency_fingerprint=frequency_fingerprint,
            verified=verified,
            model_chemistry=model_chemistry,
            frequency_model_chemistry=frequency_model_chemistry or model_chemistry,
            program=program,
            version=version,
            E=E,
            ZPE=ZPE,
            H=H,
            T=T if H is not None else None,
            geometry_hash=geometry_hash,
            job=str(job) if job is not None else None,
            configuration=str(configuration) if configuration is not None else None,
            source=source,
            date=date or datetime.datetime.now().isoformat(timespec="seconds"),
        )


class FragmentStore(GeometryEnergyMixin):
    """The fragment-energy store (see the module docstring), with the geometry
    and energy tables of :mod:`.geometry_store`.

    Parameters
    ----------
    path : str or Path, optional
        The SQLite file; default next to the ATcT cache.
    create : bool
        Create the file and schema if missing (the writer); a reader passes
        False and gets an empty store when there is no file.
    """

    def __init__(self, path=None, *, create=False):
        self.path = fragment_store_path(path)
        self._db = None
        self._readonly = not create
        if create:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._connect(create=True)

    # -- connection ---------------------------------------------------------

    def _connect(self, create=False):
        if self._db is not None:
            return self._db
        if create:
            db = sqlite3.connect(self.path, timeout=30)
        else:
            if not self.path.exists():
                raise FileNotFoundError(f"No fragment store at {self.path}")
            db = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True, timeout=30)
        db.row_factory = sqlite3.Row
        # Local-disk discipline: no memory mapping, wait for locks, the
        # rollback journal (never WAL on a network filesystem)
        db.execute("PRAGMA mmap_size = 0")
        db.execute("PRAGMA busy_timeout = 30000")
        if create:
            db.execute("PRAGMA journal_mode = DELETE")
            db.executescript(_SCHEMA + GEOMETRY_SCHEMA)
            if db.execute("PRAGMA user_version").fetchone()[0] < SCHEMA_VERSION:
                db.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
            db.commit()
        self._db = db
        return db

    def close(self):
        if self._db is not None:
            self._db.close()
            self._db = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __len__(self):
        try:
            db = self._connect()
        except FileNotFoundError:
            return 0
        return db.execute("SELECT COUNT(*) FROM fragment").fetchone()[0]

    # -- writing --------------------------------------------------------------

    def put(self, record):
        """Insert a record unless its key exists.

        Returns "inserted", "duplicate" (same key, same values) or "logged"
        (same key, values differing by more than 0.1 kJ/mol: the attempt went
        to ``duplicate_log``, the stored values stay).
        """
        if self._readonly:
            raise PermissionError("The fragment store was opened read-only.")
        db = self._connect(create=True)
        columns = [
            "species",
            "charge",
            "multiplicity",
            "energy_fingerprint",
            "frequency_fingerprint",
            "verified",
            "model_chemistry",
            "frequency_model_chemistry",
            "program",
            "version",
            "E",
            "ZPE",
            "H",
            "T",
            "geometry_hash",
            "job",
            "configuration",
            "source",
            "date",
        ]
        values = [record.get(c) for c in columns]
        values[5] = 1 if record.get("verified", True) else 0
        cursor = db.execute(
            f"INSERT OR IGNORE INTO fragment ({', '.join(columns)}) VALUES "
            f"({', '.join('?' * len(columns))})",
            values,
        )
        if cursor.rowcount == 1:
            db.commit()
            return "inserted"
        existing = db.execute(
            "SELECT * FROM fragment WHERE species = ? AND charge = ? AND "
            "multiplicity = ? AND energy_fingerprint = ? AND frequency_fingerprint = ?",
            record.key if isinstance(record, FragmentRecord) else values[:5],
        ).fetchone()
        difference = 0.0
        for field in ("E", "ZPE", "H"):
            a, b = existing[field], record.get(field)
            if a is None or b is None:
                continue
            difference = max(difference, abs(float(a) - float(b)))
        if difference <= DUPLICATE_TOLERANCE:
            return "duplicate"
        db.execute(
            "INSERT INTO duplicate_log (fragment, E, ZPE, H, geometry_hash, job, "
            "configuration, source, date, difference) VALUES (?,?,?,?,?,?,?,?,?,?)",
            (
                existing["id"],
                record.get("E"),
                record.get("ZPE"),
                record.get("H"),
                record.get("geometry_hash"),
                record.get("job"),
                record.get("configuration"),
                record.get("source"),
                record.get("date"),
                difference,
            ),
        )
        db.commit()
        logger.warning(
            f"The fragment store already has {record['species']} at this level with "
            f"values differing by {difference:.3f} kJ/mol (job "
            f"{record.get('job')}); kept the stored ones and logged the attempt."
        )
        return "logged"

    _VALUE_COLUMNS = (
        "verified",
        "model_chemistry",
        "frequency_model_chemistry",
        "program",
        "version",
        "E",
        "ZPE",
        "H",
        "T",
        "geometry_hash",
        "job",
        "configuration",
        "source",
        "date",
    )

    def supersede(self, record, reason, force=False):
        """Replace the stored values for the record's key, keeping the old ones.

        The old values go to the ``superseded`` table with today's date and
        `reason`; the record keeps its key and id. A key not yet in the store
        is simply inserted. A value the new record lacks (None) keeps the
        stored one, so a corrected energy does not erase a stored ZPE or H.

        With `force`, values that differ by less than the duplicate tolerance
        (0.1 kJ/mol) are replaced too, so that a correction applied to a whole
        level of theory reaches every record, however small; identical values
        are still left alone.

        Returns "inserted", "unchanged" (same values within 0.1 kJ/mol, or
        identical with `force`: nothing to replace) or "superseded".
        """
        if self._readonly:
            raise PermissionError("The fragment store was opened read-only.")
        if not reason:
            raise ValueError("Superseding a record needs a reason.")
        db = self._connect(create=True)
        if not isinstance(record, FragmentRecord):
            # The same key normalization as FragmentRecord.make (NULL never
            # matches in SQL, and UNIQUE treats NULLs as distinct)
            record = FragmentRecord(record)
            if not record.get("energy_fingerprint"):
                raise ValueError(
                    "A plain dict needs its energy fingerprint; build records with "
                    "FragmentRecord.make()."
                )
            record["frequency_fingerprint"] = record.get("frequency_fingerprint") or ""
            record["charge"] = int(record.get("charge") or 0)
            record["multiplicity"] = int(record.get("multiplicity") or 1)
        key = record.key
        existing = db.execute(
            "SELECT * FROM fragment WHERE species = ? AND charge = ? AND "
            "multiplicity = ? AND energy_fingerprint = ? AND frequency_fingerprint = ?",
            key,
        ).fetchone()
        if existing is None:
            return self.put(record)
        difference = 0.0
        for field in ("E", "ZPE", "H"):
            a, b = existing[field], record.get(field)
            if b is None:
                continue  # the stored value is kept
            if a is None:
                difference = float("inf")  # the new record adds a value
            else:
                difference = max(difference, abs(float(a) - float(b)))
        if difference == 0.0 or (difference <= DUPLICATE_TOLERANCE and not force):
            return "unchanged"
        columns = self._VALUE_COLUMNS
        now = datetime.datetime.now().isoformat(timespec="seconds")
        db.execute(
            f"INSERT INTO superseded (fragment, {', '.join(columns)}, superseded_on, "
            f"reason) VALUES ({', '.join('?' * (len(columns) + 3))})",
            [existing["id"]] + [existing[c] for c in columns] + [now, reason],
        )
        values = [record.get(c) for c in columns]
        values[0] = 1 if record.get("verified", True) else 0
        for i, c in enumerate(columns):
            if c in ("E", "ZPE", "H", "T") and values[i] is None:
                values[i] = existing[c]
        assignments = ", ".join(f"{c} = ?" for c in columns)
        db.execute(
            f"UPDATE fragment SET {assignments} WHERE id = ?",
            values + [existing["id"]],
        )
        db.commit()
        logger.info(
            f"Superseded {record['species']} (record {existing['id']}): {reason}"
        )
        return "superseded"

    def superseded(self, species=None):
        """The superseded values, newest first, each with the species and the
        values that replaced them; all, or those of one species."""
        try:
            db = self._connect()
        except FileNotFoundError:
            return []
        try:
            sql = (
                "SELECT s.*, f.species, f.charge, f.multiplicity, f.E AS current_E, "
                "f.ZPE AS current_ZPE, f.H AS current_H FROM superseded s "
                "JOIN fragment f ON f.id = s.fragment"
            )
            args = ()
            if species is not None:
                sql += " WHERE f.species = ?"
                args = (species,)
            rows = db.execute(sql + " ORDER BY s.id DESC", args).fetchall()
        except sqlite3.OperationalError:  # a store from before the table existed
            return []
        return [dict(r) for r in rows]

    # -- reading --------------------------------------------------------------

    def get(
        self,
        species,
        *,
        charge=0,
        multiplicity=1,
        energy_fingerprint=None,
        frequency_fingerprint=None,
        model_chemistry=None,
        version=None,
    ):
        """The best record for a species at a level, or None.

        With fingerprints, the verified record with those fingerprints (the
        frequency fingerprint may be None to accept any, including none).
        Without, or when there is no verified record, an unverified record of
        the same model chemistry (and code version if given); the caller sees
        ``record["verified"]`` and reports it.

        Any error reading the store gives None and a log message.
        """
        try:
            db = self._connect()
        except FileNotFoundError:
            return None
        except sqlite3.Error as e:
            logger.warning(f"The fragment store {self.path} could not be read: {e}")
            return None
        try:
            where = "species = ? AND charge = ? AND multiplicity = ?"
            args = [species, int(charge), int(multiplicity)]
            if energy_fingerprint:
                sql = where + " AND energy_fingerprint = ? AND verified = 1"
                args_v = args + [energy_fingerprint]
                if frequency_fingerprint is not None:
                    sql += " AND frequency_fingerprint = ?"
                    args_v.append(frequency_fingerprint)
                row = db.execute(
                    f"SELECT * FROM fragment WHERE {sql} ORDER BY id LIMIT 1", args_v
                ).fetchone()
                if row is not None:
                    return self._record(row)
            if model_chemistry:
                sql = where + " AND verified = 0 AND model_chemistry = ?"
                args_u = args + [model_chemistry]
                if version:
                    sql += " AND version = ?"
                    args_u.append(version)
                row = db.execute(
                    f"SELECT * FROM fragment WHERE {sql} ORDER BY id LIMIT 1", args_u
                ).fetchone()
                if row is not None:
                    return self._record(row)
            return None
        except sqlite3.Error as e:
            logger.warning(f"The fragment store {self.path} could not be read: {e}")
            return None

    def records(self, species=None):
        """All records, or those of one species, newest first."""
        try:
            db = self._connect()
        except FileNotFoundError:
            return []
        if species is None:
            rows = db.execute("SELECT * FROM fragment ORDER BY id DESC").fetchall()
        else:
            rows = db.execute(
                "SELECT * FROM fragment WHERE species = ? ORDER BY id DESC", (species,)
            ).fetchall()
        return [self._record(r) for r in rows]

    def duplicates(self):
        """The duplicate log, with the stored record's species and level."""
        try:
            db = self._connect()
        except FileNotFoundError:
            return []
        rows = db.execute(
            "SELECT d.*, f.species, f.model_chemistry, f.E AS stored_E, "
            "f.ZPE AS stored_ZPE, f.H AS stored_H FROM duplicate_log d "
            "JOIN fragment f ON f.id = d.fragment ORDER BY d.id DESC"
        ).fetchall()
        return [dict(r) for r in rows]

    @staticmethod
    def _record(row):
        record = FragmentRecord(dict(row))
        record["verified"] = bool(record["verified"])
        return record


# -- harvesting --------------------------------------------------------------


def _to_kj(value, units):
    factor = _TO_KJ.get(units)
    if factor is None:
        try:
            from seamm_util import Q_

            return Q_(float(value), units).m_as("kJ/mol")
        except Exception:
            raise ValueError(f"Cannot convert energy units '{units}' to kJ/mol")
    return float(value) * factor


def _property(properties, pattern):
    """(name, value in kJ/mol) of the one property matching the glob, or None."""
    if not pattern:
        return None
    values = properties.get(pattern=pattern, match="glob")
    if not values:
        return None
    if len(values) > 1:
        raise ValueError(f"Several properties match '{pattern}': {sorted(values)}")
    (name, value), *_ = values.items()
    if isinstance(value, dict):
        value = value["value"]
    units = properties.units(name)
    return name, _to_kj(value, units)


def _text_property(properties, pattern):
    if not pattern:
        return None
    values = properties.get(pattern=pattern, match="glob", types=["str"])
    if not values:
        return None
    (name, value), *_ = values.items()
    if isinstance(value, dict):
        value = value["value"]
    return str(value)


def load_species_map(path):
    """A species map from a CSV file: {(job, system name): row}.

    The columns are ``job`` (the job id, or empty for any job), ``system name``
    and ``species`` (the store key: the canonical SMILES, with the stereo of the
    structure), and optionally ``configuration name`` (the configuration to use,
    default the system's current one), ``charge``, ``multiplicity`` and
    ``conformer label``, and the thermal corrections ``E (kJ/mol)``, ``ZPE
    (kJ/mol)``, ``H298 (kJ/mol)`` and ``symmetry number`` (see
    :func:`harvest_job`). Other columns are ignored.
    """
    import csv

    out = {}
    with open(Path(path).expanduser(), newline="") as fd:
        for row in csv.DictReader(fd):
            name = (row.get("system name") or "").strip()
            species = (row.get("species") or "").strip()
            if not name or not species:
                continue
            job = (row.get("job") or "").strip() or None
            key = (job, name)
            if key in out:
                raise ValueError(f"{path}: '{name}' of job {job} is listed twice")
            out[key] = {k: (v or "").strip() for k, v in row.items() if k}
    return out


def harvest_job(
    job_directory,
    store,
    *,
    energy="total energy#*#*",
    zpe="ZPE#*#*",
    enthalpy="H(298.15 K)#*#*",
    fingerprint="settings fingerprint#*#*",
    frequency_fingerprint="frequency fingerprint#*#*",
    name_prefix="CBH:",
    all_systems=False,
    database="seamm.db",
    supersede=None,
    force=False,
    species_map=None,
    tables="all",
):
    """Promote the species energies of one job into the store.

    Every system whose name starts with `name_prefix` (``CBH:<SMILES>``: the
    SMILES is the species identity), or every system when `all_systems`
    (the name, less the prefix if present, canonicalized, is the species, as
    for systems that Read Structure named by SMILES), contributes one record
    per energy
    property found on its current configuration, with the ZPE and enthalpy
    of the same level when present, the fingerprints if the QC step stored
    them, and the provenance (job directory, configuration id, geometry
    hash, date). The model chemistry is the third ``#`` field of the energy
    property's name.

    With a `species_map` ({(job, system name): row}, see
    :func:`load_species_map`) the systems are found by their names in the map
    instead of by the prefix: the map gives each one's species (the store key,
    with the stereo of its structure), and optionally its configuration,
    charge, multiplicity and conformer label. Systems the map does not list are
    skipped. `tables` is "all" (the fragment table too) or "geometry" (only the
    geometry and energy tables, e.g. for several conformers of a species, which
    the fragment table cannot hold).

    The species map may also carry the thermal corrections of the structures
    whose jobs did not store them as properties: ``E (kJ/mol)``, ``ZPE
    (kJ/mol)`` (unscaled), ``H298 (kJ/mol)`` (H at 298.15 K) and ``symmetry
    number``. They are used only for an energy with no ZPE or enthalpy of its
    own level, as H - E added to the job's energy, and only when the map's E
    agrees with the job's within 0.1 kJ/mol (otherwise the system is reported
    and gets no geometry row).

    With `supersede` (the reason, e.g. "lowest conformer from B2PLYP
    geometries"), a record whose key is already stored with different values
    replaces them, the old values kept in the ``superseded`` table, instead of
    being logged as a differing duplicate.
    With `force` as well, values within the duplicate tolerance are replaced too
    (see :meth:`FragmentStore.supersede`).

    The same properties also fill the geometry and energy tables
    (:mod:`.geometry_store`): an energy property with a ZPE or enthalpy of the
    same level makes a geometry row at that level (the configuration's
    coordinates, its energy, ZPE and H, and the frequencies, symmetry number
    and moments of inertia when stored as ``frequencies``, ``symmetry
    number`` and ``moments of inertia`` properties of the same level), and
    every energy property makes an energy row on the configuration's geometry.
    An energy whose ``calculation`` record names a different geometry is
    refused (the pairing guard) and reported in "skipped".

    Returns {"inserted": n, "duplicate": n, "logged": n, "superseded": n,
    "skipped": [reasons], "geometry": {outcome: n}, "energy": {outcome: n}},
    the first four counting the fragment rows.
    """
    job_directory = Path(job_directory).expanduser().resolve()
    path = job_directory / database
    if not path.exists():
        raise FileNotFoundError(f"No {database} in {job_directory}")
    try:
        from molsystem import SystemDB
    except ImportError:  # pragma: no cover
        raise ImportError(
            "Harvesting a job needs molsystem (pip install molsystem), which reads "
            "the job's database."
        )
    job_id = None
    try:
        text = (job_directory / "job_data.json").read_text()
        if text.startswith("!"):  # SEAMM's header line, "!MolSSI job_data 1.0"
            text = text.split("\n", 1)[1] if "\n" in text else ""
        job_id = json.loads(text).get("job id")
    except (OSError, ValueError, AttributeError):
        pass
    job = f"{job_id}:{job_directory}" if job_id is not None else str(job_directory)

    counts = {
        "inserted": 0,
        "duplicate": 0,
        "logged": 0,
        "superseded": 0,
        "skipped": [],
        "geometry": {},
        "energy": {},
    }

    def tally(table, outcome):
        counts[table][outcome] = counts[table].get(outcome, 0) + 1

    if tables not in ("all", "geometry"):
        raise ValueError(f"tables must be 'all' or 'geometry', not '{tables}'")

    def thermal_from_map(mapped, E):
        """(ZPE, H, symmetry number) from the map's row, or None; raises
        ValueError if the map's energy is not the job's."""
        if mapped is None:
            return None
        try:
            zpe = float(mapped["ZPE (kJ/mol)"])
            h = float(mapped["H298 (kJ/mol)"])
            e = float(mapped["E (kJ/mol)"])
        except (KeyError, TypeError, ValueError):
            return None
        if abs(e - E) > DUPLICATE_TOLERANCE:
            raise ValueError(
                f"the species map's E ({e:.3f} kJ/mol) is not the job's ({E:.3f})"
            )
        try:
            symmetry = int(mapped["symmetry number"])
        except (KeyError, TypeError, ValueError):
            symmetry = None
        return zpe, E + (h - e), symmetry

    db = SystemDB(filename=f"file:{path}?mode=ro")
    try:
        for system in db.get_systems():
            mapped = None
            if species_map is not None:
                mapped = species_map.get((str(job_id), system.name)) or (
                    species_map.get((None, system.name))
                )
                if mapped is None:
                    counts["skipped"].append(f"{system.name}: not in the species map")
                    continue
                species = mapped["species"]
            elif system.name.startswith(name_prefix):
                species = system.name[len(name_prefix) :].strip()
            elif all_systems:
                species = system.name.strip()
            else:
                continue
            if mapped is None:
                try:  # the store is keyed on RDKit's canonical SMILES
                    from .atct import canonical_smiles

                    species = canonical_smiles(species)
                except Exception:
                    pass
            configuration = system.configuration
            if mapped is not None and mapped.get("configuration name"):
                try:
                    cid = system.get_configuration_id(mapped["configuration name"])
                    configuration = system.get_configuration(cid)
                except Exception:
                    counts["skipped"].append(
                        f"{system.name}: no configuration named "
                        f"'{mapped['configuration name']}'"
                    )
                    continue
            charge = configuration.charge
            multiplicity = configuration.spin_multiplicity or 1
            conformer = ""
            if mapped is not None:
                if mapped.get("charge"):
                    charge = int(mapped["charge"])
                if mapped.get("multiplicity"):
                    multiplicity = int(mapped["multiplicity"])
                conformer = mapped.get("conformer label") or ""
            properties = configuration.properties
            try:
                energies = properties.get(pattern=energy, match="glob")
            except Exception as e:  # pragma: no cover
                counts["skipped"].append(f"{species}: {e}")
                continue
            if not energies:
                counts["skipped"].append(f"{species}: no property matching '{energy}'")
                continue
            symbols = configuration.atoms.symbols
            coordinates = configuration.atoms.get_coordinates(fractionals=False)
            ghash = geometry_hash(symbols, coordinates)
            try:
                atomic_numbers = configuration.atoms.atomic_numbers
            except Exception:  # pragma: no cover
                atomic_numbers = None
            fingerprint_hash = hash_geometry(symbols, coordinates, atomic_numbers)
            for name, value in energies.items():
                if isinstance(value, dict):
                    value = value["value"]
                fields = name.split("#")
                model = fields[2] if len(fields) > 2 else None
                suffix = "#" + "#".join(fields[1:]) if len(fields) > 1 else ""
                E = _to_kj(value, properties.units(name))
                zpe_name = zpe.split("#")[0] + suffix if suffix else zpe
                h_name = enthalpy.split("#")[0] + suffix if suffix else enthalpy
                ZPE = _property(properties, zpe_name)
                H = _property(properties, h_name)
                fp = _text_property(
                    properties, fingerprint.split("#")[0] + suffix if suffix else None
                )
                mapped_symmetry = None
                from_map = False
                if ZPE is None and H is None:
                    try:
                        found = thermal_from_map(mapped, E)
                    except ValueError as e:
                        counts["skipped"].append(f"{system.name}: {e}")
                        found = None
                    if found is not None:
                        ZPE = (zpe_name, found[0])
                        H = (h_name, found[1])
                        mapped_symmetry = found[2]
                        from_map = True
                ffp = _text_property(
                    properties,
                    (frequency_fingerprint.split("#")[0] + suffix if suffix else None),
                )
                record = FragmentRecord.make(
                    species,
                    charge=charge,
                    multiplicity=multiplicity,
                    E=E,
                    ZPE=ZPE[1] if ZPE else None,
                    H=H[1] if H else None,
                    energy_fingerprint=fp,
                    frequency_fingerprint=ffp,
                    model_chemistry=model,
                    geometry_hash=ghash,
                    job=job,
                    configuration=configuration.id,
                    source=name,
                )
                if tables == "geometry":
                    pass
                elif supersede:
                    outcome = store.supersede(record, supersede, force=force)
                    counts["duplicate" if outcome == "unchanged" else outcome] += 1
                else:
                    counts[store.put(record)] += 1

                # The geometry and energy tables
                def extra(base, kind="json"):
                    if not suffix:
                        return None
                    full = base + suffix
                    try:
                        if not properties.exists(full):
                            return None
                        value = properties.get(full)
                    except Exception:
                        return None
                    if isinstance(value, dict) and full in value:
                        value = value[full]  # {name: {"value": ...}}
                    if isinstance(value, dict) and "value" in value:
                        value = value["value"]
                    if kind == "json" and isinstance(value, str):
                        try:
                            value = json.loads(value)
                        except ValueError:
                            return None
                    return value

                calculation = extra("calculation")
                calc_geometry = (
                    calculation.get("geometry fingerprint")
                    if isinstance(calculation, dict)
                    else None
                )
                program = version = None
                if isinstance(calculation, dict):
                    program = calculation.get("program")
                    version = calculation.get("version")
                if ZPE or H:
                    symmetry = extra("symmetry number", kind="int")
                    if symmetry is None:
                        symmetry = mapped_symmetry
                    try:
                        geometry = GeometryRecord.make(
                            species,
                            symbols,
                            coordinates,
                            charge=charge,
                            multiplicity=multiplicity,
                            conformer=conformer,
                            # Thermal corrections from the map were computed with
                            # the energy's own calculation
                            frequency_fingerprint=ffp or (fp if from_map else None),
                            optimization_fingerprint=fp,
                            model_chemistry=model,
                            program=program,
                            version=version,
                            E=E,
                            ZPE=ZPE[1] if ZPE else None,
                            H=H[1] if H else None,
                            frequencies=extra("frequencies"),
                            symmetry_number=(
                                int(symmetry) if symmetry is not None else None
                            ),
                            moments=extra("moments of inertia"),
                            atomic_numbers=atomic_numbers,
                            geometry_hash=fingerprint_hash,
                            job=job,
                            configuration=configuration.id,
                            source=name,
                        )
                        if supersede:
                            outcome, _ = store.supersede_geometry(
                                geometry, supersede, force=force
                            )
                        else:
                            outcome, _ = store.put_geometry(geometry)
                        tally("geometry", outcome)
                    except PairingError as e:
                        counts["skipped"].append(f"{species}: {e}")
                try:
                    energy_row = EnergyRecord.make(
                        fingerprint_hash,
                        E,
                        energy_fingerprint=fp,
                        calculation_geometry=calc_geometry,
                        model_chemistry=model,
                        program=program,
                        version=version,
                        job=job,
                        configuration=configuration.id,
                        source=name,
                    )
                except PairingError as e:
                    counts["skipped"].append(f"{species} ({name}): {e}")
                    continue
                if supersede:
                    outcome, _ = store.supersede_energy(
                        energy_row, supersede, force=force
                    )
                else:
                    outcome, _ = store.put_energy(energy_row)
                tally("energy", outcome)
    finally:
        db.close()
    return counts
