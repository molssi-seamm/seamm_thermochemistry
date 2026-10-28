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
overwritten, which catches a different conformer or SCF solution.

Reads are best effort: a missing, locked or corrupt store gives None and a
log message, and the caller computes the species.
"""

import datetime
import hashlib
import json
import logging
from pathlib import Path
import sqlite3

logger = logging.getLogger(__name__)

__all__ = ["FragmentStore", "FragmentRecord", "fragment_store_path", "harvest_job"]

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


class FragmentStore:
    """The fragment-energy store (see the module docstring).

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
            db.executescript(_SCHEMA)
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

    Returns {"inserted": n, "duplicate": n, "logged": n, "skipped": [reasons]}.
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
        with open(job_directory / "job_data.json") as fd:
            job_id = json.load(fd).get("job id")
    except (OSError, ValueError):
        pass
    job = f"{job_id}:{job_directory}" if job_id is not None else str(job_directory)

    counts = {"inserted": 0, "duplicate": 0, "logged": 0, "skipped": []}
    db = SystemDB(filename=f"file:{path}?mode=ro")
    try:
        for system in db.get_systems():
            if system.name.startswith(name_prefix):
                species = system.name[len(name_prefix) :].strip()
            elif all_systems:
                species = system.name.strip()
            else:
                continue
            try:  # the store is keyed on RDKit's canonical SMILES
                from .atct import canonical_smiles

                species = canonical_smiles(species)
            except Exception:
                pass
            configuration = system.configuration
            properties = configuration.properties
            try:
                energies = properties.get(pattern=energy, match="glob")
            except Exception as e:  # pragma: no cover
                counts["skipped"].append(f"{species}: {e}")
                continue
            if not energies:
                counts["skipped"].append(f"{species}: no property matching '{energy}'")
                continue
            ghash = geometry_hash(
                configuration.atoms.symbols,
                configuration.atoms.get_coordinates(fractionals=False),
            )
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
                ffp = _text_property(
                    properties,
                    (frequency_fingerprint.split("#")[0] + suffix if suffix else None),
                )
                record = FragmentRecord.make(
                    species,
                    charge=configuration.charge,
                    multiplicity=configuration.spin_multiplicity or 1,
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
                counts[store.put(record)] += 1
    finally:
        db.close()
    return counts
