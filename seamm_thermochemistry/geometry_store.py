# -*- coding: utf-8 -*-

"""Geometries and energies in the fragment store: a species' energy at one
level of theory on a geometry from another.

The ``fragment`` table holds E, ZPE and H of a species at one level of theory,
or at an energy level and a frequency level, with no geometry. Thermochemistry
often wants the energy at a costlier level (e.g. revDSD-PBEP86 or
CCSD(T)-F12) on the geometry and frequencies of a cheaper one (e.g.
r2SCAN-3c), and the same geometry serves several energy levels. Two tables hold
that:

``geometry``
    One row per conformer of a species at a geometry/frequency level. The key
    is (species, charge, multiplicity, frequency fingerprint, conformer
    label). The row holds the coordinates and their hash
    (:func:`seamm_util.fingerprint.geometry_fingerprint`), the energy, the
    unscaled ZPE and H(T) at that level, and when known the harmonic
    frequencies, the number of imaginary modes, the rotational symmetry number
    and the moments of inertia, so that H(T) at other temperatures (and later
    S and G) can be recomputed. ZPEs are stored unscaled; factors are applied
    at use.

``energy``
    One row per single-point energy on a geometry, keyed by the geometry's
    hash and the calculation's settings fingerprint. The hash, not the row id,
    is the key, so copies of the store on other hosts keep working; the row id
    of the geometry is kept as a convenience. Other results of the same
    calculation (e.g. the MP2 energy of a CCSD(T) run, or the (T) correction)
    go in ``byproducts``.

A composite energy (a recipe of several components) is never stored, only its
components; it is combined at use, as ZPE factors are applied at use, so a
corrected recipe or factor needs no recomputation.

Guards, applied when rows are written:

* An energy whose calculation recorded the geometry it ran on (its geometry
  fingerprint) must be on the geometry of the row it is attached to; otherwise
  it is refused.
* A geometry row is optimized and has its frequencies at one level: when the
  optimization's model chemistry is given it must be the frequency level's.
  (The settings fingerprints of the two may differ legitimately, e.g. a
  Hessian by finite differences converges the SCF more tightly, so the model
  chemistry and program are compared, and both fingerprints are kept.)
* A geometry with imaginary frequencies is stored, so that its provenance is
  kept, but :meth:`geometries` does not return it for reuse unless asked:
  modes up to ``IMAGINARY_TOLERANCE`` (cm⁻¹) are taken as numerical noise.

Writes follow the ``fragment`` table's rules: insert-or-ignore, a differing
duplicate (by more than 0.1 kJ/mol) logged in ``row_history`` and never
overwritten, and an explicit ``supersede`` with a reason that keeps the old
values in ``row_history``.
"""

import datetime
import json
import logging
import sqlite3

logger = logging.getLogger(__name__)

__all__ = [
    "GeometryRecord",
    "EnergyRecord",
    "GeometryEnergyMixin",
    "PairingError",
    "IMAGINARY_TOLERANCE",
    "SCHEMA_VERSION",
]

#: The store's schema version (``PRAGMA user_version``): 2 adds geometry,
#: energy and row_history; 0 or 1 is the fragment table alone.
SCHEMA_VERSION = 2
#: Imaginary modes at or below this (cm⁻¹, magnitude) are numerical noise
IMAGINARY_TOLERANCE = 20.0
DUPLICATE_TOLERANCE = 0.1  # kJ/mol, as for the fragment table

GEOMETRY_SCHEMA = """
CREATE TABLE IF NOT EXISTS geometry (
    id INTEGER PRIMARY KEY,
    species TEXT NOT NULL,
    charge INTEGER NOT NULL DEFAULT 0,
    multiplicity INTEGER NOT NULL DEFAULT 1,
    frequency_fingerprint TEXT NOT NULL,
    conformer TEXT NOT NULL DEFAULT '',
    verified INTEGER NOT NULL DEFAULT 1,
    geometry_hash TEXT NOT NULL,
    optimization_fingerprint TEXT,
    model_chemistry TEXT,
    program TEXT,
    version TEXT,
    symbols TEXT NOT NULL,
    coordinates TEXT NOT NULL,
    E REAL,
    ZPE REAL,
    H REAL,
    T REAL,
    frequencies TEXT,
    n_imaginary INTEGER,
    max_imaginary REAL,
    symmetry_number INTEGER,
    moments TEXT,
    conformer_provenance TEXT,
    job TEXT,
    configuration TEXT,
    source TEXT,
    date TEXT,
    UNIQUE (species, charge, multiplicity, frequency_fingerprint, conformer)
);
CREATE INDEX IF NOT EXISTS geometry_species ON geometry (species, charge, multiplicity);
CREATE INDEX IF NOT EXISTS geometry_hash_index ON geometry (geometry_hash);
CREATE TABLE IF NOT EXISTS energy (
    id INTEGER PRIMARY KEY,
    geometry_hash TEXT NOT NULL,
    energy_fingerprint TEXT NOT NULL,
    geometry INTEGER REFERENCES geometry (id),
    verified INTEGER NOT NULL DEFAULT 1,
    model_chemistry TEXT,
    program TEXT,
    version TEXT,
    E REAL NOT NULL,
    byproducts TEXT,
    job TEXT,
    configuration TEXT,
    source TEXT,
    date TEXT,
    UNIQUE (geometry_hash, energy_fingerprint)
);
CREATE TABLE IF NOT EXISTS row_history (
    id INTEGER PRIMARY KEY,
    table_name TEXT NOT NULL,
    row INTEGER NOT NULL,
    kind TEXT NOT NULL,
    old TEXT,
    new TEXT,
    difference REAL,
    reason TEXT,
    date TEXT
);
"""

_JSON_COLUMNS = ("symbols", "coordinates", "frequencies", "moments", "byproducts")


class PairingError(ValueError):
    """An energy that is not on the geometry it would be attached to, or a
    geometry optimized at one level with frequencies from another."""


def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def _unverified(level, version):
    return f"unverified:{level}@{version or '?'}"


class GeometryRecord(dict):
    """A row of the ``geometry`` table (a dict; JSON columns as Python values)."""

    @property
    def key(self):
        return (
            self["species"],
            self["charge"],
            self["multiplicity"],
            self["frequency_fingerprint"],
            self["conformer"],
        )

    @staticmethod
    def make(
        species,
        symbols,
        coordinates,
        *,
        charge=0,
        multiplicity=1,
        conformer="",
        frequency_fingerprint=None,
        optimization_fingerprint=None,
        model_chemistry=None,
        optimization_model_chemistry=None,
        program=None,
        version=None,
        E=None,
        ZPE=None,
        H=None,
        T=298.15,
        frequencies=None,
        symmetry_number=None,
        moments=None,
        conformer_provenance=None,
        atomic_numbers=None,
        geometry_hash=None,
        job=None,
        configuration=None,
        source=None,
        date=None,
    ):
        """Build a geometry row.

        Without a frequency fingerprint the row is keyed on the model chemistry
        and code version and flagged unverified. The geometry hash is computed
        from the atomic numbers (or symbols) and coordinates unless given.
        ``optimization_model_chemistry``, when given, must be the frequency
        level's (``model_chemistry``).
        """
        if (
            optimization_model_chemistry
            and model_chemistry
            and optimization_model_chemistry != model_chemistry
        ):
            raise PairingError(
                f"{species}: optimized at {optimization_model_chemistry} but the "
                f"frequencies are at {model_chemistry}. A geometry row has its "
                "geometry and frequencies at one level."
            )
        verified = bool(frequency_fingerprint)
        if not frequency_fingerprint:
            if not model_chemistry:
                raise ValueError(
                    "A geometry needs a frequency fingerprint or a model chemistry."
                )
            frequency_fingerprint = _unverified(model_chemistry, version)
        symbols = [str(s) for s in symbols]
        coordinates = [[float(x) for x in xyz] for xyz in coordinates]
        if geometry_hash is None:
            geometry_hash = hash_geometry(symbols, coordinates, atomic_numbers)
        n_imaginary = max_imaginary = None
        if frequencies is not None:
            frequencies = [float(f) for f in frequencies]
            imaginary = [-f for f in frequencies if f < 0]
            n_imaginary = sum(1 for f in imaginary if f > IMAGINARY_TOLERANCE)
            max_imaginary = max(imaginary) if imaginary else 0.0
        return GeometryRecord(
            species=species,
            charge=int(charge),
            multiplicity=int(multiplicity),
            frequency_fingerprint=frequency_fingerprint,
            conformer=conformer or "",
            verified=verified,
            geometry_hash=geometry_hash,
            optimization_fingerprint=optimization_fingerprint,
            model_chemistry=model_chemistry,
            program=program,
            version=version,
            symbols=symbols,
            coordinates=coordinates,
            E=E,
            ZPE=ZPE,
            H=H,
            T=T if H is not None else None,
            frequencies=frequencies,
            n_imaginary=n_imaginary,
            max_imaginary=max_imaginary,
            symmetry_number=symmetry_number,
            moments=[float(m) for m in moments] if moments is not None else None,
            conformer_provenance=conformer_provenance,
            job=str(job) if job is not None else None,
            configuration=str(configuration) if configuration is not None else None,
            source=source,
            date=date or _now(),
        )


class EnergyRecord(dict):
    """A row of the ``energy`` table."""

    @property
    def key(self):
        return (self["geometry_hash"], self["energy_fingerprint"])

    @staticmethod
    def make(
        geometry_hash,
        E,
        *,
        energy_fingerprint=None,
        calculation_geometry=None,
        model_chemistry=None,
        program=None,
        version=None,
        byproducts=None,
        job=None,
        configuration=None,
        source=None,
        date=None,
    ):
        """Build an energy row for the geometry with hash `geometry_hash`.

        `calculation_geometry` is the geometry fingerprint the calculation
        recorded, if any; it must equal `geometry_hash` (the pairing guard).
        Without an energy fingerprint the row is keyed on the model chemistry
        and code version and flagged unverified.
        """
        if calculation_geometry and calculation_geometry != geometry_hash:
            raise PairingError(
                "The energy was computed on a different geometry "
                f"({calculation_geometry}) from the one it would be stored with "
                f"({geometry_hash})."
            )
        verified = bool(energy_fingerprint)
        if not energy_fingerprint:
            if not model_chemistry:
                raise ValueError(
                    "An energy needs an energy fingerprint or a model chemistry."
                )
            energy_fingerprint = _unverified(model_chemistry, version)
        return EnergyRecord(
            geometry_hash=geometry_hash,
            energy_fingerprint=energy_fingerprint,
            verified=verified,
            model_chemistry=model_chemistry,
            program=program,
            version=version,
            E=float(E),
            byproducts=byproducts,
            job=str(job) if job is not None else None,
            configuration=str(configuration) if configuration is not None else None,
            source=source,
            date=date or _now(),
        )


def hash_geometry(symbols, coordinates, atomic_numbers=None):
    """The geometry's fingerprint, as QC steps record it
    (:func:`seamm_util.fingerprint.geometry_fingerprint`, 4 decimals in Å)."""
    from seamm_util.fingerprint import geometry_fingerprint

    if atomic_numbers is None:
        from seamm_util import element_data

        atomic_numbers = [
            element_data[s.rstrip(":").strip().capitalize()]["atomic number"]
            for s in symbols
        ]
    return geometry_fingerprint(atomic_numbers, coordinates)


def _encode(record, columns):
    values = []
    for c in columns:
        v = record.get(c)
        if c in _JSON_COLUMNS and v is not None:
            v = json.dumps(v)
        elif c == "verified":
            v = 1 if v is None or v else 0
        values.append(v)
    return values


def _decode(row, cls):
    record = cls(dict(row))
    for c in _JSON_COLUMNS:
        if c in record and record[c] is not None:
            record[c] = json.loads(record[c])
    if "verified" in record:
        record["verified"] = bool(record["verified"])
    return record


_GEOMETRY_COLUMNS = (
    "species",
    "charge",
    "multiplicity",
    "frequency_fingerprint",
    "conformer",
    "verified",
    "geometry_hash",
    "optimization_fingerprint",
    "model_chemistry",
    "program",
    "version",
    "symbols",
    "coordinates",
    "E",
    "ZPE",
    "H",
    "T",
    "frequencies",
    "n_imaginary",
    "max_imaginary",
    "symmetry_number",
    "moments",
    "conformer_provenance",
    "job",
    "configuration",
    "source",
    "date",
)
_GEOMETRY_KEY = 5
_GEOMETRY_VALUES = ("E", "ZPE", "H")

_ENERGY_COLUMNS = (
    "geometry_hash",
    "energy_fingerprint",
    "geometry",
    "verified",
    "model_chemistry",
    "program",
    "version",
    "E",
    "byproducts",
    "job",
    "configuration",
    "source",
    "date",
)
_ENERGY_VALUES = ("E",)


class GeometryEnergyMixin:
    """The geometry and energy tables of :class:`FragmentStore` (see the
    module docstring). The store provides ``_connect``, ``_readonly`` and
    ``path``."""

    # -- writing ----------------------------------------------------------------

    def _writer(self):
        if self._readonly:
            raise PermissionError("The fragment store was opened read-only.")
        return self._connect(create=True)

    def _history(self, db, table, row, kind, old, new, difference=None, reason=None):
        db.execute(
            "INSERT INTO row_history (table_name, row, kind, old, new, difference, "
            "reason, date) VALUES (?,?,?,?,?,?,?,?)",
            (
                table,
                row,
                kind,
                json.dumps(old, default=str),
                json.dumps(new, default=str),
                difference,
                reason,
                _now(),
            ),
        )

    @staticmethod
    def _difference(existing, record, fields):
        difference = 0.0
        for field in fields:
            a, b = existing[field], record.get(field)
            if a is None or b is None:
                continue
            difference = max(difference, abs(float(a) - float(b)))
        return difference

    def _put(self, table, columns, nkey, fields, record):
        db = self._writer()
        values = _encode(record, columns)
        cursor = db.execute(
            f"INSERT OR IGNORE INTO {table} ({', '.join(columns)}) VALUES "
            f"({', '.join('?' * len(columns))})",
            values,
        )
        if cursor.rowcount == 1:
            db.commit()
            return "inserted", cursor.lastrowid
        where = " AND ".join(f"{c} = ?" for c in columns[:nkey])
        existing = db.execute(
            f"SELECT * FROM {table} WHERE {where}", values[:nkey]
        ).fetchone()
        difference = self._difference(existing, record, fields)
        if difference <= DUPLICATE_TOLERANCE:
            return "duplicate", existing["id"]
        self._history(
            db,
            table,
            existing["id"],
            "duplicate",
            dict(existing),
            dict(record),
            difference=difference,
        )
        db.commit()
        logger.warning(
            f"The fragment store's {table} table already has this row with values "
            f"differing by {difference:.3f} kJ/mol (job {record.get('job')}); kept "
            "the stored ones and logged the attempt."
        )
        return "logged", existing["id"]

    def put_geometry(self, record):
        """Insert a geometry row unless its key exists.

        Returns ("inserted" | "duplicate" | "logged", row id), as
        :meth:`FragmentStore.put` does for the fragment table.
        """
        return self._put(
            "geometry", _GEOMETRY_COLUMNS, _GEOMETRY_KEY, _GEOMETRY_VALUES, record
        )

    def put_energy(self, record):
        """Insert an energy row unless its key exists.

        The row is linked to the geometry row with the same hash when there is
        one. Returns ("inserted" | "duplicate" | "logged", row id).
        """
        record = EnergyRecord(record)
        if record.get("geometry") is None:
            db = self._writer()
            row = db.execute(
                "SELECT id FROM geometry WHERE geometry_hash = ? ORDER BY id LIMIT 1",
                (record["geometry_hash"],),
            ).fetchone()
            record["geometry"] = row["id"] if row is not None else None
        return self._put("energy", _ENERGY_COLUMNS, 2, _ENERGY_VALUES, record)

    def _supersede(self, table, columns, nkey, fields, record, reason, force, put):
        if not reason:
            raise ValueError("Superseding a row needs a reason.")
        db = self._writer()
        values = _encode(record, columns)
        where = " AND ".join(f"{c} = ?" for c in columns[:nkey])
        existing = db.execute(
            f"SELECT * FROM {table} WHERE {where}", values[:nkey]
        ).fetchone()
        if existing is None:
            return put(record)
        difference = 0.0
        for field in fields:
            a, b = existing[field], record.get(field)
            if b is None:
                continue
            if a is None:
                difference = float("inf")
            else:
                difference = max(difference, abs(float(a) - float(b)))
        if difference == 0.0 or (difference <= DUPLICATE_TOLERANCE and not force):
            return "unchanged", existing["id"]
        self._history(
            db,
            table,
            existing["id"],
            "superseded",
            dict(existing),
            dict(record),
            difference=difference,
            reason=reason,
        )
        value_columns = columns[nkey:]
        new = dict(zip(columns, values))
        for c in fields:
            if new[c] is None:
                new[c] = existing[c]
        if table == "energy" and new.get("geometry") is None:
            new["geometry"] = existing["geometry"]
        db.execute(
            f"UPDATE {table} SET {', '.join(f'{c} = ?' for c in value_columns)} "
            "WHERE id = ?",
            [new[c] for c in value_columns] + [existing["id"]],
        )
        db.commit()
        logger.info(f"Superseded {table} row {existing['id']}: {reason}")
        return "superseded", existing["id"]

    def supersede_geometry(self, record, reason, force=False):
        """Replace a geometry row's values, keeping the old ones in
        ``row_history``. A geometry with new coordinates gets a new hash, so
        its energy rows (keyed on the old hash) no longer attach to it; they
        stay in the store for the old geometry. Returns (outcome, row id)."""
        return self._supersede(
            "geometry",
            _GEOMETRY_COLUMNS,
            _GEOMETRY_KEY,
            _GEOMETRY_VALUES,
            record,
            reason,
            force,
            self.put_geometry,
        )

    def supersede_energy(self, record, reason, force=False):
        """Replace an energy row's value, keeping the old one in
        ``row_history``. Returns (outcome, row id)."""
        return self._supersede(
            "energy",
            _ENERGY_COLUMNS,
            2,
            _ENERGY_VALUES,
            EnergyRecord(record),
            reason,
            force,
            self.put_energy,
        )

    # -- reading ----------------------------------------------------------------

    def _reader(self):
        try:
            return self._connect()
        except FileNotFoundError:
            return None
        except sqlite3.Error as e:
            logger.warning(f"The fragment store {self.path} could not be read: {e}")
            return None

    def n_geometries(self):
        """The number of geometry rows (0 for an old or missing store)."""
        db = self._reader()
        if db is None:
            return 0
        try:
            return db.execute("SELECT COUNT(*) FROM geometry").fetchone()[0]
        except sqlite3.Error:
            return 0

    def schema_version(self):
        """The store's ``PRAGMA user_version`` (0 for a new or old store)."""
        db = self._reader()
        if db is None:
            return 0
        return db.execute("PRAGMA user_version").fetchone()[0]

    def geometries(
        self,
        species,
        *,
        charge=0,
        multiplicity=1,
        frequency_fingerprint=None,
        model_chemistry=None,
        version=None,
        include_imaginary=False,
    ):
        """The geometry rows of a species at a frequency level, by conformer.

        With a frequency fingerprint, the verified rows with it; without one,
        or when there are none, the unverified rows of the model chemistry
        (and code version if given). Rows with imaginary modes beyond the
        tolerance are left out unless `include_imaginary`. Any error reading
        the store gives [] and a log message.
        """
        db = self._reader()
        if db is None:
            return []
        where = "species = ? AND charge = ? AND multiplicity = ?"
        args = [species, int(charge), int(multiplicity)]
        if not include_imaginary:
            where += " AND (n_imaginary IS NULL OR n_imaginary = 0)"
        try:
            rows = []
            if frequency_fingerprint:
                rows = db.execute(
                    f"SELECT * FROM geometry WHERE {where} AND verified = 1 AND "
                    "frequency_fingerprint = ? ORDER BY conformer, id",
                    args + [frequency_fingerprint],
                ).fetchall()
            if not rows and model_chemistry:
                sql = f"{where} AND verified = 0 AND model_chemistry = ?"
                args_u = args + [model_chemistry]
                if version:
                    sql += " AND version = ?"
                    args_u.append(version)
                rows = db.execute(
                    f"SELECT * FROM geometry WHERE {sql} ORDER BY conformer, id",
                    args_u,
                ).fetchall()
        except sqlite3.OperationalError:  # a store from before the table existed
            return []
        except sqlite3.Error as e:
            logger.warning(f"The fragment store {self.path} could not be read: {e}")
            return []
        return [_decode(r, GeometryRecord) for r in rows]

    def energy(
        self,
        geometry_hash,
        *,
        energy_fingerprint=None,
        model_chemistry=None,
        version=None,
    ):
        """The energy row on a geometry at a level, or None.

        As :meth:`geometries`: the verified row with the fingerprint, else an
        unverified row of the model chemistry (and code version if given).
        """
        db = self._reader()
        if db is None:
            return None
        try:
            if energy_fingerprint:
                row = db.execute(
                    "SELECT * FROM energy WHERE geometry_hash = ? AND "
                    "energy_fingerprint = ? AND verified = 1",
                    (geometry_hash, energy_fingerprint),
                ).fetchone()
                if row is not None:
                    return _decode(row, EnergyRecord)
            if model_chemistry:
                sql = "geometry_hash = ? AND verified = 0 AND model_chemistry = ?"
                args = [geometry_hash, model_chemistry]
                if version:
                    sql += " AND version = ?"
                    args.append(version)
                row = db.execute(
                    f"SELECT * FROM energy WHERE {sql} ORDER BY id LIMIT 1", args
                ).fetchone()
                if row is not None:
                    return _decode(row, EnergyRecord)
        except sqlite3.OperationalError:
            return None
        except sqlite3.Error as e:
            logger.warning(f"The fragment store {self.path} could not be read: {e}")
            return None
        return None

    def energies(self, geometry_hash):
        """Every energy row on a geometry."""
        db = self._reader()
        if db is None:
            return []
        try:
            rows = db.execute(
                "SELECT * FROM energy WHERE geometry_hash = ? ORDER BY id",
                (geometry_hash,),
            ).fetchall()
        except sqlite3.Error:
            return []
        return [_decode(r, EnergyRecord) for r in rows]

    def species_at(
        self,
        species,
        *,
        charge=0,
        multiplicity=1,
        frequency_fingerprint=None,
        frequency_model_chemistry=None,
        energy_fingerprint=None,
        energy_model_chemistry=None,
    ):
        """The species at an energy level on a frequency level's geometry.

        Returns a list, one entry per conformer that has both, of
        {"geometry": GeometryRecord, "energy": EnergyRecord, "E": kJ/mol,
        "ZPE": unscaled kJ/mol, "H - E": kJ/mol at the geometry's T}. The
        energy level may be the frequency level itself (its E is then in the
        geometry row, and an energy row is used when there is one).
        """
        out = []
        for g in self.geometries(
            species,
            charge=charge,
            multiplicity=multiplicity,
            frequency_fingerprint=frequency_fingerprint,
            model_chemistry=frequency_model_chemistry,
        ):
            e = self.energy(
                g["geometry_hash"],
                energy_fingerprint=energy_fingerprint,
                model_chemistry=energy_model_chemistry,
            )
            if e is None:
                same = (
                    energy_fingerprint
                    and energy_fingerprint == g.get("optimization_fingerprint")
                ) or (
                    not energy_fingerprint
                    and energy_model_chemistry
                    and energy_model_chemistry == g.get("model_chemistry")
                )
                if not same or g["E"] is None:
                    continue
                E = g["E"]
            else:
                E = e["E"]
            out.append(
                {
                    "geometry": g,
                    "energy": e,
                    "E": E,
                    "ZPE": g["ZPE"],
                    "H - E": (
                        g["H"] - g["E"]
                        if g["H"] is not None and g["E"] is not None
                        else None
                    ),
                }
            )
        return out

    def history(self, table=None):
        """The row history (duplicates and superseded values), newest first."""
        db = self._reader()
        if db is None:
            return []
        try:
            sql = "SELECT * FROM row_history"
            args = ()
            if table is not None:
                sql += " WHERE table_name = ?"
                args = (table,)
            rows = db.execute(sql + " ORDER BY id DESC", args).fetchall()
        except sqlite3.Error:
            return []
        out = []
        for r in rows:
            d = dict(r)
            for c in ("old", "new"):
                d[c] = json.loads(d[c]) if d[c] else None
            out.append(d)
        return out
