# -*- coding: utf-8 -*-

"""The supplementary reference library: enthalpies of formation that ATcT
lacks, from evaluations, measurements, literature and computed values.

``references-official.db`` is published on Zenodo, installed by the installer
and only ever read. ``references.db`` is the local layer: user references and
computed references not yet published. Both sit beside the ATcT cache in
``<...>/Parameters/thermochemistry/``. The importers that build a library are
maintainer tools; users only install the published file.

The file keeps every value a source gives (``measurement``), the one chosen
per species and source (``reference``), and how derived values were derived
(``derivation``), so every choice can be audited and redone. Uncertainties are
stored as the source states them, with their level; the expanded (95 %)
uncertainty is computed from them by :func:`expanded_uncertainty`, so the
policy can change without re-importing anything.

Rules built in (the reference-data pilot, 2026-11-12):

* Burcat's database is never a source: a source named Burcat, or a value
  citing it, is refused.
* No P species from the online JANAF tables (their PH3 and PH2 sit ~17.5
  kJ/mol off the printed edition); P species come from the printed edition.
* Estimates are stored but never chosen unless accepted, which is recorded
  with who and when.
* A later source replaces an earlier one only when its uncertainty is below a
  third of the earlier one's and they agree within the earlier one's
  (:func:`later_source_wins`); both values are kept.
* Every derived value has its derivation, and :meth:`ReferenceLibrary.check`
  recomputes it from the current values of its partners.
* Ions are stored in ATcT's stationary-electron convention, converted on
  entry from the thermal-electron one.
"""

import datetime
import json
import logging
import math
import re
import sqlite3
from pathlib import Path

logger = logging.getLogger(__name__)

#: The schema version (``PRAGMA user_version``) of this code
SCHEMA_VERSION = 1

#: (5/2) R T at 298.15 K, kJ/mol: the thermal enthalpy of a free electron, the
#: difference between the two ion conventions at 298.15 K per unit of charge
ELECTRON_H298 = 2.5 * 8.314462618e-3 * 298.15

#: The floor (95 %) for a measured value given without any uncertainty, kJ/mol
UNSTATED_FLOOR = 4.0

#: The fewest species shared with ATcT for a source's offset to be used
OFFSET_MIN_SHARED = 5

SOURCE_KINDS = ("evaluated", "measured", "literature", "derived", "computed", "user")
STATUSES = ("measured", "review", "estimate", "derived", "held", "computed", "user")
UNC_LEVELS = ("2 sigma", "95 %", "1 sigma", "unstated", "none", "estimate", "mixed")
ION_CONVENTIONS = ("stationary-electron", "thermal-electron")

#: Sources whose values for species containing these elements are refused
EXCLUDED_ELEMENTS = {"JANAF-online": {"P"}}

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT
);
CREATE TABLE IF NOT EXISTS source (
    id INTEGER PRIMARY KEY,
    name TEXT UNIQUE NOT NULL,
    kind TEXT NOT NULL,
    rank INTEGER NOT NULL,
    citation TEXT,
    licence TEXT,
    licence_checked TEXT,
    retrieved TEXT,
    offset_kj REAL,
    offset_n INTEGER,
    offset_atct TEXT,
    offset_full_kj REAL,
    offset_full_n INTEGER,
    offset_excluded TEXT,
    notes TEXT
);
CREATE TABLE IF NOT EXISTS species (
    id INTEGER PRIMARY KEY,
    smiles TEXT NOT NULL,
    charge INTEGER NOT NULL,
    multiplicity INTEGER NOT NULL,
    formula TEXT NOT NULL,
    name TEXT,
    cas TEXT,
    inchi TEXT,
    inchi_retry TEXT,
    notes TEXT,
    UNIQUE (smiles, charge, multiplicity)
);
CREATE TABLE IF NOT EXISTS measurement (
    id INTEGER PRIMARY KEY,
    species INTEGER NOT NULL REFERENCES species (id),
    source INTEGER NOT NULL REFERENCES source (id),
    dfh_298k REAL,
    dfh_0k REAL,
    uncertainty REAL,
    unc_level TEXT NOT NULL,
    status TEXT NOT NULL,
    year INTEGER,
    h298_h0 REAL,
    standard_state TEXT,
    ion_convention TEXT,
    converted TEXT,
    citation TEXT NOT NULL,
    source_ref TEXT,
    method TEXT,
    comment TEXT,
    retrieved TEXT,
    notes TEXT,
    added TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS measurement_species ON measurement (species, source);
CREATE TABLE IF NOT EXISTS reference (
    id INTEGER PRIMARY KEY,
    species INTEGER NOT NULL REFERENCES species (id),
    source INTEGER NOT NULL REFERENCES source (id),
    measurement INTEGER NOT NULL REFERENCES measurement (id),
    state TEXT NOT NULL,
    rule TEXT NOT NULL,
    spread REAL,
    flag TEXT,
    overrides INTEGER REFERENCES reference (id),
    override_reason TEXT,
    decided_by TEXT,
    decided_on TEXT,
    decision_reason TEXT,
    UNIQUE (species, source)
);
CREATE TABLE IF NOT EXISTS calibration (
    id INTEGER PRIMARY KEY,
    source INTEGER NOT NULL REFERENCES source (id),
    smiles TEXT NOT NULL,
    dfh_298k REAL NOT NULL,
    atct_version TEXT NOT NULL,
    atct_id TEXT,
    atct_298k REAL NOT NULL,
    citation TEXT,
    UNIQUE (source, smiles, atct_version)
);
CREATE TABLE IF NOT EXISTS derivation (
    id INTEGER PRIMARY KEY,
    measurement INTEGER UNIQUE NOT NULL REFERENCES measurement (id),
    reaction TEXT NOT NULL,
    drh_298k REAL NOT NULL,
    drh_uncertainty REAL,
    citation TEXT,
    notes TEXT
);
"""


class ReferenceError(ValueError):
    """A value refused by the library's rules."""


class NewerReferenceSchemaError(PermissionError):
    """The library's schema is newer than this code's, so it is not written."""


_warned_newer = set()


def references_directory(base=None):
    """The directory of the reference libraries: beside the ATcT cache."""
    from .atct import atct_directory

    return atct_directory(base).parent


def official_references_path(path=None):
    """The published library, ``references-official.db``."""
    if path is not None:
        return Path(path).expanduser()
    return references_directory() / "references-official.db"


def local_references_path(path=None):
    """The local layer, ``references.db``."""
    if path is not None:
        return Path(path).expanduser()
    return references_directory() / "references.db"


def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def _mol(smiles):
    from .smiles import mol_from_smiles

    mol = mol_from_smiles(smiles)
    if mol is None:
        raise ReferenceError(f"RDKit cannot read the SMILES {smiles!r}")
    return mol


def hill_formula(smiles):
    """The Hill formula of a SMILES, without its charge."""
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    mol = _mol(smiles)
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:  # pragma: no cover
        pass
    return re.sub(r"[+-]\d*$", "", CalcMolFormula(mol))


def elements(smiles):
    """The set of element symbols in a SMILES (with implicit H)."""
    mol = _mol(smiles)
    result = {atom.GetSymbol() for atom in mol.GetAtoms()}
    try:
        mol.UpdatePropertyCache(strict=False)
        if any(atom.GetTotalNumHs() for atom in mol.GetAtoms()):
            result.add("H")
    except Exception:  # pragma: no cover
        pass
    return result


def _formula_key(formula):
    """A formula as {element: count}, to compare formulas written differently."""
    counts = {}
    for symbol, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        counts[symbol] = counts.get(symbol, 0) + (int(n) if n else 1)
    return counts


def expanded_uncertainty(uncertainty, unc_level, status="measured"):
    """The 95 % (k = 2) uncertainty from a stated one and its level.

    Science, 2026-11-12: 95 % and 2σ as stated; 1σ and unstated levels
    doubled; a measured value with no uncertainty gets a 4 kJ/mol floor; an
    estimate without one has none (None). "mixed" (W4-17: atoms at 2σ plus
    twice the method's 1σ) is already at 95 %.
    """
    if uncertainty is None:
        if status == "estimate":
            return None
        return UNSTATED_FLOOR
    if unc_level in ("2 sigma", "95 %", "mixed"):
        return float(uncertainty)
    if unc_level in ("1 sigma", "unstated", "estimate"):
        return 2.0 * float(uncertainty)
    if unc_level == "none":
        return UNSTATED_FLOOR
    raise ReferenceError(f"Unknown uncertainty level {unc_level!r}")


def later_source_wins(earlier, later):
    """Whether a later source's value replaces an earlier one's.

    Each argument is (value, expanded uncertainty). The later one wins only
    when its uncertainty is below a third of the earlier one's **and** the two
    agree within the earlier one's uncertainty.
    """
    (v1, u1), (v2, u2) = earlier, later
    if u1 is None or u2 is None:
        return False
    return u2 < u1 / 3.0 and abs(v2 - v1) <= u1


def select_measurement(
    measurements, choosable=("measured", "review", "derived", "user")
):
    """Choose one value among a source's measurements of one species.

    `measurements` are dicts with at least ``dfh_298k`` (or ``dfh_0k`` when
    that is all a source gives), ``uncertainty``,
    ``status``, ``year`` and ``id``. The rule (thermo's pilot, science
    2026-11-12): a review value; else the most recent with a stated
    uncertainty (ties to the smaller uncertainty); else the most recent.
    Values whose ``method`` is "N/A" (derived from another value) come last.
    Only values whose status is in `choosable` are chosen: never estimates
    or held values, and computed ones only from a computed source. The spread flag is
    raised when another measurement with an uncertainty differs from the
    chosen one by more than their combined uncertainties; the species is held
    when that exceeds three times the combined uncertainties and there is no
    review value.

    Returns (chosen, rule, spread, flag, state), or None if nothing can be
    chosen. `state` is "chosen" or "held".
    """

    def value(m):
        return m["dfh_298k"] if m.get("dfh_298k") is not None else m.get("dfh_0k")

    usable = [
        m for m in measurements if value(m) is not None and m.get("status") in choosable
    ]
    if not usable:
        return None
    reviews = [m for m in usable if m["status"] == "review"]
    # A value the WebBook marks "N/A" is derived from another on the page
    # (e.g. a condensed-phase value plus a sublimation enthalpy): last resort
    direct = [m for m in usable if m.get("method") != "N/A"] or usable
    with_u = [m for m in direct if m.get("uncertainty") is not None]

    def recent(m):
        return m.get("year") or 0

    if reviews:
        chosen = max(reviews, key=recent)
        rule = "review"
    elif with_u:
        chosen = max(with_u, key=lambda m: (recent(m), -m["uncertainty"]))
        rule = "most recent with an uncertainty"
    else:
        chosen = max(direct, key=recent)
        rule = "most recent"
    if len(usable) == 1:
        rule = "single value"
    if chosen.get("dfh_298k") is None:
        rule += "; 0 K only: used only when no source has a 298.15 K value"

    spread = max(abs(value(m) - value(chosen)) for m in usable)
    flag = None
    state = "chosen"
    worst = 0.0
    outlier = None
    if chosen.get("uncertainty") is not None:
        for m in with_u:
            if m is chosen:
                continue
            combined = math.hypot(chosen["uncertainty"], m["uncertainty"])
            if combined <= 0:
                continue
            ratio = abs(value(m) - value(chosen)) / combined
            if ratio > worst:
                worst, outlier = ratio, m
    if worst > 1.0:
        flag = (
            f"spread: {value(outlier):.2f} "
            f"({outlier.get('citation') or 'measurement ' + str(outlier['id'])}) "
            f"differs by {worst:.1f}x the combined uncertainties"
        )
        if worst > 3.0 and not reviews:
            state = "held"
    return chosen, rule, spread, flag, state


class ReferenceLibrary:
    """One reference-library file (a layer).

    Parameters
    ----------
    path : str or Path
        The SQLite file.
    create : bool
        Create the file and schema if missing (a writer).
    layer : str
        For a new file: "official" or "local", recorded in ``meta``.
    immutable : bool
        Open read-only without locking (the installed official file).
    """

    def __init__(self, path, *, create=False, layer="local", immutable=False):
        self.path = Path(path).expanduser()
        self._db = None
        self._create = create
        self._immutable = immutable and not create
        self._new_layer = layer
        #: optional callable(smiles, charge, multiplicity) -> bool: whether
        #: ATcT has a species, so a derivation partner it has adds no offset
        self.in_atct = None
        if create:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._connect()

    # -- connection ---------------------------------------------------------

    def _connect(self):
        if self._db is not None:
            return self._db
        if self._create:
            db = sqlite3.connect(self.path, timeout=30)
        else:
            if not self.path.exists():
                raise FileNotFoundError(f"No reference library at {self.path}")
            flags = "mode=ro&immutable=1" if self._immutable else "mode=ro"
            db = sqlite3.connect(f"file:{self.path}?{flags}", uri=True, timeout=30)
        db.row_factory = sqlite3.Row
        db.execute("PRAGMA mmap_size = 0")
        db.execute("PRAGMA busy_timeout = 30000")
        version = db.execute("PRAGMA user_version").fetchone()[0]
        if version > SCHEMA_VERSION:
            text = (
                f"The reference library {self.path} is at schema version "
                f"{version}, newer than this seamm-thermochemistry understands "
                f"({SCHEMA_VERSION})"
            )
            if self._create:
                db.close()
                raise NewerReferenceSchemaError(
                    text + ", so it is not written to. Update seamm-thermochemistry."
                )
            if str(self.path) not in _warned_newer:
                _warned_newer.add(str(self.path))
                logger.warning(text + ": it is read as the older schema.")
        if self._create:
            db.execute("PRAGMA journal_mode = DELETE")
            db.execute("PRAGMA foreign_keys = ON")
            db.executescript(_SCHEMA)
            if version < SCHEMA_VERSION:
                db.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
            db.execute(
                "INSERT OR IGNORE INTO meta (key, value) VALUES ('layer', ?)",
                (self._new_layer,),
            )
            db.execute(
                "INSERT OR IGNORE INTO meta (key, value) VALUES ('created', ?)",
                (_now(),),
            )
            db.commit()
        self._db = db
        return db

    def _writer(self):
        if not self._create:
            raise PermissionError(f"{self.path} was opened read-only")
        return self._connect()

    def close(self):
        if self._db is not None:
            self._db.close()
            self._db = None

    @property
    def meta(self):
        db = self._connect()
        return {r["key"]: r["value"] for r in db.execute("SELECT * FROM meta")}

    @property
    def layer(self):
        return self.meta.get("layer", "local")

    @property
    def version(self):
        return self.meta.get("version")

    def set_meta(self, key, value):
        db = self._writer()
        db.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)", (key, value)
        )
        db.commit()

    # -- sources --------------------------------------------------------------

    def add_source(
        self,
        name,
        kind,
        rank,
        *,
        citation=None,
        licence=None,
        licence_checked=None,
        retrieved=None,
        notes=None,
    ):
        """Add a source, or return the id of the one with this name."""
        if "burcat" in name.lower():
            raise ReferenceError(
                f"{name}: Burcat's database is not a source of the library "
                "(Paul, 2026-10-29). Cite the evaluation it transcribes instead."
            )
        if kind not in SOURCE_KINDS:
            raise ReferenceError(f"Unknown source kind {kind!r}: one of {SOURCE_KINDS}")
        db = self._writer()
        row = db.execute("SELECT id FROM source WHERE name = ?", (name,)).fetchone()
        if row is not None:
            return row["id"]
        cursor = db.execute(
            "INSERT INTO source (name, kind, rank, citation, licence, "
            "licence_checked, retrieved, notes) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (name, kind, rank, citation, licence, licence_checked, retrieved, notes),
        )
        db.commit()
        return cursor.lastrowid

    def source(self, name_or_id):
        db = self._connect()
        column = "id" if isinstance(name_or_id, int) else "name"
        row = db.execute(
            f"SELECT * FROM source WHERE {column} = ?", (name_or_id,)
        ).fetchone()
        if row is None:
            raise KeyError(f"No source {name_or_id!r} in {self.path}")
        return dict(row)

    def sources(self):
        db = self._connect()
        return [dict(r) for r in db.execute("SELECT * FROM source ORDER BY rank")]

    def set_offset(
        self, source, offset_kj, n, atct_version, *, full=None, excluded=None
    ):
        """Record a source's offset against ATcT: the rms of (source - ATcT)
        over the `n` species both have. Fewer than OFFSET_MIN_SHARED species
        record no offset (None), only the count. `full` is (rms, n) over every
        shared species when some were `excluded` (a recorded reason)."""
        db = self._writer()
        sid = self.source(source)["id"]
        value = offset_kj if n >= OFFSET_MIN_SHARED else None
        full_kj, full_n = full if full is not None else (offset_kj, n)
        db.execute(
            "UPDATE source SET offset_kj = ?, offset_n = ?, offset_atct = ?, "
            "offset_full_kj = ?, offset_full_n = ?, offset_excluded = ? "
            "WHERE id = ?",
            (value, n, atct_version, full_kj, full_n, excluded, sid),
        )
        db.commit()

    def add_calibration(self, path, *, source_names=None):
        """Read calibration values: a source's values for species ATcT has,
        used only to measure the source's offset, never as references.

        The CSV has smiles, source, dfH298_kJmol, ATcT_version, ATcT_id,
        ATcT_dfH298_kJmol and citation. `source_names` maps the file's source
        names to the library's. Returns the number of rows read."""
        import csv

        names = {"Pedley-TRC (NIST WebBook)": "WebBook-SRD69"}
        names.update(source_names or {})
        db = self._writer()
        n = 0
        with open(path, newline="", encoding="utf-8") as fd:
            for row in csv.DictReader(fd):
                name = names.get(row["source"], row["source"])
                sid = self.source(name)["id"]
                db.execute(
                    "INSERT OR REPLACE INTO calibration (source, smiles, dfh_298k, "
                    "atct_version, atct_id, atct_298k, citation) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (
                        sid,
                        row["smiles"],
                        float(row["dfH298_kJmol"]),
                        row["ATcT_version"],
                        row.get("ATcT_id"),
                        float(row["ATcT_dfH298_kJmol"]),
                        row.get("citation"),
                    ),
                )
                n += 1
        db.commit()
        return n

    # -- species --------------------------------------------------------------

    def add_species(
        self,
        smiles,
        charge=0,
        multiplicity=1,
        *,
        formula=None,
        name=None,
        cas=None,
        inchi=None,
        inchi_retry=None,
        notes=None,
    ):
        """Add a species, or return the id of the existing one.

        The SMILES is canonicalized as the CBH code keys species. A given
        formula must match the SMILES', and a CAS number already used by a
        species of another constitution is refused.
        """
        from .atct import canonical_smiles

        can = canonical_smiles(smiles)
        if can is None:
            raise ReferenceError(f"RDKit cannot read the SMILES {smiles!r}")
        hill = hill_formula(can)
        if formula is not None and _formula_key(formula) != _formula_key(hill):
            raise ReferenceError(
                f"{smiles}: the formula {formula} does not match the SMILES ({hill})"
            )
        db = self._writer()
        if cas:
            for row in db.execute(
                "SELECT smiles FROM species WHERE cas = ? AND smiles != ?", (cas, can)
            ):
                if canonical_smiles(row["smiles"], stereo=False) != canonical_smiles(
                    can, stereo=False
                ):
                    raise ReferenceError(
                        f"{smiles}: CAS {cas} already belongs to {row['smiles']}"
                    )
        row = db.execute(
            "SELECT id FROM species WHERE smiles = ? AND charge = ? "
            "AND multiplicity = ?",
            (can, charge, multiplicity),
        ).fetchone()
        if row is not None:
            return row["id"]
        cursor = db.execute(
            "INSERT INTO species (smiles, charge, multiplicity, formula, name, cas, "
            "inchi, inchi_retry, notes) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (can, charge, multiplicity, hill, name, cas, inchi, inchi_retry, notes),
        )
        db.commit()
        return cursor.lastrowid

    def species(self, smiles, charge=0, multiplicity=None):
        """The species rows matching a SMILES (any multiplicity if None)."""
        from .atct import canonical_smiles

        can = canonical_smiles(smiles)
        db = self._connect()
        sql = "SELECT * FROM species WHERE smiles = ? AND charge = ?"
        args = [can, charge]
        if multiplicity is not None:
            sql += " AND multiplicity = ?"
            args.append(multiplicity)
        return [dict(r) for r in db.execute(sql, args)]

    def _species_row(self, species_id):
        row = (
            self._connect()
            .execute("SELECT * FROM species WHERE id = ?", (species_id,))
            .fetchone()
        )
        if row is None:
            raise KeyError(f"No species {species_id} in {self.path}")
        return dict(row)

    # -- measurements -----------------------------------------------------------

    def add_measurement(
        self,
        species_id,
        source,
        *,
        citation,
        dfh_298k=None,
        dfh_0k=None,
        uncertainty=None,
        unc_level="unstated",
        status="measured",
        year=None,
        h298_h0=None,
        standard_state=None,
        ion_convention=None,
        source_ref=None,
        method=None,
        comment=None,
        retrieved=None,
        notes=None,
    ):
        """Add one value a source gives for a species. Returns its id."""
        species = self._species_row(species_id)
        src = self.source(source)
        if dfh_298k is None and dfh_0k is None:
            raise ReferenceError(f"{species['smiles']}: no enthalpy of formation given")
        if status not in STATUSES:
            raise ReferenceError(f"Unknown status {status!r}: one of {STATUSES}")
        if unc_level not in UNC_LEVELS:
            raise ReferenceError(
                f"Unknown uncertainty level {unc_level!r}: one of {UNC_LEVELS}"
            )
        if not citation:
            raise ReferenceError(f"{species['smiles']}: every value needs a citation")
        if "burcat" in citation.lower():
            raise ReferenceError(
                f"{species['smiles']}: a value citing Burcat's database is refused; "
                "cite the evaluation it transcribes instead."
            )
        excluded = EXCLUDED_ELEMENTS.get(src["name"], set()) & elements(
            species["smiles"]
        )
        if excluded:
            raise ReferenceError(
                f"{species['smiles']}: values for species with "
                f"{', '.join(sorted(excluded))} are not taken from {src['name']} "
                "(the online JANAF P tables are off the printed edition); use the "
                "printed edition."
            )
        converted = None
        if species["charge"] != 0:
            if ion_convention not in ION_CONVENTIONS:
                raise ReferenceError(
                    f"{species['smiles']}: an ion needs its electron convention, "
                    f"one of {ION_CONVENTIONS}"
                )
            if ion_convention == "thermal-electron" and dfh_298k is not None:
                shift = species["charge"] * ELECTRON_H298
                converted = (
                    f"thermal-electron {dfh_298k:.3f} -> stationary-electron "
                    f"{dfh_298k - shift:.3f} kJ/mol at 298.15 K"
                )
                dfh_298k -= shift
        db = self._writer()
        cursor = db.execute(
            "INSERT INTO measurement (species, source, dfh_298k, dfh_0k, "
            "uncertainty, unc_level, status, year, h298_h0, standard_state, "
            "ion_convention, converted, citation, source_ref, method, comment, "
            "retrieved, notes, added) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                species_id,
                src["id"],
                dfh_298k,
                dfh_0k,
                uncertainty,
                unc_level,
                status,
                year,
                h298_h0,
                standard_state,
                ion_convention,
                converted,
                citation,
                source_ref,
                method,
                comment,
                retrieved,
                notes,
                _now(),
            ),
        )
        db.commit()
        return cursor.lastrowid

    def measurements(self, species_id, source=None):
        db = self._connect()
        sql = (
            "SELECT m.*, s.name AS source_name FROM measurement m "
            "JOIN source s ON s.id = m.source WHERE m.species = ?"
        )
        args = [species_id]
        if source is not None:
            sql += " AND s.name = ?"
            args.append(source)
        return [dict(r) for r in db.execute(sql + " ORDER BY m.id", args)]

    def _measurement(self, measurement_id):
        row = (
            self._connect()
            .execute("SELECT * FROM measurement WHERE id = ?", (measurement_id,))
            .fetchone()
        )
        if row is None:
            raise KeyError(f"No measurement {measurement_id} in {self.path}")
        return dict(row)

    # -- choices ----------------------------------------------------------------

    def choose(self, species_id, source):
        """Choose the value of a source for a species by
        :func:`select_measurement`, and record it. Returns the reference id,
        or None if nothing can be chosen."""
        choosable = ("measured", "review", "derived", "user")
        if self.source(source)["kind"] in ("computed", "literature"):
            choosable += ("computed",)
        result = select_measurement(self.measurements(species_id, source), choosable)
        if result is None:
            return None
        chosen, rule, spread, flag, state = result
        return self._set_reference(
            species_id, source, chosen["id"], rule, spread, flag, state
        )

    def _set_reference(
        self, species_id, source, measurement_id, rule, spread, flag, state
    ):
        db = self._writer()
        sid = self.source(source)["id"]
        db.execute(
            "INSERT INTO reference (species, source, measurement, state, rule, "
            "spread, flag) VALUES (?, ?, ?, ?, ?, ?, ?) "
            "ON CONFLICT (species, source) DO UPDATE SET measurement = "
            "excluded.measurement, state = excluded.state, rule = excluded.rule, "
            "spread = excluded.spread, flag = excluded.flag, decided_by = NULL, "
            "decided_on = NULL, decision_reason = NULL",
            (species_id, sid, measurement_id, state, rule, spread, flag),
        )
        db.commit()
        return db.execute(
            "SELECT id FROM reference WHERE species = ? AND source = ?",
            (species_id, sid),
        ).fetchone()["id"]

    def accept(self, measurement_id, *, by, reason):
        """A person accepts a value the rules would not choose: an estimate,
        or a species held for its spread. Recorded with who, when and why."""
        if not by or not reason:
            raise ReferenceError("An accept needs who accepted it and why")
        m = self._measurement(measurement_id)
        if m["dfh_298k"] is None and m["dfh_0k"] is None:  # pragma: no cover
            raise ReferenceError("Nothing to accept")
        source = self.source(m["source"])["name"]
        ref_id = self._set_reference(
            m["species"], source, measurement_id, "accepted", None, None, "chosen"
        )
        db = self._writer()
        db.execute(
            "UPDATE reference SET decided_by = ?, decided_on = ?, "
            "decision_reason = ? WHERE id = ?",
            (by, _now(), reason, ref_id),
        )
        db.commit()
        return ref_id

    def accept_value(self, smiles, charge, source, dfh_298k, *, by, reason):
        """Accept a source's value for a species, found by its value at
        298.15 K (within 0.01 kJ/mol), so that an accept can be kept in a
        file and reapplied whenever the library is rebuilt."""
        matches = []
        for row in self.species(smiles, charge):
            for m in self.measurements(row["id"], source):
                if m["dfh_298k"] is not None and abs(m["dfh_298k"] - dfh_298k) <= 0.01:
                    matches.append(m)
        if len(matches) != 1:
            raise ReferenceError(
                f"{smiles} ({source}): {len(matches)} values of {dfh_298k} kJ/mol "
                "to accept, not one"
            )
        return self.accept(matches[0]["id"], by=by, reason=reason)

    def hold(self, smiles, charge, source, *, by, reason):
        """A person holds a source's value for a species (not used), with
        why: e.g. two modern measurements in conflict, awaiting a computed
        arbitration. Recorded with who and when.

        `source` may name several sources separated by "/" (e.g.
        "JANAF-1998 / NBS-82": the two disagree, so neither is used). Each
        must be a source in the library; a source with no value for the
        species, or whose values its own file already holds, is passed over.
        Returns the number of values held here."""
        if not by or not reason:
            raise ReferenceError("A hold needs who decided it and why")
        names = [name.strip() for name in source.split("/") if name.strip()]
        sids = [self.source(name)["id"] for name in names]
        rows = self.species(smiles, charge)
        db = self._writer()
        n = 0
        found = 0
        for sid in sids:
            for row in rows:
                k = db.execute(
                    "UPDATE reference SET state = 'held', decided_by = ?, "
                    "decided_on = ?, decision_reason = ? "
                    "WHERE species = ? AND source = ?",
                    (by, _now(), reason, row["id"], sid),
                ).rowcount
                n += k
                if k == 0:
                    # Already held by its source file: nothing was chosen
                    k = db.execute(
                        "SELECT COUNT(*) FROM measurement WHERE species = ? "
                        "AND source = ?",
                        (row["id"], sid),
                    ).fetchone()[0]
                found += min(k, 1)
        db.commit()
        if found == 0 or (len(names) == 1 and n > 1):
            raise ReferenceError(f"{smiles} ({source}): {n} values to hold, not one")
        return n

    def accept_file(self, path):
        """Apply the decisions kept in a CSV: smiles, charge, source,
        dfh_298k, by, reason, and optionally action ("accept", the default,
        or "hold", which needs no value). Returns the reference ids of the
        accepts."""
        import csv

        result = []
        with open(path, newline="", encoding="utf-8") as fd:
            for row in csv.DictReader(fd):
                if (row.get("action") or "accept").strip() == "hold":
                    self.hold(
                        row["smiles"],
                        int(row.get("charge") or 0),
                        row["source"],
                        by=row["by"],
                        reason=row["reason"],
                    )
                    continue
                result.append(
                    self.accept_value(
                        row["smiles"],
                        int(row.get("charge") or 0),
                        row["source"],
                        float(row["dfh_298k"]),
                        by=row["by"],
                        reason=row["reason"],
                    )
                )
        return result

    def references(self, species_id, *, include_held=False):
        """The chosen values of every source for a species, best rank first,
        with their expanded uncertainty."""
        db = self._connect()
        rows = db.execute(
            "SELECT r.id AS reference_id, r.state, r.rule, r.spread, r.flag, "
            "r.overrides, r.override_reason, r.decided_by, r.decided_on, "
            "r.decision_reason, m.*, s.name AS source_name, s.kind AS source_kind, "
            "s.rank AS source_rank, s.offset_kj FROM reference r "
            "JOIN measurement m ON m.id = r.measurement "
            "JOIN source s ON s.id = r.source "
            "WHERE r.species = ? ORDER BY s.rank, s.id",
            (species_id,),
        ).fetchall()
        result = []
        for row in rows:
            row = dict(row)
            if row["state"] != "chosen" and not include_held:
                continue
            if row["status"] == "estimate" and row["decided_by"] is None:
                continue
            row["expanded_uncertainty"] = expanded_uncertainty(
                row["uncertainty"], row["unc_level"], row["status"]
            )
            row["offset_from"] = row["source_name"] if row["offset_kj"] else None
            if row["offset_kj"] is None:
                inherited = self.inherited_offset(row["id"])
                if inherited is not None:
                    row["offset_kj"], row["offset_from"] = inherited
            row["total_uncertainty"] = (
                None
                if row["expanded_uncertainty"] is None
                else math.hypot(row["expanded_uncertainty"], row["offset_kj"] or 0.0)
            )
            result.append(row)
        return result

    def inherited_offset(self, measurement_id):
        """The offset a value derived from others inherits (science,
        2026-11-12): the partners' source offsets propagated through the
        derivation, sqrt(Σ (ν_i offset_i)²). A partner ATcT has (see
        :attr:`in_atct`) is resolved from ATcT and adds nothing. Returns
        (offset, "from ...") or
        None when there is no derivation or no partner source has an offset."""
        db = self._connect()
        d = db.execute(
            "SELECT reaction FROM derivation WHERE measurement = ?", (measurement_id,)
        ).fetchone()
        if d is None:
            return None
        total = 0.0
        names = []
        for p in json.loads(d["reaction"]):
            if p.get("target") or p.get("element"):
                continue
            if self.in_atct is not None and self.in_atct(
                p["smiles"], p.get("charge", 0), p.get("multiplicity")
            ):
                continue  # resolved from ATcT: no offset
            try:
                offset = self.source(p["source"])["offset_kj"]
            except KeyError:
                offset = None
            if offset:
                total += (p["coefficient"] * offset) ** 2
                names.append(p["source"])
        if not names:
            return None
        return math.sqrt(total), "inherited from " + ", ".join(sorted(set(names)))

    def apply_overrides(self):
        """Mark, for every species, which source's value wins by
        :func:`later_source_wins`, going down the ranks. Both values stay;
        the winner records the one it overrides and why."""
        db = self._writer()
        db.execute("UPDATE reference SET overrides = NULL, override_reason = NULL")
        changed = 0
        for (species_id,) in db.execute("SELECT DISTINCT species FROM reference"):
            refs = [
                r
                for r in self.references(species_id)
                if r["dfh_298k"] is not None and r["expanded_uncertainty"] is not None
            ]
            if not refs:
                continue
            best = refs[0]
            for later in refs[1:]:
                if later_source_wins(
                    (best["dfh_298k"], best["expanded_uncertainty"]),
                    (later["dfh_298k"], later["expanded_uncertainty"]),
                ):
                    reason = (
                        f"{later['source_name']} ± {later['expanded_uncertainty']:.2f} "
                        f"is below a third of {best['source_name']}'s ± "
                        f"{best['expanded_uncertainty']:.2f} and agrees within it "
                        f"({later['dfh_298k'] - best['dfh_298k']:+.2f} kJ/mol)"
                    )
                    db.execute(
                        "UPDATE reference SET overrides = ?, override_reason = ? "
                        "WHERE id = ?",
                        (best["reference_id"], reason, later["reference_id"]),
                    )
                    changed += 1
                    best = later
        db.commit()
        return changed

    def best(self, species_id):
        """The value to use for a species: the best-ranked chosen value,
        unless a later source overrides it. A value at 0 K only (W4-17) is
        used only when no source gives one at 298.15 K. None if there is
        none."""
        refs = self.references(species_id)
        if not refs:
            return None
        at_298 = [r for r in refs if r["dfh_298k"] is not None]
        refs = at_298 or refs
        by_id = {r["reference_id"]: r for r in refs}
        best = refs[0]
        while True:
            later = [r for r in refs if r["overrides"] == best["reference_id"]]
            if not later:
                return best
            best = by_id[later[0]["reference_id"]]

    # -- derivations --------------------------------------------------------------

    def add_derivation(
        self,
        measurement_id,
        reaction,
        drh_298k,
        *,
        drh_uncertainty=None,
        citation=None,
        notes=None,
    ):
        """Record how a value was derived from others: every derived value
        has one, and so may a value a source derived (e.g. ΔfH from a D0 and
        the atoms), so that its offset is inherited and it can be checked.

        `reaction` is a list of dicts, one per species, each with its
        ``coefficient`` (products positive): the derived species has
        ``"target": True``; a partner gives ``smiles`` (with ``charge`` and
        ``multiplicity``) and the ``source`` its value comes from, or
        ``"element": True`` for an element in its reference state (ΔfH = 0).
        """
        self._measurement(measurement_id)
        targets = [p for p in reaction if p.get("target")]
        if len(targets) != 1:
            raise ReferenceError("A derivation needs exactly one target")
        for p in reaction:
            if "coefficient" not in p:
                raise ReferenceError(f"No coefficient in {p}")
            if not (p.get("target") or p.get("element") or p.get("source")):
                raise ReferenceError(f"A partner needs its source: {p}")
        db = self._writer()
        cursor = db.execute(
            "INSERT INTO derivation (measurement, reaction, drh_298k, "
            "drh_uncertainty, citation, notes) VALUES (?, ?, ?, ?, ?, ?)",
            (
                measurement_id,
                json.dumps(reaction),
                drh_298k,
                drh_uncertainty,
                citation,
                notes,
            ),
        )
        db.commit()
        return cursor.lastrowid

    def _partner_value(self, partner):
        for row in self.species(
            partner["smiles"], partner.get("charge", 0), partner.get("multiplicity")
        ):
            for ref in self.references(row["id"]):
                if ref["source_name"] == partner["source"]:
                    return ref["dfh_298k"]
        return None

    def check(self):
        """Check every value with a derivation against its partners' current
        values, and that every derived value has a derivation.

        Returns a list of problems (dicts with ``measurement``, ``smiles``
        and ``problem``); empty when all is well."""
        db = self._connect()
        problems = []
        for m in db.execute(
            "SELECT m.*, sp.smiles FROM measurement m "
            "JOIN species sp ON sp.id = m.species WHERE m.status = 'derived' "
            "OR m.id IN (SELECT measurement FROM derivation)"
        ):
            m = dict(m)
            d = db.execute(
                "SELECT * FROM derivation WHERE measurement = ?", (m["id"],)
            ).fetchone()
            if d is None:
                problems.append(
                    {
                        "measurement": m["id"],
                        "smiles": m["smiles"],
                        "problem": "no derivation",
                    }
                )
                continue
            reaction = json.loads(d["reaction"])
            total = d["drh_298k"]
            target = None
            missing = []
            for p in reaction:
                if p.get("target"):
                    target = p
                    continue
                if p.get("element"):
                    continue
                value = self._partner_value(p)
                if value is None:
                    missing.append(f"{p['smiles']} ({p['source']})")
                    continue
                total -= p["coefficient"] * value
            if missing:
                problems.append(
                    {
                        "measurement": m["id"],
                        "smiles": m["smiles"],
                        "problem": "no value for " + ", ".join(missing),
                    }
                )
                continue
            recomputed = total / target["coefficient"]
            difference = recomputed - m["dfh_298k"]
            allowed = m["uncertainty"] if m["uncertainty"] is not None else 0.0
            if abs(difference) > max(allowed, 0.01):
                problems.append(
                    {
                        "measurement": m["id"],
                        "smiles": m["smiles"],
                        "problem": (
                            f"stored {m['dfh_298k']:.2f}, its derivation gives "
                            f"{recomputed:.2f} kJ/mol ({difference:+.2f}, more than "
                            f"its uncertainty {allowed:.2f})"
                        ),
                        "recomputed": recomputed,
                    }
                )
        return problems


class LayeredReferences:
    """The local and official reference libraries, read as one.

    For each species the values of both layers are merged, best rank first;
    within a rank the local layer comes first. Every value carries its
    ``layer``.
    """

    def __init__(self, local_path=None, official_path=None):
        self.layers = []
        official = official_references_path(official_path)
        local = local_references_path(local_path)
        if local.exists():
            library = ReferenceLibrary(local)
            if library.layer == "official":
                raise ValueError(
                    f"{local} is an official reference library, not a local one."
                )
            self.layers.append(("local", library))
        if official.exists():
            self.layers.append(("official", ReferenceLibrary(official, immutable=True)))

    def set_in_atct(self, function):
        """Tell every layer how to know whether ATcT has a species."""
        for _, library in self.layers:
            library.in_atct = function

    @property
    def versions(self):
        """{layer: (version, path)} of the libraries read."""
        return {name: (lib.version, str(lib.path)) for name, lib in self.layers}

    def candidates(self, smiles, charge=0, multiplicity=None, *, include_held=False):
        """Every usable value for a species, from both layers (and the held
        ones too, if asked, never preferred)."""
        result = []
        for name, library in self.layers:
            for row in library.species(smiles, charge, multiplicity):
                best = library.best(row["id"])
                for ref in library.references(row["id"], include_held=include_held):
                    ref["layer"] = name
                    ref["preferred"] = (
                        best is not None and ref["reference_id"] == best["reference_id"]
                    )
                    ref["smiles"] = row["smiles"]
                    ref["formula"] = row["formula"]
                    result.append(ref)
        order = {"local": 0, "official": 1}
        result.sort(key=lambda r: (r["source_rank"], order[r["layer"]]))
        return result

    def lookup(self, smiles, charge=0, multiplicity=None):
        """The value to use: each layer's preferred value (overrides applied
        within the layer), the best-ranked of them, the local layer winning a
        tie. None if neither layer has the species."""
        for ref in self.candidates(smiles, charge, multiplicity):
            if ref["preferred"]:
                return ref
        return None

    def close(self):
        for _, library in self.layers:
            library.close()


def compute_offsets(library, atct=None, *, exclude=None):
    """Measure each source's offset against ATcT and record it.

    The differences source - ATcT come from two places: the library's own
    values for species ATcT also has (with `atct`), and the calibration
    values (:meth:`ReferenceLibrary.add_calibration`), which are used for
    nothing else. A source's offset is the rms of its differences over at
    least OFFSET_MIN_SHARED species (science, 2026-11-12), else none. It is
    added in quadrature to the source's uncertainties.

    Caveat (science): the calibration species are well studied, so the
    offset is a *lower bound* for the obscure species the library serves.
    It is not a measure of their reliability; the stated uncertainties and
    the spread rule are.

    `exclude` is {source: (set of SMILES, reason)}: species left out of the
    offset for a recorded reason (e.g. CODATA key values known to be
    superseded); the offset over every species is kept beside it.

    Returns {source: (rms, n, [(smiles, diff)])}, over the species used.
    """
    from .atct import ATcTError

    db = library._connect()
    result = {}
    for source in library.sources():
        diffs = []
        versions = set()
        if atct is None:
            rows = []
        else:
            rows = db.execute(
                "SELECT sp.* FROM reference r JOIN species sp ON sp.id = r.species "
                "WHERE r.source = ? AND r.state = 'chosen'",
                (source["id"],),
            ).fetchall()
        for row in rows:
            mine = [
                r
                for r in library.references(row["id"])
                if r["source_name"] == source["name"]
            ]
            if not mine or mine[0]["dfh_298k"] is None:
                continue
            try:
                entry = atct.lookup(
                    row["smiles"],
                    charge=row["charge"],
                    multiplicity=row["multiplicity"],
                )
            except (ATcTError, KeyError):
                continue
            diffs.append((row["smiles"], mine[0]["dfh_298k"] - entry["dfH0_298K"]))
            versions.add(atct.version)
        seen = {smiles for smiles, _ in diffs}
        for row in db.execute(
            "SELECT * FROM calibration WHERE source = ?", (source["id"],)
        ):
            if row["smiles"] in seen:
                continue
            versions.add(row["atct_version"])
            diffs.append((row["smiles"], row["dfh_298k"] - row["atct_298k"]))

        def rms_of(values):
            return (
                math.sqrt(sum(d * d for _, d in values) / len(values))
                if values
                else None
            )

        full = (rms_of(diffs), len(diffs))
        excluded = None
        if exclude and source["name"] in exclude:
            from .atct import canonical_smiles

            smiles_set, reason = exclude[source["name"]]
            drop = {canonical_smiles(x) for x in smiles_set}
            kept = [(x, d) for x, d in diffs if canonical_smiles(x) not in drop]
            if len(kept) != len(diffs):
                names = [x for x, _ in diffs if canonical_smiles(x) in drop]
                excluded = f"{reason} ({', '.join(names)})"
            diffs = kept
        n = len(diffs)
        rms = rms_of(diffs)
        library.set_offset(
            source["name"],
            rms,
            n,
            ", ".join(sorted(versions)),
            full=full if excluded else None,
            excluded=excluded,
        )
        result[source["name"]] = (rms, n, diffs)
    return result
