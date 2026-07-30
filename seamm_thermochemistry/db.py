"""SQLite-backed store for atomic reference energies.

Every SEAMM code that reports a formation-referenced energy (Gaussian, Psi4,
ORCA, VASP, ...) needs the same two pieces of data per element: an
experimental reference (heat of formation, entropy, standard-state
description) and one or more *computed* atomic reference energies (one per
code/method/settings combination). Today each plugin carries its own copy of
both as a wide, mostly-empty CSV. This module replaces that with a single,
shared, relational store -- one `element` row per element, one `atom_energy`
row per (element, code, method, ref_type, settings) combination -- with room
for the provenance a shared reference dataset needs (who computed it, when,
with what).

Two reference conventions are both first-class (`ref_type`):

- ``"atom"``: energy of the isolated, gas-phase atom. This is the
  Gaussian/Psi4/ORCA convention, and the one to use for a physically
  anchored, cross-code-comparable energy of formation.
- ``"element_phase"``: energy per atom of the element's standard-state phase
  (bulk metal, graphite, O2(g), ...). This is VASP's existing convention
  (``element_energies.csv``) -- no experimental anchor needed, but the
  result is referenced to the computed elemental phases, not real formation
  energies. It also serves as a fallback for elements (e.g. Mn) where the
  free atom is a poor DFT reference.

See ``formation.py`` for the arithmetic that consumes this store.
"""

import configparser
import csv
from contextlib import contextmanager
import sqlite3
from pathlib import Path

__all__ = ["ThermoDB", "DEFAULT_DB_PATH"]

_SCHEMA_VERSION = 1

_SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS schema_info (
    key   TEXT PRIMARY KEY,
    value TEXT
);

CREATE TABLE IF NOT EXISTS element (
    atomic_number            INTEGER PRIMARY KEY,
    symbol                   TEXT NOT NULL UNIQUE,
    multiplicity             INTEGER,
    term_symbol               TEXT,
    standard_state            TEXT,      -- e.g. "1/2 H2(g)", "C(s,gr)"
    dfH0_0K                   REAL,      -- experimental DfH(X, 0 K), kJ/mol
    dfH0_298K                 REAL,      -- experimental DfH(X, 298 K, gas), kJ/mol
    dfH0_298K_stderr          REAL,
    h298_minus_h0_atom        REAL,      -- kJ/mol, gas atom H(298K) - H(0K)
    h298_minus_h0_std_state   REAL,      -- kJ/mol, standard-state phase, per atom
    s298_gas                  REAL,      -- J/(mol K), gas atom
    s298_gas_stderr           REAL,
    reference                 TEXT,      -- source URL
    reference_note             TEXT      -- e.g. "JANAF"
);

CREATE TABLE IF NOT EXISTS atom_energy (
    id                 INTEGER PRIMARY KEY AUTOINCREMENT,
    atomic_number      INTEGER NOT NULL REFERENCES element(atomic_number),
    code               TEXT NOT NULL,               -- "gaussian", "psi4", "vasp", ...
    method             TEXT NOT NULL,                -- e.g. "CBS-QB3", "PBE-D3BJ"
    ref_type           TEXT NOT NULL DEFAULT 'atom'  -- "atom" | "element_phase"
                           CHECK (ref_type IN ('atom', 'element_phase')),
    settings           TEXT NOT NULL DEFAULT '',     -- e.g. "encut=700eV"
    energy             REAL NOT NULL,                -- kJ/mol, canonical unit
    correction         REAL,                          -- optional correction, kJ/mol
    spin_multiplicity  INTEGER,
    source             TEXT,                          -- file / job path
    computed_date      TEXT,                          -- ISO 8601, e.g. "2026-07-24"
    code_version       TEXT,
    notes              TEXT,
    UNIQUE (atomic_number, code, method, ref_type, settings)
);

CREATE INDEX IF NOT EXISTS idx_atom_energy_lookup
    ON atom_energy (code, method, ref_type, settings);
"""

_BUNDLED_DB_PATH = Path(__file__).parent / "data" / "thermochemistry.db"
_SEAMM_INI_PATH = Path("~/.seamm.d/seamm.ini").expanduser()


def _resolve_default_db_path():
    """The installer-managed database path from seamm.ini's
    [thermochemistry] section, if set, else the bundled package path.

    Reads seamm.ini directly with the stdlib configparser rather than
    importing seamm_installer, so this core module keeps its only real
    dependency (seamm_util). seamm_installer (the `installer` extra) is
    only needed to *populate* database-path in the first place -- see
    installer.py -- not to read it back here.
    """
    if _SEAMM_INI_PATH.exists():
        parser = configparser.ConfigParser(interpolation=None)
        parser.read(_SEAMM_INI_PATH)
        if parser.has_section("thermochemistry"):
            value = parser.get("thermochemistry", "database-path", fallback="")
            if value:
                return Path(value).expanduser().resolve()
    return _BUNDLED_DB_PATH


DEFAULT_DB_PATH = _resolve_default_db_path()

# kJ/mol per unit of energy, for the handful of units this data actually
# shows up in. Kept local (rather than pulling in seamm_util.Q_ for every
# call) since this is a tight, fixed set; formation.py uses Q_ for the more
# open-ended unit handling it needs.
_TO_KJ_PER_MOL = {
    "kJ/mol": 1.0,
    "kcal/mol": 4.184,
    "eV": 96.485332,
    "hartree": 2625.499639,
    "E_h": 2625.499639,
}


def _convert(value, from_units, to_kJ=True):
    try:
        factor = _TO_KJ_PER_MOL[from_units]
    except KeyError:
        raise ValueError(
            f"Unknown units {from_units!r}; known: {sorted(_TO_KJ_PER_MOL)}"
        )
    return value * factor if to_kJ else value / factor


class ThermoDB:
    """Helper around the shared atomic reference-energy SQLite database.

    Parameters
    ----------
    path : str or Path, optional
        Database file. Defaults to ``DEFAULT_DB_PATH``: the installer-
        managed location registered in ``~/.seamm.d/seamm.ini``'s
        ``[thermochemistry]`` section (``database-path``, normally
        ``~/SEAMM/Parameters/thermochemistry/thermochemistry.db`` -- see
        ``installer.py``) if the installer has been run, else the bundled
        ``data/thermochemistry.db`` shipped with the package (a stale
        prototype snapshot, useful only before the installer has run).
    read_only : bool
        Open without creating/migrating the schema, and without permission
        to write. Use for a package-installed, already-built database.
    """

    def __init__(self, path=None, read_only=False):
        self.path = Path(path) if path is not None else DEFAULT_DB_PATH
        if read_only:
            if not self.path.exists():
                raise FileNotFoundError(self.path)
            self._db = sqlite3.connect(
                f"file:{self.path}?mode=ro", uri=True, timeout=10.0
            )
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._db = sqlite3.connect(str(self.path), timeout=10.0)
        self._db.row_factory = sqlite3.Row
        self._db.execute("PRAGMA foreign_keys = ON")
        self._autocommit = True
        if not read_only:
            self._create_schema()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def close(self):
        self._db.close()

    def _maybe_commit(self):
        if self._autocommit:
            self._db.commit()

    @contextmanager
    def batch(self):
        """Defer commits until the end of the block, for bulk imports.

        `add_element`/`add_atom_energy` normally commit immediately --
        the right default for interactive use, but one fsync per row makes
        a bulk import of hundreds of thousands of rows (e.g. the full
        Gaussian/Psi4 composite-method grids) prohibitively slow. Wrap
        those calls in ``with db.batch():`` to commit once at the end
        instead.
        """
        previous = self._autocommit
        self._autocommit = False
        try:
            yield self
            self._db.commit()
        finally:
            self._autocommit = previous

    def _create_schema(self):
        self._db.executescript(_SCHEMA_SQL)
        self._db.execute(
            "INSERT OR IGNORE INTO schema_info (key, value) VALUES ('version', ?)",
            (str(_SCHEMA_VERSION),),
        )
        self._db.commit()

    # ------------------------------------------------------------------
    # element (experimental reference) table
    # ------------------------------------------------------------------

    def add_element(
        self,
        atomic_number,
        symbol,
        *,
        multiplicity=None,
        term_symbol=None,
        standard_state=None,
        dfH0_0K=None,
        dfH0_298K=None,
        dfH0_298K_stderr=None,
        h298_minus_h0_atom=None,
        h298_minus_h0_std_state=None,
        s298_gas=None,
        s298_gas_stderr=None,
        reference=None,
        reference_note=None,
    ):
        """Insert or update the experimental reference row for one element."""
        self._db.execute(
            """
            INSERT INTO element (
                atomic_number, symbol, multiplicity, term_symbol, standard_state,
                dfH0_0K, dfH0_298K, dfH0_298K_stderr, h298_minus_h0_atom,
                h298_minus_h0_std_state, s298_gas, s298_gas_stderr,
                reference, reference_note
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(atomic_number) DO UPDATE SET
                symbol=excluded.symbol,
                multiplicity=excluded.multiplicity,
                term_symbol=excluded.term_symbol,
                standard_state=excluded.standard_state,
                dfH0_0K=excluded.dfH0_0K,
                dfH0_298K=excluded.dfH0_298K,
                dfH0_298K_stderr=excluded.dfH0_298K_stderr,
                h298_minus_h0_atom=excluded.h298_minus_h0_atom,
                h298_minus_h0_std_state=excluded.h298_minus_h0_std_state,
                s298_gas=excluded.s298_gas,
                s298_gas_stderr=excluded.s298_gas_stderr,
                reference=excluded.reference,
                reference_note=excluded.reference_note
            """,
            (
                atomic_number,
                symbol,
                multiplicity,
                term_symbol,
                standard_state,
                dfH0_0K,
                dfH0_298K,
                dfH0_298K_stderr,
                h298_minus_h0_atom,
                h298_minus_h0_std_state,
                s298_gas,
                s298_gas_stderr,
                reference,
                reference_note,
            ),
        )
        self._maybe_commit()

    def get_element(self, symbol=None, atomic_number=None):
        """Return the experimental reference row for one element as a dict, or None."""
        if symbol is not None:
            cur = self._db.execute("SELECT * FROM element WHERE symbol = ?", (symbol,))
        elif atomic_number is not None:
            cur = self._db.execute(
                "SELECT * FROM element WHERE atomic_number = ?", (atomic_number,)
            )
        else:
            raise ValueError("Give either symbol or atomic_number")
        row = cur.fetchone()
        return dict(row) if row is not None else None

    def dfH0(self, symbol, *, at_0K=True):
        """Experimental atomic heat of formation, kJ/mol.

        Parameters
        ----------
        at_0K : bool
            If True (default, and the recommended anchor for an "energy of
            formation" -- see the design doc), return DfH(X, 0K). If False,
            return the 298 K gas-phase value used by the existing
            enthalpy-of-formation code paths.
        """
        el = self.get_element(symbol=symbol)
        if el is None:
            return None
        return el["dfH0_0K"] if at_0K else el["dfH0_298K"]

    def elements(self):
        """Return all element rows, ordered by atomic number."""
        cur = self._db.execute("SELECT * FROM element ORDER BY atomic_number")
        return [dict(row) for row in cur.fetchall()]

    # ------------------------------------------------------------------
    # atom_energy (computed reference) table
    # ------------------------------------------------------------------

    def add_atom_energy(
        self,
        symbol,
        code,
        method,
        energy,
        *,
        ref_type="atom",
        settings="",
        units="kJ/mol",
        correction=None,
        spin_multiplicity=None,
        source=None,
        computed_date=None,
        code_version=None,
        notes=None,
    ):
        """Insert or update one computed atomic reference energy.

        Parameters
        ----------
        symbol : str
            Element symbol, e.g. "O". Must already exist via `add_element`.
        code : str
            Originating code: "gaussian", "psi4", "orca", "vasp", ...
        method : str
            Method/functional label, e.g. "CBS-QB3", "PBE-D3BJ". Where the
            basis is folded into the label (as in the legacy Gaussian/Psi4
            tables, e.g. "CCD-FC/6-31++G(2d,2p)"), pass the full label.
        energy : float
            The atom's energy, in `units` (converted to kJ/mol on storage).
        ref_type : {"atom", "element_phase"}
            "atom" = isolated gas-phase atom. "element_phase" = energy per
            atom of the standard-state phase (VASP's existing convention).
        settings : str
            Free-form disambiguator that is part of the identity key, e.g.
            "encut=700eV". Defaults to "" (not NULL, so the UNIQUE
            constraint dedups consistently).
        units : str
            Units of `energy` and `correction`: one of "kJ/mol", "kcal/mol",
            "eV", "hartree"/"E_h".
        correction : float, optional
            Additive correction (e.g. Gaussian's "<method> correction"
            column), same units as `energy`.
        """
        atno = self._require_atomic_number(symbol)
        energy_kJ = _convert(energy, units)
        correction_kJ = None if correction is None else _convert(correction, units)
        self._db.execute(
            """
            INSERT INTO atom_energy (
                atomic_number, code, method, ref_type, settings, energy,
                correction, spin_multiplicity, source, computed_date,
                code_version, notes
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(atomic_number, code, method, ref_type, settings) DO UPDATE SET
                energy=excluded.energy,
                correction=excluded.correction,
                spin_multiplicity=excluded.spin_multiplicity,
                source=excluded.source,
                computed_date=excluded.computed_date,
                code_version=excluded.code_version,
                notes=excluded.notes
            """,
            (
                atno,
                code,
                method,
                ref_type,
                settings,
                energy_kJ,
                correction_kJ,
                spin_multiplicity,
                source,
                computed_date,
                code_version,
                notes,
            ),
        )
        self._maybe_commit()

    def get_atom_energy(
        self, symbol, code, method, *, ref_type="atom", settings="", units="kJ/mol"
    ):
        """Return one atomic reference energy (including its correction), or None."""
        atno = self._require_atomic_number(symbol)
        cur = self._db.execute(
            """
            SELECT energy, correction FROM atom_energy
            WHERE atomic_number=? AND code=? AND method=? AND ref_type=? AND settings=?
            """,
            (atno, code, method, ref_type, settings),
        )
        row = cur.fetchone()
        if row is None:
            return None
        value = row["energy"] + (row["correction"] or 0.0)
        return _convert(value, units, to_kJ=False)

    def get_reference_energies(
        self, code, method, *, ref_type="atom", settings="", units="kJ/mol"
    ):
        """Return ``{symbol: energy}`` for every element tabulated for this
        (code, method, ref_type, settings) combination.

        This is the lookup table `formation.formation_energy` needs.
        """
        cur = self._db.execute(
            """
            SELECT e.symbol, a.energy, a.correction FROM atom_energy a
            JOIN element e ON e.atomic_number = a.atomic_number
            WHERE a.code=? AND a.method=? AND a.ref_type=? AND a.settings=?
            """,
            (code, method, ref_type, settings),
        )
        result = {}
        for row in cur.fetchall():
            value = row["energy"] + (row["correction"] or 0.0)
            result[row["symbol"]] = _convert(value, units, to_kJ=False)
        return result

    def missing(self, code, method, symbols, *, ref_type="atom", settings=""):
        """Return the subset of `symbols` with no tabulated (code, method) energy."""
        have = self.get_reference_energies(
            code, method, ref_type=ref_type, settings=settings
        )
        return [s for s in symbols if s not in have]

    def list_methods(self, code=None):
        """Return distinct (code, method, ref_type, settings) tuples present."""
        if code is None:
            cur = self._db.execute(
                "SELECT DISTINCT code, method, ref_type, settings FROM atom_energy "
                "ORDER BY code, method, ref_type, settings"
            )
        else:
            cur = self._db.execute(
                "SELECT DISTINCT code, method, ref_type, settings FROM atom_energy "
                "WHERE code=? ORDER BY method, ref_type, settings",
                (code,),
            )
        return [tuple(row) for row in cur.fetchall()]

    def coverage(self, code, method, *, ref_type="atom", settings=""):
        """Return (n_elements_present, max_atomic_number) for one (code, method)."""
        cur = self._db.execute(
            """
            SELECT COUNT(*) AS n, MAX(atomic_number) AS max_z FROM atom_energy
            WHERE code=? AND method=? AND ref_type=? AND settings=?
            """,
            (code, method, ref_type, settings),
        )
        row = cur.fetchone()
        n = row["n"]
        return (n, row["max_z"] if n else None)

    # ------------------------------------------------------------------
    # review / export -- keeps the machine-generated store honest against
    # a human-reviewable snapshot (see the design doc's storage tradeoff)
    # ------------------------------------------------------------------

    def dump_elements_csv(self, path):
        """Write the `element` table to CSV, ordered by atomic number."""
        rows = self.elements()
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh)
            if rows:
                writer.writerow(rows[0].keys())
                for row in rows:
                    writer.writerow(list(row.values()))

    def dump_atom_energies_csv(self, path):
        """Write the `atom_energy` table (joined to element symbols) to CSV."""
        cur = self._db.execute("""
            SELECT e.symbol, a.* FROM atom_energy a
            JOIN element e ON e.atomic_number = a.atomic_number
            ORDER BY a.code, a.method, a.ref_type, a.settings, e.atomic_number
            """)
        rows = cur.fetchall()
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh)
            if rows:
                writer.writerow(rows[0].keys())
                for row in rows:
                    writer.writerow(list(row))

    # ------------------------------------------------------------------
    def _require_atomic_number(self, symbol):
        el = self.get_element(symbol=symbol)
        if el is None:
            raise KeyError(
                f"Unknown element symbol {symbol!r}; add it with add_element() first"
            )
        return el["atomic_number"]
