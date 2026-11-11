# -*- coding: utf-8 -*-

"""Tools for the two-layer fragment store (see :mod:`.layered_store`):

* :func:`migrate_store` brings an older store up to date (the ``meta`` table, the
  model chemistries stored bare with their program), after a dated backup;
* :func:`export_local` writes the local rows that the official store lacks, or
  has with other values, to a contribution file for the next official version;
* :func:`merge_stores` adds the rows of other stores (copies from other
  installations, contributions, campaign stores) to a store;
* :func:`build_official` builds a new official store from the previous one and
  contributions;
* :func:`layer_conflicts` lists the keys whose values differ between the local
  and the official store.

All writes are insert-or-ignore with provenance: a row whose key is already
there with values differing by more than 0.1 kJ/mol is logged, never
overwritten, and superseded values travel as history.
"""

import datetime
import hashlib
import logging
from pathlib import Path
import shutil
import sqlite3

from .fragment_store import FragmentRecord, FragmentStore
from .geometry_store import (
    _ENERGY_COLUMNS,
    _GEOMETRY_COLUMNS,
    _GEOMETRY_KEY,
    DUPLICATE_TOLERANCE,
    SCHEMA_VERSION,
    EnergyRecord,
    GeometryRecord,
    _decode,
    split_model_chemistry,
)

logger = logging.getLogger(__name__)

__all__ = [
    "migrate_store",
    "export_local",
    "merge_stores",
    "build_official",
    "layer_conflicts",
]

#: The key columns of each table
KEYS = {
    "fragment": (
        "species",
        "charge",
        "multiplicity",
        "energy_fingerprint",
        "frequency_fingerprint",
    ),
    "geometry": _GEOMETRY_COLUMNS[:_GEOMETRY_KEY],
    "energy": _ENERGY_COLUMNS[:2],
}
VALUES = {
    "fragment": ("E", "ZPE", "H"),
    "geometry": ("E", "ZPE", "H"),
    "energy": ("E",),
}


def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fd:
        for chunk in iter(lambda: fd.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _rows(store, table):
    """All rows of a table as records (JSON columns decoded); [] if none."""
    try:
        db = store._connect()
    except FileNotFoundError:
        return []
    try:
        rows = db.execute(f"SELECT * FROM {table} ORDER BY id").fetchall()
    except sqlite3.OperationalError:  # a store from before the table existed
        return []
    if table == "fragment":
        out = []
        for r in rows:
            record = FragmentRecord(dict(r))
            record["verified"] = bool(record["verified"])
            out.append(record)
        return out
    return [
        _decode(r, GeometryRecord if table == "geometry" else EnergyRecord)
        for r in rows
    ]


def _row_by_key(store, table, record):
    """The row of `table` with `record`'s key, or None."""
    columns = KEYS[table]
    where = " AND ".join(f"{c} = ?" for c in columns)
    try:
        row = (
            store._connect()
            .execute(
                f"SELECT * FROM {table} WHERE {where}", [record[c] for c in columns]
            )
            .fetchone()
        )
    except (FileNotFoundError, sqlite3.OperationalError):
        return None
    if row is None:
        return None
    return dict(row)


def _difference(a, b, fields):
    d = 0.0
    for f in fields:
        x, y = a.get(f), b.get(f)
        if x is None or y is None:
            continue
        d = max(d, abs(float(x) - float(y)))
    return d


def _put(target, table, record):
    """Insert-or-ignore one row into `target`; returns the outcome."""
    record = dict(record)
    record.pop("id", None)
    if table == "fragment":
        return target.put(FragmentRecord(record))
    if table == "geometry":
        return target.put_geometry(GeometryRecord(record))[0]
    record["geometry"] = None  # relinked to the target's geometry row
    return target.put_energy(EnergyRecord(record))[0]


def _copy_history(source, target, table, record, label):
    """Copy the superseded values of `record`'s row in `source` to the same key's
    row in `target`, as history (the date and source in the reason)."""
    old_values = source.superseded_values(table, record["id"])
    if not old_values:
        return 0
    row = _row_by_key(target, table, record)
    if row is None:
        return 0
    db = target._writer()
    n = 0
    for old in old_values:
        reason = f"{old.get('reason')} (superseded {old.get('date')} in {label})"
        if table == "fragment":
            db.execute(
                "INSERT INTO superseded (fragment, E, ZPE, H, superseded_on, reason) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (row["id"], old.get("E"), old.get("ZPE"), old.get("H"), _now(), reason),
            )
        else:
            target._history(
                db, table, row["id"], "superseded", old, None, reason=reason
            )
        n += 1
    db.commit()
    return n


# -- migrate -------------------------------------------------------------------


def migrate_store(path=None, *, backup=True):
    """Bring a store up to date: the ``meta`` table and schema 3, and the model
    chemistries stored bare with the program filled in (``ORCA:CMP@R2SCAN-3C``
    -> ``CMP@R2SCAN-3C`` + ``ORCA``). A dated copy of the file is kept first.
    Running it again changes nothing (and makes no copy).

    Keys do not change: verified rows are keyed by fingerprints, and the keys of
    unverified rows keep the prefixed spelling. Returns {"backup": path or None,
    "changed": {table: n}, "schema": (old, new)}.
    """
    store = FragmentStore(path)
    if not store.path.exists():
        raise FileNotFoundError(f"No fragment store at {store.path}")
    old_schema = store.schema_version()
    plan = {}
    for table, columns in (
        ("fragment", ("model_chemistry", "frequency_model_chemistry")),
        ("geometry", ("model_chemistry",)),
        ("energy", ("model_chemistry",)),
    ):
        try:
            rows = (
                store._connect()
                .execute(f"SELECT id, program, {', '.join(columns)} FROM {table}")
                .fetchall()
            )
        except sqlite3.OperationalError:
            continue
        for row in rows:
            program = row["program"]
            new = {}
            for c in columns:
                model, prefix = split_model_chemistry(row[c], program)
                if model != row[c]:
                    new[c] = model
                program = program or prefix
            if program != row["program"]:
                new["program"] = program
            if new:
                plan.setdefault(table, []).append((row["id"], new))
    needs_meta = old_schema < SCHEMA_VERSION or "layer" not in store.meta()
    store.close()
    result = {"backup": None, "changed": {}, "schema": (old_schema, old_schema)}
    if not plan and not needs_meta:
        return result
    if backup:
        stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        copy = store.path.with_name(f"{store.path.name}.pre-migrate-{stamp}")
        shutil.copy2(store.path, copy)
        result["backup"] = copy
    with FragmentStore(store.path, create=True) as writer:  # adds meta, schema 3
        db = writer._connect(create=True)
        for table, changes in plan.items():
            for row_id, new in changes:
                assignments = ", ".join(f"{c} = ?" for c in new)
                db.execute(
                    f"UPDATE {table} SET {assignments} WHERE id = ?",
                    list(new.values()) + [row_id],
                )
            result["changed"][table] = len(changes)
        db.commit()
        result["schema"] = (old_schema, writer.schema_version())
    return result


# -- merge, export, official -----------------------------------------------------


def merge_stores(target, sources, *, label=None, allow_official=False):
    """Add the rows of `sources` (paths or stores) to `target` (a path or a
    writable store): insert-or-ignore, differing duplicates logged, superseded
    values copied as history. An official target is refused unless
    `allow_official` (for :func:`build_official`).

    Returns {table: {outcome: n}} plus "history": n.
    """
    own = not isinstance(target, FragmentStore)
    store = FragmentStore(target, create=True) if own else target
    try:
        if store.layer == "official" and not allow_official:
            raise ValueError(
                f"{store.path} is an official store; merge into a local store, or "
                "build a new official one."
            )
        counts = {"fragment": {}, "geometry": {}, "energy": {}, "history": 0}
        for source in sources:
            src = source if isinstance(source, FragmentStore) else FragmentStore(source)
            name = label or str(src.path)
            try:
                for table in ("fragment", "geometry", "energy"):
                    for record in _rows(src, table):
                        outcome = _put(store, table, record)
                        tally = counts[table]
                        tally[outcome] = tally.get(outcome, 0) + 1
                        counts["history"] += _copy_history(
                            src, store, table, record, name
                        )
                store.set_meta(**{f"merged {_now()}": name})
            finally:
                if src is not source:
                    src.close()
        return counts
    finally:
        if own:
            store.close()


def export_local(out, local=None, official=None, *, label=None):
    """Write the local store's rows that the official store lacks, or has with
    values differing by more than 0.1 kJ/mol (flagged in the history as
    "differs from official", with the official values), to the contribution file
    `out`, with their provenance and superseded values. `label` is free text
    recorded in the file; nothing about the host or user is added.

    Returns {table: n rows exported}, plus "differing": n.
    """
    from .fragment_store import official_store_path

    out = Path(out).expanduser()
    if out.exists():
        raise FileExistsError(f"{out} exists; export to a new file.")
    local_store = FragmentStore(local)
    official_path = official_store_path(official)
    official_store = (
        FragmentStore(official_path, immutable=True) if official_path.exists() else None
    )
    counts = {"fragment": 0, "geometry": 0, "energy": 0, "differing": 0}
    try:
        with FragmentStore(out, create=True, layer="contribution") as target:
            target.set_meta(
                label=label or "",
                exported=_now(),
                official_version=(
                    official_store.meta().get("version") if official_store else ""
                ),
            )
            for table in ("fragment", "geometry", "energy"):
                for record in _rows(local_store, table):
                    theirs = (
                        _row_by_key(official_store, table, record)
                        if official_store
                        else None
                    )
                    if theirs is not None and (
                        _difference(record, theirs, VALUES[table])
                        <= DUPLICATE_TOLERANCE
                    ):
                        continue
                    _put(target, table, record)
                    counts[table] += 1
                    _copy_history(local_store, target, table, record, "the local store")
                    if theirs is not None:
                        counts["differing"] += 1
                        row = _row_by_key(target, table, record)
                        db = target._writer()
                        target._history(
                            db,
                            table,
                            row["id"],
                            "differs from official",
                            {k: theirs.get(k) for k in VALUES[table]},
                            {k: record.get(k) for k in VALUES[table]},
                            difference=_difference(record, theirs, VALUES[table]),
                        )
                        db.commit()
    finally:
        local_store.close()
        if official_store is not None:
            official_store.close()
    return counts


def build_official(out, version, sources, *, previous=None, source=None):
    """Build an official store: the previous official one (if given) plus the
    contributions `sources`, marked official with `version` (and `source`, e.g.
    its Zenodo DOI), compacted, with a ``.sha256`` file beside it.

    A contribution whose values differ from a row already there stops the build,
    and nothing is written: the maintainer resolves each with ``supersede`` and a
    reason first. Returns (path, sha256, counts).
    """
    out = Path(out).expanduser()
    if out.exists():
        raise FileExistsError(f"{out} exists; build to a new file.")
    if previous is not None:
        with FragmentStore(previous, immutable=True) as old:
            if old.layer != "official":
                raise ValueError(f"{previous} is not an official store.")
        shutil.copy2(Path(previous).expanduser(), out)
    try:
        with FragmentStore(out, create=True, layer="official") as store:
            counts = merge_stores(store, sources, allow_official=True)
            logged = sum(counts[t].get("logged", 0) for t in VALUES)
            if logged:
                conflicts = [
                    f"{h['table_name']} row {h['row']}: {h['difference']:.2f} kJ/mol"
                    for h in store.history()
                    if h["kind"] == "duplicate"
                ]
                conflicts += [
                    f"fragment {d['species']}: {d['difference']:.2f} kJ/mol"
                    for d in store.duplicates()
                ]
                raise ValueError(
                    f"{logged} contributed values differ from those already in the "
                    "store; resolve them with 'supersede' first:\n  "
                    + "\n  ".join(conflicts[:50])
                )
            store.set_meta(
                layer="official",
                version=version,
                created=_now(),
                source=source or "",
                schema=SCHEMA_VERSION,
            )
            store._connect(create=True).execute("VACUUM")
    except Exception:
        out.unlink(missing_ok=True)
        raise
    digest = _sha256(out)
    out.with_name(out.name + ".sha256").write_text(f"{digest}  {out.name}\n")
    return out, digest, counts


def layer_conflicts(local=None, official=None):
    """The keys in both the local and the official store whose values differ by
    more than 0.1 kJ/mol: [{"table", "key", "local", "official", "difference",
    "superseded"}], "superseded" when the local value is one the official store
    superseded (the official value is then used)."""
    from .fragment_store import official_store_path

    official_path = official_store_path(official)
    if not official_path.exists():
        return []
    out = []
    with (
        FragmentStore(local) as mine,
        FragmentStore(official_path, immutable=True) as theirs,
    ):
        for table in ("fragment", "geometry", "energy"):
            for record in _rows(mine, table):
                other = _row_by_key(theirs, table, record)
                if other is None:
                    continue
                d = _difference(record, other, VALUES[table])
                if d <= DUPLICATE_TOLERANCE:
                    continue
                superseded = any(
                    _difference(record, old, VALUES[table]) <= DUPLICATE_TOLERANCE
                    for old in theirs.superseded_values(table, other["id"])
                )
                out.append(
                    {
                        "table": table,
                        "key": tuple(record[c] for c in KEYS[table]),
                        "local": {k: record.get(k) for k in VALUES[table]},
                        "official": {k: other.get(k) for k in VALUES[table]},
                        "difference": d,
                        "superseded": superseded,
                        "job": record.get("job"),
                    }
                )
    return out
