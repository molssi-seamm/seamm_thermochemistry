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
import json
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
    is_unusable,
    split_model_chemistry,
)

logger = logging.getLogger(__name__)

__all__ = [
    "ground_state_multiplicity",
    "sanity_check",
    "load_exclusions",
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


def _package_version():
    from . import __version__

    return __version__


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
    """Bring a store up to date: the ``meta`` and ``writeback_log`` tables
    (schema 4), and the model chemistries stored bare with the program filled in
    (``ORCA:CMP@R2SCAN-3C``
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
    with FragmentStore(
        store.path, create=True
    ) as writer:  # adds meta, writeback_log: schema 4
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
                    if is_unusable(record.get("model_chemistry")):
                        continue  # quarantined rows are not contributed
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


def ground_state_multiplicity(species, charge=0, atct=None):
    """(multiplicity, source) of a species' ground state where it is known, or
    (None, None).

    Singlet for the species of :data:`.cbh.SINGLET_GROUND_STATES` (silylenes,
    dihalocarbenes and BH, which a SMILES writes as triplets), triplet for those of
    :data:`.cbh.TRIPLET_GROUND_STATES` (O2, S2, SO, OCCO, cyclic C3, which a
    SMILES writes closed-shell); else from ATcT when
    it lists the species in labelled spin states: the plain entry is the ground
    state, and the labelled entry with the same enthalpy of formation names its
    multiplicity (CH2: plain = triplet; SiH2: plain = singlet). Unlike
    :func:`.cbh.ground_state_multiplicity` this does not fall back on the radical
    count, which is not knowledge of the ground state (O2 would be a singlet).
    """
    from .atct import canonical_smiles
    from .cbh import SINGLET_GROUND_STATES

    try:
        from .cbh import TRIPLET_GROUND_STATES
    except ImportError:  # pragma: no cover
        TRIPLET_GROUND_STATES = ()

    try:
        can = canonical_smiles(species, stereo=False) or species
    except Exception:
        can = species
    for states, multiplicity, name in (
        (SINGLET_GROUND_STATES, 1, "singlet ground state (cbh.SINGLET_GROUND_STATES)"),
        (TRIPLET_GROUND_STATES, 3, "triplet ground state (cbh.TRIPLET_GROUND_STATES)"),
    ):
        for smiles in states:
            try:
                if canonical_smiles(smiles, stereo=False) == can:
                    return multiplicity, name
            except Exception:
                continue
    if atct is None:
        return None, None
    try:
        candidates = [e for e in atct.candidates(can, charge=charge) if e.is_gas]
    except Exception:
        return None, None
    plain = [e for e in candidates if e.is_plain]
    labelled = [e for e in candidates if e.multiplicity is not None]
    if len(plain) != 1 or not labelled:
        return None, None
    value = plain[0].get("dfH0_298K")
    for e in labelled:
        if (
            value is not None
            and e.get("dfH0_298K") is not None
            and (abs(float(e["dfH0_298K"]) - float(value)) < 1e-6)
        ):
            return e.multiplicity, f"ATcT {atct.version} ({e['id']})"
    return None, None


def sanity_check(store, *, atct=None):
    """The problems that keep rows out of an official store, by species state:
    {(species, charge, multiplicity): [reason, ...]}.

    * a geometry with an imaginary mode beyond 20i (a saddle point);
    * a geometry whose bonds, perceived again from its coordinates, are not its
      species', whose charge or multiplicity do not fit its SMILES, or with atoms
      too close together (:func:`.fragment_store.check_structure`);
    * a multiplicity that is not the species' ground state, where that is known
      (:func:`ground_state_multiplicity`).

    Ions and open-shell species are checked like the rest, not excluded.
    """
    from .fragment_store import check_structure
    from .geometry_store import IMAGINARY_TOLERANCE

    problems = {}

    def add(key, reason):
        problems.setdefault(key, [])
        if reason not in problems[key]:
            problems[key].append(reason)

    states = set()
    for g in _rows(store, "geometry"):
        key = (g["species"], g["charge"], g["multiplicity"])
        states.add(key)
        label = f"conformer '{g['conformer']}'" if g["conformer"] else "geometry"
        if g.get("n_imaginary"):
            add(
                key,
                f"{label}: {g['n_imaginary']} imaginary mode(s), the largest "
                f"{g['max_imaginary']:.0f}i cm^-1 (> {IMAGINARY_TOLERANCE:.0f}i)",
            )
        problem = check_structure(
            g["symbols"], g["coordinates"], g["species"], g["charge"], g["multiplicity"]
        )
        if problem:
            add(key, f"{label}: {problem}")
    for f in _rows(store, "fragment"):
        states.add((f["species"], f["charge"], f["multiplicity"]))
    for species, charge, multiplicity in sorted(states, key=str):
        ground, source = ground_state_multiplicity(species, charge, atct)
        if ground is not None and ground != multiplicity:
            add(
                (species, charge, multiplicity),
                f"multiplicity {multiplicity} is not the ground state's ({ground}, "
                f"{source})",
            )
    return problems


def _select(store, frequency_fingerprints=None, energy_fingerprints=None):
    """Delete the rows not reachable at the given levels: geometry rows of other
    frequency fingerprints, energy rows of other energy fingerprints or not on a
    remaining geometry, and fragment rows of neither. Returns {table: n deleted}."""
    db = store._connect(create=True)
    deleted = {}
    freq = list(frequency_fingerprints or [])
    energy = list(energy_fingerprints or [])
    if freq:
        marks = ", ".join("?" * len(freq))
        deleted["geometry"] = db.execute(
            f"DELETE FROM geometry WHERE frequency_fingerprint NOT IN ({marks})", freq
        ).rowcount
    sql = (
        "DELETE FROM energy WHERE geometry_hash NOT IN "
        "(SELECT geometry_hash FROM geometry)"
    )
    args = []
    if energy:
        marks = ", ".join("?" * len(energy))
        sql += f" OR energy_fingerprint NOT IN ({marks})"
        args = energy
    deleted["energy"] = db.execute(sql, args).rowcount
    if freq or energy:
        allowed_e = freq + energy
        marks_e = ", ".join("?" * len(allowed_e))
        sql = f"DELETE FROM fragment WHERE energy_fingerprint NOT IN ({marks_e})"
        args = list(allowed_e)
        if freq:
            marks_f = ", ".join("?" * len(freq))
            sql += (
                f" OR (frequency_fingerprint != '' AND frequency_fingerprint NOT IN "
                f"({marks_f}))"
            )
            args += freq
        deleted["fragment"] = db.execute(sql, args).rowcount
    db.commit()
    return deleted


def _drop_unusable(store):
    """Delete the rows quarantined by hand (model chemistry "UNUSABLE: <reason>")
    and return [{table, species or geometry, reason}] for the report."""
    db = store._connect(create=True)
    dropped = []
    for table, who in (
        ("geometry", "species || ' ' || conformer"),
        ("energy", "geometry_hash"),
        ("fragment", "species"),
    ):
        rows = db.execute(
            f"SELECT id, {who}, model_chemistry FROM {table} "
            "WHERE UPPER(TRIM(model_chemistry)) LIKE 'UNUSABLE:%'"
        ).fetchall()
        for row in rows:
            dropped.append(
                {
                    "table": table,
                    "row": row[1],
                    "reason": row[2].split(":", 1)[1].strip(),
                }
            )
        db.execute(
            f"DELETE FROM {table} WHERE UPPER(TRIM(model_chemistry)) LIKE 'UNUSABLE:%'"
        )
    db.commit()
    return dropped


def _exclude(store, excluded):
    """Delete the excluded rows: a whole species state, or one conformer of it
    (the energies on a deleted geometry, and the fragment-table records made from
    it, go too). `excluded` is {(species, charge, multiplicity, conformer or
    None): reason}; returns {key: rows deleted}."""
    from .fragment_store import geometry_hash as fragment_geometry_hash

    db = store._connect(create=True)
    removed = {}
    for key in excluded:
        species, charge, multiplicity, conformer = key
        where = "species = ? AND charge = ? AND multiplicity = ?"
        args = [species, int(charge), int(multiplicity)]
        if conformer is not None:
            where_g = where + " AND conformer = ?"
            args_g = args + [conformer]
        else:
            where_g, args_g = where, args
        geometries = db.execute(
            f"SELECT geometry_hash, symbols, coordinates FROM geometry WHERE {where_g}",
            args_g,
        ).fetchall()
        n = db.execute(f"DELETE FROM geometry WHERE {where_g}", args_g).rowcount
        if conformer is None:
            n += db.execute(f"DELETE FROM fragment WHERE {where}", args).rowcount
        else:
            # The fragment-table records made from this conformer's geometry
            for g in geometries:
                old_hash = fragment_geometry_hash(
                    json.loads(g["symbols"]), json.loads(g["coordinates"])
                )
                n += db.execute(
                    f"DELETE FROM fragment WHERE {where} AND geometry_hash = ?",
                    args + [old_hash],
                ).rowcount
        for g in geometries:
            still = db.execute(
                "SELECT COUNT(*) FROM geometry WHERE geometry_hash = ?",
                (g["geometry_hash"],),
            ).fetchone()[0]
            if not still:
                n += db.execute(
                    "DELETE FROM energy WHERE geometry_hash = ?", (g["geometry_hash"],)
                ).rowcount
        removed[key] = n
    db.commit()
    return removed


def load_exclusions(path):
    """Exclusions from a CSV file with the columns species, charge, multiplicity,
    reason and, optionally, conformer (to leave out one conformer instead of the
    whole species state): {(canonical species, charge, multiplicity, conformer or
    None): reason}."""
    import csv

    from .atct import canonical_smiles

    out = {}
    with open(Path(path).expanduser(), newline="") as fd:
        for row in csv.DictReader(fd):
            species = (row.get("species") or "").strip()
            reason = (row.get("reason") or "").strip()
            if not species:
                continue
            if not reason:
                raise ValueError(f"{path}: excluding {species} needs a reason")
            try:
                species = canonical_smiles(species) or species
            except Exception:
                pass
            key = (
                species,
                int(row.get("charge") or 0),
                int(row.get("multiplicity") or 1),
                (row.get("conformer") or "").strip() or None,
            )
            out[key] = reason
    return out


def build_official(
    out,
    version,
    sources,
    *,
    previous=None,
    source=None,
    frequency_fingerprints=None,
    energy_fingerprints=None,
    exclude=None,
    atct=None,
    check=True,
):
    """Build an official store: the previous official one (if given) plus the
    contributions `sources`, marked official with `version` (and `source`, e.g.
    its Zenodo DOI), compacted, with a ``.sha256`` file and a ``.report.json``
    beside it.

    `frequency_fingerprints` and `energy_fingerprints`, when given, keep only the
    rows reachable at those levels (e.g. the presets' r2SCAN-3c geometries and
    the revDSD energies on them). `exclude` ({(species, charge, multiplicity):
    reason}, or with a fourth element, a conformer label, to leave out one
    conformer; see :func:`load_exclusions`) leaves rows out, with the reason in
    the report.

    Nothing is written when a contribution's values differ from a row already
    there (resolve it with ``supersede`` and a reason first), or when the sanity
    pass (:func:`sanity_check`; `atct` for the multiplicities) finds a problem in
    a species state that is not excluded. Returns (path, sha256, counts).
    """
    out = Path(out).expanduser()
    if out.exists():
        raise FileExistsError(f"{out} exists; build to a new file.")
    # Older three-element keys exclude the whole species state
    exclude = {
        (tuple(k) + (None,) if len(k) == 3 else tuple(k)): v
        for k, v in (exclude or {}).items()
    }
    if previous is not None:
        with FragmentStore(previous, immutable=True) as old:
            if old.layer != "official":
                raise ValueError(f"{previous} is not an official store.")
        shutil.copy2(Path(previous).expanduser(), out)
    # Rows quarantined by hand ("UNUSABLE: ...") are dropped first and listed
    report = {
        "version": version,
        "built": _now(),
        "previous": str(previous) if previous else None,
        "sources": [str(s) for s in sources],
        "frequency fingerprints": list(frequency_fingerprints or []),
        "energy fingerprints": list(energy_fingerprints or []),
        "atct": getattr(atct, "version", None),
        "generator": f"seamm-thermochemistry {_package_version()}",
    }
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
            report["quarantined"] = _drop_unusable(store)
            report["deselected"] = _select(
                store, frequency_fingerprints, energy_fingerprints
            )
            removed = _exclude(store, exclude)
            report["excluded"] = [
                {
                    "species": k[0],
                    "charge": k[1],
                    "multiplicity": k[2],
                    "conformer": k[3],
                    "reason": reason,
                    "rows": removed.get(k, 0),
                }
                for k, reason in exclude.items()
            ]
            if check:
                problems = sanity_check(store, atct=atct)
                if problems:
                    lines = [
                        f"{s} (charge {q}, multiplicity {m}): " + "; ".join(r)
                        for (s, q, m), r in sorted(problems.items(), key=str)
                    ]
                    raise ValueError(
                        f"{len(problems)} species states fail the sanity checks; fix "
                        "them, or exclude them with a reason (--exclude):\n  "
                        + "\n  ".join(lines[:100])
                    )
            db = store._connect(create=True)
            report["rows"] = {
                t: db.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
                for t in ("fragment", "geometry", "energy")
            }
            report["species"] = db.execute(
                "SELECT COUNT(*) FROM (SELECT DISTINCT species, charge, "
                "multiplicity FROM geometry)"
            ).fetchone()[0]
            report["checked"] = bool(check)
            store.set_meta(
                layer="official",
                version=version,
                created=_now(),
                source=source or "",
                schema=SCHEMA_VERSION,
            )
            db.execute("VACUUM")
    except Exception:
        out.unlink(missing_ok=True)
        raise
    digest = _sha256(out)
    report["sha256"] = digest
    out.with_name(out.name + ".sha256").write_text(f"{digest}  {out.name}\n")
    out.with_name(out.name + ".report.json").write_text(json.dumps(report, indent=1))
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
