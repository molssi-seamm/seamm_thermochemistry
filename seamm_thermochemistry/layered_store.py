# -*- coding: utf-8 -*-

"""The fragment store in two layers: the user's local store and the official
one, looked up local first.

``fragments.db`` (the local layer) is the store write-back, harvest, merge and
supersede write to. ``fragments-official.db`` (the official layer) is published
on Zenodo, installed and replaced whole by the installer, and only ever read,
without locking (it may sit on a network filesystem). See the P5 design note,
``docs/developer_guide/campaigns/2026-11-10``, §4.

The rules, for rows of the same key in both layers:

R1
    The local row is used, unless R2 applies.
R2
    The official row has been superseded and the local one equals a value it
    superseded (within 0.1 kJ/mol): the official row is used, and noted.
R3
    Otherwise, values differing by more than 0.1 kJ/mol are noted as a
    conflict, and the local row is used.
R4
    Conformers come from both layers; the lowest E + ZPE wins, a tie within
    0.1 kJ/mol going to the local layer (:meth:`LayeredStore.lowest`).
    Verified rows in either layer come before unverified ones in either.
R5
    Energies pair with geometries only on the same geometry hash: looked up in
    the geometry's own layer first, then the other. Every row returned carries
    its ``layer``, so a caller can say when a geometry and its energy come from
    different layers.

Notes from R2 and R3 accumulate in :attr:`LayeredStore.notes` for the caller to
report.
"""

import logging

from .fragment_store import FragmentStore, official_store_path
from .geometry_store import DUPLICATE_TOLERANCE, _as_list, same_model_chemistry

logger = logging.getLogger(__name__)

__all__ = ["LayeredStore"]

_VALUES = {
    "fragment": ("E", "ZPE", "H"),
    "geometry": ("E", "ZPE", "H"),
    "energy": ("E",),
}


def _difference(a, b, fields):
    """The largest difference of the fields both rows have (0 if none)."""
    d = 0.0
    for f in fields:
        x, y = a.get(f), b.get(f)
        if x is None or y is None:
            continue
        d = max(d, abs(float(x) - float(y)))
    return d


class LayeredStore:
    """The local and official fragment stores, read as one (see the module
    docstring).

    Parameters
    ----------
    local_path : str or Path, optional
        The local store; by default ``fragments.db`` next to the ATcT cache.
    official_path : str or Path, optional
        The official store; by default ``fragments-official.db`` beside it. A
        missing file means no official layer.
    create_local : bool
        Open the local store for writing (creating it if missing).
    """

    def __init__(self, local_path=None, official_path=None, *, create_local=False):
        self.notes = []
        self.local = FragmentStore(local_path, create=create_local)
        if self.local.path.exists() and self.local.layer == "official":
            raise ValueError(
                f"{self.local.path} is an official store, not a local one: the "
                "official store is never written."
            )
        self.official = None
        self.official_version = None
        path = official_store_path(official_path)
        if path.exists() and path.resolve() != self.local.path.resolve():
            official = FragmentStore(path, immutable=True)
            if official.layer != "official":
                self._note(
                    f"{path} is not an official store (its layer is "
                    f"'{official.layer}'); it is not used."
                )
                official.close()
            else:
                self.official = official
                self.official_version = official.meta().get("version")

    # -- housekeeping ------------------------------------------------------

    @property
    def path(self):
        return self.local.path

    def close(self):
        self.local.close()
        if self.official is not None:
            self.official.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __len__(self):
        return sum(len(store) for _, store in self._layers())

    def n_geometries(self):
        return sum(store.n_geometries() for _, store in self._layers())

    def schema_version(self):
        return self.local.schema_version()

    def _layers(self):
        yield "local", self.local
        if self.official is not None:
            yield "official", self.official

    def _note(self, text):
        if text not in self.notes:
            self.notes.append(text)
            logger.info(text)

    def layer_name(self, layer):
        """'local', or 'official <version>'."""
        if layer == "official" and self.official_version:
            return f"official {self.official_version}"
        return layer

    # -- R1-R3 -------------------------------------------------------------

    def _choose(self, table, what, local_row, official_row):
        """The row to use for one key found in both layers."""
        fields = _VALUES[table]
        if _difference(local_row, official_row, fields) <= DUPLICATE_TOLERANCE:
            return local_row  # R1
        for old in self.official.superseded_values(table, official_row["id"]):
            if _difference(local_row, old, fields) <= DUPLICATE_TOLERANCE:
                self._note(  # R2
                    f"The local {what} is the value the "
                    f"{self.layer_name('official')} store superseded on "
                    f"{old.get('date')} ({old.get('reason')}); the official value "
                    "is used."
                )
                return official_row
        self._note(  # R3
            f"The local and {self.layer_name('official')} stores differ for {what} "
            f"by {_difference(local_row, official_row, fields):.2f} kJ/mol; the "
            "local value is used."
        )
        return local_row

    # -- reading -----------------------------------------------------------

    def geometries(self, species, **kwargs):
        """The geometry rows of a species at a frequency level from both layers
        (see :meth:`FragmentStore.geometries` for the arguments), one per key,
        each with its ``layer``."""
        found = []
        for name, store in self._layers():
            for row in store.geometries(species, **kwargs):
                row["layer"] = name
                found.append(row)
        if any(r["verified"] for r in found):
            found = [r for r in found if r["verified"]]
        chosen = {}
        for row in found:  # local rows first
            key = row.key
            if key not in chosen:
                chosen[key] = row
            elif row["layer"] == "official" and chosen[key]["layer"] == "local":
                label = f"{species}" + (
                    f" (conformer '{row['conformer']}')" if row["conformer"] else ""
                )
                chosen[key] = self._choose(
                    "geometry", f"geometry of {label}", chosen[key], row
                )
        return sorted(
            chosen.values(), key=lambda r: (r["conformer"], r["layer"] != "local")
        )

    @staticmethod
    def lowest(rows):
        """The lowest conformer (E + ZPE) of geometry rows, a tie within 0.1
        kJ/mol going to the local layer (R4); None if none has E and ZPE."""
        usable = [r for r in rows if r.get("E") is not None and r.get("ZPE")]
        if not usable:
            return None

        def total(r):
            return r["E"] + r["ZPE"]

        best = min(usable, key=total)
        local = [r for r in usable if r.get("layer", "local") == "local"]
        if best.get("layer", "local") != "local" and local:
            candidate = min(local, key=total)
            if total(candidate) - total(best) <= DUPLICATE_TOLERANCE:
                return candidate
        return best

    def energy(self, geometry_hash, *, prefer="local", **kwargs):
        """The energy row on a geometry at a level (see
        :meth:`FragmentStore.energy`), from the layer `prefer` (the geometry's)
        first, then the other (R5); with its ``layer``."""
        found = {}
        for name, store in self._layers():
            row = store.energy(geometry_hash, **kwargs)
            if row is not None:
                row["layer"] = name
                found[name] = row
        if not found:
            return None
        if len(found) == 1:
            return next(iter(found.values()))
        local, official = found["local"], found["official"]
        if local["verified"] != official["verified"]:
            return local if local["verified"] else official
        if local["energy_fingerprint"] == official["energy_fingerprint"]:
            return self._choose(
                "energy",
                f"energy on geometry {geometry_hash[:15]}...",
                local,
                official,
            )
        return found[prefer] if prefer in found else local

    def energies(self, geometry_hash):
        out = []
        for name, store in self._layers():
            for row in store.energies(geometry_hash):
                row["layer"] = name
                out.append(row)
        return out

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
        """As :meth:`FragmentStore.species_at`, over both layers: each entry
        also has "layers" (the geometry's layer, the energy's)."""
        energy_fingerprints = _as_list(energy_fingerprint)
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
                prefer=g["layer"],
                energy_fingerprint=energy_fingerprint,
                model_chemistry=energy_model_chemistry,
            )
            if e is None:
                same = (
                    energy_fingerprints
                    and g.get("optimization_fingerprint") in energy_fingerprints
                ) or (
                    not energy_fingerprints
                    and energy_model_chemistry
                    and g.get("model_chemistry")
                    and same_model_chemistry(
                        energy_model_chemistry,
                        g["model_chemistry"],
                        None,
                        g.get("program"),
                    )
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
                    "layers": (g["layer"], e["layer"] if e else g["layer"]),
                }
            )
        return out

    def get(self, species, **kwargs):
        """The best fragment-table record of a species at a level (see
        :meth:`FragmentStore.get`) over both layers, with its ``layer``."""
        found = {}
        for name, store in self._layers():
            record = store.get(species, **kwargs)
            if record is not None:
                record["layer"] = name
                found[name] = record
        if not found:
            return None
        if len(found) == 1:
            return next(iter(found.values()))
        local, official = found["local"], found["official"]
        if local["verified"] != official["verified"]:
            return local if local["verified"] else official
        if local.key == official.key:
            return self._choose("fragment", f"{species}", local, official)
        return local

    # -- writing: the local layer only -------------------------------------

    def put(self, record):
        return self.local.put(record)

    def put_geometry(self, record):
        return self.local.put_geometry(record)

    def put_geometry_conformer(self, record, rmsd=None):
        return self.local.put_geometry_conformer(record, rmsd)

    def put_energy(self, record):
        return self.local.put_energy(record)

    def supersede(self, record, reason, force=False):
        return self.local.supersede(record, reason, force)

    def history(self, table=None):
        return self.local.history(table)

    def records(self, species=None):
        return self.local.records(species)

    def duplicates(self):
        return self.local.duplicates()
