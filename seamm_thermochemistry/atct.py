# -*- coding: utf-8 -*-

"""A local copy of the Active Thermochemical Tables (ATcT) and the rules for
choosing an entry.

ATcT (https://atct.anl.gov) publishes enthalpies of formation at 0 K and 298.15
K for a few thousand species, derived from one thermochemical network (TN), so
the values are mutually consistent and correlated. This module keeps a complete
copy of one version of the tables on disk, builds an index from canonical
SMILES to the ATcT entries, chooses the right entry for a molecule, and
propagates the network covariance through a reaction.

Layout of the cache, one directory per TN version::

    <base>/ATcT/<version>/
        atct_all.json       the /api/v1/all/ download, unmodified
        metadata.json       version, source URL, download date, SHA-256
        index.json          canonical SMILES -> ATcT IDs, built locally
        covariance.json     covariance blocks fetched so far, by version

where ``<base>`` is the installation's ``Parameters/thermochemistry`` directory
by default. Nothing here is distributed with SEAMM: each installation downloads
its own copy, and the copy is never modified, so results can be reproduced.

Conventions that matter (Ruscic 2004, 2005; Ruscic & Bross 2025):

* The listed uncertainties, and the covariance matrix, are 95 % confidence
  limits, valid at 298.15 K only. There is no 0 K uncertainty.
* A plain ``(g)`` entry is the species in its ground state: at 0 K the lowest
  level of the lowest conformer; at 298.15 K a thermal average over conformers.
  Conformers, stereoisomers, spin states and nuclear-spin isomers are
  separate entries, labelled in the formula's state, e.g. ``(g, trans)``,
  ``(g, triplet)``, ``(g, para)``.
* Ions use the stationary-electron convention.
* Never mix TN versions in one reaction.

RDKit is needed only to build the index and to canonicalize query SMILES; it
is imported lazily, so the rest of the package does not depend on it.
"""

from datetime import datetime, timezone
import hashlib
import json
import logging
import math
from pathlib import Path
import re
import urllib.error
import urllib.parse
import urllib.request

from . import __version__

__all__ = [
    "ATcT",
    "ATcTEntry",
    "ATcTError",
    "AmbiguousEntry",
    "MissingEntry",
    "atct_directory",
    "canonical_smiles",
    "download_all",
    "normalize_entry",
    "reaction_uncertainty",
]

logger = logging.getLogger(__name__)

API_BASE = "https://atct.anl.gov/api/v1"
USER_AGENT = (
    f"SEAMM-thermochemistry/{__version__} (+https://molssi-seamm.github.io; "
    "contact seamm@molssi.org)"
)

# State descriptors that name nuclear-spin isomers; never selected.
_NUCLEAR_SPIN = {"ortho", "para", "normal"}
# ... spin states, keyed by multiplicity
_SPIN_STATES = {
    1: "singlet",
    2: "doublet",
    3: "triplet",
    4: "quartet",
    5: "quintet",
    6: "sextet",
    7: "septet",
    8: "octet",
}
_SPIN_NAMES = {v: k for k, v in _SPIN_STATES.items()}


class ATcTError(Exception):
    """Base class for errors from the ATcT cache."""


class MissingEntry(ATcTError, KeyError):
    """No ATcT entry matches."""


class AmbiguousEntry(ATcTError):
    """More than one ATcT entry matches; the candidates are listed."""

    def __init__(self, message, candidates):
        super().__init__(message)
        self.candidates = candidates


# ---------------------------------------------------------------------------
# Entries
# ---------------------------------------------------------------------------


def _key(name):
    """A key name reduced to lower-case letters, digits and underscores, so
    that the /all/ names (``Delta_Hf_0K``) and the per-species names
    (``∆fH_0K``) compare equal after the prefixes are stripped."""
    s = re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()
    s = s.replace("delta_hf", "fh").replace("delta_hf298", "fh_298")
    return s


_KEY_MAP = {
    "atct_tn_version": "version",
    "atct_id": "id",
    "atctid": "id",
    "name": "name",
    "formula": "formula",
    "fh_0k": "dfH0_0K",
    "fh_298k": "dfH0_298K",
    "fh_298k_uncertainty": "uncertainty_298K",
    "fh298k_uncertainty": "uncertainty_298K",
    "unit": "units",
    "units": "units",
    "mass": "mass",
    "mass_uncertainty": "mass_uncertainty",
    "smiles": "smiles",
    "casrn": "casrn",
    "charge": "charge",
    "inchi": "inchi",
    "inchi_key": "inchikey",
    "inchikey": "inchikey",
    "xyz": "xyz",
}

_STATE_RE = re.compile(r"^\s*(?P<formula>.*?)\s*\((?P<state>[^)]*)\)\s*$")


def _number(value):
    """A float from the API's string numbers; None for missing/'exact'."""
    if value is None or value == "":
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text.lower() in ("exact", "n/a", "none"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def normalize_entry(raw):
    """One ATcT record, from ``/all/`` or a per-species route, as a dict with
    fixed keys: ``version, id, name, formula, phase, descriptor, charge,
    dfH0_0K, dfH0_298K, uncertainty_298K (None if exact), exact, units,
    smiles, casrn, inchi, inchikey, mass``.

    The formula column carries the state, e.g. ``CH4  (g)`` or
    ``C2H2Cl2  (g, trans)``. It is split into ``formula``, ``phase`` (``g``,
    ``cr``, ``l``, ``cr,l``, ``aq``, ...) and ``descriptor`` (``trans``,
    ``triplet``, ``ortho``, ``syn``, ... or ``""``).
    """
    entry = {}
    for name, value in raw.items():
        key = _KEY_MAP.get(_key(name))
        if key is not None:
            entry[key] = value

    formula = str(entry.get("formula", "")).strip()
    match = _STATE_RE.match(formula)
    if match:
        formula = match.group("formula").strip()
        state = match.group("state")
        parts = [p.strip() for p in state.split(",")]
        phase = parts[0]
        descriptor = ", ".join(parts[1:])
        # "(cr,l)" is one phase label
        if descriptor and parts[1] in ("l", "cr", "g", "s") and len(parts) == 2:
            phase = f"{parts[0]},{parts[1]}"
            descriptor = ""
    else:
        phase = ""
        descriptor = ""

    uncertainty = entry.get("uncertainty_298K")
    exact = isinstance(uncertainty, str) and uncertainty.strip().lower() == "exact"

    charge = entry.get("charge")
    try:
        charge = int(float(charge)) if charge not in (None, "") else 0
    except (TypeError, ValueError):
        charge = 0

    return {
        "version": str(entry.get("version", "")),
        "id": str(entry.get("id", "")),
        "name": str(entry.get("name", "")),
        "formula": formula,
        "phase": phase,
        "descriptor": descriptor,
        "charge": charge,
        "dfH0_0K": _number(entry.get("dfH0_0K")),
        "dfH0_298K": _number(entry.get("dfH0_298K")),
        "uncertainty_298K": 0.0 if exact else _number(uncertainty),
        "exact": exact,
        "units": str(entry.get("units", "kJ/mol")),
        "smiles": str(entry.get("smiles", "")).strip(),
        "casrn": str(entry.get("casrn", "")).strip(),
        "inchi": entry.get("inchi"),
        "inchikey": entry.get("inchikey"),
        "mass": _number(entry.get("mass")),
    }


class ATcTEntry(dict):
    """A normalized ATcT entry (a dict) with a few conveniences."""

    @property
    def casrn_group(self):
        """The CAS number part of the ID: conformers and states of one species
        share it (``74-82-8*0``, ``74-82-8*1`` ...)."""
        return self["id"].split("*")[0]

    @property
    def index(self):
        """The ``*n`` index of the ID: 0 is the plain entry."""
        try:
            return int(self["id"].split("*")[1])
        except (IndexError, ValueError):
            return 0

    @property
    def is_gas(self):
        return self["phase"].startswith("g")

    @property
    def is_plain(self):
        """A plain entry, no conformer, state or isomer label."""
        return self["descriptor"] == ""

    @property
    def multiplicity(self):
        """The spin multiplicity named in the descriptor, or None."""
        for word in self["descriptor"].replace(",", " ").split():
            if word.lower() in _SPIN_NAMES:
                return _SPIN_NAMES[word.lower()]
        return None

    @property
    def is_nuclear_spin_isomer(self):
        words = {w.strip().lower() for w in self["descriptor"].split(",")}
        return bool(words & _NUCLEAR_SPIN)

    def label(self):
        state = self["phase"]
        if self["descriptor"]:
            state += ", " + self["descriptor"]
        unc = "exact" if self["exact"] else f"± {self['uncertainty_298K']}"
        return (
            f"{self['id']}  {self['name']}  {self['formula']} ({state})  "
            f"{self['smiles']}  ΔfH°(298.15 K) = {self['dfH0_298K']} {unc} "
            f"{self['units']}"
        )


# ---------------------------------------------------------------------------
# SMILES
# ---------------------------------------------------------------------------


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit import RDLogger
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "RDKit is needed to index and look up ATcT entries by SMILES: "
            "install seamm-thermochemistry[cbh] or rdkit."
        ) from e
    RDLogger.DisableLog("rdApp.*")
    return Chem


def canonical_smiles(smiles, *, stereo=True, keep_single_center=False):
    """RDKit canonical SMILES for matching.

    Double-bond stereochemistry is kept (ATcT lists cis and trans isomers as
    separate species). Tetrahedral stereocentres are dropped when the molecule
    has only one, since enantiomers have the same enthalpy of formation and
    ATcT has no enantiomer entries; with two or more, the relative
    configuration is kept. Returns None if RDKit cannot parse the SMILES.
    """
    Chem = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if not stereo:
        # Drop the stereo first, so that a hydrogen that only defined it
        # ([H]/N=C\O) is removed too
        Chem.RemoveStereochemistry(mol)
        mol = Chem.RemoveHs(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=False)
    if not keep_single_center:
        centers = Chem.FindMolChiralCenters(
            mol, includeUnassigned=False, useLegacyImplementation=False
        )
        if len(centers) == 1:
            atom = mol.GetAtomWithIdx(centers[0][0])
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _natural_mass(smiles):
    """The average molecular mass at natural abundance, or None."""
    Chem = _rdkit()
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Descriptors.MolWt(mol)


# ---------------------------------------------------------------------------
# The cache on disk
# ---------------------------------------------------------------------------


def atct_directory(base=None):
    """The directory holding the ATcT versions.

    `base` overrides everything. Otherwise the installation's
    ``Parameters/thermochemistry/ATcT``, found as ``db.py`` finds the
    reference database: the current root first, then the directory of the
    database registered in seamm.ini, then ``~/SEAMM/Parameters/thermochemistry``.
    """
    if base is not None:
        return Path(base).expanduser()
    from .db import DEFAULT_DB_PATH, _BUNDLED_DB_PATH

    try:
        from seamm_util import current_root

        own = current_root() / "Parameters" / "thermochemistry"
        if own.exists():
            return own / "ATcT"
    except ImportError:  # pragma: no cover
        pass
    if DEFAULT_DB_PATH != _BUNDLED_DB_PATH:
        return Path(DEFAULT_DB_PATH).parent / "ATcT"
    return Path("~/SEAMM/Parameters/thermochemistry/ATcT").expanduser()


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fd:
        for chunk in iter(lambda: fd.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _get(url, *, user_agent=USER_AGENT, timeout=60):
    """GET a URL and return the decoded JSON. Raises ATcTError with the HTTP
    status for 403 (blocked), 429 (rate limited) and other failures."""
    request = urllib.request.Request(
        url, headers={"User-Agent": user_agent, "Accept": "application/json"}
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = response.read()
    except urllib.error.HTTPError as e:
        if e.code == 403:
            raise ATcTError(
                f"ATcT refused the request ({url}): HTTP 403. The site blocks "
                "some clients; try again later or download the data in a browser."
            ) from None
        if e.code == 429:
            raise ATcTError(
                "ATcT rate limit reached (100 requests per minute per address); "
                "wait a minute and try again."
            ) from None
        raise ATcTError(f"ATcT request failed ({url}): HTTP {e.code}") from None
    except urllib.error.URLError as e:
        raise ATcTError(f"Cannot reach ATcT ({url}): {e.reason}") from None
    try:
        return json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        raise ATcTError(f"ATcT returned something other than JSON from {url}")


def download_all(base=None, *, user_agent=USER_AGENT, source=None):
    """Download ``/all/`` into ``<base>/<version>/atct_all.json`` and write its
    metadata. An existing version is never overwritten. `source` may be a local
    file (a copy saved from a browser) instead of the API.

    Returns the version directory.
    """
    base = atct_directory(base)
    url = f"{API_BASE}/all/"
    if source is not None:
        source = Path(source).expanduser()
        with open(source, encoding="utf-8") as fd:
            data = json.load(fd)
        origin = str(source)
    else:
        data = _get(url, user_agent=user_agent)
        origin = url
    if not isinstance(data, list) or not data:
        raise ATcTError("The ATcT download is not a list of species")

    versions = {str(e.get("ATcT_TN_Version", "")) for e in data}
    if len(versions) != 1:
        raise ATcTError(f"The ATcT download mixes TN versions: {sorted(versions)}")
    version = versions.pop()

    directory = base / version
    target = directory / "atct_all.json"
    if target.exists():
        raise ATcTError(f"ATcT version {version} is already in {directory}")
    directory.mkdir(parents=True, exist_ok=True)
    partial = directory / "atct_all.json.download"
    with open(partial, "w", encoding="utf-8") as fd:
        json.dump(data, fd, ensure_ascii=False, indent=1)
    partial.replace(target)
    metadata = {
        "version": version,
        "source": origin,
        "downloaded": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": _sha256(target),
        "n_entries": len(data),
        "citation": (
            "B. Ruscic and D. H. Bross, Active Thermochemical Tables (ATcT) values "
            f"based on ver. {version} of the Thermochemical Network, Argonne "
            "National Laboratory, Lemont, Illinois, USA; available at ATcT.anl.gov"
        ),
    }
    with open(directory / "metadata.json", "w", encoding="utf-8") as fd:
        json.dump(metadata, fd, indent=1)
    return directory


class ATcT:
    """One version of the ATcT tables, from the local cache.

    Parameters
    ----------
    version : str, optional
        The TN version, e.g. "1.220". The newest cached version by default.
    base : path-like, optional
        The cache directory; see `atct_directory`.
    """

    def __init__(self, version=None, base=None):
        self.base = atct_directory(base)
        versions = self.available(self.base)
        if not versions:
            raise FileNotFoundError(
                f"No ATcT data in {self.base}. Run "
                "'seamm-thermochemistry-installer atct' or "
                "seamm_thermochemistry.atct.download_all()."
            )
        if version is None:
            version = versions[-1]
        elif version not in versions:
            raise FileNotFoundError(
                f"ATcT version {version} is not in {self.base}; have {versions}"
            )
        self.version = version
        self.directory = self.base / version
        with open(self.directory / "metadata.json", encoding="utf-8") as fd:
            self.metadata = json.load(fd)
        with open(self.directory / "atct_all.json", encoding="utf-8") as fd:
            raw = json.load(fd)
        self.entries = {}
        for record in raw:
            entry = ATcTEntry(normalize_entry(record))
            self.entries[entry["id"]] = entry
        self._index = None
        self._covariance = None

    # -- versions ----------------------------------------------------------

    @staticmethod
    def available(base=None):
        """The cached TN versions, oldest first."""
        base = atct_directory(base)
        if not base.exists():
            return []
        versions = [p.name for p in base.iterdir() if (p / "atct_all.json").exists()]
        return sorted(versions, key=lambda v: [int(x) for x in re.findall(r"\d+", v)])

    @property
    def citation(self):
        return self.metadata.get("citation", "")

    # -- the index ---------------------------------------------------------

    @property
    def index(self):
        """{"stereo": {canonical SMILES: [ids]}, "flat": {SMILES without
        stereo: [ids]}, "inchikey": {key: [ids]}, "unparsed": [ids]}, built
        with RDKit on first use and kept in ``index.json``."""
        if self._index is None:
            path = self.directory / "index.json"
            if path.exists():
                with open(path, encoding="utf-8") as fd:
                    self._index = json.load(fd)
            else:
                self._index = self.build_index()
                try:
                    with open(path, "w", encoding="utf-8") as fd:
                        json.dump(self._index, fd, indent=1)
                except OSError:  # a read-only copy
                    pass
        return self._index

    def build_index(self):
        """Canonicalize every entry's SMILES with RDKit."""
        import rdkit

        stereo = {}
        flat = {}
        inchikey = {}
        unparsed = []
        for entry in self.entries.values():
            smiles = entry["smiles"]
            key = entry.get("inchikey")
            if key:
                inchikey.setdefault(key, []).append(entry["id"])
            if not smiles:
                unparsed.append(entry["id"])
                continue
            can = canonical_smiles(smiles)
            if can is None:
                unparsed.append(entry["id"])
                continue
            stereo.setdefault(can, []).append(entry["id"])
            flat.setdefault(canonical_smiles(smiles, stereo=False), []).append(
                entry["id"]
            )
        return {
            "rdkit": rdkit.__version__,
            "stereo": stereo,
            "flat": flat,
            "inchikey": inchikey,
            "unparsed": unparsed,
        }

    # -- lookup -------------------------------------------------------------

    def candidates(self, smiles, *, charge=None, gas_only=True):
        """All entries whose canonical SMILES matches, as ATcTEntry objects.

        If `smiles` carries no double-bond or tetrahedral stereo but matching
        entries do, those entries are returned too (they come from the
        stereo-free index), so that the caller can see them.
        """
        can = canonical_smiles(smiles)
        if can is None:
            raise ATcTError(f"RDKit cannot parse the SMILES '{smiles}'")
        ids = list(self.index["stereo"].get(can, []))
        flat = canonical_smiles(smiles, stereo=False)
        if can == flat:
            # No stereo in the query: include the stereo entries too
            for _id in self.index["flat"].get(flat, []):
                if _id not in ids:
                    ids.append(_id)
        elif not ids:
            # Stereo in the query that the tables do not have (e.g. the E/Z of
            # an N-H imine): fall back to entries written without stereo
            for _id in self.index["flat"].get(flat, []):
                entry_smiles = self.entries[_id]["smiles"]
                if canonical_smiles(entry_smiles) == canonical_smiles(
                    entry_smiles, stereo=False
                ):
                    ids.append(_id)
        result = [self.entries[_id] for _id in ids]
        if gas_only:
            result = [e for e in result if e.is_gas]
        if charge is not None:
            result = [e for e in result if e["charge"] == charge]
        return result

    def lookup(
        self,
        smiles=None,
        *,
        atct_id=None,
        charge=0,
        multiplicity=None,
        conformer=None,
    ):
        """The ATcT entry for a molecule, by the selection rules.

        Parameters
        ----------
        smiles : str
            The molecule, with double-bond stereo if it has any.
        atct_id : str, optional
            An explicit ATcT ID, which overrides every rule.
        charge : int
            The net charge; ATcT lists ions separately.
        multiplicity : int, optional
            The spin multiplicity. A labelled spin-state entry is used when it
            matches; otherwise the plain entry, assumed to be the ground state.
        conformer : str, optional
            A conformer or isomer label to select (``syn``, ``anti``, ``cis``,
            ``trans``, ...) instead of the plain entry.

        Returns
        -------
        ATcTEntry

        Raises
        ------
        MissingEntry, AmbiguousEntry
        """
        if atct_id is not None:
            if atct_id not in self.entries:
                raise MissingEntry(f"No ATcT entry with ID {atct_id}")
            return self.entries[atct_id]
        if smiles is None:
            raise ValueError("Give a SMILES or an ATcT ID")

        candidates = self.candidates(smiles, charge=charge)
        if not candidates:
            raise MissingEntry(
                f"No gas-phase ATcT entry (ver. {self.version}) for '{smiles}' "
                f"with charge {charge}"
            )

        # Nuclear-spin isomers are never wanted
        candidates = [e for e in candidates if not e.is_nuclear_spin_isomer]

        # A requested conformer or isomer label
        if conformer is not None:
            wanted = conformer.strip().lower()
            chosen = [
                e
                for e in candidates
                if wanted in {w.strip().lower() for w in e["descriptor"].split(",")}
            ]
            if len(chosen) == 1:
                return chosen[0]
            if not chosen:
                raise MissingEntry(
                    f"No ATcT entry for '{smiles}' labelled '{conformer}'; "
                    + "candidates:\n  "
                    + "\n  ".join(e.label() for e in candidates)
                )
            raise AmbiguousEntry(
                f"Several ATcT entries for '{smiles}' labelled '{conformer}':\n  "
                + "\n  ".join(e.label() for e in chosen),
                chosen,
            )

        # Spin state: a labelled entry that matches the multiplicity wins
        if multiplicity is not None:
            spin = [e for e in candidates if e.multiplicity == multiplicity]
            if len(spin) == 1:
                return spin[0]
            if len(spin) > 1:
                raise AmbiguousEntry(
                    f"Several ATcT {_SPIN_STATES.get(multiplicity, '')} entries for "
                    f"'{smiles}':\n  " + "\n  ".join(e.label() for e in spin),
                    spin,
                )
            # Otherwise drop the other labelled spin states
            candidates = [e for e in candidates if e.multiplicity is None]

        # The plain entries. Several can remain for stereoisomers without a
        # plain mixture entry (cis/trans with a stereo-free query), or for
        # positional isomers that canonicalize alike (they don't).
        plain = [e for e in candidates if e.is_plain]
        if len(plain) == 1:
            return plain[0]
        if len(plain) > 1:
            # Prefer the lowest *0 entry of a single CAS group
            groups = {e.casrn_group for e in plain}
            if len(groups) == 1:
                plain.sort(key=lambda e: e.index)
                return plain[0]
            # Isotopologues (H2, HD, D2) share a SMILES: keep the natural
            # abundance one, whose mass matches the query's
            natural = _natural_mass(smiles)
            if natural is not None:
                close = [
                    e
                    for e in plain
                    if e["mass"] is not None and abs(e["mass"] - natural) < 0.3
                ]
                if len(close) == 1:
                    return close[0]
                if close:
                    plain = close
            raise AmbiguousEntry(
                f"Several plain ATcT entries match '{smiles}' (different species "
                "with the same connectivity, e.g. stereoisomers): give the stereo "
                "in the SMILES or an ATcT ID.\n  "
                + "\n  ".join(e.label() for e in plain),
                plain,
            )

        # No plain entry: only labelled ones (e.g. (g, cis) and (g, trans) for
        # 2-butene). A single match is fine; otherwise ask for the stereo.
        if len(candidates) == 1:
            return candidates[0]
        raise AmbiguousEntry(
            f"ATcT has only labelled entries for '{smiles}'; give the stereo or "
            "a conformer label, or an ATcT ID:\n  "
            + "\n  ".join(e.label() for e in candidates),
            candidates,
        )

    # -- covariance --------------------------------------------------------

    @property
    def covariance_cache(self):
        if self._covariance is None:
            path = self.directory / "covariance.json"
            if path.exists():
                with open(path, encoding="utf-8") as fd:
                    self._covariance = json.load(fd)
            else:
                self._covariance = {"version": self.version, "blocks": []}
        return self._covariance

    def _save_covariance(self):
        path = self.directory / "covariance.json"
        try:
            with open(path, "w", encoding="utf-8") as fd:
                json.dump(self.covariance_cache, fd, indent=1)
        except OSError:  # a read-only copy
            pass

    def covariance(self, ids, *, fetch=True, user_agent=USER_AGENT):
        """The covariance matrix (kJ²/mol², 95 % scale, 298.15 K) for the
        given ATcT IDs, in the order given.

        A cached block that contains every ID is used; otherwise the block is
        fetched from ``/covariance/species/`` (one request for all IDs),
        checked against the listed uncertainties, and cached.

        Returns None if the block is not cached and `fetch` is False or the
        request fails; callers then fall back to independent uncertainties.
        """
        ids = list(ids)
        for block in self.covariance_cache["blocks"]:
            if all(_id in block["ids"] for _id in ids):
                return self._extract(block, ids)
        if not fetch:
            return None
        try:
            block = self.fetch_covariance(ids, user_agent=user_agent)
        except ATcTError as e:
            logger.warning(str(e))
            return None
        return self._extract(block, ids)

    def fetch_covariance(self, ids, *, user_agent=USER_AGENT):
        """Fetch and cache the covariance block for `ids`."""
        query = urllib.parse.urlencode({"atctids": ",".join(ids)})
        data = _get(f"{API_BASE}/covariance/species/?{query}", user_agent=user_agent)
        block = self.parse_covariance(data)
        if block.get("version") not in (None, "", self.version):
            raise ATcTError(
                f"The covariance block is for TN version {block['version']}, "
                f"the cache is {self.version}"
            )
        self.covariance_cache["blocks"].append(block)
        self._save_covariance()
        return block

    def parse_covariance(self, data):
        """The API's covariance response as {"ids": [...], "matrix": [[...]]}.

        The rows come back in the API's own order (``species_atctids``), not
        the request order. The square root of each diagonal element must equal
        the listed 298.15 K uncertainty, which confirms that the matrix is on
        the 95 % scale; a mismatch is an error.
        """
        ids = data.get("species_atctids")
        if ids is None:
            ids = [s["atctid"] for s in data.get("species", [])]
        matrix = data.get("matrix")
        if not ids or matrix is None or len(matrix) != len(ids):
            raise ATcTError("Malformed covariance response from ATcT")
        for i, _id in enumerate(ids):
            entry = self.entries.get(_id)
            if entry is None or entry["uncertainty_298K"] is None:
                continue
            listed = entry["uncertainty_298K"]
            diagonal = math.sqrt(max(matrix[i][i], 0.0))
            if abs(diagonal - listed) > 0.002 + 0.01 * listed:
                raise ATcTError(
                    f"The covariance diagonal for {_id} gives an uncertainty of "
                    f"{diagonal:.4f} but the table lists {listed}; the matrix "
                    "is not on the expected 95 % scale"
                )
        return {
            "version": str(data.get("version", "")) or None,
            "units": data.get("units", "kJ/mol"),
            "ids": list(ids),
            "matrix": [[float(x) for x in row] for row in matrix],
        }

    @staticmethod
    def _extract(block, ids):
        pos = {_id: i for i, _id in enumerate(block["ids"])}
        return [[block["matrix"][pos[a]][pos[b]] for b in ids] for a in ids]


def reaction_uncertainty(coefficients, uncertainties, covariance=None):
    """The 95 % uncertainty of Σ νᵢ ΔfHᵢ at 298.15 K.

    Parameters
    ----------
    coefficients : [float]
        The stoichiometric coefficients ν of the reference species.
    uncertainties : [float]
        Their listed 95 % uncertainties, same order.
    covariance : [[float]], optional
        Their covariance matrix, same order. If None the uncertainties are
        combined in quadrature, which ignores the correlations.

    Returns
    -------
    (float, bool)
        The uncertainty, and whether the covariance was used.
    """
    c = [float(v) for v in coefficients]
    if covariance is None:
        return math.sqrt(sum((ci * ui) ** 2 for ci, ui in zip(c, uncertainties))), False
    n = len(c)
    total = 0.0
    for i in range(n):
        for j in range(n):
            total += c[i] * covariance[i][j] * c[j]
    return math.sqrt(max(total, 0.0)), True
