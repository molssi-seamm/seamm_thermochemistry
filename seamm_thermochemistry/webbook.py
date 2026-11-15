# -*- coding: utf-8 -*-

"""The NIST Chemistry WebBook (SRD 69) as a source of the reference library.

A maintainer tool, in two stages so that the network part is separate and
can be repeated:

``fetch``
    Looks each species up by InChI (gas-phase thermochemistry, SI units),
    saves the page in a cache directory, waits between requests and logs
    each one (URL, date, sha256). A stereo InChI that the WebBook does not
    find is retried without its stereo layers, and the retry is recorded.
``import``
    Offline and deterministic: parses every gas-phase enthalpy of formation
    on the cached pages into the library, with the original citation the
    WebBook gives, then chooses one value per species by
    :func:`.references.select_measurement`.

Rows the WebBook marks as evaluations (class "cal", e.g. JANAF reviews) are
kept but held: the value is taken from the evaluation's own source.

Two cases the plain InChI lookup does not settle are followed by ``fetch``:

* An InChI that matches several WebBook entries gives a list of them. Each is
  fetched, and the one whose own InChI is the one asked for and that has
  gas-phase enthalpies of formation is used; if not exactly one does, the
  species has no page.
* A page that gives only an average of several values ("AVG", with a link to
  the individual data points) has that second page fetched beside it
  (``<page>.hfg.html``). ``import`` reads the individual values from it; the
  average, an evaluation, is kept but held.

Every page fetched is logged with its URL, which ``import`` records as where
each value was retrieved from.
"""

import csv
import datetime
import hashlib
import html
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

from .references import ReferenceError

SOURCE = "WebBook-SRD69"
SOURCE_RANK = 50
SOURCE_CITATION = (
    "NIST Chemistry WebBook, NIST Standard Reference Database 69, eds. P. J. "
    "Linstrom and W. G. Mallard, https://webbook.nist.gov (values cited to the "
    "original measurements)"
)
URL = "https://webbook.nist.gov/cgi/cbook.cgi"
USER_AGENT = "seamm-thermochemistry reference-library builder"

#: The index of a cache directory: one row per page
INDEX = "index.csv"
INDEX_FIELDS = ("file", "smiles", "charge", "multiplicity", "inchi", "inchi_retry")


def _clean(text):
    return html.unescape(
        re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", text or ""))
    ).strip()


def parse_value(text):
    """(value, uncertainty) from a WebBook cell such as ``132.9 ± 0.6``; the
    uncertainty is None when not given, and both None for no number."""
    m = re.match(r"\s*([-+]?\d+\.?\d*)\s*(?:±\s*(\d+\.?\d*))?", text or "")
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2)) if m.group(2) else None


def year_of(reference):
    """The last four-digit year in a reference such as "Roux et al., 2008"."""
    years = re.findall(r"(1[89]\d\d|20\d\d)", reference or "")
    return int(years[-1]) if years else None


def parse_page(text):
    """The gas-phase enthalpies of formation on a WebBook page.

    Returns {"title", "found", "cas", "formula", "rows"}, each row a dict
    with ``class`` (exp or cal), ``value``, ``uncertainty``, ``units``,
    ``method``, ``reference`` (short), ``citation`` (full), ``year`` and
    ``comment``.
    """
    m = re.search(r"<title>(.*?)</title>", text, re.S)
    title = _clean(m.group(1)) if m else ""
    result = {
        "title": title,
        "found": bool(title)
        and "not found" not in title.lower()
        and "not found" not in text[:3000].lower()
        and "Information from the InChI" not in title
        and title != "Search Results",
        "rows": [],
        "hits": [],
        "individual": None,
    }
    m = re.search(r"(InChI=1S?/[^\s<\"]+)", text)
    result["inchi"] = html.unescape(m.group(1)) if m else None
    if title == "Search Results":
        # Several entries match: their pages, in the WebBook's order
        for href in re.findall(r'<li><a href="(/cgi/cbook\.cgi\?ID=[^"]+)"', text):
            result["hits"].append(_absolute(href))
    m = re.search(r"CAS Registry Number:</strong>\s*([0-9-]+)", text)
    result["cas"] = m.group(1) if m else None
    m = re.search(r"Formula</a>:</strong>(.*?)</li>", text, re.S)
    result["formula"] = (
        re.sub(r"<sup>.*?</sup>", "", m.group(1))
        .replace("<sub>", "")
        .replace("</sub>", "")
        .strip()
        if m
        else None
    )
    refs = {}
    for rid, body in re.findall(
        r'<span id="(ref-\d+)">(.*?)(?=<span id="ref-|</p>\s*(?:<p>\s*<span '
        r'id="ref-|</div|<hr|<h2))',
        text,
        re.S,
    ):
        body = re.sub(r'<span class="Z3988".*?</span>', "", body, flags=re.S)
        refs[rid] = _clean(body).replace(" [ all data ]", "")
    m = re.search(r'<h2 id="Thermo-Gas">(.*?)(<h2 |<hr)', text, re.S)
    if m:
        for cls, row in re.findall(
            r'<tr class="(exp|cal)">(.*?)</tr>', m.group(1), re.S
        ):
            cells = re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)
            if len(cells) < 5:
                continue
            quantity = _clean(cells[0]).replace(" ", "")
            if not quantity.startswith("ΔfH°gas"):
                continue
            link = re.search(
                r'href="([^"]*Type=HFG[^"]*)"[^>]*>\s*Individual\s+data\s+points',
                row,
            )
            if link:
                result["individual"] = _absolute(link.group(1))
            value, uncertainty = parse_value(_clean(cells[1]))
            link = re.search(r'href="#(ref-\d+)"', cells[4])
            short = _clean(cells[4])
            full = refs.get(link.group(1), "") if link else ""
            result["rows"].append(
                {
                    "class": cls,
                    "value": value,
                    "uncertainty": uncertainty,
                    "units": _clean(cells[2]),
                    "method": _clean(cells[3]),
                    "reference": short,
                    "citation": full or short,
                    "year": year_of(short),
                    "comment": _clean(cells[5]) if len(cells) > 5 else "",
                }
            )
    # The individual data points (Type=HFG): value, method, reference, comment
    m = re.search(
        r'<table class="data" aria-label="Enthalpy of formation of gas[^"]*">'
        r"(.*?)</table>",
        text,
        re.S,
    )
    if m:
        for cls, row in re.findall(
            r'<tr class="(exp|cal)">(.*?)</tr>', m.group(1), re.S
        ):
            cells = re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)
            if len(cells) < 3:
                continue
            value, uncertainty = parse_value(_clean(cells[0]))
            link = re.search(r'href="#(ref-\d+)"', cells[2])
            short = _clean(cells[2])
            full = refs.get(link.group(1), "") if link else ""
            result["rows"].append(
                {
                    "class": cls,
                    "value": value,
                    "uncertainty": uncertainty,
                    "units": "kJ/mol",
                    "method": _clean(cells[1]),
                    "reference": short,
                    "citation": full or short,
                    "year": year_of(short),
                    "comment": _clean(cells[3]) if len(cells) > 3 else "",
                }
            )
    return result


def _absolute(href):
    """A WebBook link as an absolute URL."""
    href = html.unescape(href)
    if href.startswith("/"):
        return "https://webbook.nist.gov" + href
    return href


def enthalpy_rows(page):
    """The rows of a parsed page with a gas-phase enthalpy of formation."""
    return [
        r for r in page["rows"] if r["value"] is not None and r["units"] == "kJ/mol"
    ]


def companion(path):
    """The file of a page's individual data points, beside it."""
    path = Path(path)
    return path.with_name(path.stem + ".hfg.html")


def logged_urls(directory):
    """{file: URL} from a cache's fetch log (the last fetch of each file)."""
    log = Path(directory) / "fetch_log.csv"
    if not log.exists():
        return {}
    with open(log, newline="", encoding="utf-8") as fd:
        return {r["file"]: r["url"] for r in csv.DictReader(fd) if r.get("file")}


def non_stereo_inchi(inchi):
    """An InChI without its stereo layers (/t, /m, /s)."""
    return re.sub(r"/[tms][^/]*", "", inchi)


def page_url(inchi):
    return (
        URL + "?" + urllib.parse.urlencode({"InChI": inchi, "Units": "SI", "Mask": "1"})
    )


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_index(directory):
    path = Path(directory) / INDEX
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fd:
        return list(csv.DictReader(fd))


def fetch(species, directory, *, delay=3.0, user_agent=USER_AGENT, opener=None):
    """Fetch WebBook pages into a cache directory.

    `species` is a list of (smiles, charge, multiplicity). Pages already in
    the cache are not fetched again. Returns the index rows written.
    """
    from rdkit import Chem

    from .smiles import mol_from_smiles

    directory = Path(directory).expanduser()
    directory.mkdir(parents=True, exist_ok=True)
    opener = opener or urllib.request.urlopen
    index = {
        (r["smiles"], r["charge"], r["multiplicity"]) for r in read_index(directory)
    }
    log = directory / "fetch_log.csv"
    new_log = not log.exists()
    written = []
    with (
        open(directory / INDEX, "a", newline="", encoding="utf-8") as fi,
        open(log, "a", newline="", encoding="utf-8") as fl,
    ):
        iw = csv.DictWriter(fi, INDEX_FIELDS)
        lw = csv.writer(fl)
        if fi.tell() == 0:
            iw.writeheader()
        if new_log:
            lw.writerow(("date", "url", "file", "sha256", "purpose"))

        def get(url, path, purpose):
            """Fetch a page into the cache, unless it is there, and log it."""
            if path.exists():
                return
            request = urllib.request.Request(url, headers={"User-Agent": user_agent})
            with opener(request, timeout=60) as response:
                path.write_bytes(response.read())
            lw.writerow(
                (
                    datetime.date.today().isoformat(),
                    url,
                    path.name,
                    _sha256(path),
                    purpose,
                )
            )
            time.sleep(delay)

        for smiles, charge, multiplicity in species:
            if (smiles, str(charge), str(multiplicity)) in index:
                continue
            inchi = Chem.MolToInchi(mol_from_smiles(smiles))
            tries = [inchi]
            if non_stereo_inchi(inchi) != inchi:
                tries.append(non_stereo_inchi(inchi))
            row = None
            for attempt in tries:
                url = page_url(attempt)
                name = hashlib.sha1(attempt.encode()).hexdigest()[:12] + ".html"
                path = directory / name
                get(url, path, f"WebBook gas thermochemistry {smiles}")
                page = parse_page(path.read_text(errors="ignore"))
                if not page["found"] and page["hits"]:
                    # Several entries: the one that is this species, with data
                    hits, name, page = page["hits"], None, None
                    for hit in hits:
                        hit_name = hashlib.sha1(hit.encode()).hexdigest()[:12] + ".html"
                        get(hit, directory / hit_name, f"WebBook entry for {smiles}")
                        candidate = parse_page(
                            (directory / hit_name).read_text(errors="ignore")
                        )
                        if (
                            candidate["found"]
                            and candidate["inchi"] == attempt
                            and (enthalpy_rows(candidate) or candidate["individual"])
                        ):
                            if name is not None:
                                name, page = None, None  # ambiguous
                                break
                            name, page = hit_name, candidate
                    if page is None:
                        continue
                    path = directory / name
                if page["found"]:
                    if page["individual"]:
                        get(
                            page["individual"],
                            companion(path),
                            f"WebBook individual data points {smiles}",
                        )
                    row = {
                        "file": name,
                        "smiles": smiles,
                        "charge": charge,
                        "multiplicity": multiplicity,
                        "inchi": inchi,
                        "inchi_retry": attempt if attempt != inchi else "",
                    }
                    break
            if row is None:
                row = {
                    "file": "",
                    "smiles": smiles,
                    "charge": charge,
                    "multiplicity": multiplicity,
                    "inchi": inchi,
                    "inchi_retry": "",
                }
            iw.writerow(row)
            written.append(row)
    return written


def import_webbook(library, directory):
    """Read the cached WebBook pages listed in a cache's index into `library`.

    Returns a list of (smiles, outcome, detail): outcome is "added" (with
    the number of values), "no page", "no gas-phase value" or "refused".
    """
    directory = Path(directory).expanduser()
    library.add_source(SOURCE, "measured", SOURCE_RANK, citation=SOURCE_CITATION)
    report = []
    urls = logged_urls(directory)
    for entry in read_index(directory):
        smiles = entry["smiles"]
        if not entry.get("file"):
            report.append((smiles, "no page", ""))
            continue
        path = directory / entry["file"]
        page = parse_page(path.read_text(errors="ignore"))
        retry = entry.get("inchi_retry") or None
        url = urls.get(path.name) or page_url(retry or entry.get("inchi") or "")
        rows = [(r, f"{url} sha256 {_sha256(path)}") for r in enthalpy_rows(page)]
        individual = companion(path)
        if individual.exists():
            points = parse_page(individual.read_text(errors="ignore"))
            where = (
                f"{urls.get(individual.name) or page['individual'] or ''} "
                f"sha256 {_sha256(individual)}"
            ).strip()
            rows += [(r, where) for r in enthalpy_rows(points)]
        if not page["found"] or not rows:
            report.append((smiles, "no gas-phase value", page["title"]))
            continue
        try:
            species = library.add_species(
                smiles,
                int(entry.get("charge") or 0),
                int(entry.get("multiplicity") or 1),
                formula=page["formula"] or None,
                name=page["title"],
                cas=page["cas"],
                inchi=entry.get("inchi") or None,
                inchi_retry=retry,
                notes=(
                    "found on the non-stereo InChI; the stereoisomer is not verified"
                    if retry
                    else None
                ),
            )
        except ReferenceError as e:
            report.append((smiles, "refused", str(e)))
            continue
        for r, retrieved in rows:
            if r["class"] == "cal":
                status = "held"
                note = (
                    "an evaluation shown by the WebBook: taken from the "
                    "evaluation's own source instead"
                )
            else:
                status = "review" if r["method"] == "Review" else "measured"
                note = (
                    "derived from another value on the page"
                    if r["method"] == "N/A"
                    else None
                )
            try:
                library.add_measurement(
                    species,
                    SOURCE,
                    citation=r["citation"],
                    dfh_298k=r["value"],
                    uncertainty=r["uncertainty"],
                    unc_level="unstated",
                    status=status,
                    year=r["year"],
                    source_ref=page["title"],
                    method=r["method"],
                    comment=r["comment"] or None,
                    retrieved=retrieved,
                    notes=note,
                )
            except ReferenceError as e:
                report.append((smiles, "refused", str(e)))
                break
        else:
            library.choose(species, SOURCE)
            report.append((smiles, "added", str(len(rows))))
    library.apply_overrides()
    return report
