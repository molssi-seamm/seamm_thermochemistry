#!/usr/bin/env python3
"""Compute ORCA atomic reference energies robustly, by multi-start SCF.

Why: a single SCF start can land an open-shell atom in the wrong state. Fe with
revDSD from ORCA's PModel guess converges to 3d7 4s1, 120 kJ/mol above the
3d6 4s2 solution, with no warning. The reference energy of an atom for a method
is defined here as

    the LOWEST-energy SCF solution of the method at the experimental
    ground-state spin multiplicity

which is well defined even where the method's ordering of states disagrees
with experiment. (That is where the atomization approach to formation energies
is breaking down anyway; for careful work, reaction energies to well-known
species are the better tool.)

For each element, method and basis this runs several SCF starts:

* ``pbe0``: orbitals of PBE0/def2-SV(P) from a damped PModel guess (the
  protocol of the original SEAMM atom-energy flowchart, which converges well);
* ORCA's own ``PModel``, ``HCore``, ``Hueckel`` and ``PAtom`` guesses -- except
  that from La (Z = 57) on, where Hueckel and PAtom are unavailable, a
  ``smear`` start replaces them: PBE/def2-SV(P) orbitals from a Fermi-smeared
  SCF;
* then, for each distinct state found in ANY basis (identified by the Mulliken
  s/p/d/f charge and spin populations), that state's orbitals read into every
  other basis where it has not been found yet. So all bases of an element end
  up sharing the same lowest state.

It keeps the converged solution with the lowest SCF energy (the SCF part is
variational; the MP2 part of a double hybrid is not, so the SCF energy is what
decides). <S**2> does not disqualify a solution, since spin contamination is
normal in unrestricted DFT; more than 20% above S(S+1) is flagged.

Settings that were found to matter:

* ``NoCOSX``: ORCA 6.1.1's default RIJCOSX mis-builds the virtual orbitals of
  some lone atoms (Na, Na+ and Mg are ~5 kJ/mol off in the MP2 part; Li and H
  are off at DEFGRID3).
* ``ConvCheckMode 0``: otherwise ORCA can stop on the energy change alone with
  loose orbitals, which leaves 0.01-0.1 kJ/mol noise in the MP2 part.

Outputs, in ``--out``:

* ``runs.csv``: every SCF run (energies, <S**2>, populations, status);
* ``Results.csv``: the chosen energies, in the wide format that
  ``seamm_thermochemistry.importers.import_orca_atom_results`` reads;
* ``choices.csv``: which start won for each entry, how many distinct states
  were seen, and any flags (open-shell singlet, <S**2>, no convergence).

It is resumable: finished runs are parsed (also if their ``orca.out`` was
compressed to ``orca.out.gz``), not rerun, so a new version of the script can
extend a finished run with just its new starts. ``--rechoose`` rebuilds
``choices.csv`` and ``Results.csv`` from an existing ``runs.csv`` without
running anything (e.g. after changing the selection rule).

Example::

    orca_atom_multistart.py --elements 26 --methods REVDSD-PBEP86-D4/2021 \\
        --bases def2-SVP def2-TZVPPD --workers 4 --out ~/atoms_fe
"""

import argparse
import concurrent.futures
import configparser
import csv
import gzip
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import time

KJ_PER_EH = 2625.499639

DEFAULT_BASES = [
    "def2-SV(P)",
    "def2-SVP",
    "def2-TZVP",
    "def2-TZVP(-f)",
    "def2-QZVPP",
    "def2-SVPD",
    "def2-TZVPD",
    "def2-TZVPPD",
    "def2-QZVPD",
    "def2-QZVPPD",
    "ma-def2-SVP",
    "ma-def2-mSVP",
    "ma-def2-TZVP",
    "ma-def2-TZVP(-f)",
    "ma-def2-TZVPP",
    "ma-def2-QZVPP",
]
DIRECT_GUESSES = ["PModel", "HCore", "Hueckel", "PAtom"]
# ORCA's Hueckel and PAtom guesses need an extended-Hueckel minimal basis, which
# stops at Z = 56 ("Atomic number (57) too high"). From La on those starts are
# skipped and a smeared seed (below) is added instead.
EHT_MAX_Z = 56

# The damping block of the original SEAMM atom-energy flowchart.
SCF_BLOCK = """%scf
  MaxIter 500
  ConvCheckMode 0
  Shift Shift 0.3 ErrOff 0.05 end
  DIISBfac 1.1
{guess}end
"""

# <S**2> above S(S+1) by more than this fraction is flagged (not rejected: the
# lowest solution is kept regardless, since spin contamination is a normal
# feature of unrestricted DFT, e.g. Sc's 4s2 pair polarized by its 3d electron).
S2_REL_TOL = 0.20


def safe(name):
    return re.sub(r"[^A-Za-z0-9.+-]", "_", name)


def db_method(method):
    """The ThermoDB spelling of an ORCA keyword ('/' is reserved)."""
    return method.replace("/", "_")


def open_shell_singlet(term):
    """True for an experimental term with S=0 but L>0 (e.g. Ce 1G4): no single
    determinant describes it, so any SCF value is questionable."""
    m = re.match(r"(\d+)([A-Z])", term or "")
    return bool(m) and m.group(1) == "1" and m.group(2) != "S"


# --------------------------------------------------------------------------
# One ORCA run
# --------------------------------------------------------------------------
def write_input(path, keywords, symbol, mult, guess_block, moinp=None):
    lines = [f"! {keywords}"]
    if moinp:
        lines.append(f'%moinp "{moinp}"')
    lines.append(SCF_BLOCK.format(guess=guess_block).rstrip())
    lines += [f"* xyz 0 {mult}", f"{symbol} 0 0 0", "*", ""]
    path.write_text("\n".join(lines))


def read_output(out_path):
    """The text of orca.out, or of orca.out.gz if the output was compressed."""
    if out_path.exists():
        return out_path.read_text(errors="replace")
    gz = out_path.with_name(out_path.name + ".gz")
    if gz.exists():
        with gzip.open(gz, "rt", errors="replace") as fh:
            return fh.read()
    return ""


def parse(out_path):
    """Energies, <S**2>, convergence and Mulliken l-populations from orca.out
    (or orca.out.gz)."""
    text = read_output(out_path)
    r = {"status": "missing"}
    if not text:
        return r
    if "ORCA TERMINATED NORMALLY" not in text:
        r["status"] = "error"
    elif "SCF NOT CONVERGED" in text or "SCF CONVERGED AFTER" not in text:
        r["status"] = "not converged"
    else:
        r["status"] = "ok"
    m = re.findall(r"^Total Energy\s+:\s+(-?\d+\.\d+) Eh", text, re.M)
    if m:
        r["scf"] = float(m[-1])
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    if m:
        r["total"] = float(m[-1])
    m = re.findall(r"Expectation value of <S\*\*2>\s+:\s+(-?\d+\.\d+)", text)
    if m:
        r["s2"] = float(m[-1])
    # Open shell: "...CHARGES AND SPIN POPULATIONS" with CHARGE and SPIN
    # sub-blocks; closed shell: "...CHARGES" alone.
    blk = re.split(r"MULLIKEN REDUCED ORBITAL CHARGES(?: AND SPIN POPULATIONS)?", text)
    if len(blk) > 1:
        body = blk[-1].split("LOEWDIN")[0]
        charge, _, spin = body.partition("\nSPIN")
        for label, part in (("q", charge), ("spin", spin)):
            for lval in "spdfg":
                mm = re.search(rf"\b{lval} :\s+(-?\d+\.\d+)", part)
                if mm:
                    r[f"{label}_{lval}"] = float(mm.group(1))
    return r


def fingerprint(r):
    """A basis-independent label for the electronic state: rounded l-shell
    charge and spin populations."""
    parts = []
    for label in ("q", "spin"):
        for lval in "spdf":
            v = r.get(f"{label}_{lval}")
            # Only occupied shells: a zero p/d/f population just says whether
            # the basis has such functions, which is not a property of the state.
            if v is not None and round(v * 2) != 0:
                parts.append(f"{label}{lval}{round(v * 2) / 2:g}")
    return " ".join(parts) or None


def run_orca(orca, workdir, timeout):
    """Run ORCA in workdir unless a finished output is already there."""
    out = workdir / "orca.out"
    done = parse(out)
    if done["status"] != "missing":
        return done
    t0 = time.perf_counter()
    try:
        with out.open("w") as fh:
            subprocess.run(
                [orca, "orca.inp"],
                cwd=workdir,
                stdout=fh,
                stderr=subprocess.STDOUT,
                timeout=timeout,
            )
    except subprocess.TimeoutExpired:
        pass
    r = parse(out)
    r["seconds"] = round(time.perf_counter() - t0, 1)
    # Keep orca.inp/orca.out (the record) and orca.gbw (needed to share states
    # between bases, removed later unless --keep-gbw).
    for pattern in (
        "*.tmp*",
        "*.densities*",
        "*.bibtex",
        "*.property.txt",
        "*.bas*",
        "guess.gbw",
    ):
        for junk in workdir.glob(pattern):
            junk.unlink(missing_ok=True)
    return r


# --------------------------------------------------------------------------
# Choosing the solution
# --------------------------------------------------------------------------
def converged(rec):
    return rec.get("status") == "ok" and rec.get("scf") not in (None, "")


def select(runs, mult, term):
    """The chosen run among `runs` (one element, method and basis) and its
    flags: the lowest SCF energy of any converged run. <S**2> does not
    disqualify a run; contamination above S2_REL_TOL is flagged."""
    cands = [r for r in runs if converged(r)]
    if not cands:
        return None, "no converged solution"
    best = min(cands, key=lambda r: float(r["scf"]))
    flags = []
    if open_shell_singlet(term):
        flags.append(f"open-shell singlet ({term})")
    S = (mult - 1) / 2
    s2 = best.get("s2")
    if S > 0 and s2 not in (None, ""):
        pure = S * (S + 1)
        excess = (float(s2) - pure) / pure
        if abs(excess) > S2_REL_TOL:
            flags.append(f"<S**2> {float(s2):.3f} vs {pure:.3f} ({100 * excess:+.0f}%)")
    return best, "; ".join(flags)


# --------------------------------------------------------------------------
# The search for one element and method
# --------------------------------------------------------------------------
class ElementSearch:
    def __init__(self, args, Z, symbol, mult, term, method):
        self.args, self.Z, self.symbol, self.mult = args, Z, symbol, mult
        self.term, self.method = term, method
        self.root = Path(args.out) / f"{Z:03d}_{symbol}"
        self.keywords = f"{method} {{basis}} AutoAux {args.grid} TIGHTSCF NoCOSX" + (
            f" {args.extra}" if args.extra else ""
        )
        self.runs = []  # dicts: basis, start, result, dir

    def pbe0_seed(self):
        """PBE0/def2-SV(P) orbitals from a damped PModel guess, once per
        element (shared by all methods)."""
        d = self.root / "seed_pbe0"
        d.mkdir(parents=True, exist_ok=True)
        if not (d / "orca.inp").exists():
            write_input(
                d / "orca.inp",
                "PBE0 def2-SV(P) AutoAux TIGHTSCF NoCOSX SlowConv",
                self.symbol,
                self.mult,
                "  Guess PModel\n",
            )
        r = run_orca(self.args.orca, d, self.args.timeout)
        return d / "orca.gbw" if r.get("status") == "ok" else None

    def smear_seed(self):
        """PBE/def2-SV(P) orbitals from a Fermi-smeared SCF (SmearTemp 5000 K),
        which lets near-degenerate orbitals share electrons before settling --
        a start that explores configurations differently from the others. Used
        from La on, where Hueckel and PAtom are unavailable."""
        d = self.root / "seed_smear"
        d.mkdir(parents=True, exist_ok=True)
        if not (d / "orca.inp").exists():
            write_input(
                d / "orca.inp",
                "PBE def2-SV(P) AutoAux TIGHTSCF NoCOSX SlowConv",
                self.symbol,
                self.mult,
                "  SmearTemp 5000\n  Guess PModel\n",
            )
        r = run_orca(self.args.orca, d, self.args.timeout)
        return d / "orca.gbw" if r.get("status") == "ok" else None

    def one(self, basis, start, guess_block, moinp_src=None):
        d = self.root / safe(self.method) / safe(basis) / safe(start)
        d.mkdir(parents=True, exist_ok=True)
        if not (d / "orca.inp").exists():
            moinp = None
            if moinp_src is not None:
                shutil.copy(moinp_src, d / "guess.gbw")
                moinp = "guess.gbw"
            write_input(
                d / "orca.inp",
                self.keywords.format(basis=basis),
                self.symbol,
                self.mult,
                guess_block,
                moinp,
            )
        r = run_orca(self.args.orca, d, self.args.timeout)
        r["fingerprint"] = fingerprint(r)
        rec = {"basis": basis, "start": start, "dir": d, **r}
        self.runs.append(rec)
        return rec

    def search(self, pool):
        seeds = {"pbe0": self.pbe0_seed()}
        guesses = DIRECT_GUESSES
        if self.Z > EHT_MAX_Z:
            seeds["smear"] = self.smear_seed()
            guesses = [g for g in DIRECT_GUESSES if g not in ("Hueckel", "PAtom")]
        # Round 1: every direct start in every basis.
        jobs = []
        for basis in self.args.bases:
            for name, seed in seeds.items():
                if seed is not None:
                    jobs.append((basis, name, "  Guess MORead\n", seed))
            for g in guesses:
                jobs.append((basis, g, f"  Guess {g}\n", None))
        list(pool.map(lambda j: self.one(*j), jobs))
        # Round 2: share every state found anywhere with every basis.
        states = {}
        for rec in self.runs:
            fp = rec.get("fingerprint")
            if converged(rec) and fp and (rec["dir"] / "orca.gbw").exists():
                best = states.get(fp)
                if best is None or rec["scf"] < best["scf"]:
                    states[fp] = rec
        jobs = []
        for basis in self.args.bases:
            seen = {
                r.get("fingerprint")
                for r in self.runs
                if r["basis"] == basis and converged(r)
            }
            for fp, src in states.items():
                if fp not in seen:
                    tag = f"state_{safe(src['basis'])}_{safe(src['start'])}"
                    gbw = src["dir"] / "orca.gbw"
                    jobs.append((basis, tag, "  Guess MORead\n", gbw))
        list(pool.map(lambda j: self.one(*j), jobs))
        return len(states)

    def choose(self, basis):
        return select(
            [r for r in self.runs if r["basis"] == basis], self.mult, self.term
        )


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def element_table(db_path):
    c = sqlite3.connect(db_path)
    return {
        Z: (sym, mult, term)
        for Z, sym, mult, term in c.execute(
            "SELECT atomic_number, symbol, multiplicity, term_symbol FROM element"
        )
    }


def parse_elements(text, table):
    out = []
    for part in text.replace(",", " ").split():
        if "-" in part:
            lo, hi = part.split("-")
            out += range(int(lo), int(hi) + 1)
        elif part.isdigit():
            out.append(int(part))
        else:
            out += [Z for Z, (sym, _, _) in table.items() if sym == part]
    return [Z for Z in out if Z in table]


def default_orca():
    ini = configparser.ConfigParser()
    ini.read(Path("~/SEAMM/orca.ini").expanduser())
    for section in ini.sections():
        code = ini.get(section, "code", fallback="")
        if code:
            return code
    return shutil.which("orca")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--elements", default=None, help="e.g. '1-36', '26 Kr'")
    p.add_argument(
        "--methods",
        nargs="+",
        default=["REVDSD-PBEP86-D4/2021", "PWLDA", "VWN", "VWN3"],
    )
    p.add_argument("--bases", nargs="+", default=DEFAULT_BASES)
    p.add_argument("--grid", default="DEFGRID3")
    p.add_argument("--extra", default="", help="extra '!' keywords")
    p.add_argument("--orca", default=default_orca())
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--timeout", type=float, default=4 * 3600, help="s per run")
    p.add_argument("--out", required=True)
    p.add_argument(
        "--rechoose",
        action="store_true",
        help="only rebuild choices.csv/Results.csv from the runs.csv in --out",
    )
    p.add_argument("--db", default=None, help="ThermoDB for the element table")
    p.add_argument("--keep-gbw", action="store_true", help="keep all orbital files")
    args = p.parse_args()

    if args.db is None:
        from seamm_thermochemistry import DEFAULT_DB_PATH

        args.db = DEFAULT_DB_PATH
    table = element_table(args.db)
    out = Path(args.out).expanduser()
    if args.rechoose:
        return rechoose(out, table)
    if args.elements is None:
        p.error("--elements is required unless --rechoose is given")
    Zs = parse_elements(args.elements, table)
    out.mkdir(parents=True, exist_ok=True)
    args.out = str(out)

    results = {}  # Z -> {column: value}
    all_runs, choices = [], []

    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        for Z in Zs:
            sym, mult, term = table[Z]
            for method in args.methods:
                es = ElementSearch(args, Z, sym, mult, term, method)
                n_states = es.search(pool)
                for rec in es.runs:
                    all_runs.append(
                        {
                            "Z": Z,
                            "El": sym,
                            "method": method,
                            **{k: rec.get(k) for k in RUN_FIELDS[3:]},
                        }
                    )
                for basis in args.bases:
                    best, flags = es.choose(basis)
                    record(
                        results, choices, table, Z, method, basis, best, flags, n_states
                    )
                if not args.keep_gbw:
                    for rec in es.runs:
                        (rec["dir"] / "orca.gbw").unlink(missing_ok=True)
                print(f"{Z:3d} {sym:2s} {method}: {n_states} state(s)", flush=True)
            # Rewrite the outputs after every element, so partial runs are usable.
            write_csv(out / "runs.csv", RUN_FIELDS, all_runs)
            write_outputs(out, choices, results)


RUN_FIELDS = [
    "Z",
    "El",
    "method",
    "basis",
    "start",
    "status",
    "scf",
    "total",
    "s2",
    "fingerprint",
    "seconds",
    "dir",
]
CHOICE_FIELDS = [
    "Z",
    "El",
    "Multiplicity",
    "Term",
    "method",
    "basis",
    "start",
    "energy_kJ",
    "scf",
    "s2",
    "n_states",
    "flags",
]


def record(results, choices, table, Z, method, basis, best, flags, n_states):
    """Add one chosen entry to its Results.csv row and to choices.csv."""
    sym, mult, term = table[Z]
    row = results.setdefault(
        Z, {"Atomic Number": Z, "Element": sym, "Multiplicity": mult}
    )
    col = f"DFT@{db_method(method)}/{basis}"
    energy = round(float(best["total"]) * KJ_PER_EH, 3) if best else ""
    if best is not None:
        row[f"E {col} (kJ/mol)"] = energy
        row[f"S^2 {col}"] = best.get("s2", "")
    choices.append(
        {
            "Z": Z,
            "El": sym,
            "Multiplicity": mult,
            "Term": term,
            "method": method,
            "basis": basis,
            "start": best["start"] if best else "",
            "energy_kJ": energy,
            "scf": best["scf"] if best else "",
            "s2": best.get("s2", "") if best else "",
            "n_states": n_states,
            "flags": flags,
        }
    )


def write_outputs(out, choices, results):
    write_csv(out / "choices.csv", CHOICE_FIELDS, choices)
    fixed = ["Atomic Number", "Element", "Multiplicity"]
    cols = fixed + sorted(
        {k for r in results.values() for k in r} - set(fixed),
        key=lambda k: (k.split(" ", 2)[1], k.split(" ")[0]),
    )
    write_csv(out / "Results.csv", cols, [results[z] for z in sorted(results)])


def rechoose(out, table):
    """Rebuild choices.csv and Results.csv from runs.csv with `select`."""
    with open(out / "runs.csv", newline="") as fh:
        runs = list(csv.DictReader(fh))
    groups = {}
    for r in runs:
        groups.setdefault((int(r["Z"]), r["method"]), []).append(r)
    results, choices = {}, []
    for (Z, method), recs in sorted(groups.items()):
        _, mult, term = table[Z]
        n_states = len(
            {r["fingerprint"] for r in recs if converged(r) and r["fingerprint"]}
        )
        for basis in dict.fromkeys(r["basis"] for r in recs):
            best, flags = select([r for r in recs if r["basis"] == basis], mult, term)
            record(results, choices, table, Z, method, basis, best, flags, n_states)
    write_outputs(out, choices, results)
    print(f"rechose {len(choices)} entries from {len(runs)} runs in {out}")


def write_csv(path, fields, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
