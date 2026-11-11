2026-11-10 -- Fragment store P5: write-back, in-job reuse, two layers, ATcT download
====================================================================================

Status: design note, written 2026-11-10 (design). **Nothing coded.** Paul's
decisions 1-5 below are settled. Science reviewed the note the same evening, and
its points are folded in:

- the per-species expected fingerprint (§2a);
- the structure guard and distinct conformers (§3);
- R5 (§4);
- the uncertainty-method column (§5);
- the official sanity pass (§6);
- the validation additions.

**Decided 2026-11-10 (Paul): D1-D8 accepted as recommended** (see `D-items
for Paul`_). Implementation can start with P5a.

**P5a implemented 2026-11-11** (branch ``p5a-store-reuse`` in seamm_thermochemistry,
orca_step, seamm_exec and thermochemistry_step; not yet released). Found while
doing it: the step kept the fingerprints of the *last* species it computed as
the ones to look up with, so after a boron fragment every later lookup in the
iteration missed; it now keeps the target's. And an energy-level single point
found on a configuration in the job is reused only when its record names the
configuration's current geometry.

**P5a released 2026-11-11** (orca-step, seamm-exec, seamm-thermochemistry and
thermochemistry-step 2026.11.11). **P5b implemented 2026-11-11** (branch
``p5b-writeback`` in seamm_thermochemistry, thermochemistry_step and
seamm_jobserver): write-back, the structure guard, ``auto:`` conformers, the
``job_finished`` hook and the marker, as in §3. One change from §3: a geometry
that is a contact too close to read as a bond is caught by the contact check,
but ligands closer still (the old PF₆⁻) read as bonded and are refused by the
connectivity check first.

**P5b released 2026-11-11** (seamm-jobserver 2026.11.11, seamm-thermochemistry and
thermochemistry-step 2026.11.11.1). **P5c implemented 2026-11-11** (branch
``p5c-layers``): the meta table (schema 3), ``LayeredStore`` with R1-R5,
``migrate``, ``merge``, ``export``, ``official build`` and ``duplicates
--layers``, and the step reading both layers. Recording the official file's
sha256 in the installation, and ``installer check`` verifying it, moved to P5d
with the rest of the installer; the official build's sanity pass is P5e.
``migrate`` was run on a copy of ChemAI's store (2026-11-11 06:14): 100
fragment, 50 geometry and 639 energy rows normalized, all keys unchanged, a
second run a no-op.

**P5c released 2026-11-11** (seamm-thermochemistry and thermochemistry-step
2026.11.11.2). **P5e code implemented 2026-11-11** (branch ``p5e-official``): the
sanity pass, selection by fingerprint, exclusions with reasons and the report.
The multiplicity check uses ``cbh.SINGLET_GROUND_STATES`` (design's, from the CBH
review: SiH2, SiCl2, SiF2, CF2, CCl2, BH) and ATcT's labelled spin states; ATcT alone cannot catch BH (only a
plain entry) or SiCl2 (absent). Writing it showed that the P5b structure guard
refused singlet BH and carbenes (RDKit gives them two radical electrons); it now
checks only the parity, and it held fluorines on one atom (CF2, CF3) to the 2.2 Å
meant for fluorines on different atoms (now 2.0 Å). Building v1 waits for the DEFGRID3 harvest into ChemAI's
store: frequency key sha256:6b46c8c4..., energy keys sha256:1e16835c... and
sha256:7ca1cbbc... (boron), the [CH3-] saddle from 5684 excluded.

**Official v1 built 2026-11-11** (``~/Work/SEAMM/official_store_v1/``, from a
read-only copy of ChemAI's settled store, sha256 2303508002b9..., migrated):
458 species states, 597 geometries, 1078 energies (r2SCAN-3c, revDSD 1e16 and
7ca1), 912 fragment records. The first build refused 11 species states; four were
false positives of the bond perception (the diagonals of azete, oxetene and the
N2O3 cage; one long N-N bond), which led to perceiving bonds at several
covalent-radius factors. Six exclusions remain, with reasons in
``exclude-v1.csv`` and the report: the [CH3-] saddle; C3 (ring opened); C1=CCC=1
(became bicyclic); and three dissociated conformers, two of them the lowest of
their species (N2O3 as NO2 + NO, Cl2O2 as two ClO), which lookups would otherwise
have used. To publish: rebuild with the released seamm-thermochemistry, then Paul
uploads to Zenodo.

**Released 2026-11-11:** seamm-thermochemistry 2026.11.11.4 (design's CBH fixes)
and 2026.11.11.5 (P5e: official-build checks, structure-guard fixes). ChemAI's
store was then quarantined by hand (Paul's OK): the six problem rows, with
``UNUSABLE: <reason>`` as model chemistry; 2026.11.11.6 makes the tools honour
that, and a dry run on a copy of the quarantined store gives the same v1 (458
species states; the quarantine covers all six exclusions). A ``quarantine``
column (schema 4) would be cleaner than the convention; not needed for v1.

Numbering: these phases are the fragment store's P5 (P5a-P5e). The 11-08
presets plan also had a P5 (long-bond and %[(T)] flags) and a P6 (pilot). Those
keep their names and have nothing to do with this note.

Why
---

The store (P1-P4, seamm-thermochemistry 2026.11.8.x) is read by the step but
filled only by hand with ``seamm-thermochemistry-fragments harvest``. ChemAI job
5684 (the 1694 set, the default preset, a **serial** Loop, ``parallel: 'no'``)
shows what that costs:

- By 18:53 on 11-10, 188 iterations had run 675 species calculations, but only
  ~170 distinct species. That is about 500 repeats (~75 %): H₂ 114 times, CH₄ 74,
  HF 58, C₂H₆ 49, CH₃F 39, HCl 35, ...
- Every repeat already existed in the job's own database as ``CBH:<SMILES>``
  with its fingerprints. The step never looks there: ``_compute_species`` finds
  the system, then optimizes it and computes its Hessian again.
- When Science harvested 5684 into the ChemAI store by hand (18:35),
  iterations 178 and later got CH₄, H₂, C₂H₆ and C₃H₈ from the store ("store;
  energy store", same E). Harvesting during a job does give reuse within the
  job, so write-back would have done the same automatically.
- The 598 seed rows were not reused, for two reasons:

  - **A real grid difference.** The sub-step rows were run with
    ``! R2SCAN-3C TIGHTSCF`` (DEFGRID2). The engine adds ``TIGHTSCF DEFGRID3``
    (``orca_step/data/orca_mdi.py`` ``ACCURACY_KEYWORDS``). The fingerprint guard
    is correct. Paul's decision on the grid is "(a) now, then (b)": the
    DEFGRID3 recompute is under way.
  - **Different names for one model chemistry.** The sub-step rows say
    ``CMP@R2SCAN-3C``, but the step asks for ``ORCA:CMP@R2SCAN-3C``. This
    matters only for the unverified fallback, but it will bite every
    name-based lookup.

- **A third miss, by elements.** The step looks up every species'
  energy-level energy with the **target's** settings fingerprint.
  ``_energy_level_energies`` passes the target's record fingerprint to
  ``store.energy``. But ORCA chooses exact exchange (``NoCOSX``) per calculation
  by its own elements: a one-center job, or any atom of Na, Mg, Zn, B or P
  (orca_step#44). So the stored revDSD energies of boron species have
  fingerprint ``sha256:7ca1cbbc…``, and those of ordinary species have
  ``sha256:1e16835c…``.

  - CBH fragments are sub-structures of the target, so the miss is for the
    fragments **without** those elements under a target **with** them. A B or P
    target looks up its CH₄, H₂ and HF with ``7ca1…`` and misses. They are then
    computed with COSX and stored as ``1e16…``, so they miss again in every
    such job.
  - P matters most: PF₆⁻, phosphates and the P-F references.
  - The geometry key is not split: r2SCAN-3c has no exact exchange, so COSX is
    moot and is dropped from its fingerprint (``sha256:6b46…``). A hybrid
    geometry level would split the same way.

The repeats are mostly small molecules, so the saving in core-hours is less
than 75 %. Validation run V1 will measure it.

Paul's decisions (settled)
--------------------------

1. **Write-back.** The step adds each species it computes to the store: the
   geometry and frequency row, and an energy row per level. The existing guards
   apply: imaginary modes above 20i are rejected, the pairing and same-level
   guards hold, and duplicates agree within 0.1 kJ/mol. Jobs on remote clusters
   are harvested when they come back, by a JobServer stage-out hook or the CLI.
   A merge tool combines copies of the store on different installations.
2. **In-job reuse.** "design should check that a job's fragments are reused".
   Each iteration checks the store as it runs, so write-back gives reuse within
   a job. Duplicates computed by parallel iterations are acceptable. An
   optional "in progress" marker per species and fingerprint would let a
   parallel iteration wait instead.
3. **Two layers.** "keep two databases, our 'official' one plus the user's
   local one, and preferentially get from the user's".

   - The official layer is on Zenodo: read-only and versioned. Its sha256 is
     recorded in the installation, and the installer replaces it whole, as it
     does ``thermochemistry.db``.
   - The local layer is writable, filled only by write-back and harvest, and
     never overwritten.
   - Lookup goes to local first, then official.

   Science's refinements:

   - An official "superseded" mark beats an identical old local row (logged).
   - The same species and fingerprint differing by more than 0.1 kJ/mol is
     logged as a duplicate.
   - The lowest conformer is chosen across both layers, but energies always come
     from the same layer as their geometry.

   A CLI exports the local rows that are not in official, with their provenance,
   for the next official version. The single ``fragments.db`` now on each
   machine becomes its local layer.
4. **ATcT.** "ATcT is open access, so we should automatically download what is
   needed so the user doesn't have to do anything."

   - The step or the installer downloads the ATcT version when it is missing.
     The never-modified-copy rule holds, and the metadata and SHA-256 are kept.
   - The covariance default becomes "fetch if needed", falling back to
     quadrature when offline.
   - One installer run sets up ATcT, ``thermochemistry.db`` and the seed store.
5. **Model-chemistry names.** Lookup, write-back and harvest all use the bare
   model plus the program column (``CMP@R2SCAN-3C`` + ``ORCA``).

Design
------

1. Model-chemistry names (P5a)
''''''''''''''''''''''''''''''

``seamm_thermochemistry.fragment_store.split_model_chemistry(text, program=None)``
returns ``(model, program)``:

- The text splits at its first ``:`` only when that colon comes before the first
  ``@`` and the prefix has no ``@``. So ``ORCA:CMP@R2SCAN-3C`` gives
  (``CMP@R2SCAN-3C``, ``ORCA``), and ``ORCA:DFT@B3LYP/bse:def2-SVP`` gives
  (``DFT@B3LYP/bse:def2-SVP``, ``ORCA``).
- ``QC@CCSD(T)-F12D_RI/cc-pVTZ-F12`` is left whole.
- Programs are compared case-insensitively and stored as the calculation record
  writes them (``ORCA``).
- An explicit ``program`` that contradicts the prefix is logged, and the explicit
  one is kept.

It is applied everywhere a model chemistry enters or queries the store:

- ``FragmentRecord.make``, ``GeometryRecord.make`` and ``EnergyRecord.make``
  store the bare model and fill ``program``.
- In ``harvest_job``, the property name's model field is split. The program
  comes from the ``calculation`` record, or else from the prefix. The
  property's second field is a step name (``Thermochemistry``), never a
  program.
- ``get``, ``geometries``, ``energy`` and ``species_at`` normalize the model
  they are given and match ``model_chemistry = ? AND (program IS NULL OR
  program = ?)``. Verified lookups by fingerprint are unaffected.
- Unverified keys become ``unverified:<program or ?>:<model>@<version or ?>``,
  one spelling for one level.

2. In-job reuse (P5a, with P5b)
'''''''''''''''''''''''''''''''

Before computing a species, the step tries three sources in order: the job's
database, the store (local, then official, see §4), and finally a computation.

**The job's database** (new, in ``cbh_approach``: ``_from_job(system_db,
can)``):

- It looks at every system named ``CBH:<can>``; a snapshot merge can leave two.
  For each one it checks the current configuration first, then the others.
- **Geometry level.** It takes the properties ``total energy``, ``ZPE`` and
  ``H(298.15 K)`` at ``#Thermochemistry#<model>``, but only when the stored
  ``frequency fingerprint`` (split levels) or ``settings fingerprint`` (one
  level) equals the fingerprint of this run. The configuration must also have
  no imaginary mode above 20i.
- **Energy level.** It takes ``total energy#Thermochemistry#<level>`` when the
  ``settings fingerprint#Thermochemistry#<level>`` is one of **this species'**
  expected fingerprints (§2a), not the target's.
- It reuses the existing configuration. No new "from the fragment store"
  configuration is created, as ``_configuration_from_geometry`` does today.
- A fingerprint must match, so results from another code version or other
  settings are never reused silently. The subflowchart route stores no
  fingerprints, so it gets no in-job reuse from the job's database; write-back
  covers it.
- The fingerprints are known only after the step's first engine calculation,
  which is always the target's. This is the same constraint the store lookup
  has now.

This alone fixes 5684-type serial loops, and it needs no store at all.

2a. Per-species expected fingerprints (P5a; Science 11-10)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Each species is looked up, in the job's database and in the store, with the
fingerprints that ORCA would produce **for that species**: NoCOSX or not, chosen
by its own elements. This applies to the energy level and to the geometry
(frequency) level. Only the elements change the settings this way, because
geometry, charge and multiplicity are removed before hashing
(``orca_step.fingerprint``).

**Interface.** The code plug-in predicts; the step never re-implements ORCA's
rules.

- **orca_step**: ``batch.expected_fingerprints(model_chemistry, symbols, *,
  kind="energy" | "hessian", version, options=None)``.

  - It builds the ``!`` line and blocks with the same helpers that
    ``get_task`` and the MDI engine use (``engine_helpers``,
    ``exact_exchange_method``), and hashes them with
    ``orca_calculation_record``.
  - It returns a small set: the fingerprint as written, and the one with COSX
    and NoCOSX dropped. Whether COSX is moot is known only from ORCA's output
    (``cosx_is_moot``), so the stored fingerprint is always one of the two.
  - ``version`` comes from the target's calculation record in the same run (one
    ORCA per run).

- **seamm_exec**: ``Evaluator.expected_fingerprints(configuration, kind)``
  calls the code's function when it has one. For codes without it, it returns
  ``None``, and the step falls back to today's rule (the target's fingerprint).

- **thermochemistry_step**: ``_energy_level_energies`` and ``_from_job`` (and
  the store lookups) pass each species' set. The store's ``energy`` and
  ``geometries`` accept a set of fingerprints.

**Self-check.** After computing a species, the step compares the recorded
fingerprint with the predicted set. A miss means the prediction has drifted from
ORCA's rules, for example after an orca_step change. The step prints a warning
naming the species and both fingerprints. The species is still used and written
back under its recorded fingerprint, so a drift costs reuse, never correctness.

**The store, per iteration** (exists). With write-back (§3), each species is in
the local layer seconds after any iteration computes it. That covers what the
job's database cannot:

- iterations of a **parallel** Loop, which each get a snapshot ``seamm.db``
  holding only their own structures;
- later jobs.

**Parallel duplicates.** The first wave of a parallel Loop (the iterations
running at once) computes the shared fragments concurrently, and later
iterations read them. The waste is bounded by iterations-at-once × shared
fragments. A claim ("in progress") table could remove it; see D2.

3. Write-back (P5b)
'''''''''''''''''''

**One code path.** The per-configuration body of ``harvest_job`` becomes
``harvest_configuration(store, configuration, species, *, charge,
multiplicity, conformer, job, ...)``. ``harvest_job`` loops over systems and
calls it. The step calls it right after each species it **computes**: the
target, every fragment, and the "species only" mode (46e2c3f). Write-back and a
later harvest of the same job therefore give the same rows, and the later
harvest finds only duplicates. Species that came from the store or from the
job's database are not written again.

**Guards:**

- The pairing guard and the same-level guard are as now.
- Writes are insert-or-ignore, and differences over 0.1 kJ/mol are logged, never
  overwritten.
- Write-back **refuses** a geometry with an imaginary mode above 20i, together
  with the energy rows on it, and says so in the step output. Paul: "rejected".
- The explicit ``harvest`` CLI keeps today's behaviour: such a geometry is stored
  for provenance and never returned for reuse.
- **Structure guard** (Science 11-10). Before a computed species is written,
  its bonds are perceived again from the optimized coordinates, not taken from
  the configuration's bonds, which come from the SMILES.

  - The canonical SMILES without stereo, the charge and the multiplicity must
    equal the key, also compared without stereo.
  - Unusual centres (hypervalent, metal, high coordination) must have a
    shortest ligand-ligand distance above the molsystem threshold for atoms too
    close, as for PF₆⁻.
  - Otherwise the species is refused and the step says why. A fragment that
    rearranged during its optimization would otherwise become an ``auto:``
    "lower conformer" (D5) and win every later lookup. Hypervalent and
    zwitterionic fragments are the likely cases; the 1694 set has 16
    hypervalent hydrides.
  - The job keeps using what it computed; only the store refuses it.
  - Perception that fails (no connectivity found) also refuses.
  - The same check runs in ``harvest_configuration``, so the CLI applies it too.
    ``--no-structure-check`` overrides it, with the reason recorded in
    ``source``.

**When the step writes:**

- The parameter ``fragment store`` gets three choices:

  - "use and add what is computed" (the new default);
  - "use only" (today's "use if available");
  - "ignore".

- The local layer is opened writable only if it is on a **local filesystem of
  the machine running the step**. ``seamm_thermochemistry.fragment_store.
  is_local_filesystem(path)`` reads the mount table: ``/proc/mounts`` on Linux,
  ``getmntinfo`` through ``os.statvfs``/``stat -f`` on macOS.

  - zfs, apfs, ext4, xfs and btrfs count as local.
  - nfs*, gpfs, lustre, beegfs, ceph, cifs/smbfs and fuse.sshfs, and anything
    unknown, count as not local.
  - ``[thermochemistry] fragment-store-writable = auto | yes | no`` in seamm.ini
    overrides the check.

- This keeps the reason for the 10-28 rule (SQLite locking is unreliable on
  GPFS/NFS; the 5398 SIGBUS crashes) and drops its letter ("a single writer;
  never from a job"). See D1.
- Several iterations on one host may write at once. Each write is one short
  transaction per species, using the rollback journal and the existing 30 s
  busy timeout.
- A write that still fails is logged and skipped. It never fails the species,
  and the harvest picks it up later.
- When the step does not write, the properties stay on the configurations, as
  now (they are the harvest's input). The step prints "not added to the store
  (<fstype>); harvested when the job returns" and writes the marker
  ``<job>/harvest_fragments.json`` (species, levels, why).

**Conformers (D5).** The geometry key has ``conformer = ''``, so today a lower
conformer found later is only logged. The proposal:

- Write-back stores a computed geometry as a new conformer,
  ``auto:<geometry hash[:12]>``, when both of these hold:

  - its E + ZPE is lower than every stored conformer's by more than
    0.1 kJ/mol;
  - it is **structurally distinct** from each of them: a symmetry-aware Kabsch
    RMSD above 0.05 Å (RDKit ``GetBestRMS`` over the graph).

- Lower alone is not enough: finite-difference ZPE noise on the same minimum
  can reach about 0.1 kJ/mol (Science).
- It must also pass the structure guard above.
- ``_geometry_from_store`` already takes the lowest E + ZPE.
- A repeat of the same minimum, or a higher conformer, goes to
  ``duplicate_log`` / ``row_history`` as now. Nothing is overwritten.

**Remote jobs: the JobServer hook (D6).**

- A generic hook in seamm_jobserver: after ``_finalize_job_status`` (that is,
  after stage-out, so the files are local), the JobServer runs the entry points
  in ``org.molssi.seamm.jobserver.job_finished``.
- They run in a background thread with a time limit, and are called as
  ``hook(job_id=, wdir=, state=, logger=)``.
- Exceptions are logged and never change the job's status. Nothing
  thermochemistry-specific goes into the JobServer.
- The hook must never block the poll loop: compare the 28 GB stage-in that
  stalled ChemAI for 3.5 h.
- seamm-thermochemistry registers ``harvest``, which:

  - acts only when ``harvest_fragments.json`` exists, in any terminal state,
    since finished iterations are valid even if the job failed later;
  - copies ``seamm.db`` with SQLite's backup API first (the job database is
    WAL), as Science did by hand for 5684;
  - harvests into the host's local layer, and appends the counts to the job's
    ``job.out``.

Jobs run without a JobServer (``run_flowchart`` by hand on ARC) are harvested
with the existing CLI on the host that holds the local layer.

4. Two layers (P5c)
'''''''''''''''''''

**Files**, in ``<root>/Parameters/thermochemistry/``:

``fragments.db``
    The local layer. The name is unchanged, so today's files become local layers
    where they are. It is the only file write-back, harvest, merge and
    supersede write to.

``fragments-official.db``
    The official layer. The installer installs it and replaces it whole. It is
    opened ``mode=ro&immutable=1`` (no locking), so clusters can read it in place
    on GPFS.

Both files get a ``meta`` table with these keys: ``layer`` (local | official),
``version`` (official: CalVer), ``created``, ``schema`` and ``source`` (the
Zenodo DOI). The official file's sha256 goes in seamm.ini ``[thermochemistry]
fragment-store-sha256``, next to ``database-path``, and ``installer check``
verifies it. The schema becomes version 3 (the ``meta`` table, plus ``claim``
if D2 is (b)). Readers of schema-2 files keep working.

**``LayeredStore(local, official)``** (new ``layered_store.py``):

- It has the reader API the step uses: ``geometries``, ``energy``,
  ``species_at`` and ``get``.
- Every row it returns carries ``layer``.
- Writes go to local only.
- It refuses a file whose ``meta.layer`` is wrong for its role.

**Lookup rules:**

R1
    The same key in both layers: the **local** row is used, except under R2.
R2
    The official row for the key has been superseded (its ``superseded`` or
    ``row_history`` kind ``superseded``), and the local row equals one of the
    superseded old values within 0.1 kJ/mol. Then the **official** current row
    is used, and the step logs "local <species> is the value official <version>
    superseded on <date>: <reason>".
R3
    The same key, values differing by more than 0.1 kJ/mol, and not R2: local is
    used and the conflict is logged in the step output. The new CLI
    ``duplicates --layers`` lists such conflicts across the two files.
R4
    Conformers: the candidates come from both layers, after the imaginary
    filter. The lowest E + ZPE wins, and a tie within 0.1 kJ/mol goes to local.
    Verified rows in either layer come before unverified rows in either.
R5
    Energies pair with geometries only on the **same geometry hash**, never
    across hashes. They are looked up in the geometry's own layer first, then
    in the other layer. The energy fingerprint must still match (§2a). This is
    Science's answer of 11-10, and it lets a user's revDSD single point on an
    official r2SCAN-3c geometry be reused. Science attached two conditions:

    (i) The returned energy row carries its own ``layer``, and the step output
    shows any mix, e.g. "geometry official 2026.12.1, energy local".

    (ii) If a later official version supersedes that geometry with new
    coordinates, its hash changes and the local energy is simply not found. A
    test covers this.

**Export, merge and build (CLI ``seamm-thermochemistry-fragments``):**

``export --out contrib.db [--label TEXT]``
    Writes the local rows whose key is not in official, and rows that differ
    from official by more than 0.1 (flagged), with every provenance column and
    their ``row_history``. ``--label`` is free text the user picks. No host or
    user names are added automatically.
``merge SOURCE.db ... [--store TARGET]``
    Insert-or-ignore of every fragment, geometry and energy row into the target
    (by default the local layer), with provenance kept and ``source`` noting the
    origin label. Differing duplicates are logged, and supersede history is
    copied as history rows. It refuses an official file as the target. It is
    used for copies on several installations (Mac and ChemAI), for folding a
    campaign store into a local layer, and for building the next official.
``official build --version V --previous OLD.db --from A.db B.db ... --out NEW.db``
    Starts from the previous official and merges the contributions. A differing
    duplicate stops the build until the maintainer resolves it with
    ``supersede`` and a reason. Then it writes ``meta``, runs VACUUM and writes
    a ``.sha256`` sidecar. Paul publishes the file on Zenodo, as he does for
    ``thermochemistry.db`` (a new concept record).
``migrate [--store FILE]``
    Makes a dated backup, adds ``meta``, normalizes model names (§1) and fills
    in ``program``. It reports any key collision, and running it twice changes
    nothing.

5. ATcT and the installer (P5d)
'''''''''''''''''''''''''''''''

**One installer run.** ``seamm-thermochemistry-installer install`` and
``update``, which seamm-manager runs, set up:

1. ``thermochemistry.db`` (as now);
2. ``fragments-official.db`` (new Zenodo concept id; checksum verified,
   sha256 recorded);
3. the newest ATcT version, if none is cached (``download_all``).

``check`` reports all three. The installer never touches ``fragments.db``.

**ATcT download:**

- From ``atct.anl.gov``, with SEAMM's own User-Agent and no spoofing.
- ``metadata.json`` already keeps the source, date, sha256 and the number of
  entries. An existing version is never overwritten.
- ``check`` recomputes the sha256 against ``metadata.json``, to catch a
  modified copy.
- "latest" stays "the newest **cached** version", so results do not change
  silently. A newer TN version arrives only through an explicit ``update``, and
  old versions are kept.

**At run time:**

- If the requested version, or any version for "latest", is missing, the step
  downloads it when the machine has a network. Without one, the error names the
  installer command.
- Two iterations may race to download. ``download_all`` must then treat
  "already there" as success: today it raises ``ATcTError``.

**Access:** In October, ``/api/v1/all/`` answered 403 to non-browser clients.
On 2026-11-10, one request to the API root with SEAMM's User-Agent answered 200;
``/all/`` was not tried. Hence D4.

**Covariances:**

- The default becomes "fetch if needed". When the fetch fails (offline, 403 or
  429), the independent uncertainties are combined in quadrature. The step says
  so as a warning, not an error.
- The results also record how each uncertainty was combined (Science 11-10),
  so that a Loop's table cannot silently mix the two:

  - a column "uncertainty method" in the ``CBH`` table and the results, with the
    values "covariance" or "quadrature (no covariance: <reason>)";
  - the matching result property, so it reaches the results CSV and the
    dashboard.
- ``covariance.json`` is written atomically and merged on write: re-read it,
  take the union of the blocks, then ``os.replace``. Parallel iterations write
  it.
- A read-only cache keeps the blocks in memory, as now.
- ``atct_all.json`` is never modified. ``covariance.json`` is SEAMM's own cache
  and only grows.

6. Official v1 (P5e)
''''''''''''''''''''

Official v1 is built with ``official build`` from the ChemAI local layer after
the DEFGRID3 recompute. Its contents are D3. Paul publishes it on Zenodo, and
the installer's concept id is filled in by a seamm-thermochemistry release.

**Sanity pass before publishing** (Science 11-10; part of ``official build``,
which refuses to write the file while any check fails, unless a row is
explicitly excluded with a reason):

- **Imaginary modes:** no geometry row with a mode above 20i.
- **Multiplicity:** checked against the ATcT species where ATcT has the
  species. BH and SiCl₂ are singlets, but the 1694 set had them as triplets, so
  no triplet rows for them may get in.
- **Connectivity:** the structure guard of §3, run over every row: bonds
  perceived again, SMILES without stereo, charge and multiplicity.
- **Unusual centres:** the shortest ligand-ligand distance.
- **Ions and doublets** from 5684 belong in v1: they are checked, not excluded.

The report (rows checked, refused, why) is published with the file.

D-items for Paul
----------------

D1 -- Write-back from a running job.
    The step writes the local layer directly when it is on a local disk of the
    machine running the step, and otherwise leaves the species to the harvest.
    This relaxes the 10-28 rule "no store writes from jobs; one writer" to "no
    SQLite writes over network filesystems; several writer processes on one
    host are fine". **Recommended: yes.**
D2 -- "In progress" marker.
    \(a) None. (b) A ``claim`` table: an iteration that finds a fresh claim polls
    for up to a time limit, then computes; stale claims are those whose process
    is dead or whose heartbeat is too old.
    **Recommended: (a).** The job's database covers serial loops. Parallel
    duplicates are limited to the first wave, and claims bring stale-lock
    handling across processes and hosts. Revisit if V2 shows a large waste.
D3 -- Official v1 contents.
    \(a) Only rows the presets can reach: the engine (DEFGRID3) r2SCAN-3c
    geometries with frequencies, and the revDSD energies on them, from the
    DEFGRID3 recompute plus 5684. That is about 600 species.
    \(b) Also the DEFGRID2 sub-step rows and the campaign stores (F12, B2PLYP,
    PBE0, GFN2).
    **Recommended: (a).** Rows no preset can reach only add size and confusion.
    \(b) can follow in v2 if wanted.
D4 -- ATcT source.
    \(a) ``atct.anl.gov`` only. When it refuses, a clear message and the existing
    ``download --from FILE``. (b) As (a), plus a SEAMM mirror on Zenodo of the
    unmodified ``atct_all.json`` as a fallback. That needs ANL's permission, a
    question for the Ruscic/Bross contact that is already open.
    **Recommended: (a) now, and ask about (b).**
D5 -- A lower conformer found by write-back.
    \(a) Added as a new conformer row (``auto:<hash>``), and the lowest is used.
    \(b) First wins, and the lower one is only logged (today).
    **Recommended: (a)**, with the structure guard and the 0.05 Å distinctness
    test of §3. It matches "lowest conformer across layers" within a layer, and
    nothing is overwritten.
D6 -- Harvesting remote jobs.
    \(a) A generic JobServer ``job_finished`` hook plus the CLI. (b) The CLI only.
    **Recommended: (a).** Without it, results from ARC queues reach no store
    unless someone remembers.
D7 -- Clusters.
    ARC's ``Parameters`` is on GPFS. ARC gets the official layer only, read in
    place (immutable), and no local layer. ARC results reach ChemAI's local
    layer through the hook or a harvest or merge on ChemAI.
    **Recommended: yes.**
D8 -- Names.
    ``fragments.db`` is local and ``fragments-official.db`` is official.
    ``fragment store`` gets the choices "use and add what is computed"
    (default), "use only" and "ignore". ``fragment store file`` names the local
    layer. **Recommended: yes.**

**Paul's decision (11-10): all eight as recommended.** D1 yes, D2 (a), D3
(a), D4 (a) and ask ANL about (b), D5 (a) with the structure guard, D6 (a),
D7 yes, D8 yes.

**Science's view (11-10):** D2 (a), D3 (a), D4 (a) now and ask ANL about (b),
D5 (a) with the structure guard, D6 (a), D7 yes. D1 and D8 are engineering
calls, and Science has no objection to them. Paul decides.

Science's answer on R5
----------------------

Energies may pair across layers on the same geometry hash: the geometry's own
layer first, then the other. They never pair across hashes, and the energy
fingerprint must match. Science added two conditions, the visible layer mix and
a test for a superseded geometry (now in R5 in §4).

File-level changes
------------------

seamm_thermochemistry:

- ``fragment_store.py``:

  - ``split_model_chemistry``, applied in the three ``make``\ s, ``get`` and
    ``harvest_job`` (P5a);
  - ``harvest_configuration`` factored out of ``harvest_job``;
  - ``is_local_filesystem``;
  - the imaginary refusal for write-back;
  - ``check_structure`` (bonds perceived again, SMILES without stereo, charge
    and multiplicity, unusual centres; P5b), which is also used by
    ``official build`` (P5e);
  - the RMSD distinctness test for ``auto:`` conformers (P5b);
  - the ``meta`` table, schema 3, the ``immutable`` read-only open and
    ``migrate`` (P5c).

- ``geometry_store.py``:

  - normalized ``geometries``, ``energy`` and ``species_at``, which also accept
    a set of fingerprints (P5a);
  - the ``auto:`` conformer path (P5b);
  - the ``layer`` field on returned rows (P5c).

- ``layered_store.py`` (new): ``LayeredStore`` with rules R1-R5 (P5c).
- ``fragment_cli.py``: ``export``, ``merge``, ``official build``, ``migrate``
  and ``duplicates --layers`` (P5c).
- ``hooks.py`` (new): ``harvest_finished_job``; ``setup.py`` gets the entry
  point in ``org.molssi.seamm.jobserver.job_finished`` (P5b).
- ``atct.py``: "already there" counts as success, and ``_save_covariance``
  is atomic with merge-on-write (P5d).
- ``installer.py``: the official store (Zenodo concept id, sha256 in seamm.ini),
  ATcT on install, and ``check`` covering all three (P5d).
- ``data/configuration.txt``: ``fragment-store-writable = auto`` (P5b).

thermochemistry_step:

- ``cbh_approach.py``:

  - ``_from_job``, tried before ``_geometry_from_store`` and
    ``_energies_from_store`` in ``run_cbh``, and in ``_energy_level_energies``
    for the energy level (P5a);
  - normalized model names in the store calls (P5a);
  - each species' expected fingerprints from
    ``Evaluator.expected_fingerprints``, used in ``_from_job``,
    ``_geometry_from_store``, ``_energies_from_store`` and
    ``_energy_level_energies`` instead of the target's, plus the self-check
    warning (P5a);
  - the "uncertainty method" column and property (P5d);
  - layers named in the output when an energy and its geometry come from
    different layers (P5c);
  - write-back after ``_compute_species`` and after the energy level's computed
    single points, the marker file, and ``_fragment_store`` opening writable
    when allowed (P5b);
  - ``LayeredStore`` in ``_fragment_store`` (P5c);
  - the run-time ATcT download where ``ATcT(...)`` is built (P5d).

- ``thermochemistry_parameters.py``: the ``fragment store`` choices (P5b), and
  the ``covariances`` default "fetch if needed" (P5d).
- ``tk_thermochemistry.py``: the new choices, nothing else (P5b).
- The documentation of the step's store parameters (P5b, P5c).

seamm_jobserver:

- ``jobserver.py``: run the ``job_finished`` entry points in a background
  thread after ``_finalize_job_status`` (P5b).

orca_step (P5a):

- ``batch.py``: ``expected_fingerprints(model_chemistry, symbols, *, kind,
  version, options=None)``, built with ``engine_helpers`` and
  ``exact_exchange_method`` (the MDI engine's rule) and hashed with
  ``fingerprint.orca_calculation_record``. It returns the as-written and the
  COSX-moot variants.

seamm_exec (P5a):

- ``evaluator.py``: ``Evaluator.expected_fingerprints(configuration, kind)``,
  which calls the code's batch-contract function, or returns None for codes
  without one.

Tests
-----

Unit tests need no network: Zenodo and ATcT are mocked.

P5a
    - ``split_model_chemistry`` on ``ORCA:CMP@R2SCAN-3C``, a bare name plus a
      program, ``ORCA:DFT@B3LYP/bse:def2-SVP``, ``QC@CCSD(T)-F12D_RI/...`` and a
      contradicting program.
    - Prefixed and bare names find the same rows. A harvest of a sub-step job
      and of an engine job gives the same model and program.
    - ``_from_job``: a ``CBH:C`` with matching fingerprints is reused (with
      ``_compute_species`` mocked); a mismatched fingerprint is computed; two
      ``CBH:C`` systems work; the energy level is reused; an imaginary mode
      above 20i is not reused.
    - ``expected_fingerprints`` (orca_step): for CH₄, BF₃, PF₆⁻, Na⁺ and a lone
      atom, at revDSD and at r2SCAN-3c, the predicted set contains the
      fingerprint recorded by a real ORCA run. The recorded fingerprints come
      from fixtures cut from 5602 and 5684 outputs, with no ORCA in the tests.
      Ordinary species give ``1e16…``, boron and phosphorus species give
      ``7ca1…``, and r2SCAN-3c gives the same fingerprint for all.
    - In a P target, the CH₄, H₂ and HF fragments are looked up with ``1e16…``
      and hit the stored rows; the P fragments are looked up with ``7ca1…``.
    - A code without ``expected_fingerprints`` falls back to the target's
      fingerprint.
    - The self-check warning appears when the prediction is made to miss.
P5b
    - Golden: write-back rows equal ``harvest_job`` rows for the same job.
    - An imaginary geometry is refused, with its energies.
    - Structure guard: a fragment whose optimized geometry has different
      connectivity (a fixture of a rearranged hypervalent hydride) is refused.
      So are a wrong multiplicity (BH as a triplet against the key's singlet)
      and an unusual centre with ligands too close (the old bad PF₆⁻ geometry).
      The good octahedral PF₆⁻ passes.
    - A lower conformer that is structurally distinct gets ``auto:``. The same
      minimum with E + ZPE lower by 0.12 kJ/mol and RMSD below 0.05 Å goes to
      the duplicate log. A higher conformer is logged.
    - ``is_local_filesystem`` against monkeypatched mount tables.
    - Eight processes writing 50 species each to one temporary store give every
      row and no error.
    - A store locked past the timeout is logged, and the species still succeeds.
    - The marker is written only when the step did not write.
    - The JobServer hook is discovered, runs off the poll thread, its
      exceptions are isolated, and a second harvest finds only duplicates.
P5c
    - R1-R5 with two temporary stores, including R2 (official supersede beats an
      identical old local row) and R3 (logged).
    - The lowest conformer across layers.
    - R5: a local energy on an official geometry (same hash) is used, and the
      output names both layers. After a newer official supersedes that
      geometry with new coordinates, the local energy is not found and the
      species is computed (Science's condition (ii)).
    - The official file is read-only (a write raises) and opened immutable.
    - ``export`` gives only rows that are not in official.
    - Merge round trip: export, then merge into an empty store, gives the
      source's rows.
    - ``migrate`` on a schema-2 fixture copied from ChemAI's layout: backup,
      normalization, idempotence.
P5d
    - Install with mocks gives three files and the recorded sha256s. ``check``
      flags a modified official file.
    - Run-time download when missing. A download race counts as success.
    - A failed covariance fetch gives quadrature, the warning, and "quadrature
      (no covariance: …)" in the uncertainty-method column; a cached block gives
      "covariance".
    - ``covariance.json`` merge with two writers.
P5e
    - The ``official build`` sanity pass refuses a store with an imaginary row,
      a BH triplet, a rearranged structure or a bad unusual centre, and builds
      once those rows are excluded with reasons.

**Validation runs** (cluster). Science approved them on 11-10: about 10-15
core-h, within the per-item limit.

V1
    A serial 20-compound slice of the 1694 set, default preset, ChemAI queue.

    - It includes at least one **P or B target** (e.g. a phosphate or a
      fluoroborane). Its ordinary fragments (CH₄, H₂, HF) must be hits on rows
      that the ordinary targets computed, and its P/B fragments must be looked
      up with the NoCOSX fingerprint.
    - Acceptance: each distinct fragment is computed once, and results equal
      5684's within 0.01 kJ/mol. The fingerprint self-check gives no warning,
      in particular on the ordinary-fragment rows of the P/B target (Science).
    - It also measures the core-hours saved.
    - About 5 core-h.
V2
    The same slice as a parallel Loop, 4 iterations at once. Store hits after
    the first wave. The report gives the duplicate count per fragment in the
    first wave (Science), which settles D2. About 5-7 core-h.
V3
    Two compounds through an ssh-transport queue to TinkerCliffs. The hook
    harvests into ChemAI's local layer when the job returns. About 1 core-h.

Total about 10-15 core-h.

Migration of the existing stores
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 35 30

   * - File (2026-11-10)
     - Contents
     - Becomes
   * - ChemAI ``/home/seamm/SEAMM/Parameters/thermochemistry/fragments.db`` (local
       ZFS)
     - Schema 2. fragment 100 rows (``ORCA:``-prefixed names, ``program`` NULL);
       geometry 648 (598 ``CMP@R2SCAN-3C`` from 5602-5604, plus 50
       ``ORCA:CMP@R2SCAN-3C`` from 5684); energy 1321; no history. Backup:
       ``fragments.db.pre-5684-harvest-20261110``.
     - ChemAI's local layer, in place, after ``migrate``. It is the source of
       official v1 once the DEFGRID3 recompute is in.
   * - Mac ``~/SEAMM/Parameters/thermochemistry/fragments.db`` (11-09, the
       598-row P1 store, which is also
       ``presets/revdsd_store/fragments_p1.db``)
     - Schema 2, the seed before 5684.
     - The Mac's local layer, after ``migrate``. Optionally ``merge`` ChemAI's
       copy into it.
   * - ChemAI ``/D3/psaxe/ATcT/{fragments, f12_pbe0geom, gfn2_campaign,
       pbe0_campaign, r2scan3c_atct, r2scan3c_campaign}.db``, with read-only
       copies on ARC ``/projects/seamm/psaxe/ATcT``
     - Campaign stores (ladders, the benchmark); some are schema 1.
     - They stay named stores (``fragment store file``, ``ladder --store``).
       ``migrate`` normalizes them when run on them. No automatic merge.
   * - ARC ``/projects/seamm/SEAMM/Parameters/thermochemistry`` (GPFS)
     - No ``fragments.db``.
     - Official layer only (D7).

``migrate`` is safe on ChemAI's store: every row there is verified, so it
changes the ``model_chemistry`` and ``program`` columns but no key. It will
still report any collision.

Phases and estimate
-------------------

.. list-table::
   :header-rows: 1
   :widths: 10 40 12 18 20

   * - Phase
     - Content
     - Dev time
     - Releases
     - Needs
   * - P5a
     - Model-chemistry names; in-job reuse from the job's database; per-species
       expected fingerprints and the self-check
     - 1.5-2 days
     - seamm-thermochemistry, thermochemistry-step, orca-step, seamm-exec
     - No D-item. Can start on Paul's go.
   * - P5b
     - Write-back, ``harvest_configuration``, the filesystem rule, the
       structure guard, distinct conformers, the JobServer hook, the marker
     - 2 days
     - seamm-thermochemistry, thermochemistry-step, seamm-jobserver
     - D1, D2, D5, D6
   * - P5c
     - Two layers, R1-R5, export, merge, official build, migrate
     - 1.5-2 days
     - seamm-thermochemistry, thermochemistry-step
     - D7, D8 (R5 answered by Science)
   * - P5d
     - Installer (official store, ATcT), run-time ATcT download, covariance
       default and cache, the uncertainty-method column
     - 1-1.5 days
     - seamm-thermochemistry, thermochemistry-step
     - D4
   * - P5e
     - Official v1 build with the sanity pass, Zenodo record, concept id
       release
     - 0.5-1 day, plus Paul's Zenodo upload
     - seamm-thermochemistry
     - D3; the DEFGRID3 recompute (already approved, about 90 core-h, not
       counted here)
   * - V1-V3
     - Validation runs
     - 0.5 day elapsed
     - none
     - Approved by Science, 10-15 core-h

**Total: 7.5-9 working days** (5.5-6.5 before Science's review; the additions
are the per-species fingerprints, the structure guard, the uncertainty column
and the sanity pass).

- About 10 PRs: seamm-thermochemistry 4 (P5a, b, c, d+e), thermochemistry-step
  3, orca-step 1, seamm-exec 1, seamm-jobserver 1.
- P5a's orca-step and seamm-exec PRs merge before its thermochemistry-step PR.
- Order: P5a → P5b → P5c → P5e. P5d can go alongside P5b or P5c.
- Installs: ChemAI and the Mac after P5a+P5b (V1-V3 need them), and ChemAI,
  ARC and the Mac after P5c+P5d, as Paul batches updates.
- ChemAI is updated only while no thermochemistry job is running there.

Not in this note
----------------

- A dedupe pre-pass across a Loop's targets (requirement 3 of 10-28) is a
  separate design.
- The step raises for all-trivial targets that are in ATcT (Science, 11-10
  17:50). This is a separate bug.
- S and G from the stored frequencies.
- Fingerprints for codes other than ORCA over MDI. Write-back works for any
  code that stores them, and is unverified otherwise.
