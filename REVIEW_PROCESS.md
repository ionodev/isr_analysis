# Review process for changes to `main`

Nothing is merged into `main` until it has passed the review below. The aim
is that the code on `main` is free of known bugs and gives correct results.
The process was agreed on 2026-09-30. It applies to everyone who changes the
code, including Claude sessions, supervised or autonomous.

## 1. Which gate?

Every change goes through one of two gates, depending on what it does to the
products: the ACFs, system temperatures, densities, temperatures and
velocities the pipeline writes.

| | Gate A: output-neutral | Gate B: product-changing |
|---|---|---|
| examples | new scripts, tests, refactoring, speed-ups, new options that are off by default | bug fixes that change numbers, new defaults, calibration changes |
| decided by | bit-identical benchmark outputs, and the reviewer's judgement of what the benchmark does not cover | any output differs, or a change the benchmark cannot see |
| approval | none, merged once the gate is passed | Henrik; Juha too for calibration or physics |

A change that is meant to be neutral but does not give bit-identical outputs
is either a bug or a gate-B change. The reverse does not hold: bit-identical
outputs on the benchmark are necessary for gate A, but not sufficient. The
benchmark does not cover:
- `run_analysis.py` itself (how the configuration maps to arguments);
- how the theory tables are generated (`isr_spec.py`, `il_interp.py` table
  building: the benchmark reads the existing tables);
- mode 800, options that are off by default, and data outside the benchmark.

A change to any of these must be shown neutral in another way that the
reviewer accepts, for example a targeted test. Otherwise it is gate B.
Changes to how the tables are generated are always gate B.

Changes to the review tooling itself (`review/regression.py`,
`review/benchmark.json`, this file) need Henrik's approval, like gate B.

Changes to documentation only (Markdown, LaTeX) need no gate. A `.py` change
that is meant to touch comments only still goes through gate A, because
someone has to confirm that it does.

## 2. Gate A: output-neutral changes

1. **Tests.** `python3 -m pytest -q` passes. New code comes with unit tests
   (`test_*.py`) where it can be tested without the raw data.
2. **Benchmark**, whenever `git diff --name-only main...<branch>` lists a
   `.py` file other than `test_*.py` and the files under `review/`. First
   merge or rebase `main` into the branch, so that only its own changes are
   compared. Then run **main's copy** of the tool, so that a branch cannot
   judge itself with a changed tool:

   ```bash
   python3 ~/isr_project/isr_analysis/review/regression.py <branch>
   ```

   The tool refuses to run if its `review/regression.py` or
   `review/benchmark.json` is not identical to `main`'s copy. The flag
   `--allow-other-tool` overrides this. It is only for testing a change to
   the tool itself, and such a change needs Henrik's approval. When `main`
   has no tool yet (the branch that introduced it), the tool warns and runs
   its own copy. The path above is `main`'s copy only while that
   checkout is on `main`.

   It must exit with status 0 ("Every output is bit-identical"). The record
   gives the tested commit. If commits are added after the run, it is
   repeated, unless those commits change documentation only.
3. **Independent code review.** A reviewer who did not write the change reads
   the diff (`git diff main...<branch>`) against the checklist in section 5.
   For Claude this means a fresh session or a subagent with no part in the
   work, for example `/code-review` at high effort. Findings are fixed, or
   answered in the record, and the reviewer checks the fixes.
4. **Review record.** Write `review/records/<YYYY-MM-DD>-<branch>.md` from
   `review/records/TEMPLATE.md` (with any `/` in the branch name replaced by
   `-`) and commit it on the branch.
5. **Merge.** Right before merging, check that `main` has not moved since the
   benchmark: `git merge-base --is-ancestor main <tested commit>`. If it has,
   merge `main` into the branch and repeat the benchmark. Then merge into
   `main` with a merge commit that names the record, and push `main` and the
   branch.

## 3. Gate B: product-changing changes

1. **Tests**, as in gate A.
2. **Benchmark, quantified.** Run main's copy of the tool on the branch, as
   in gate A, and keep the report. It shows which outputs change and by how much: in
   standard deviations for the ACFs, relative for the rest, and lost or
   recovered values. Every change in it must be explained by the fix. A
   change nobody can explain is a bug until shown otherwise.
3. **Independent verification of the scientific claim.** The claim that the
   fix is right (for example "the metadata are 8.6 s late", or "the
   injection window catches the spike") is rechecked with separately written
   scripts, not the author's. For Claude, a fresh session or subagent does
   this. Memo 30's last section is an example: it found two overstated
   numbers.
4. **Independent code review**, as in gate A.
5. **Memo.** A memo says what changed, why, by how much, and what was not
   checked.
6. **Review record**, as in gate A, with the benchmark report and the
   verification attached.
7. **Pull request and approval.** Open a pull request within the fork
   (section 6), with the record's content. Henrik approves the merge, and
   Juha too for calibration or physics, in the pull request or in person.
   Until then the branch is pushed but not merged, and it is listed in
   TODO.tex (Q8).
8. **Merge**, as in gate A. Products made before the merge are marked
   out of date in TODO.tex.

## 4. The benchmark

`review/benchmark.json` fixes the data. Each period was chosen because a
memo found something at it:
- zenith-l: a quiet period, a cycle change with wrong-antenna pulses, a
  period made only of them, power-line impulses, a long bright satellite
  pass, and a quiet night;
- misa-l: lost lags at the clutter gates, a cycle change, power-line
  impulses, and the start of the recording;
- a 300 s `fit_lpi` stretch, and 13 consecutive mode-300 range–Doppler
  periods per channel.

`review/regression.py` runs these stages with the code of `main` and of the
branch:
- `outlier_lpi.lpi_files`, with the arguments `run_analysis.py` passes for
  `config/millstone_eclipse2024.json`;
- `fit_lpi.fit_lpifiles`;
- `avg_range_doppler_spec` for mode 300, on 13 consecutive periods per
  channel;
- `fit_lp.fit_spectra` on its output. This uses a 60 s averaging window
  instead of `run_analysis.py`'s 600 s, which would need 61 range–Doppler
  periods.

It then compares every dataset and attribute of every output file, byte for
byte. Each commit gets a detached worktree under
`~/isr_project/regression/worktrees/`, removed after a successful run. Each
job runs in its own memory-limited scope. A commit's outputs are cached under
`~/isr_project/regression/runs/<sha>-<key>/`. The key is a hash of the tool,
the benchmark, the installed Python packages and the table files, so a
change to any of them starts a fresh run. Two runs of the same commit wait
for each other. Reports go to `~/isr_project/regression/reports/`.

```bash
python3 ~/isr_project/isr_analysis/review/regression.py <branch> [--base main] [--jobs 16]
```

Exit status:
- 0: every output is bit-identical;
- 1: some output differs;
- 2: the run is not valid. A job failed, an expected output is missing on
  either side, nothing was compared, the branch does not contain `main`, or
  the tool failed. The logs are under the run directory.

One commit takes about 19 minutes with 16 jobs (measured on 2026-09-30:
1129 s for `main`, 1113 s for a branch), so a review with a fresh `main`
takes about 40 minutes, and about 20 when `main`'s outputs are cached.

Add a period to the benchmark when a new problem is found. This changes the
key, so every commit is run afresh.

## 5. Review checklist

The reviewer checks at least:

- **Correctness:** indexing and off-by-one errors (slices, range gates, lag
  indices), units (µs or samples, Hz or kHz, K), signs (Doppler, delays),
  NaN and empty-array handling, integer division, dtype precision
  (complex64 roundoff, Memo 17).
- **Pulse selection and timing:** antenna and power tests, mode tables
  (`radar_timing.py`), window boundaries such as `last_echo`, `noise0` and the
  diode switch-on (Memos 27, 28, 30).
- **Defaults:** does a default change? If so, it is gate B.
- **Tests:** do they test the change, and would they fail without it?
- **Claims:** does every number in the commit message or memo come from a
  run? Is anything stated that was not checked?
- **Parallel runs:** memory isolation (`systemd-run --scope -p MemoryMax
  -p MemorySwapMax=0`), no stray processes.
- **Commit identity:** `ionodev <44322493+ionodev@users.noreply.github.com>`,
  no private address, no attribution trailers.

## 6. Pull requests

The review record in the branch is the durable record of the review. Gate-B
changes also get a pull request, so that Henrik and Juha can read and approve
them on GitHub. Gate-A changes may have one. Pull requests are opened within
the fork, branch into `main`, never against `jvierine/isr_analysis`. The
repository is a fork, and `gh` would otherwise offer the parent as the
target, so always pass `--repo`:

```bash
gh pr create --repo ionodev/isr_analysis --base main --head <branch> \
  --title "<title>" --body-file review/records/<record>.md
```

The description follows `.github/pull_request_template.md`. After approval,
merge with a merge commit, either locally (then push) or with
`gh pr merge <number> --merge --repo ionodev/isr_analysis`.

## 7. Autonomous runs

Autonomous Claude runs (see `~/isr_project/AUTONOMOUS_RUNS.md`) follow the
same gates. They may merge gate-A changes themselves, with a subagent as the
independent reviewer. They push gate-B branches without merging them, with
the benchmark report and the verification ready in the record, and they may
open the pull request.
