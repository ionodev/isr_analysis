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
| decided by | the benchmark: bit-identical outputs | the benchmark: any output differs |
| approval | none, merged once the gate is passed | Henrik; Juha too for calibration or physics |

If in doubt, the benchmark decides. A change that is meant to be neutral but
does not give bit-identical outputs is either a bug or a gate-B change.

Changes to documentation only (Markdown, LaTeX, comments) need no gate.

## 2. Gate A: output-neutral changes

1. **Tests.** `python3 -m pytest -q` passes. New code comes with unit tests
   (`test_*.py`) where it can be tested without the raw data.
2. **Benchmark**, if the change touches pipeline code (anything imported by
   `run_analysis.py`, `outlier_lpi.py`, `fit_lpi.py`,
   `avg_range_doppler_spec.py`, `fit_lp.py` or the modules they use):
   `python3 review/regression.py <branch>` must exit with status 0 ("Every
   output is bit-identical").
3. **Independent code review.** A reviewer who did not write the change reads
   the diff (`git diff main...<branch>`) against the checklist in section 5.
   For Claude this means a fresh session or a subagent with no part in the
   work, for example `/code-review` at high effort. Findings are fixed, or
   answered in the record, and the reviewer checks the fixes.
4. **Review record.** Write `review/records/<YYYY-MM-DD>-<branch>.md` from
   `review/records/TEMPLATE.md` and commit it on the branch.
5. **Merge.** Merge into `main` with a merge commit that names the record,
   then push `main` and the branch.

## 3. Gate B: product-changing changes

1. **Tests**, as in gate A.
2. **Benchmark, quantified.** Run `python3 review/regression.py <branch>`
   and keep the report. It shows which outputs change and by how much: in
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
7. **Approval.** Henrik approves the merge, and Juha too for calibration or
   physics. Until then the branch is pushed but not merged, and it is listed
   in TODO.tex (Q8).
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
- a 300 s `fit_lpi` stretch and a mode-300 range–Doppler period per channel.

`review/regression.py` runs the pipeline's default stages on these periods
with the code of `main` and of the branch. Each commit gets a detached
worktree under `~/isr_project/regression/worktrees/`. Each job runs in its
own memory-limited scope. A commit's outputs are cached under
`~/isr_project/regression/runs/<sha>/`, so `main` is computed only once per
commit. Reports go to `~/isr_project/regression/reports/`.

```bash
python3 review/regression.py <branch> [--base main] [--jobs 16]
```

Exit status: 0 bit-identical, 1 some output differs, 2 a job failed or an
output is missing (see the logs under the run directory). One commit takes
about half an hour with 16 jobs.

Add a period to the benchmark when a new problem is found. After changing
the benchmark, delete the cached runs (`~/isr_project/regression/runs/`).

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

The review record in the branch is the record of the review, so a pull
request is optional. When one is used, it is opened within the fork
(`ionodev/isr_analysis`, branch into `main`), never against
`jvierine/isr_analysis`. The description follows
`.github/pull_request_template.md`.

## 7. Autonomous runs

Autonomous Claude runs (see `~/isr_project/AUTONOMOUS_RUNS.md`) follow the
same gates. They may merge gate-A changes themselves, with a subagent as the
independent reviewer. They push gate-B branches without merging them, with
the benchmark report and the verification ready in the record.
