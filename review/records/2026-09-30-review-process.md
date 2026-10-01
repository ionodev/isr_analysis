# Review record: `review-process`

- **Gate:** A (output-neutral: no pipeline code changes). It changes the
  review tooling, so it needs Henrik's approval under REVIEW_PROCESS.md,
  section 1.
- **Author:** Claude (supervised session with Henrik, 2026-09-30).
- **Reviewer:** an independent Claude subagent with no part in writing the
  change. Four rounds, on commits b8b9188 and d587f1e, then 3f7e96b,
  9ca6e26 and 612031e.
- **Commits:** `b8b9188..612031e`, based on `main` at `241190d`.

## What the change does

It introduces the review process for changes to `main`: `REVIEW_PROCESS.md`
(gate A for output-neutral changes, gate B for product-changing ones); the
regression benchmark `review/benchmark.json` and `review/regression.py`; this
record template; and a pull request template. The benchmark runs the
pipeline's LPI, fit_lpi, mode-300 range–Doppler and fit_lp stages on fixed
periods, each chosen because an earlier memo found something at it. It does
so with the code of `main` and of the branch, and compares every output byte
for byte.

## Tests

- `python3 -m pytest -q`: 23 passed, with main's interpolation tables linked into the worktree's
  `data/`. Without them, 7 tests skip: a fresh worktree has no tables.
- The comparison logic was tested on synthetic HDF5 files for:
  - identical files;
  - value, dtype, NaN-payload and signed-zero changes;
  - lost values;
  - nested groups;
  - file, group and dataset attributes;
  - variable-length strings against fixed bytes;
  - variable-length numeric data;
  - compound types with variable-length fields;
  - a mismatched variance shape;
  - files present on one side only;
  - empty directories.
- The worker's check that every loaded pipeline module comes from the tested
  commit's worktree was tested with an empty job.

## Benchmark

This branch introduces the tool, so it could not use `main`'s copy. The
tool ran its own copy (bootstrap) and said so.

- On 9ca6e26: exit 0. 79 of 79 output files bit-identical (base 1129 s,
  branch 1113 s). This branch changes no pipeline code, so it is also a test
  that two independent runs of the same pipeline code give byte-identical
  results. They do.
- On 612031e, the tested commit: exit 0. 79 of 79 output files
  bit-identical, environment key 287c81c2822b (base 1216 s, branch
  1986 s, while another job shared the machine). Report:
  `~/isr_project/regression/reports/review-process-612031e_vs_main-241190d.md`.
  Later commits on the branch change documentation only.

## Code review

Round 1 found nine problems. The most serious:
- the 300 s fit stretch had no data;
- missing outputs counted as bit-identical;
- the cache was keyed by commit only;
- the comparison was not byte-exact;
- fit_lp was not run;
- there was no ancestry check, crashes gave exit 1, and nothing locked a run.

Round 2 found two more that made every run invalid:
- the module-origin check failed on the tool itself;
- fit_lp can never produce output from one file.

It also found smaller ones: an output missing only on the branch counted as
invalid, and the key lacked the system libraries.

Round 3 found:
- benchmark.json was not checked against the base's copy;
- compound datasets with variable-length fields were compared by pointer;
- stale directories were never cleaned up;
- the documentation had leftovers.

All were fixed in 3f7e96b, 9ca6e26 and 612031e. Round 4 verified the
round-3 fixes. It found nothing that blocks gate A.

Accepted and not fixed:
- In `clean_stale()`, DONE is checked before the lock is taken. A run that
  finishes in that gap could have its completed directory removed. The only
  cost is a recompute. It is to be fixed in a follow-up change: check DONE
  again after taking the lock.
- HDF5 object and region references compare as different even when equal.
  This errs on the safe side, and the outputs contain none.
- UTF-8 and ASCII strings with the same bytes compare as identical.
  Cosmetic.

The review also found two existing pipeline bugs, now in TODO.tex:
- fit_lpi never fits the last group of LPI files (branch
  `fit-lpi-last-group`);
- fit_lp never uses the first range–Doppler file of a run.

## Added on 2026-10-01: the scientific review

At Henrik's request, the process now also requires a scientific review of
results, for gate B and for every memo that reports results. It is done by
an independent agent briefed as an experienced space physicist and ISR
expert: `.claude/agents/isr-science-reviewer.md`, read-only (section 5 of
REVIEW_PROCESS.md). These commits change documentation and the agent
definition only, so neither the tool, the benchmark nor the tested commit
changes. The tiering of section 6a (Sonnet 5.5 for gate A code review, Opus
5.5 for gate B and the scientific review) was added the same day, also at
Henrik's request.

Trial run on Memo 37 (the zenith beam's pointing), on Opus 5.5, with the
agent's brief:
- verdict "sound with changes";
- two errors found:
  - the round Gaussian beam model is contradicted by a one-sided northern
    sidelobe in the data, which biases the centre north and explains the
    memo's unexplained drift with the SNR cut (refit without the sidelobe
    detections: 40.53 degrees, model-free centroid 40.60-40.62, against the
    memo's 40.57-40.80);
  - the memo's picture of detections through MISA's sidelobes is not
    supported by the echo strengths;
- consequence: the aperture efficiency is 0.22-0.27, not 0.20-0.24;
- it reproduced the memo's fits before judging them;
- it gave five questions for the supervisor.

The memo's correction is separate work.

Also on 2026-10-01, at Henrik's request, the gates were restructured. Every
code change goes through gate A, which has the only code review (Sonnet 5.5).
A change that alters the products then also goes through gate B: one
independent verification and scientific review by `isr-science-reviewer`
(Opus 5.5), plus approval. Gate B no longer has a second code review. The
trial above did the verification and the scientific review together, which
is the form gate B now takes.

Later the same day, again to save usage and at Henrik's request, gate B's
review was split. First comes a scientific review with spot checks: it
recomputes the one or two key numbers from existing result files, with no
new data runs, and lists any claims it could not confirm. A full independent
verification follows only when that review flags something or leaves a key
claim unconfirmed, and always for changes to the calibration or the
absolute density scale.

## Decision

Approved by Henrik on 2026-10-01, and merged by him on GitHub as pull request #1
(merge commit 51b9a0c, containing the reviewed head 6e998aa unchanged).
