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

## Decision

Waiting for Henrik's approval, as a change to the review tooling.
