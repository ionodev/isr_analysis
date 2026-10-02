# Review record: `changelog`

- **Gates:** A only (documentation; no pipeline code). The REVIEW_PROCESS.md
  change is a change to the review tooling, so it needs Henrik's approval
  (REVIEW_PROCESS.md, section 1).
- **Author:** Claude (Opus 5.5) session for Henrik, 2026-10-02, at Henrik's
  request ("we should start writing changelogs").
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent),
  2026-10-02, two rounds.
- **Commits:** `8c986cf..1d3c408`, based on `main` at `8ce21f2`.

## What the change does

Adds `CHANGELOG.md`: every change to `main` since the fork's work began on
2026-09-16, by date under Fixed / Added / Changed / Process, with the memo,
the commit and the effect on the products. It also lists the seven fix
branches waiting for gate B, with their approval status and merge order.
REVIEW_PROCESS.md: gate A step 4 adds the entry on the branch;
documentation-only merges add one too; gate B step 7 moves the entry to the
merge date in a commit on the branch; step 5 needs no new benchmark when
`main` moved only in Markdown or LaTeX files. The record and pull-request
templates ask for the entry.

## Tests

Not applicable: no code touched.

## Benchmark

Not applicable: no pipeline code touched (`git diff --name-only main...changelog`:
CHANGELOG.md, REVIEW_PROCESS.md, .github/pull_request_template.md,
review/records/).

## Scientific review (gate B, and memos with results)

Not applicable.

## Full verification (when flagged, and always for calibration)

Not needed.

## Code review

First round (8c986cf):
- **Every hash, date and memo number checked against the repository:**
  correct. No first-parent commit since 2026-09-16 is missing.
- **Fixed in a9f4886:**
  - the merge dependencies were understated;
  - the raw_reader figure was wrong (Memo 19: 25-30 times in one process,
    up to about 40 times with 4-16 readers), and its `fast_read` default
    was not mentioned;
  - entries that change numbers lacked the effect on the products;
  - the site's T_sys log was overstated, now called preliminary;
  - the process text: what to cite before a merge, documentation-only
    merges, the gate B move as a commit, a conflict note, a Changelog
    section in the record template;
  - minor: `refine_reader`, Memos 12 and 14-16, Memo 7 for 70429cd, an
    entry for this change.

Second round (a9f4886): fixes confirmed. Four minor points, fixed in
1d3c408:
- the median's scale;
- "the last commit of the change itself" instead of "the branch's tip";
- the step 5 exemption: Markdown or LaTeX only, checked against the diff
  from the tested commit to the branch tip;
- two more effect figures (0330e30, 95aba03). The third suggested figure,
  for 61c79aa, was not found in its commit message and is left out.

## Changelog

The 2026-10-02 Process entry in `CHANGELOG.md`.

## Decision

- Gate A: passed on 2026-10-02 (code review answered; no tests or benchmark
  apply).
- Waiting for Henrik's approval, as a change to the review tooling.
