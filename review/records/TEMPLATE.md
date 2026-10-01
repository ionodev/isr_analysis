# Review record: `<branch>`

Copy to `review/records/<YYYY-MM-DD>-<branch>.md`, fill it in, and commit it on
the branch before merging. See REVIEW_PROCESS.md.

- **Gates:** A only (output-neutral) / A and B (product-changing)
- **Author:** who wrote the change (person, or Claude session and date)
- **Reviewer:** who reviewed it (must not be the author; for Claude, a fresh
  session or subagent)
- **Commits:** `<first>..<last>`, based on `main` at `<sha>`

## What the change does

One paragraph: what, and why.

## Tests

`python3 -m pytest -q`: N passed. New tests: which, and what they check.

## Benchmark

`python3 review/regression.py <branch>`: exit status N. Report:
`~/isr_project/regression/reports/<name>.md`.

For gate A: "Every output is bit-identical", or "not applicable: no pipeline
code touched".

For gate B: which outputs change and by how much, and why each change is
expected.

## Verification and scientific review (gate B, and memos with results)

Reviewer (`isr-science-reviewer`): what it recomputed and whether it agreed,
its verdict, its findings and what was done about each, and its questions
for the supervisor.

## Code review

Findings, and what was done about each: fixed in `<sha>`, or answered.

## Decision

- Gate A: passed on `<date>`, merged in `<merge sha>`.
- Gate B: approved by `<name>` on `<date>`, merged in `<merge sha>`; or
  waiting for approval.
