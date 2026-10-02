# Review record: `antenna-switch-gap`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/antenna-switch-gap.md`
- **Commits:** `54eb8ea`, `8c15038` (with a merge of main 8206c84: `8fc50c1`); based on `main` at `8206c84`. Main is now `8ce21f2`, which differs only in `.claude/agents/`; merge main before the pull request.

## What the change does

The antenna metadata change antenna at the event that opens the next record cycle, but the radar changes at the end of the previous one, a median of 8.6 s earlier (2.3-20.8 s; Memo 30), so a channel's inversion took in pulses sent on the other antenna or not at all. `get_antenna_select` now returns 0 (unknown) from `switch_guard_s` = 2.5 s before each closing event until the opening event, and every caller rejects 0. The guard was chosen on this recording.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 29 passed. New: `test_millstone_radar_state.py` (a normal cycle, a recording that starts inside a cycle, the last event, the clamp on a short cycle, a change without a closing event, and that the values stay steps in {-1, 0, 1}); the start-of-recording test fails on `54eb8ea`.

## Benchmark

`8c15038`: exit 1, 72 of 79 files bit-identical; report `~/isr_project/regression/reports/8c15038-8c15038_vs_8206c84-8206c84.md`. Before the fix (`8fc50c1`) it was 70 of 79: misa-l's first period and its fit were lost, because the first metadata event (+52.3 s) closes a cycle and its 0 was held back to the recording start. The remaining 7 files are the intended changes at antenna switches: one zenith-l 10-s period (`lpi-1712594080`) is no longer written, as none of its pulses passes the antenna test any more (Memo 30's periods made only of wrong pulses); misa-l `lpi-1712592900` loses 1183 values (pulses at the end of its cycle removed; P_tx changes by 12.6 %), and its fit `pp-1712592750` loses all 108 gates on this branch alone, which `lpi-unconstrained-gates` recovers (see the integration below).

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Not done yet. Gate B waits for the supervisor's decision on the fix branches (TODO.tex, Q8), and is an L task (Henrik's go-ahead first).

## Full verification (when flagged, and always for calibration)

Not decided yet: the scientific review says whether it is needed.

## Code review

- **Real bug: the start of the recording counted as unknown** (60.6 s of MISA pulses dropped): fixed in `8c15038` (the step functions hold the first recorded value before the first event).
- **No test:** fixed in `8c15038`.
- **A missing `cycle_name` field** would have silently guarded on an arbitrary basis or crashed: `8c15038` falls back to the metadata as recorded, with a message.
- **Switches not preceded by a closing event** were left unguarded silently: `8c15038` reports how many.
- **The docstring's "8.2-20.8 s"** was not Memo 30's: corrected to 2.3-20.8 s (median 8.6 s).
- **A trailing closing event** would hold its antenna to the end: not present in this recording; left.
- **`tsys_harvest.py`'s comment** now documents 0.
- **The satellite catalogue reads the raw metadata itself** and keeps the 8-s lag: outside this branch, noted.
- **Interplay with `tx-delay-own-antenna`:** the delay is measured at the recording start; with the start fixed, both measure MISA's first cycle again.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: waiting for the supervisor's decision on the fix branches (Q8), then the scientific review, then Henrik's approval (and Juha's, for calibration or physics).
