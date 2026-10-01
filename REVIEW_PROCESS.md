# Review process for changes to `main`

Nothing is merged into `main` until it has passed the review below. The aim
is that the code on `main` is free of known bugs and gives correct results.
The process was agreed on 2026-09-30. It applies to everyone who changes the
code, including Claude sessions, supervised or autonomous.

## 1. The two gates

Every code change goes through **gate A**: tests, the benchmark, an
independent code review, and a review record. The benchmark then decides
whether the change also needs **gate B**:

- **The outputs are bit-identical** (the change is *output-neutral*): gate A
  is enough, and the change is merged without further approval. Examples:
  new scripts, tests, refactoring, speed-ups, new options that are off by
  default.
- **Some output differs** (the change is *product-changing*): gate B follows.
  Gate B is a scientific review of the new results, a full independent
  verification when that review calls for it, and Henrik's approval, and
  Juha's too for calibration or physics.
  Examples: bug fixes that change numbers, new defaults, calibration changes.

Gate B does not repeat the code review. It checks what gate A cannot: whether
the changed results are right.

"Products" means the ACFs, system temperatures, densities, temperatures and
velocities the pipeline writes.

A change that is meant to be neutral but does not give bit-identical outputs
is either a bug or a product-changing change. The reverse does not hold:
bit-identical outputs on the benchmark are necessary, but not sufficient. The
benchmark does not cover:
- `run_analysis.py` itself (how the configuration maps to arguments);
- how the theory tables are generated (`isr_spec.py`, `il_interp.py` table
  building: the benchmark reads the existing tables);
- mode 800, options that are off by default, and data outside the benchmark.

A change to any of these must be shown neutral in another way that the code
reviewer accepts, for example a targeted test. Otherwise it goes through
gate B. Changes to how the tables are generated always do.

Changes to the review tooling itself (`review/regression.py`,
`review/benchmark.json`, `.claude/agents/`, this file) need Henrik's
approval.

Changes to documentation only (Markdown, LaTeX) need no gate. A `.py` change
that is meant to touch comments only still goes through gate A, because
someone has to confirm that it does.

## 2. Gate A: every code change

1. **Tests.** `python3 -m pytest -q` passes, with nothing skipped. A fresh
   worktree has no interpolation tables, so link them in first:
   `ln -s ~/isr_project/isr_analysis/data/ion_line_interpolate_*.h5 data/`
   (remove the links afterwards). New code comes with unit tests
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
   its own copy. The path above is `main`'s copy only while that checkout is
   on `main`.

   Exit status 0 ("Every output is bit-identical") means the change is
   output-neutral. Status 1 means it changes products, so gate B follows;
   keep the report. Status 2 means the run is not valid; fix it and repeat.
   The record gives the tested commit. If commits are added after the run,
   it is repeated, unless those commits change documentation only.
3. **Independent code review.** A reviewer who did not write the change reads
   the diff (`git diff main...<branch>`) against the checklist in section 6.
   For Claude this means the `isr-code-reviewer` agent on Sonnet 5.5
   (section 6a), which has no part in the work. Findings are fixed, or
   answered in the record, and the reviewer checks the fixes. This is the
   only code review, for both kinds of change.
4. **Review record.** Write `review/records/<YYYY-MM-DD>-<branch>.md` from
   `review/records/TEMPLATE.md` (with any `/` in the branch name replaced by
   `-`) and commit it on the branch.
5. **Merge, if output-neutral.** Right before merging, check that `main` has
   not moved since the benchmark: `git merge-base --is-ancestor main <tested
   commit>`. If it has, merge `main` into the branch and repeat the
   benchmark. Then merge into `main` with a merge commit that names the
   record, and push `main` and the branch. A product-changing change goes on
   to gate B instead.

## 3. Gate B: in addition, for changes that alter the products

1. **The benchmark report, explained.** It shows which outputs change and by
   how much: in standard deviations for the ACFs, relative for the rest, and
   lost or recovered values. Every change in it must be explained by the
   fix. A change nobody can explain is a bug until shown otherwise.
2. **Memo.** A memo says what changed, why, by how much, and what was not
   checked.
3. **Scientific review**, by the independent `isr-science-reviewer` agent on
   Opus 5.5 (sections 5, 6a). Are the changed results physically sound, and
   does the conclusion follow from the evidence? It includes **spot checks**:
   the reviewer recomputes the one or two numbers the conclusion rests on,
   from the result files that already exist, without new runs on the raw
   data. It lists every key claim it could not confirm, and says whether a
   full verification is needed.
4. **Full independent verification, only when needed.** The key claims are
   recomputed with separately written scripts, not the author's, including
   new runs on the data where necessary. Memo 30's last section is an
   example: it found two overstated numbers. This is done:
   - when the scientific review flags a result as questionable, or leaves a
     key claim unconfirmed;
   - **always** for a change to the calibration or the absolute density scale
     (noise injection, system temperature, the conversion to density, the
     magic constant). An error there silently shifts every product.

   The same reviewer continues with it, so it keeps its context.

   Findings of steps 3 and 4 are fixed in the memo or the code, or answered
   in the record.
5. **Review record**, completed with the benchmark report, the scientific
   review, and the verification if one was done.
6. **Pull request and approval.** Open a pull request within the fork
   (section 7), with the record's content. Henrik approves the merge, and
   Juha too for calibration or physics, in the pull request or in person.
   Until then the branch is pushed but not merged, and it is listed in
   TODO.tex (Q8).
7. **Merge**, as in gate A. Products made before the merge are marked out of
   date in TODO.tex.

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

One commit takes 19 to 33 minutes with 16 jobs, depending on what else is
running (measured on 2026-09-30: 1113 to 1986 s). A review with a fresh `main`
therefore takes 40 to 60 minutes, and about half that when `main`'s outputs
are cached.

Add a period to the benchmark when a new problem is found. This changes the
key, so every commit is run afresh.

## 5. Scientific review

Correct code can still give wrong science, and an author tends to read their
own results kindly. So results get a scientific review by an independent
reviewer who did not do the work. This is required for gate B, and for every
memo that reports results (not only those tied to a code change), before the
memo goes to Juha.

For Claude, the reviewer is the agent `isr-science-reviewer`
(`.claude/agents/isr-science-reviewer.md`, on Opus 5.5). It is briefed as an experienced
space physicist and ISR expert, it is read-only, and it may read the
literature on the web. It checks:
- physical plausibility, and the eclipse response against earlier eclipse
  studies;
- the assumptions of ISR theory, and fits on table edges;
- calibration;
- whether the uncertainties match the scatter;
- whether the claims follow from the evidence;
- alternative explanations;
- comparison with independent sources (IRI, ionosondes, GNSS TEC, published
  results).

It reports a verdict, findings with evidence, the checks it made, and
questions for the supervisor. To run it, ask Claude to "use the
isr-science-reviewer agent on Memo N" (or on a branch's record). The
agents in `.claude/agents/` are found by sessions started in the repository
or below it, and by every session when they are linked into
`~/.claude/agents/`. A new `agents` directory is noticed only by sessions
started after it was created.

The review includes spot checks of the key numbers. A full independent
verification follows when the review calls for it, and always for results on
the calibration or the absolute density scale (section 3, step 4).

This review is an AI's second look. It can share blind spots with the agent
that did the work, and it does not replace Juha's judgement. Its purpose is
to catch what can be caught before results reach him.

## 6. Code review checklist

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

## 6a. Which model does which review

The reviews are done by AI agents, and they cost usage. The best model is
kept for the steps where a missed error would reach the results.

| Step | Who | Model |
|---|---|---|
| tests, benchmark, bit-identity, commit identity | scripts | none |
| code review (gate A, every change) | `isr-code-reviewer` | Claude Sonnet 5.5 (`claude-sonnet-5-5`) |
| scientific review with spot checks (gate B, and memos with results) | `isr-science-reviewer` | Claude Opus 5.5 (`claude-opus-5-5`) |
| full verification (when flagged, and always for calibration) | `isr-science-reviewer`, continued | Claude Opus 5.5 |

The agents' models are fixed by full ID in `.claude/agents/`, so they do not
change when the short names `sonnet` and `opus` move to newer versions.
Change them there, deliberately, when a newer model is adopted.

To keep the usage down:
- continue the same reviewer for follow-up rounds rather than starting a
  fresh one, since it keeps its context;
- give it the diff, the record and the memo it needs, not the whole project.

Sonnet 5.5 costs half as much per token as Opus 5.5 (2 and 10 dollars per
million input and output tokens, against 4 and 20; Anthropic's list prices,
September 2026).

## 7. Pull requests

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

## 8. Autonomous runs

Autonomous Claude runs (see `~/isr_project/AUTONOMOUS_RUNS.md`) follow the
same gates. They may merge output-neutral changes themselves once gate A is
passed, with the `isr-code-reviewer` agent as the independent reviewer. They push gate-B branches without merging them, with
the benchmark report and the verification ready in the record, and they may
open the pull request.

They also run the scientific review (section 5) on every memo they write
that reports results, with the `isr-science-reviewer` agent. They fix what
it finds, or list it in the report as open.
