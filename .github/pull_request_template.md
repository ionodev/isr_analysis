<!-- Within ionodev/isr_analysis only (branch into main); never against jvierine/isr_analysis. See REVIEW_PROCESS.md. -->

**Gates:** A only (output-neutral) / A and B (product-changing)

**What and why:**

**Tests:** `python3 -m pytest -q`:

**Benchmark:** `python3 review/regression.py <branch>`: exit status, and for
gate B the changes and their explanation.

**Verification and scientific review (gate B):** what was recomputed, the verdict, the findings, and what was done.

**Code review:** reviewer, findings, and what was done.

**Review record:** `review/records/<YYYY-MM-DD>-<branch>.md`

**Memo (gate B):**
