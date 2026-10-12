<!-- Within ionodev/isr_analysis only (branch into main); never against jvierine/isr_analysis. See REVIEW_PROCESS.md. -->

**Gates:** A only (output-neutral) / A and B (product-changing)

**What and why:**

**Tests:** `python3 -m pytest -q`:

**Benchmark:** `python3 review/regression.py <branch>`: exit status, and for
gate B the changes and their explanation.

**For Henrik's review (gate B):** what was wrong, what changes in the products and by how much, the evidence, and what was not checked. Short.

**Scientific review (only if Henrik asked for one):** verdict, findings, and what was done.

**Code review:** reviewer, findings, and what was done.

**Review record:** `review/records/<YYYY-MM-DD>-<branch>.md`

**Changelog:** the entry added to `CHANGELOG.md`

**Memo (gate B):**
