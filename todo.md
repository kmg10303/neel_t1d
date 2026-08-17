To do for 8/17 — prep punch list for 1-hour client session:

Blocker (your action, 2 min, do this before the session):
- OpenGWAS token expired 2026-08-15 (confirmed dead via `gwas.gwasinfo` — returns "Invalid token").
  Renew at https://api.opengwas.io and paste into .ieugwaspy.json. Everything live (tophits,
  associations) is blocked until this happens. Nothing below Phase 3 depends on it.

Phase 1 (0–25 min) — stats engine. No API needed, works entirely offline:
- Add MR-Egger regression (weighted, with intercept) as a function next to the IVW/WLS block in
  mr_t1d.py — operates on the already-harmonized `merged` df (beta_exp, beta_out, se_out).
- Add weighted median estimator (median of per-SNP Wald ratios, weighted by inverse variance).
- Add Cochran's Q heterogeneity stat on the real merged data.
- Restore the refutation/falsification block (negative control outcome, Sargan overID,
  leave-one-out) from archived_outputs/causality.py into causal.py — it got stripped out when
  causal.py was split off and PROJECT_STATUS.md still describes it as done. causal.py:133-141 is
  currently a dead "Real estimate" stub where this belongs.
- These three estimators are the standard IVW/Egger/weighted-median sensitivity trio —
  PROJECT_STATUS.md calls this "the highest-value remaining work" and it's fully buildable/testable
  without the API using the hardcoded INSTRUMENTS calibration data already in causal.py.

Phase 2 (25–40 min) — harden mr_t1d.py so a token expiry doesn't block work a third time:
- Fix the harmonization flip bug (~line 131): `swap` mask is built pre-filter but applied
  post-filter with `.loc`. Recompute `swap` after `merged = merged[same | swap].copy()`.
- Add palindromic SNP handling (A/T, C/G alleles) — currently silently dropped; resolve by EAF or
  exclude explicitly and print a count.
- Remove the dead `snps=HLA_SNPs` default arg on `mr_analysis` (line 76) — it's overwritten
  unconditionally on line 100 area, so the default never does anything.
- Replace the bare `except Exception` with something that surfaces the real error, not just a
  one-line print.
- Batch instrument/outcome queries instead of `return snp_list[:63]` truncation — the 63 is a
  workaround for the API's N(id)*N(variant)<=64 limit and is currently silently discarding
  instruments.
- Add CSV caching: write the harmonized `merged` df to outputs/ on every successful run. This is
  the fix for "no cached data" biting the project twice already per PROJECT_STATUS.md.

Phase 3 (40–60 min, contingent on the token being renewed) — get a real, defensible IL-6 number:
- Verify what prot-a-1538 actually is: `gwas.gwasinfo(["prot-a-1538"])`. The write-up shouldn't say
  "IL-6" until this is confirmed on the record.
- Re-run `tophits(clump=1, kb=10000, r2=0.001)` for the IL-6 exposure. Archived outputs show 4
  different rsIDs returning identical betas to 6 sig figs — they're proxies in near-perfect LD, not
  independent instruments, which means the current headline 3.483 / p=0.038 understates its own
  uncertainty. This re-run is the fix.
- Run the newly clumped instruments through the Phase 1 sensitivity stack (Egger/weighted
  median/Cochran's Q) for the first *real* (non-simulated) sensitivity analysis in the project.
- If time remains: reconcile DoWhy IV estimate (3.254) vs hand-rolled 2SLS (3.349), ~3% apart per
  PROJECT_STATUS.md.

End-of-hour deliverable if this order holds: a sensitivity toolkit that doesn't depend on API
uptime, a pipeline that won't lose work to the next token expiry, and — if the token comes back in
time — a corrected, LD-independent IL-6 estimate with real (not simulated) MR-Egger/weighted
median/Cochran's Q behind it.

To do for 8/8:
- Methodology is slightly different now. This is the final version. Can do a quick write-up on the methodology this week.
- Using GAI, attempt to add the real estimate section of the causal graph. 


To do for 5/31:
- Methodology up to this point (keep 1 paragraph)
- Editing intro and data source section
- python dowhy --> reading through the intro documentation


Implementation Plan:

Task 1: Add more models. 
- WLS w/ different intercepts.
- Median, weighted median (medians of our Wald ratios)

Task 2: Fix the scatter plot.

Task 3: Turn our models into causal models.
- use python dowhy package





**For Neel to research over the next few weeks**

Write-up:
- How have past researchers shown that HLA-DQ1 is a factor for risk of T1D?
- Is HLA-DQ1 this a good candidate for an MR study?
    - Too complex, lots of confounding results.
    - Protein could have a lot of different exposures. 
        - rs9272346
- Data sources we can write out. Look into the Sun BB IL2B study. 
- Expected outcomes. Do we already know that IL2B might increase risk for t1D? If we do, we probably expect our result to show this. 
    - As evidence for your expectations, you should use the DAG.

- Create a final hypothesis. ("Increase in IL2B proteins increases risk for T1D"... "Increase in IL2B proteins decrease/increases X/Y which, in turn, increases risk for T1D").

How to run function:

- get API key:
    - https://api.opengwas.io/
    - Sign in with github or microsoft
    - Create and copy your token
    - Create file called .ieugwaspy.json
    - In that file paste in this:
        - {"base_url": "https://api.opengwas.io/api/", "jwt": [INSERT TOKEN]}

. venv/bin/activate
python mr_t1d.py

Code:
- Add more SNPs from the dataset we found in mr_t1d.py

