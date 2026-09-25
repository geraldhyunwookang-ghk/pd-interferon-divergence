# Deviations from the pre-analysis plan

The plan (`PRE_ANALYSIS_PLAN_v2.0.md`) was frozen by the commit that added it to this repository. Its version history (v1.0 → v2.0, all made before any expression data were loaded) is in the plan's Appendix B.

Every change made **after** that commit is recorded below, with:
- the date,
- what changed,
- why, and
- whether it was made **before or after** the relevant result was seen.

| Date | Change | Reason | Before or after seeing the relevant result? |
|---|---|---|---|
| 25 Sep 2026 | S-6 depth matching uses 20 fixed random draws (seeds 20260926–20260945) | The plan fixed the depth (T = 201) but not the number of draws; 20 keeps memory within a 16 GB laptop. S-1 keeps its pre-registered 100 | Before — no expression data loaded |
