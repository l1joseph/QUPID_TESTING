# Draft caption + Results text — external-tools runtime benchmark

Ready-to-paste draft prose for `fig_runtime_external.{png,pdf}` (qupid vs.
MatchIt, R Matching, and CEM on the AGP IBD cohort). Not yet wired into the
manuscript — drop in wherever you decide (main figure, supplementary, or a
paragraph in the "Qupid scales to large background pools" Results subsection).
Numbers are from the latest barnacle run (`benchmark_real/agp_external_results.tsv`,
median of 3; 159 cases / 4,400 controls; sex + age_cat + bmi_cat).

---

## Draft figure caption

**Figure SX. Qupid amortizes multi-matching where established tools do not.**
Wall-clock time to produce *k* matched sets on the AGP IBD cohort (159 cases,
4,400 controls; exact matching on sex, age category, and BMI category), for
qupid and three established host-covariate matching tools: MatchIt
(`method="nearest"` with exact strata), the R `Matching` package
(`Match`, 1:1 exact), and CEM (coarsened exact matching). Both axes are log
scaled; points are the median of 3 runs. Qupid builds the case-to-valid-control
bipartite graph once and re-randomizes the matching per iteration, so its
runtime is near-flat in *k* (0.36 s at *k*=1 to 1.22 s at *k*=100). The three
external tools have no native multi-matching capability: producing *k* matched
sets requires *k* independent invocations, so each scales linearly in *k* and
crosses above qupid by *k* ≈ 3–12. At *k*=100 qupid is 3.8× faster than MatchIt
(4.59 s), 2.9× faster than R Matching (3.48 s), and 4.1× faster than CEM
(4.94 s). **Note on CEM:** CEM produces a *single, deterministic,
order-invariant* coarsened-exact stratification — there is no way to obtain *k*
distinct matchings from it (reshuffling the input leaves its output unchanged).
Its curve therefore reflects the cost of naively re-running CEM *k* times, shown
for reference; in practice a researcher would run it once. The comparison is
thus between qupid's native multi-matching workflow and *k* repeated
single-matching calls from each external tool.

---

## Draft Results sentence/paragraph

To place qupid's scaling in the context of established software, we benchmarked
it head-to-head against three widely used host-covariate matching tools —
MatchIt, the R `Matching` package, and CEM — on the same AGP IBD configuration
(Fig. SX). Because none of these tools generates multiple distinct matchings
natively, obtaining *k* matched sets requires *k* separate invocations; each
therefore scales linearly in *k*, whereas qupid's runtime is near-flat because
the bipartite graph is constructed once. By *k* = 100, qupid (1.2 s) was
2.9–4.1× faster than the external tools (R Matching 3.5 s, MatchIt 4.6 s, CEM
4.9 s), and the gap widens with *k*. An important asymmetry underlies this
comparison: MatchIt and R Matching are deterministic in their standard
single-matching configuration but *could* be re-seeded to yield distinct
matchings (at the same per-call cost), whereas CEM is *fundamentally
order-invariant* — it yields one coarsened-exact stratification and cannot
produce *k* distinct matchings at all. Its linear curve is shown to illustrate
the cost of repeated invocation, not a recommended workflow.

---

## Methods note (if needed)

All four tools matched on the same three categorical covariates (sex, age
category, BMI category) with exact matching, on the AGP IBD cohort. Each
external tool was timed by looping its standard single-matching call *k* times
within one R session (amortizing interpreter startup, so the measured cost is
the matching call itself, not subprocess launch). qupid was timed via
`match_by_multiple` + `create_matched_pairs(iterations=k)`. SPSS FUZZY and
miMatch were not benchmarked here (SPSS is proprietary and unavailable; miMatch
matches on microbial metabolic background rather than host covariates and
returns a single cohort). Reproduce with `benchmarking/benchmark_external_tools.py`
and `benchmarking/figures/fig_runtime_external.py`.
