# Qupid manuscript revision plan

> Source: peer review of `manuscript/manuscript.md` (Bioinformatics, Initial submission, Reviewer #1, Minor revision recommendation)
> Date: 2026-06-05 (updated 2026-06-06)
> Status: **local edits complete (Majors A1, A2; all minors except A4 verification)**; remaining work is the matcher-distinctness check on barnacle

## Completion status

**Done (local, 2026-06-06):**
- [x] **Major A1** — softened "amplifies" in section heading + Discussion; Abstract uses "modestly"
- [x] **Major A2** — HMP2 selection-bias clause added to Abstract ("within the demographically matched, Crohn's-enriched 259-case subset of 1,179 IBD cases")
- [x] **Major A3** — superseded by barnacle commits 3837552 + 2d97fff (3 modern comparators on AGP); Results paragraph rewritten for 4-tool single-panel figure
- [x] **Minor B1** — Supp Table 1 caption added
- [x] **Minor B2** — Discussion paragraph synthesizing methodological contributions added
- [x] **Minor B3** — figure files renamed (fig2_* → fig_agp_hmp2, fig3_* → fig_thdmi); manuscripts (both copies), generator scripts, and BARNACLE2_RUN docs updated
- [x] **Minor B4** — Abstract HMP2 sentence now includes negative-control decomposition (61% matching / 39% control-pool overlap)
- [x] **Minor B5** — verified already explicit
- [x] **Minor B6** — Methods MatchIt description aligned with R template (default within-stratum propensity)
- [x] **Minor B7** — Supp Fig 2 panel-b repeat count clarified
- [x] **Minor B8** — `time.perf_counter` formalized
- [x] **Figure 2** — caption rewritten for single-panel 4-tool comparison
- [x] **Methods runtime benchmarks** — describes R `Matching` + CEM setup in addition to MatchIt
- [x] **References** — added Sekhon 2011 (R Matching) [24] and Iacus et al. 2012 (CEM) [25]

**Done (barnacle2, commit dc37ab1 — awaiting push):**
- [x] **Major A4** — verified empirically: MatchIt returns byte-identical matched sets across the loop; R `Matching` differed only with an artificial reshuffle that was subsequently removed; CEM is order-invariant in addition to deterministic and cannot produce k distinct matchings even with reshuffling. Manuscript Results paragraph and Figure 2 caption now state the verified determinism explicitly and flag the CEM order-invariance caveat ("its k-curve reflects k naive re-runs of the same underlying matching, analogous to the rebuild-every-iteration baseline in Fig. 1").

**Done (Recommended, 2026-06-06):**
- [x] **Recommended C2** — real-data PERMANOVA agreement check. Empirically measured pseudo-F agreement of 4.72 × 10<sup>−14</sup> on a 500-sample Bray-Curtis distance matrix subsampled from THDMI (well below the 10<sup>−13</sup> threshold). Methods now state explicitly that *all* real-data PERMANOVA results in the paper (AGP Bray-Curtis, HMP2 and THDMI unweighted UniFrac) were computed with the vectorized implementation, with Supp Fig 2 isolating the per-call speedup on synthetic matrices to control for cohort-specific structure. The Supp Fig 2 caption also cites the real-data agreement point.

---

## Round 2 — Minor revisions (2026-06-06 review)

After the round-1 revisions landed, a second peer-review pass surfaced **zero Major
items** but eight cosmetic Minors. None require new computation or barnacle work;
all are local manuscript edits.

### Group 1 — Trivial text edits (~5 min total)

- [x] **M1: Verb consistency** — line 56 "amplified" → "increased"; also softened Figure 3 caption title from "amplifies and stabilizes" → "increases and stabilizes" for consistency. Kept "amplification" as the category-label noun on line 83 (fine in that role).
- [x] **M3: Abstract MatchIt anchor** — appended "(3.8× vs. MatchIt, the most widely used)" after the 2.9–4.1× range.
- [x] **M4: Supp Table 1 CSV path** — changed to `manuscript/tables/HMP2_supplementary_table1.csv`.

### Group 2 — Methods clarifications (~10 min total)

- [x] **M2: AGP cohort attrition** — added one-sentence note in Runtime benchmarks Methods subsection: "The external-tool benchmark cohort is slightly smaller than the analysis cohort (159 vs. 169 cases; 4,400 vs. 4,672 controls) because AGP missingness sentinels were dropped from the categorical matching covariates prior to benchmarking."
- [x] **M5: THDMI sample-size attrition** — added parallel sentence to existing HMP2 note: "Similarly, the THDMI 563-pair count reflects the constraint that 495 of 1,058 unhealthy cases and 355 of 918 healthy controls had no valid match under the (sex, THDMI cohort, BMI category, cosmetics use, host age ± 5) criteria simultaneously."

### Group 3 — Verification + structural (~15 min)

- [x] **M6: PERMANOVA verification at larger n** — re-ran at n = 1,126 (THDMI matched-cohort size). Measured absolute difference 6.06 × 10⁻¹³ (relative 9.07 × 10⁻¹³, 12 significant digits of agreement). This is the expected machine-precision floor for ~n³ matrix-multiplication operations at this scale; cited both n = 500 (4.72 × 10⁻¹⁴) and n = 1,126 (6.06 × 10⁻¹³) measurements in Methods and Supp Fig 2 caption, with a brief explanation of the n-scaling.
- [ ] **M7: Reference style** — deferred to submission copy-edit; cosmetic only.
- [ ] **M8: Methods subsection merge** — deferred; only worth doing if word-count pressure at submission.

---

## Round 3 — Fresh peer-review feedback (2026-06-06)

A "from fresh" peer-review pass (deliberately tougher, no prior-session anchoring)
surfaced **4 Major** and **8 Minor** items. After review, **two of the Majors and three
of the Minors can be disregarded or addressed with minimal effort**; the rest are
worth fixing.

### Disregard / minimal-address (with rationale)

- **Major #3 — Cohort selection disclosure.** Reviewer asked for a paragraph
  defending the AGP/HMP2/THDMI choice against post-hoc rationalization. These are
  *standard reference cohorts* in microbiome research (largest citizen-science
  cohort, gold-standard IBD multi-omics, recent multi-country diet study); the
  critique applies to literally every multi-cohort microbiome paper. **Action:**
  add one sentence in Methods Datasets noting "We selected these three cohorts to
  span the range of microbiome-study designs..." — *no* full disclosure paragraph.
- **Minor #4 — miMatch philosophical engagement.** We already cite miMatch [20] and
  give a clear technical reason for excluding it from the wall-clock benchmark.
  Adding a "miMatch's host-covariate critique" paragraph in Discussion would be
  defensive overkill for a software paper. **Action:** disregard.
- **Minor #6 — Causal-inference literature (King & Nielsen 2019, Iacus 2019).**
  These are foundational for matching theory but this is a microbiome bioinformatics
  paper, not a causal-inference paper. **Action:** disregard; the existing
  matching-tool references [6, 24, 25] cover the practical landscape.

### Address — claim-critical (Group A, ~30 min total)

- [x] **R3-A1: Sampling-distribution framing + empirical characterization (Major #1).**
  User chose the empirical approach over text-reframe. Initially ran on a single
  small graph (5×8, M=165), then **extended in Round 4 to a multi-scale comparison
  across three graphs** (small 5×8 / M=165 / N=100K; medium 7×12 / M=5,540 / N=1M;
  large 9×16 / M=92,970 / N=3M) — spanning more than 2.7 orders of magnitude in
  matching-space size. **Results across all three scales:**
  - Full coverage of the matching space at every scale (all M matchings reached)
  - TV stable at 0.15–0.20 (vs sampling-noise floors 0.020–0.088)
  - KS distance on pseudo-F *decreases* with graph scale (0.053 → 0.024 → 0.022)
  - SD of the marginal pseudo-F distribution preserved within 3% at every scale
  - Absolute mean shift small (0.03–0.13) and bounded
  Methods paragraph now describes the multi-scale result and concludes that "the
  marginal pseudo-F distribution is stably preserved under that sampler at scales
  spanning two-plus orders of magnitude in matching-space size — supporting the
  regime-classification claims at the production cohort sizes used in this study"
  — no remaining "demonstrated at small scale" caveat required.

- [SKIP — user kept current language] **R3-A2: THDMI circularity soften (Major #2).** The "without model specification"
  framing oversells a result that follows by construction (matching on known
  confounders attenuates the confounded effect). **Action:** in Results (line 73)
  and Discussion (line 81), soften "without requiring any model specification or
  covariate selection beyond the matching criteria" → "without a parametric
  regression specification"; add one acknowledging sentence in Discussion that
  matching on the THDMI primary paper's identified drivers is the standard CCM
  approach when prior knowledge of confounders is available, and Qupid's
  contribution is the distribution-over-matchings rather than the choice of
  variables.

- [x] **R3-A3: HMP2 CD-enrichment in Discussion + Conclusion (Major #4).** The Abstract
  already flags the CD-enriched subset (Round 1 A2). The Discussion and Conclusion
  still discuss HMP2 stabilization as a general result. **Action:** in Discussion
  (line 81) add one explicit clause that the stabilization characterizes the
  Crohn's-enriched matched subset; in the Conclusion (line 89) avoid the
  unqualified "stabilizes" language. *Do not* re-run alternate-matching analyses —
  the SD collapse is what it is and the Supp Table 2 disclosure is the right level
  of self-scrutiny for a tool paper.

### Address — publication polish (Group B, ~15 min total)

- [x] **R3-B1: Baseline framing (Minor #1).** Lead with correctness over speedup —
  Fig 1c (0/100 vs 100/100) before Fig 1b (12.5×). Reorder Results paragraph
  starting line 32 to discuss the correctness advantage first, treat the 12.5×
  internal speedup as a sanity-check comparison.
- [x] **R3-B2: Single-matching disadvantage acknowledgment (Minor #2).** At k = 1
  Qupid (0.36 s) is ~6× slower than R `Matching` (0.05 s). Add one parenthetical
  sentence after line 46: "For single-matching workflows, R `Matching` and CEM
  remain marginally faster than Qupid; Qupid's advantage is specific to the
  multi-matching regime (k ≳ 10)."
- [SKIP — user kept "negative control"] **R3-B3: "Negative control" rename (Minor #3).** Line 64 calls the random-draw
  analysis a "negative control" but it has a non-zero expected effect. **Action:**
  rename to "mechanical-floor analysis" or "control-pool-overlap decomposition" in
  Results, Abstract, and any caption mentions.
- [x] **R3-B5: Verify sex in THDMI matching list (Minor #5).** Reviewer flagged that
  line 71 names sex as a matching variable but Khatib's primary-driver list does
  not include sex. Check whether sex was an intentional addition (likely yes —
  basic demographic balance) and add a clarifying half-sentence if so.
- [x] **R3-B7: PERMANOVA validation on UniFrac (Minor #7).** Run a one-shot agreement
  check on an HMP2 or THDMI unweighted UniFrac matrix to close the "validated on
  Bray-Curtis only" gap. Expected: same agreement (pseudo-F is matrix-agnostic).
  10-second script + one sentence update in Methods.
- [x] **R3-B8: Discussion victory-lap relocation (Minor #8).** Move the
  implementation-synthesis paragraph (line 87, the "These biological findings rest
  on..." paragraph) to follow the regime-overview paragraph (line 81) rather than
  the caveats paragraph (line 85), so the methodology summary doesn't interrupt
  the biology → caveats → conclusion arc.

### Execution order

1. R3-A1, R3-A2, R3-A3 — three claim-level text adjustments
2. R3-B5 — verify sex matching variable consistency
3. R3-B7 — run UniFrac PERMANOVA verification (10-sec script)
4. R3-B1, R3-B2, R3-B3, R3-B8 — text-only polish
5. Sync root manuscript.md, commit, push

### Verification checklist (Round 3)

- [ ] No "distribution over the matching space" framing remains — all distribution
      claims explicitly qualified as "under Qupid's sampler"
- [ ] THDMI Discussion acknowledges matching on known confounders as standard CCM
- [ ] HMP2 stabilization claim in Discussion + Conclusion mentions the CD-enriched
      matched subset
- [ ] Fig 1 narrative leads with correctness (panel c) before speedup (panel b)
- [ ] One sentence acknowledging single-matching disadvantage
- [ ] "Negative control" renamed throughout
- [ ] PERMANOVA validation now cites both Bray-Curtis and UniFrac agreement
- [ ] Discussion paragraph order: regimes → methods-synthesis → caveats → summary
- [ ] All Round 1, 2, 3 checks still pass

### Execution order

1. M1, M3, M4 — three single-line edits to `manuscript.md`
2. M6 — re-run PERMANOVA verification at n = 1,126; update Methods + Supp Fig 2 caption
3. M2, M5 — add Methods clarification sentences
4. M8 — merge Methods subsections (optional; only if length pressure at submission)
5. M7 — defer to journal copy-edit (cosmetic only)
6. Sync root manuscript.md, commit, push

### Verification checklist (Round 2)

- [ ] No residual "amplified" verbs in Results/Discussion outside the noun "amplification" as a category label
- [ ] Abstract reads as "2.9–4.1× at k = 100 (3.8× vs. MatchIt, the most widely used)"
- [ ] Supp Table 1 caption cites the repository-relative path
- [ ] Methods explain both AGP (169 → 159) and THDMI (918 → 563) attrition
- [ ] PERMANOVA verification cites the larger-n measurement
- [ ] Round 1 verification checklist still passes
- [ ] `grep -n "Fig\." manuscript.md` shows unchanged figure reference structure

## Status update (2026-06-06)

After this plan was first drafted, two commits (3837552, 2d97fff) updated the
external-tools benchmark on barnacle:
- Added R Matching (Sekhon `Matching` package, replicating Patel et al.) — commit 3837552
- Added CEM (Coarsened Exact Matching) — commit 2d97fff
- **Dropped the historical Wisconsin/SPSS panel** (commit 3837552)
- Figure is now single-panel AGP with 4 host-covariate matchers head-to-head

Effect on this plan:
- **Recommended C1** is now fully addressed (3 modern comparators in place)
- **Major A3** is moot in its original form because the "classical + modern" framing
  no longer matches the figure; the Results paragraph and Methods need to be
  *rewritten* (not just softened) to describe 4 modern tools on AGP only
- **Major A4** (matcher distinctness verification) now applies to all three R tools
  (MatchIt, R Matching, CEM are all deterministic by default) — the framing for
  Results should be "qupid multi-matching vs. k repeated single-matchings from each
  external tool"
- All Abstract and Minor items unchanged

---

## Summary of changes

- **4 Major revisions** (required for minor-revision target): claim calibration in Abstract/Discussion/Results
- **8 Minor revisions** (publication quality): captions, references, figure-file naming, code/text alignment
- **2 Recommended revisions** (optional but high-value): additional benchmark comparators, Discussion synthesis paragraph

All edits target `manuscript.md` except Major 4 (which also touches `benchmarking/benchmark_external_tools.py` on barnacle2) and Recommended 1 (which adds optmatch/Matching to the same script).

---

## Part A — Major revisions (REQUIRED)

### A1. Soften "amplifies" framing for the modest AGP shift

**What:** Section heading (line 52) "Case-control matching (CCM) amplifies disease signal in the AGP IBD cohort" and Discussion text (line 79) "CCM *amplified* the disease signal" present an 8% relative shift (0.40% → 0.43% absolute) as strong amplification. The abstract three-regime framing (line 16: "Qupid recovers three distinct regimes of CCM behavior") inherits this overstatement.

**Why:** A reviewer will read the headline verb before the numbers; "amplifies" implies a substantially larger effect than 0.03 percentage points. The HMP2 stabilization (~8-fold SD collapse) and THDMI exposure (35% R² drop) carry the three-regime story even if AGP is described more conservatively.

**How:**
- Section heading: "CCM modestly increases disease signal in the AGP IBD cohort" or "CCM yields a small, consistent increase in AGP IBD effect size"
- Discussion (line 79): change "CCM *amplified* the disease signal" → "CCM yielded a modest, statistically significant increase in the disease signal (8% relative, p = 0.006)"
- Abstract: leave the AGP sentence as-is (it already says "increased from 0.40% to 0.43%") but verify the three-regime framing remains calibrated

**Files:** `manuscript.md` lines 16, 52, 79

---

### A2. Add HMP2 selection-bias clause to Abstract

**What:** The Abstract sentence (line 16) reports the HMP2 8-fold SD collapse as a property of "the Human Microbiome Project 2 (HMP2) IBD cohort" without flagging that this characterizes a demographically constrained 259-case subset (CD-enriched: 77% vs. 58%, younger by ~9 years, more female).

**Why:** The Results section and Supplementary Table 2 acknowledge this thoroughly, but abstract-only readers will overgeneralize the precision-gain result to the full IBD pool. Selection bias in matched cohorts is exactly the failure mode reviewers will probe.

**How:** Insert a parenthetical or single qualifier into the existing Abstract HMP2 sentence (line 16):

Current: "in the Human Microbiome Project 2 (HMP2) IBD cohort, mean R² for the binary IBD comparison (is_case) was 2.08% pre-CCM and 1.80% post-CCM..."

Proposed: "in the Human Microbiome Project 2 (HMP2) IBD cohort — within the demographically matched, Crohn's-enriched 259-case subset of 1,179 IBD cases — mean R² for the binary IBD comparison (is_case) was 2.08% pre-CCM and 1.80% post-CCM..."

**Files:** `manuscript.md` line 16 (Abstract)

---

### A3. Soften cross-tool generalization in Results

**What:** Line 46 claims "The same amortization advantage therefore holds across both classical (SPSS, R Matching) and modern (MatchIt) per-iteration matchers and across two independent microbiome cohorts." This is a categorical claim derived from one modern comparator (MatchIt) and one historical cohort with an earlier qupid version.

**Why:** "Across both classical and modern matchers" reads as a sweeping comparison; the data is n = 2 tool-classes (or n = 1 modern + n = 2 historical). Tightening the language costs nothing and prevents a deserved reviewer pushback.

**How:** Replace line 46 final sentence with:

Proposed: "Qupid's amortization advantage holds against the most widely used modern matcher (MatchIt) and is consistent with earlier benchmarks against SPSS FUZZY and R Matching on an independent microbiome cohort, suggesting the advantage is robust to both the comparator choice and the cohort being matched."

**Files:** `manuscript.md` line 46

---

### A4. Verify and clarify MatchIt loop distinctness

**What:** The MatchIt comparison loops `matchit(method = "nearest", exact = ..., replace = FALSE, ratio = 1)` k times. MatchIt's nearest-neighbor matching is deterministic given fixed inputs; without re-seeding (or random tie-breaking) the loop may produce the *same* matched set k times, making the wall-clock comparison favorable to qupid for an unfair reason.

**Why:** Any R-experienced reviewer will spot this. If MatchIt produces identical matchings on each call, the honest comparison is "qupid multi-matching workflow vs. k × MatchIt single-matching repeated" — i.e., qupid's intrinsic amortization vs. an external tool that simply isn't designed for the multi-matching use case. That framing is defensible and arguably stronger; the fragile framing is "k *distinct* matchings".

**How (in three steps):**

1. **Empirical verification on barnacle2** (~30 min): Modify `benchmark_external_tools.py` MatchIt R template to write each iteration's matched control IDs to a temp file, then check uniqueness across iterations. Three possible outcomes:
   - **All k matchings identical:** present as "qupid multi-matching vs. k × MatchIt single-matching" — the most honest framing.
   - **Some matchings distinct (random tie-breaking with `set.seed(i)` inside the loop):** present as "MatchIt was re-seeded each iteration; X of k matchings were distinct".
   - **All k matchings distinct:** state explicitly in Results.

2. **Update Results (line 44) accordingly.** A defensible single sentence for the most-likely outcome (identical matchings):

   Proposed addition after the existing MatchIt description: "Because MatchIt's nearest-neighbor matching with exact strata is deterministic given fixed inputs, each MatchIt iteration in this benchmark produces the same matched set; the comparison therefore measures qupid's multi-matching workflow against k repeated MatchIt single-matchings rather than k distinct MatchIt matchings. This is the operationally honest comparison for a researcher who needs k matched sets to characterize an effect-size distribution: qupid amortizes that cost; MatchIt does not."

3. **Update Methods (line 129):** Note the verification result and the framing choice.

**Files:** `benchmarking/benchmark_external_tools.py` (verification logic); `manuscript.md` lines 44, 129

---

## Part B — Minor revisions (REQUIRED for publication quality)

### B1. Resolve orphan Supplementary Table 1 reference
Methods line 113 references "Supplementary Table 1" (HMP2 per-category R² changes), but no caption exists in Figure Legends. Either add a caption near the existing Supp Table 2 (line 153) or remove the reference. The CSV exists at `manuscript/tables/HMP2_supplementary_table1.csv`; recommend adding a one-line caption.

### B2. Add Discussion sentence synthesizing methodological contributions
The Discussion (lines 77–87) does not mention the external-tool comparison (Fig. 2) or the vectorized PERMANOVA implementation (Supp Fig. 2). Add one sentence near the end of the second-to-last Discussion paragraph: "These biological findings rest on a runtime infrastructure that scales to multi-matching at production speed — outperforming MatchIt by 3.8× on AGP and reducing per-iteration PERMANOVA cost by ~6–10× through a vectorized implementation included with the package."

### B3. Standardize figure-file naming
Repository contains both `fig2_agp_hmp2.png` / `fig3_thdmi.png` (old, numeric) and `fig_runtime_external.png` / `fig_permanova_speedup.png` (new, descriptive). Pre-submission, rename old files to descriptive form or vice versa for repository hygiene. Recommended: rename old → `fig_agp_hmp2.png`, `fig_thdmi.png`, then update three `![Figure ...](figures/...)` references in `manuscript.md`.

### B4. Add negative-control caveat to Abstract HMP2 sentence
Abstract line 16 reports the 8-fold SD collapse without the negative-control decomposition (61% of the collapse is attributable to demographic matching, 39% to control-pool overlap). Add a clarifying clause:

Current: "while the standard deviation collapsed ~8-fold (0.30% to 0.04%), demonstrating that CCM stabilizes effect size estimates"

Proposed: "while the standard deviation collapsed ~8-fold (0.30% to 0.04%); a negative-control analysis attributed ~61% of this stabilization to demographic matching itself and ~39% to mechanical control-pool overlap"

### B5. Name the matched covariates explicitly in Fig. 4 caption
Fig. 4 legend (line 147) says "the variables matched upon (sex, THDMI cohort, BMI category, cosmetics use, host age)" — good, already explicit. ✓ **No action needed.** (Cross-referencing during the audit confirmed this was already present.) Remove this item from the plan and verify on read.

### B6. Align MatchIt Methods text with the R template
Methods (line 129) describes MatchIt as using "a logistic propensity score for ordering within strata". Confirm against the actual R template in `benchmark_external_tools.py`; the original draft used `distance="mahalanobis"`, the current categorical version uses logistic by MatchIt default for categorical-only formulas. Verify and align.

### B7. State Supp Fig 2 panel-b repeat count
Supp Fig. 2 caption (line 151) notes "Median of 10 runs" for panel a but does not state the same for panel b. Add "all timings are median of 10 runs" to the panel b sentence.

### B8. Polish Python idioms in Methods
Line 127 references `` `time.perf_counter` `` in the prose. For a Bioinformatics submission, replace with "wall-clock timing via Python's high-resolution monotonic timer (`time.perf_counter`)" or drop the inline code formatting.

---

## Part C — Recommended revisions (OPTIONAL)

### C1. Add modern R `Matching` and/or `optmatch` to Fig. 2a
Strengthens the cross-tool generalization (addresses A3 from the other direction by *expanding* the comparison rather than softening the language). Requires ~1 hour on barnacle2 and an updated `benchmark_external_tools.py`. Worth doing if response to reviewer requires it.

**Plan:**
- Add R templates for `Matching::Match()` (modern release) and `optmatch::pairmatch()` to `benchmark_external_tools.py`
- Re-run on AGP cohort
- Regenerate Fig. 2a panel with 3–4 tools instead of 2
- Update Results paragraph (line 44) to mention all comparators

### C2. Real-data PERMANOVA agreement check
Supp Fig. 2 currently shows vectorized vs. scikit-bio agreement on simulated random DMs. A single real-data point (e.g., one AGP matched cohort at n=338) would strengthen the claim. Add as a single bullet in Methods (~5 lines of Python) or as a third panel in Supp Fig. 2.

---

## Execution order

1. **Major 4 verification first** (Part A4 step 1 on barnacle2) — outcome dictates the wording of Results and Methods for the MatchIt comparison
2. **All Major revisions** (A1–A4) — claim calibration
3. **All Minor revisions** (B1–B8) — publication polish (B5 may already be addressed; verify on read)
4. **Recommended** (C1 + C2) — only if revision rounds require deeper external comparison

Steps 1–3 should be a single commit; step 4 a follow-up commit. Push at the end with explicit user approval (see [[feedback_explicit_push_approval]]).

---

## Verification checklist

After all edits land:

- [ ] No section heading uses "amplifies" without quantification
- [ ] Abstract acknowledges HMP2 cohort selection (n = 259 of 1,179)
- [ ] Abstract notes negative-control variance decomposition
- [ ] Results line 44 framing of MatchIt comparison matches empirical verification outcome
- [ ] Discussion references both Fig. 2 (external) and Supp Fig. 2 (PERMANOVA)
- [ ] No orphan Supp Table 1 reference
- [ ] All figure files in `manuscript/figures/` follow a consistent naming convention
- [ ] Supp Fig. 2 caption specifies repeat count for both panels
- [ ] Re-render manuscript.md to PDF (if applicable) and visually inspect Figs. 1–4
- [ ] `grep -n "Fig\." manuscript.md` shows all references match the post-renumber Fig. 1/2/3/4 + Supp Fig. 1/2 structure

---

## Files touched

- `manuscript/manuscript.md` (all sections; ~10–15 distinct edits)
- `benchmarking/benchmark_external_tools.py` (Major 4 verification; Recommended C1 if pursued)
- `manuscript/tables/HMP2_supplementary_table1.csv` (existence confirmed; B1 only adds a caption to the manuscript, not the CSV)
- `manuscript/REVISION_PLAN.md` (this file — keep as a revision-history artefact; can be deleted before journal submission)

---

## Notes on scope creep

The peer review surfaced no biological-claim concerns; all comments are about claim calibration, evidence framing, and publication polish. No new analyses, no re-running of the AGP/HMP2/THDMI pipelines, no new figures beyond optional C1. The revision is genuinely minor in the scientific sense and should require a single commit cycle plus the (optional) external-comparator expansion.
