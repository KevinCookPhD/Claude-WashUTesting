# Analysis map: which code section and which output answers which question

Script: `protectTrajectories_analysis.R` (section headers are the lines starting with `# <number>.`; line numbers below are
from the version at the time of writing and will drift, so search for the header text). Flowchart: `flowchart/make_flowchart.py`.
Draft captions: `figure_captions_draft.md`.

**Output naming** (files are written to `figure/` when `SAVE_FIGS` / `SAVE_TABLES` are on):
`Fig<n>_<Name>`, `Table<n>_<Name>.csv`, `Data_<Name>.csv` (edge-wise results), and `SI<nn>_Fig_<Name>` / `SI<nn>_Table_<Name>.csv`.
SI numbers are set in one place, the `SI_ID` registry in section 0, so renumbering is a single edit.

| SI | Topic | Section |
|---|---|---|
| SI01 | usable scans by PMA (QA) | 1 |
| SI02 | mean connectivity matrices (Fetal, EPT, VPT, MPT) | 4g |
| SI03 | harmonization checks (scanner difference, PMA preserved) | 2 |
| SI04 | post hoc slope contrasts on the 16 edges | 8 |
| SI05 to SI16 | controls 9a to 9k (SI09 and SI10 are both from 9e) | 9 |
| SI17, SI18 | sex analyses | 11 |
| SI19 | four-group model | 12 |
| SI20 | Table 1 group comparisons | 13 |
| SI21 | preterm scans re-run through the fetal pipeline | 14 |
| SI22 | global level removed | 15 |
| SI23 | fetal vs preterm x PMA interaction (two groups): interaction results, group difference by age, lobe/hemisphere/distance summary, motion checks (see "Section 16 output files" below) | 16 |
| SI24 | sensitivity checks for the primary additive model: joint 2-df test, splines | 16 |

## Main text

| Question | Code section | Output | Key result |
|---|---|---|---|
| Who is in the sample and how was it formed? | 1 (QA), 3a to 3e | `Fig1_Flowchart` (flowchart script); `Table1_ByScan`, `Table1_ByParticipant` (section 3d) | 179 preterm infants / 346 scans, 161 fetal infants / 192 scans; groups EPT 72, VPT 73, MPT 34 |
| Is the data harmonized properly (scanner effects removed, biology kept)? | 2 | SI03 | see section 2 comments |
| How do fetal and preterm connectivity differ across all edges? | 4 (model), 4a tallies, 4b matrices, 4c chords | `Fig2_FetalVsPreterm` (A, B); `Data_FetalVsPreterm_AllEdges` | 1,404 FDR edges: 1,118 preterm > fetal, 286 fetal > preterm |
| Which homotopic (left-right) regions differ? | 4d, 4e | `Fig2_FetalVsPreterm` (C); `Table2_HomotopicFetalVsPreterm` | |
| Which connections show PMA x prematurity? | 5 (5a stats, 5b FDR, 5c Table 3) | `Table3_PmaByPrematurity`; `Data_PmaByPrematurity_AllEdges` | 16 edges; main effects (PMA, scanner, prematurity) not significant after FDR in those edges |
| What do the trajectories look like? | 6 (panels B and C), 7 and 7b (brain, panel A) | `Fig3_PmaByPrematurity` | 13 edges in panel B (MPT steeper than EPT and VPT), 3 oddball edges in panel C |
| How do EPT, VPT and MPT slopes differ pairwise? | 8 | `SI04_Table_*PosthocSlopeContrasts*`; `Data_PmaByPrematurity_SlopesSigEdges` | MPT vs EPT and MPT vs VPT differ on 15 of 16 (both FDR-significant); VPT vs EPT on 3 uncorrected, 2 after FDR |

## Controls for the fetal-vs-preterm result

| Concern | Section | Output | Result |
|---|---|---|---|
| Different amount of data (volumes) | 9e A to C, 9k | SI09, SI16 | direction preserved for nearly all primary edges; strict 80-volume cut: all 1,404 same direction, r about 0.99. The volume numbers in the 9e comments are still marked `[finalVol]` and need refreshing with the usable-volume results |
| Motion | 9e D to F, 9j | SI10, SI15 | motion is almost collinear with group (r = -0.93, VIF 7.2) so it can't be adjusted for; within groups it doesn't drive FC |
| Different PMA ranges | 9f | SI11 | age-matched: 1,174 to 1,354 of 1,404 edges stay significant, all same direction |
| Harmonization (ComBat) | 9h | SI13 | fetal-vs-preterm fully robust unharmonized |
| Distance between ROIs | 9i | SI14 | differences unrelated to distance |
| Overall connectivity level (and pipeline) | 14, 15 | SI21, SI22 | see below |

## Controls for the PMA x prematurity result

| Concern | Section | Output | Result |
|---|---|---|---|
| Birth-group cutoffs | 9a | SI05 | 8 of 16 hold under the clinical convention |
| Continuous birth GA | 9b | SI06 | 0 significant (the effect is a step, not a gradient) |
| Earliest scans driving it | 9c | SI07 (incl. trajectories figure) | 6 of 16 hold after 32 weeks, same slope pattern |
| Within vs between infants | 9d | SI08 | omnibus not significant (underpowered), within-infant slopes match |
| Postnatal age vs PMA | 9g | SI12 | 14 of 16 hold adjusted for postnatal age |
| Unharmonized data | 9h | SI13 | 9 of 16 hold |
| Motion within preterm | 9e F | SI10 | adjusting for motion leaves results unchanged |

## Additional analyses

| Question | Section | Output | Result |
|---|---|---|---|
| Sex differences, fetal vs preterm and within preterm | 11a to 11d | SI17, SI18 | no sex main effect or sex interaction survives FDR; focused test on the 16 edges: smallest FDR p = .39 (MPT-EPT) and .82 (MPT-VPT); the effect is present in both sexes, somewhat weaker in males |
| All four groups in one model | 12 | SI19 | each preterm group differs from fetal on 965 to 1,178 edges; EPT, VPT, MPT do not differ from each other in average level |
| Are the groups comparable on Table 1 variables? | 13 | SI20 | no sex difference; clinical severity falls steeply from EPT to MPT; MPT scanned later with more volumes, slightly lower FD |
| Does the processing pipeline explain the fetal-vs-preterm difference? (six preterm scans run through the fetal pipeline) | 14 (a to g) | SI21 | fetal pipeline lowers overall connectivity (mean about -0.10; -0.14 at 1.5T) and gives a different edge pattern per scan (r 0.14 to 0.51, 19 to 74% of what reliability allows); the shift is not due to censoring, global signal, band-pass range or drift; the edge differences are not consistent across scans (mean r 0.005); the cause of the level shift is unresolved |
| Does the fetal-vs-preterm difference change with age (two-group PMA x group interaction)? | 16 | `Fig4_GroupByPma` (candidate main-text figure); SI23 tables and figure; SI24 sensitivity checks | 346 edges with an interaction (63 in the shared PMA range), 92% of the 1,404 retained in the joint test (SI24); of the 346, 318 have a falling fetal and rising preterm slope |
| Where are the PMA x cohort interactions (lobe, hemisphere)? | 16 (end; the distance parts need 9i run first for `edge_dist`, `Dmat`) | SI23 `LobeHemisphereEnrichment`, `LobePairs`, `DistanceModel` | within-lobe edges 25% vs 7% significant (OR 4.7); frontal-frontal, parietal-parietal, occipital-occipital 27 to 29%, temporal-temporal at baseline (12%), frontal-occipital 0.5%; in a joint model distance dominates (OR 0.26 per SD), within-lobe (1.6) and same-hemisphere (1.5) remain independent |
| Is the interaction larger for short connections? | 16 (end) | SI23 `SizeVsDistance`, `SizeVsDistance_Pvalues` | \|coefficient\| vs distance rho = -0.39 (-0.50 within, -0.26 between hemispheres), permutation p < 0.0002; not explained by edge strength (partial rho -0.34, \|coef\|/SD -0.36); within minus between hemisphere is -0.24 but is not a matched-distance test |
| Which group produces the distance pattern? | 16 (end) | SI23 `SlopesVsDistance` | both: fetal slopes rise toward zero with distance (rho +0.25), preterm slopes fall toward zero (rho -0.40); between hemispheres the fetal side is near zero (+0.08) |
| Does head motion explain the distance pattern? | 16 (motion check) | SI23 `MotionCheck`, `MotionVsPma`, `MotionVariants`, `PartialCorrelations` (+ `_AllEdges` files) | no: meanFD (centered within cohort) as a covariate changes the interaction coefficients by a median 0.0003 (r = 0.9995; 351 vs 346 significant) and leaves the distance rho unchanged (-0.39); with motion x group and motion x PMA terms rho is -0.38; 199 motion x group terms are significant; edge-wise motion effects are unrelated to the interaction (\|rho\| <= 0.09); partial correlation controlling for each edge's motion sensitivity is unchanged (-0.39). Raw (uncentered) meanFD gives identical interaction results because group is in the model. Group-level motion cannot be adjusted for (meanFD vs group r = -0.93); fetal meanFD is unrelated to PMA (r 0.03), preterm falls with PMA (r -0.29) |
| Is the group difference just an overall-level difference? | 15 | SI22 | preterm overall level is higher (0.050 [0.036, 0.063]); after removing each scan's level, 644 of 1,118 preterm > fetal edges (58%) and all 286 fetal > preterm edges remain significant |

## Section 16 output files (all `SI23_Table_GroupByPma_<name>.csv` unless noted; SI24 files are `SI24_Table_...`)

| Name | Contents |
|---|---|
| `Summary`, `GroupDifferenceByAge` | counts of significant interactions (all scans, shared range) and the preterm - fetal difference at 30, 33, 36, 38 wk |
| `InteractionThenGroupEffect` | hierarchical summary: interaction first, then the group main effect among edges without one (1,166 of 2,657, all scans) |
| `TopEdges`, `SignificantEdges_AllScans`, `SignificantEdges_SharedRange` | the strongest 20 edges and every significant edge |
| `InteractionTable`, `InteractionTable_SharedRange` | publication-style table (full region names, slopes, coefficient, p, FDR p) |
| `AllEdges_AllScans`, `AllEdges_SharedRange` | edge-wise results for all 3,003 edges (data files, numeric p-values) |
| `JointTest_Summary`, `SplineSensitivity`, `SplineSensitivity_AllEdges` (SI24) | joint 2-df test (92% of the 1,404 retained) and the spline sensitivity check |
| `LobeHemisphereEnrichment`, `LobePairs`, `DistanceModel` | where the interactions fall: lobe, hemisphere, distance |
| `SizeVsDistance`, `SizeVsDistance_Pvalues`, `SlopesVsDistance`, `PartialCorrelations` | interaction size vs distance, with normalization checks, permutation p-values and partial correlations |
| `MotionCheck`, `MotionCheck_AllEdges`, `MotionVsPma`, `MotionVariants`, `MotionVariants_AllEdges` | motion checks |
| figures | `Fig4_GroupByPma` (main-text candidate: A = lobe-ordered matrix of slope differences, B = fetal vs preterm slope per edge), `SI23_Fig_GroupByPma` |

Slow steps (the interaction, spline and motion models, and the three permutation analyses) are cached; set `USE_CACHE <- FALSE` in section 0 to recompute. Only the distance-based outputs need section 9i run first.

## Section 14 sub-tests (all about the six preterm scans)

| Test | Question |
|---|---|
| (a) | how closely does the fetal pipeline reproduce the standard connectivity of the same scan |
| (b), (c) | does the pipeline shift resemble the group difference (with and without the global level) |
| (d) | the same for the 1.5T scans only |
| (e) | shared volumes: is it just censoring |
| (f) | rough global signal regression |
| (g) | power in 0.08 to 0.1 Hz and re-filtering to the standard band |

## Where to change things

- Switches for every analysis (`RUN_*`, `SAVE_FIGS`, `SAVE_TABLES`, `FIG_FORMATS`): top of section 0.
- Birth-group rule: `birth_cat()` in section 0 (EPT <= 28, VPT > 28 to <= 32, MPT > 32).
- Inclusion rule (4 minutes of data: 120 volumes at TR 2 s, 80 at TR 3 s): section 1 (`good4Min`, preterm scans) and the fetal controls in section 2; the flowchart script (`QA_REASON_*`) states it.
- Clinical data for Table 1: `clin_file` in section 3d (workbook `2026-08-03_KC_fetal_preterm.xlsx`).
- Figure layout settings: Fig. 2 in 4f (`fig2_*`), Fig. 3 in sections 6 and 7b (`bn_*`, `panelB_*`), frames `fig2_boxes` and `fig3_boxes`.

## Open items

- 9e volume results in the comments still marked `[finalVol]`.
- Kidokoro cutoff wording and step order for Methods and the flowchart; source of the fetal flowchart counts (535 / 385 and 343 / 224).
- `figure_captions_draft.md` has `[verify]` flags and predates the new SI numbering.
- Section 14: source of the pipeline-level difference (nuisance regression and motion handling are the remaining candidates).
