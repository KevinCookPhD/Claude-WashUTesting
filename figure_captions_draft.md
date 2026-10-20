# Draft captions and notes (edit before use)

Items marked **[verify]** need a check against your data or protocol. Numbers come from the script output unless marked.

## Figure 1. Participant flow (flowchart)
Flow of fetal and preterm scans into the analysis. Counts are scans; participant counts give participants who lose all of
their scans at that step, so a participant with several scans can remain in the sample after losing some. Preterm:
data-quality exclusions (191 scans, 29 participants) include scans with fewer than 120 volumes (4 min at TR = 2 s),
scans that failed processing, and scans excluded for high motion; severe brain injury (Kidokoro score >= 9 **[verify cutoff]**) excluded 14 scans (9 participants).
Fetal: scans with fewer than 80 volumes (4 min at TR = 3 s) were excluded. **[verify: the fetal dataset and exclusion counts (535 / 385 and 343 / 224) come from outside the analysis project; cite the source table]**

## Figure 2. Fetal versus preterm connectivity
(A) Preterm and (B) fetal: chord diagrams and triangular matrices of the edges that differ significantly between groups
and are larger in that group (A: 1,118 edges, preterm > fetal; B: 286 edges, fetal > preterm; FDR-corrected p < .05 across
3,003 cortical edges). Matrix cells show the absolute group difference (beta, Fisher z) for each edge; regions are ordered by
lobe (frontal, temporal, parietal, occipital). (C) Estimated marginal mean connectivity (Fisher z) of homotopic
(left-right) connections in preterm and fetal scans (error bars: SE); asterisks mark regions with an FDR-corrected
significant group difference. Models: connectivity ~ PMA at scan + scanner + group with a random intercept per infant, on
LongComBat-harmonized data (scanner = batch). Group difference = preterm minus fetal.

## Figure 3. PMA x prematurity
(A) The 16 connections with a significant PMA x birth-group interaction (FDR-corrected p < .05), on a neonatal surface:
right lateral, coronal and dorsal views, with a magnified view of the boxed frontal region (the six most medial
frontal-frontal connections). Nodes are colored by lobe. (B) Thirteen connections in which connectivity (Fisher z,
scanner-adjusted) rises with PMA at scan in moderately premature infants but not in extremely or very premature infants.
(C) Three connections with a different pattern. Points are scans; lines are model-estimated trends over each group's
observed ages. Birth groups: EPT <= 28 weeks, VPT > 28 to <= 32 weeks, MPT > 32 weeks (gestational age at birth).
Model: connectivity ~ scanner + PMA x birth group with a random intercept per infant; 2-df Wald chi-square test of the
interaction. **[verify: the 13-connection description against the group slopes in your table]**

## Supplemental figure (section 9c): Fig. 3 B and C on scans after 32 weeks PMA
The same connections, order and layout as Fig. 3B-C, refit on scans after 32 weeks PMA only (301 scans from 176 infants;
45 earlier scans removed). Six of the 16 connections remain FDR-significant in this subset; the group pattern (MPT slope
steeper than EPT and VPT) is unchanged.

## Supplemental figure 1: usable scans
Usable scans per participant by PMA at scan: preterm (top) and fetal (bottom), colored by scanner (1.5 T, 3 T). Fetal scans
were all acquired at 1.5 T.

## Supplemental figure 2: mean connectivity
Mean connectivity (r, back-transformed from the mean Fisher z over scans) for each pair of regions in the fetal and
preterm groups (EPT, VPT, MPT), regions ordered by lobe.

## Table 1 notes
Volumes per scan are the usable (non-censored) volumes in the ROI time series from which connectivity was computed;
minutes of data = volumes x TR / 60 (TR = 2 s preterm, 3 s fetal). Scan-level rows are per scan; the other rows are per
infant. Mean framewise displacement (mm) was higher in fetal scans with no overlap between groups; fetal and infant
displacement may not be measured on the same footing **[verify with your pipeline]**.

## Methods sentences (from the script comments)
- Inclusion: scans needed at least 4 minutes of usable data after censoring (>= 120 volumes for preterm infants at TR 2 s;
  >= 80 volumes for fetuses at TR 3 s).
- Motion: "Mean framewise displacement was almost completely collinear with group (r = -0.93; variance inflation 7.2), so
  it could not be included as a covariate. Within the fetal scans, where motion varied most, it was unrelated to
  connectivity (1 of 3,003 edges significant after FDR) ..." (see section 9j).
- Equal volume: see the REPORT text in section 9k.
