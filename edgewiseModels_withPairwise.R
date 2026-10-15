{
  ###############################################################################
  #
  #  EDGE-WISE GROUP ANALYSES OF HARMONIZED FUNCTIONAL CONNECTIVITY
  #  (main-effect models + group-by-age interaction models + pairwise subgroup
  #   contrasts)
  #
  #  PURPOSE
  #  For every connection (edge), this script fits six linear mixed-effects
  #  models. The first three ask whether groups differ in connectivity. The
  #  second three ask whether groups differ in how connectivity CHANGES with
  #  postmenstrual age (PMA), i.e. whether developmental trajectories differ.
  #  Pairwise contrasts among the preterm subgroups (EPT, VPT, MPT) are then
  #  computed from Models 2 and 2i.
  #
  #  MAIN-EFFECT MODELS (group differences in level)
  #    Model 1   Fetal vs Preterm
  #    Model 2   Fetal vs each preterm subgroup (not nested)
  #    Model 3   Nested: Fetal vs preterm average, plus subgroup differences
  #              within Preterm
  #
  #  INTERACTION MODELS (group differences in slope with age)
  #    Model 1i  Does the FC-age slope differ between Fetal and Preterm?
  #    Model 2i  Does each preterm subgroup's FC-age slope differ from Fetal's?
  #    Model 3i  Nested: does the preterm-average slope differ from Fetal's,
  #              and do the preterm subgroups differ in slope from each other?
  #
  #  PAIRWISE SUBGROUP CONTRASTS (VPT-EPT, MPT-EPT, MPT-VPT)
  #    pair2         level differences, from Model 2
  #    pair2i_level  level differences at pma_center, from Model 2i
  #    pair2i_slope  slope differences, from Model 2i
  #
  #  All models adjust for scanner and include a random intercept for each
  #  participant to account for repeated scans.
  #
  #  IMPORTANT: HARMONIZATION MUST PROTECT THE INTERACTIONS
  #  ComBat only preserves the effects listed in its formula. If the ComBat run
  #  protected only group and PMA main effects, part of any group-by-PMA
  #  interaction may have been removed as if it were scanner variance, which
  #  would bias the interaction models toward zero. Before trusting Models
  #  1i-3i, rerun ComBat with the interactions in its formula (see the chat
  #  notes) and use that harmonized output as goodFCData.
  #
  #  INPUTS (must already exist in the workspace)
  #    goodFCData
  #      Harmonized connectivity values, one row per scan and one column per
  #      connection. Row names are the full scan IDs.
  #    goodData
  #      Scan-level demographics with columns: fullID, birthPMA, scanPMA,
  #      scanner, scanType.
  #    columnMembershipCortOnly
  #      Table with a column named Connection listing the goodFCData column
  #      names to analyze.
  #
  #  OUTPUTS
  #    res1, res2, res3     main-effect model results
  #    res1i, res2i, res3i  interaction model results
  #    pair2, pair2i_level, pair2i_slope   pairwise subgroup contrasts
  #    One row per connection: coefficients (b), p-values (p), FDR-adjusted
  #    p-values (p.Adj), and an error column for connections that failed to fit.
  #
  ###############################################################################

  library(nlme)      # lme(): linear mixed-effects models
  library(emmeans)   # emmeans(): estimated group means; emtrends(): group slopes


  ###############################################################################
  #  PART 1. SETUP (runs once)
  #
  #  Only FC changes between connections, so everything else is built once here.
  ###############################################################################

  # ---- Optimizer settings -----------------------------------------------------
  # Settings from the original script: more iterations and the "optim"
  # optimizer, used when the default has trouble converging. Interaction models
  # are slightly harder to fit, so keeping these is reasonable.
  ctrl <- lmeControl(msMaxIter = 200, tolerance = 1e-4, msTol = 1e-5, opt = "optim")


  # ---- Base data frame, one row per scan --------------------------------------
  # Built from goodFCData's row names so it is in the same row order as the
  # connectivity data. Each connection's FC values can then be dropped straight
  # in inside the loop without merging.
  base <- data.frame(longID = rownames(goodFCData), stringsAsFactors = FALSE)


  # ---- Participant ID (shortID) -----------------------------------------------
  # The participant ID is the first two "_"-separated pieces of the scan ID.
  # Scans from the same participant share a shortID, which the random intercept
  # groups on. strsplit() splits every ID at once; sapply() rejoins the first two
  # pieces for each ID.
  base$shortID <- factor(sapply(strsplit(base$longID, "_", useBytes = TRUE),
                                function(x) paste(x[1:2], collapse = "_")))


  # ---- Attach demographics ----------------------------------------------------
  # match() finds each scan's row in goodData without re-sorting, keeping `base`
  # aligned with goodFCData.
  idx <- match(base$longID, goodData$fullID)

  # Should print 0. Anything else = scans with no demographics row.
  sum(is.na(idx))

  base$BirthPMA <- goodData$birthPMA[idx]   # gestational age at birth
  base$ScanPMA  <- goodData$scanPMA[idx]    # postmenstrual age at scan
  base$scanner  <- factor(goodData$scanner[idx])

  # Fetal first = reference group, so group coefficients read "Preterm minus
  # Fetal". Change the labels if yours are spelled differently.
  base$scanType <- factor(goodData$scanType[idx], levels = c("Fetal", "Preterm"))


  # ---- Birth category ---------------------------------------------------------
  # WHO definitions: EPT < 28 wk, VPT 28 to <32 wk, MPT 32 to <37 wk.
  # right = FALSE makes intervals [a, b), so exactly 28.0 wk is VPT.
  bc <- cut(base$BirthPMA, breaks = c(-Inf, 28, 32, 37),
            labels = c("EPT", "VPT", "MPT"), right = FALSE)

  # Group comes from scanType; birth age only splits the preterm group, so a
  # fetal participant later born preterm is still counted as Fetal.
  base$BirthCat <- factor(ifelse(base$scanType == "Fetal", "Fetal", as.character(bc)),
                          levels = c("Fetal", "EPT", "VPT", "MPT"))

  # Check: Fetal scans all in "Fetal"; Preterm scans all in EPT/VPT/MPT.
  # An NA column means a preterm scan with missing or 37+ wk birth age.
  table(base$scanType, base$BirthCat, useNA = "ifany")


  # ---- Nested subgroup contrasts (Models 3 and 3i) ----------------------------
  # Sum-to-zero coding within Preterm; both columns are 0 for Fetal scans.
  #
  #               pt_EPT   pt_VPT
  #     Fetal        0        0
  #     EPT          1        0
  #     VPT          0        1
  #     MPT         -1       -1
  #
  #   scanTypePreterm = average of EPT, VPT, MPT minus Fetal
  #   pt_EPT, pt_VPT  = deviation of EPT, VPT from the preterm average
  #   MPT deviation   = -(pt_EPT + pt_VPT), derived after the loop
  base$pt_EPT <- ifelse(base$BirthCat %in% "EPT", 1, ifelse(base$BirthCat %in% "MPT", -1, 0))
  base$pt_VPT <- ifelse(base$BirthCat %in% "VPT", 1, ifelse(base$BirthCat %in% "MPT", -1, 0))


  # ---- Centered PMA (interaction models) --------------------------------------
  # WHY CENTER: once a model includes group x PMA, the group coefficient becomes
  # "the group difference at PMA = 0". PMA = 0 weeks is meaningless, so that
  # coefficient would be an extrapolation far outside the data. Subtracting a
  # centering value makes the group coefficient "the group difference at that
  # PMA" instead. Centering does not change the slopes or the interaction terms.
  #
  # WHERE TO CENTER: ideally at an age where both fetal and preterm scans exist,
  # so the group difference is estimated where both groups have data. The code
  # below finds the range of PMA shared by both groups and centers at its
  # midpoint. You can replace pma_center with a meaningful value of your own
  # choosing (e.g. a specific week) as long as it lies inside that overlap.
  pma_range   <- tapply(base$ScanPMA, base$scanType, range)
  pma_range   # PMA range for each group

  overlap_lo  <- max(pma_range$Fetal[1], pma_range$Preterm[1])
  overlap_hi  <- min(pma_range$Fetal[2], pma_range$Preterm[2])

  if (overlap_lo < overlap_hi) {
    pma_center <- (overlap_lo + overlap_hi) / 2
  } else {
    # No shared age range: any group difference in the interaction models is an
    # extrapolation for at least one group. Slopes and interaction terms are
    # still interpretable; the group difference at pma_center is not.
    warning("Fetal and preterm PMA ranges do not overlap; group differences in ",
            "the interaction models are extrapolated.")
    pma_center <- mean(base$ScanPMA)
  }
  pma_center   # the age at which interaction-model group differences are estimated

  base$ScanPMA_c <- base$ScanPMA - pma_center


  # ---- Coefficient names and reference scanner --------------------------------
  # Coefficients are pulled by name, built from the actual factor levels, so a
  # change in formula or levels cannot silently point to the wrong term.
  scanner_term <- paste0("scanner", levels(base$scanner)[2])   # e.g. "scanner3T"

  # Reference scanner for estimated means; should print the 1.5T label.
  ref_scanner <- levels(base$scanner)[1]
  ref_scanner


  # ---- Empty results tables ---------------------------------------------------
  # Pre-allocated, one row per connection. "b" = estimate, "p" = p-value,
  # "se" = standard error.
  n_conn <- nrow(columnMembershipCortOnly)
  conn   <- columnMembershipCortOnly$Connection

  # --- Main-effect models ---

  # Model 1: Fetal vs Preterm
  res1 <- data.frame(Connection = conn,
                     scanPMAb = NA, scanPMAp = NA,
                     scannerb = NA, scannerp = NA,
                     scanTypeb = NA, scanTypep = NA,
                     emmMean.Fetal = NA, emmSE.Fetal = NA,
                     emmMean.Preterm = NA, emmSE.Preterm = NA,
                     error = NA)

  # Model 2: each preterm subgroup vs Fetal
  res2 <- data.frame(Connection = conn,
                     scanPMAb = NA, scanPMAp = NA,
                     scannerb = NA, scannerp = NA,
                     EPTb = NA, EPTp = NA, VPTb = NA, VPTp = NA, MPTb = NA, MPTp = NA,
                     birthCatF = NA, birthCatp = NA,
                     error = NA)

  # Model 3: nested
  res3 <- data.frame(Connection = conn,
                     scanPMAb = NA, scanPMAp = NA,
                     scannerb = NA, scannerp = NA,
                     scanTypeb = NA, scanTypep = NA,
                     ptEPTb = NA, ptEPTp = NA, ptVPTb = NA, ptVPTp = NA,
                     subF = NA, subp = NA,
                     error = NA)

  # --- Interaction models ---
  # In these tables:
  #   group columns (scanType, EPT, ...)  = group difference AT pma_center
  #   slope columns                       = change in FC per week of PMA
  #   "x" columns                         = interaction: difference in slope
  #                                         relative to Fetal (or to the
  #                                         preterm average, for Model 3i)

  # Model 1i: Fetal vs Preterm slopes
  res1i <- data.frame(Connection = conn,
                      scannerb = NA, scannerp = NA,
                      scanTypeb = NA, scanTypep = NA,            # group difference at pma_center
                      xScanTypeb = NA, xScanTypep = NA,          # Preterm slope minus Fetal slope
                      slope.Fetal = NA, slopeSE.Fetal = NA,      # FC-age slope within each group
                      slope.Preterm = NA, slopeSE.Preterm = NA,
                      error = NA)

  # Model 2i: each subgroup's slope vs Fetal's slope
  res2i <- data.frame(Connection = conn,
                      scannerb = NA, scannerp = NA,
                      xEPTb = NA, xEPTp = NA,                    # EPT slope minus Fetal slope
                      xVPTb = NA, xVPTp = NA,                    # VPT slope minus Fetal slope
                      xMPTb = NA, xMPTp = NA,                    # MPT slope minus Fetal slope
                      xBirthCatF = NA, xBirthCatp = NA,          # omnibus: any slope differences
                      slope.Fetal = NA, slopeSE.Fetal = NA,
                      slope.EPT = NA,   slopeSE.EPT = NA,
                      slope.VPT = NA,   slopeSE.VPT = NA,
                      slope.MPT = NA,   slopeSE.MPT = NA,
                      error = NA)

  # Model 3i: nested slopes
  res3i <- data.frame(Connection = conn,
                      scannerb = NA, scannerp = NA,
                      scanTypeb = NA, scanTypep = NA,            # group difference at pma_center
                      fetalSlopeb = NA, fetalSlopep = NA,        # Fetal FC-age slope
                      xScanTypeb = NA, xScanTypep = NA,          # preterm-average slope minus Fetal
                      xptEPTb = NA, xptEPTp = NA,                # EPT slope minus preterm-average slope
                      xptVPTb = NA, xptVPTp = NA,                # VPT slope minus preterm-average slope
                      xSubF = NA, xSubp = NA,                    # omnibus: do subgroup slopes differ
                      error = NA)


  # ---- Helper functions -------------------------------------------------------
  # fit_lme(): fits one mixed model with a random intercept per participant.
  # tryCatch() returns the error instead of stopping the loop if a fit fails.
  fit_lme <- function(fixed, data) tryCatch(
    lme(fixed = fixed, random = ~ 1 | shortID, data = data, control = ctrl),
    error = function(e) e)

  # bp(): estimate ("Value") and p-value for one named term from the model's
  # coefficient table (summary(model)$tTable).
  bp <- function(tt, term) tt[term, c("Value", "p-value")]

  # lin_contrast(): estimate, SE and p-value for a linear combination of fixed
  # effects, e.g. VPT minus EPT. `weights` is a named vector such as
  # c(BirthCatVPT = 1, BirthCatEPT = -1). The df is the smallest lme
  # denominator df among the terms involved (conservative).
  lin_contrast <- function(model, weights) {
    beta <- nlme::fixef(model)
    V    <- vcov(model)
    w    <- setNames(numeric(length(beta)), names(beta))
    w[names(weights)] <- weights
    est  <- sum(w * beta)
    se   <- sqrt(drop(t(w) %*% V %*% w))
    df   <- min(model$fixDF$X[names(weights)])
    c(b = est, se = se, p = 2 * pt(-abs(est / se), df))
  }

  # The three pairwise contrasts among preterm subgroups
  pair_defs  <- list(VPTvsEPT = c(BirthCatVPT = 1, BirthCatEPT = -1),
                     MPTvsEPT = c(BirthCatMPT = 1, BirthCatEPT = -1),
                     MPTvsVPT = c(BirthCatMPT = 1, BirthCatVPT = -1))
  # Same contrasts on the interaction (slope) terms of Model 2i
  slope_defs <- lapply(pair_defs, function(w) {
    names(w) <- paste0(names(w), ":ScanPMA_c"); w })

  # Empty results table: Connection plus .b / .se / .p for each contrast
  make_pair_tbl <- function(conn, defs) {
    cols <- as.vector(t(outer(names(defs), c(".b", ".se", ".p"), paste0)))
    data.frame(Connection = conn,
               matrix(NA_real_, length(conn), length(cols), dimnames = list(NULL, cols)),
               check.names = FALSE)
  }

  # --- Pairwise subgroup contrast tables ---
  pair2        <- make_pair_tbl(conn, pair_defs)    # Model 2: level differences
  pair2i_level <- make_pair_tbl(conn, pair_defs)    # Model 2i: differences at pma_center
  pair2i_slope <- make_pair_tbl(conn, pair_defs)    # Model 2i: slope differences


  ###############################################################################
  #  PART 2. LOOP OVER CONNECTIONS
  ###############################################################################

  for (i in seq_len(n_conn)) {

    if (i %% 100 == 0) print(paste0("Iteration: ", i, " of ", n_conn))

    # Copy covariates and add this connection's FC values (same row order)
    roiAnalysisFrame    <- base
    roiAnalysisFrame$FC <- as.numeric(goodFCData[, conn[i]])


    # =========================================================================
    #  MAIN-EFFECT MODELS
    # =========================================================================

    # ---- MODEL 1: Fetal vs Preterm ------------------------------------------
    #   FC ~ ScanPMA + scanner + scanType
    #   scanType = Preterm minus Fetal, at the same age and scanner
    m1 <- fit_lme(FC ~ ScanPMA + scanner + scanType, roiAnalysisFrame)

    if (inherits(m1, "error")) {
      res1$error[i] <- conditionMessage(m1)
    } else {
      tt <- summary(m1)$tTable
      res1[i, c("scanPMAb",  "scanPMAp")]  <- bp(tt, "ScanPMA")
      res1[i, c("scannerb",  "scannerp")]  <- bp(tt, scanner_term)
      res1[i, c("scanTypeb", "scanTypep")] <- bp(tt, "scanTypePreterm")

      # Estimated group means at 1.5T (where both groups have data) and at the
      # sample-mean PMA. data = is required because the model was fit inside
      # fit_lme(), so emmeans cannot find the dataset on its own.
      emm <- summary(emmeans(m1, ~ scanType, at = list(scanner = ref_scanner),
                             data = roiAnalysisFrame))
      res1[i, c("emmMean.Fetal",   "emmSE.Fetal")]   <- emm[emm$scanType == "Fetal",   c("emmean", "SE")]
      res1[i, c("emmMean.Preterm", "emmSE.Preterm")] <- emm[emm$scanType == "Preterm", c("emmean", "SE")]
    }

    # ---- MODEL 2: each preterm subgroup vs Fetal (not nested) ---------------
    #   FC ~ ScanPMA + scanner + BirthCat
    #   Three separate comparisons with Fetal: EPT, VPT, MPT.
    m2 <- fit_lme(FC ~ ScanPMA + scanner + BirthCat, roiAnalysisFrame)

    if (inherits(m2, "error")) {
      res2$error[i] <- conditionMessage(m2)
    } else {
      tt <- summary(m2)$tTable
      res2[i, c("scanPMAb", "scanPMAp")] <- bp(tt, "ScanPMA")
      res2[i, c("scannerb", "scannerp")] <- bp(tt, scanner_term)
      res2[i, c("EPTb", "EPTp")]         <- bp(tt, "BirthCatEPT")
      res2[i, c("VPTb", "VPTp")]         <- bp(tt, "BirthCatVPT")
      res2[i, c("MPTb", "MPTp")]         <- bp(tt, "BirthCatMPT")

      # Omnibus: do the four groups differ at all? (joint F-test)
      omni <- anova(m2, Terms = "BirthCat")
      res2[i, c("birthCatF", "birthCatp")] <- omni[, c("F-value", "p-value")]

      # Pairwise preterm subgroup differences (VPT-EPT, MPT-EPT, MPT-VPT)
      for (nm in names(pair_defs))
        pair2[i, paste0(nm, c(".b", ".se", ".p"))] <- lin_contrast(m2, pair_defs[[nm]])
    }

    # ---- MODEL 3: nested ----------------------------------------------------
    #   FC ~ ScanPMA + scanner + scanType + pt_EPT + pt_VPT
    #   scanType       = preterm average minus Fetal
    #   pt_EPT, pt_VPT = subgroup deviations within Preterm
    m3 <- fit_lme(FC ~ ScanPMA + scanner + scanType + pt_EPT + pt_VPT, roiAnalysisFrame)

    if (inherits(m3, "error")) {
      res3$error[i] <- conditionMessage(m3)
    } else {
      tt <- summary(m3)$tTable
      res3[i, c("scanPMAb",  "scanPMAp")]  <- bp(tt, "ScanPMA")
      res3[i, c("scannerb",  "scannerp")]  <- bp(tt, scanner_term)
      res3[i, c("scanTypeb", "scanTypep")] <- bp(tt, "scanTypePreterm")
      res3[i, c("ptEPTb",    "ptEPTp")]    <- bp(tt, "pt_EPT")
      res3[i, c("ptVPTb",    "ptVPTp")]    <- bp(tt, "pt_VPT")

      # Omnibus: do EPT, VPT, MPT differ from one another? (joint F-test)
      omni <- anova(m3, Terms = c("pt_EPT", "pt_VPT"))
      res3[i, c("subF", "subp")] <- omni[, c("F-value", "p-value")]
    }


    # =========================================================================
    #  INTERACTION MODELS
    #
    #  These add group x PMA terms, which let each group have its own FC-age
    #  slope. The interaction coefficient is the DIFFERENCE in slope between a
    #  group and its reference. A significant interaction means connectivity
    #  develops at a different rate in that group.
    #
    #  In R formulas, A * B expands to A + B + A:B (both main effects plus the
    #  interaction). Group is written first so coefficient names come out as
    #  "group:ScanPMA_c", which is what the code below pulls.
    #
    #  ScanPMA_c (centered PMA) is used so the group coefficients are group
    #  differences at pma_center rather than at PMA = 0.
    #
    #  Interpretation caution: each group's slope is estimated only over the
    #  ages where that group was scanned. If fetal and preterm scans cover
    #  different PMA ranges, a slope difference can reflect genuinely different
    #  development OR the same nonlinear curve sampled at different ages.
    # =========================================================================

    # ---- MODEL 1i: Fetal vs Preterm slopes ----------------------------------
    #   FC ~ scanner + scanType * ScanPMA_c
    #     ScanPMA_c                  = Fetal slope (Fetal is the reference)
    #     scanTypePreterm            = Preterm minus Fetal, at pma_center
    #     scanTypePreterm:ScanPMA_c  = Preterm slope minus Fetal slope
    m1i <- fit_lme(FC ~ scanner + scanType * ScanPMA_c, roiAnalysisFrame)

    if (inherits(m1i, "error")) {
      res1i$error[i] <- conditionMessage(m1i)
    } else {
      tt <- summary(m1i)$tTable
      res1i[i, c("scannerb",   "scannerp")]   <- bp(tt, scanner_term)
      res1i[i, c("scanTypeb",  "scanTypep")]  <- bp(tt, "scanTypePreterm")
      res1i[i, c("xScanTypeb", "xScanTypep")] <- bp(tt, "scanTypePreterm:ScanPMA_c")

      # Slope within each group, with standard errors. emtrends() computes the
      # FC-age slope for each group directly from the model (Fetal slope, and
      # Fetal slope + interaction for Preterm), including correct SEs.
      trd <- summary(emtrends(m1i, ~ scanType, var = "ScanPMA_c",
                              data = roiAnalysisFrame))
      res1i[i, c("slope.Fetal",   "slopeSE.Fetal")]   <- trd[trd$scanType == "Fetal",   c("ScanPMA_c.trend", "SE")]
      res1i[i, c("slope.Preterm", "slopeSE.Preterm")] <- trd[trd$scanType == "Preterm", c("ScanPMA_c.trend", "SE")]
    }

    # ---- MODEL 2i: each subgroup's slope vs Fetal's slope -------------------
    #   FC ~ scanner + BirthCat * ScanPMA_c
    #     BirthCatEPT:ScanPMA_c  = EPT slope minus Fetal slope (same for VPT, MPT)
    m2i <- fit_lme(FC ~ scanner + BirthCat * ScanPMA_c, roiAnalysisFrame)

    if (inherits(m2i, "error")) {
      res2i$error[i] <- conditionMessage(m2i)
    } else {
      tt <- summary(m2i)$tTable
      res2i[i, c("scannerb", "scannerp")] <- bp(tt, scanner_term)
      res2i[i, c("xEPTb", "xEPTp")]       <- bp(tt, "BirthCatEPT:ScanPMA_c")
      res2i[i, c("xVPTb", "xVPTp")]       <- bp(tt, "BirthCatVPT:ScanPMA_c")
      res2i[i, c("xMPTb", "xMPTp")]       <- bp(tt, "BirthCatMPT:ScanPMA_c")

      # Omnibus: do any of the four groups differ in slope? (joint F-test of
      # all three interaction terms)
      omni <- anova(m2i, Terms = "BirthCat:ScanPMA_c")
      res2i[i, c("xBirthCatF", "xBirthCatp")] <- omni[, c("F-value", "p-value")]

      # Slope within each of the four groups
      trd <- summary(emtrends(m2i, ~ BirthCat, var = "ScanPMA_c",
                              data = roiAnalysisFrame))
      for (g in c("Fetal", "EPT", "VPT", "MPT")) {
        res2i[i, c(paste0("slope.", g), paste0("slopeSE.", g))] <-
          trd[trd$BirthCat == g, c("ScanPMA_c.trend", "SE")]
      }

      # Pairwise preterm subgroup differences: level at pma_center, and slope
      for (nm in names(pair_defs)) {
        pair2i_level[i, paste0(nm, c(".b", ".se", ".p"))] <- lin_contrast(m2i, pair_defs[[nm]])
        pair2i_slope[i, paste0(nm, c(".b", ".se", ".p"))] <- lin_contrast(m2i, slope_defs[[nm]])
      }
    }

    # ---- MODEL 3i: nested slopes --------------------------------------------
    #   FC ~ scanner + (scanType + pt_EPT + pt_VPT) * ScanPMA_c
    #   which expands to all main effects plus:
    #     scanTypePreterm:ScanPMA_c  = preterm-average slope minus Fetal slope
    #     pt_EPT:ScanPMA_c           = EPT slope minus preterm-average slope
    #     pt_VPT:ScanPMA_c           = VPT slope minus preterm-average slope
    #     (MPT's slope deviation = -(EPT + VPT deviations), derived after loop)
    #
    #   This mirrors Model 3: one question about fetal vs preterm development,
    #   one about whether preterm subgroups develop differently from each other.
    m3i <- fit_lme(FC ~ scanner + (scanType + pt_EPT + pt_VPT) * ScanPMA_c,
                   roiAnalysisFrame)

    if (inherits(m3i, "error")) {
      res3i$error[i] <- conditionMessage(m3i)
    } else {
      tt <- summary(m3i)$tTable
      res3i[i, c("scannerb",    "scannerp")]    <- bp(tt, scanner_term)
      res3i[i, c("scanTypeb",   "scanTypep")]   <- bp(tt, "scanTypePreterm")
      res3i[i, c("fetalSlopeb", "fetalSlopep")] <- bp(tt, "ScanPMA_c")
      res3i[i, c("xScanTypeb",  "xScanTypep")]  <- bp(tt, "scanTypePreterm:ScanPMA_c")
      res3i[i, c("xptEPTb",     "xptEPTp")]     <- bp(tt, "pt_EPT:ScanPMA_c")
      res3i[i, c("xptVPTb",     "xptVPTp")]     <- bp(tt, "pt_VPT:ScanPMA_c")

      # Omnibus: do the preterm subgroups differ in slope from each other?
      # Joint F-test of the two subgroup-by-PMA terms. This is the test to
      # report for "do EPT, VPT, and MPT develop at different rates".
      omni <- anova(m3i, Terms = c("pt_EPT:ScanPMA_c", "pt_VPT:ScanPMA_c"))
      res3i[i, c("xSubF", "xSubp")] <- omni[, c("F-value", "p-value")]
    }
  }


  ###############################################################################
  #  PART 3. AFTER THE LOOP
  ###############################################################################

  # ---- FDR correction ---------------------------------------------------------
  # Benjamini-Hochberg FDR across connections, separately for each effect. Use
  # the ".Adj" columns for significance. Scanner is a nuisance covariate and is
  # not corrected.

  # Model 1
  res1$scanPMAp.Adj  <- p.adjust(res1$scanPMAp,  "fdr")
  res1$scanTypep.Adj <- p.adjust(res1$scanTypep, "fdr")

  # Model 2
  res2$scanPMAp.Adj  <- p.adjust(res2$scanPMAp,  "fdr")
  res2$birthCatp.Adj <- p.adjust(res2$birthCatp, "fdr")
  res2$EPTp.Adj      <- p.adjust(res2$EPTp,      "fdr")
  res2$VPTp.Adj      <- p.adjust(res2$VPTp,      "fdr")
  res2$MPTp.Adj      <- p.adjust(res2$MPTp,      "fdr")

  # Model 3
  res3$ptMPTb        <- -(res3$ptEPTb + res3$ptVPTb)   # MPT deviation
  res3$scanPMAp.Adj  <- p.adjust(res3$scanPMAp,  "fdr")
  res3$scanTypep.Adj <- p.adjust(res3$scanTypep, "fdr")
  res3$subp.Adj      <- p.adjust(res3$subp,      "fdr")

  # Model 1i
  res1i$scanTypep.Adj  <- p.adjust(res1i$scanTypep,  "fdr")
  res1i$xScanTypep.Adj <- p.adjust(res1i$xScanTypep, "fdr")

  # Model 2i
  res2i$xBirthCatp.Adj <- p.adjust(res2i$xBirthCatp, "fdr")
  res2i$xEPTp.Adj      <- p.adjust(res2i$xEPTp,      "fdr")
  res2i$xVPTp.Adj      <- p.adjust(res2i$xVPTp,      "fdr")
  res2i$xMPTp.Adj      <- p.adjust(res2i$xMPTp,      "fdr")

  # Model 3i
  # MPT slope deviation from the preterm average (sum-to-zero coding)
  res3i$xptMPTb        <- -(res3i$xptEPTb + res3i$xptVPTb)
  # Slope within each group, as point estimates (no SEs):
  #   Fetal            = fetal slope
  #   preterm average  = fetal slope + preterm interaction
  #   each subgroup    = preterm-average slope + that subgroup's deviation
  res3i$slope.Fetal    <- res3i$fetalSlopeb
  res3i$slope.PTavg    <- res3i$fetalSlopeb + res3i$xScanTypeb
  res3i$slope.EPT      <- res3i$slope.PTavg + res3i$xptEPTb
  res3i$slope.VPT      <- res3i$slope.PTavg + res3i$xptVPTb
  res3i$slope.MPT      <- res3i$slope.PTavg + res3i$xptMPTb
  res3i$scanTypep.Adj  <- p.adjust(res3i$scanTypep,  "fdr")
  res3i$xScanTypep.Adj <- p.adjust(res3i$xScanTypep, "fdr")
  res3i$xSubp.Adj      <- p.adjust(res3i$xSubp,      "fdr")

  # Pairwise subgroup contrasts: FDR across connections, separately for each
  for (nm in names(pair_defs)) {
    for (tb in c("pair2", "pair2i_level", "pair2i_slope")) {
      d <- get(tb)
      d[[paste0(nm, ".p.Adj")]] <- p.adjust(d[[paste0(nm, ".p")]], "fdr")
      assign(tb, d)
    }
  }


  # ---- Checks -----------------------------------------------------------------
  # Failed fits per model (ideally all 0; see the error column if not)
  c(m1  = sum(!is.na(res1$error)),  m2  = sum(!is.na(res2$error)),
    m3  = sum(!is.na(res3$error)),  m1i = sum(!is.na(res1i$error)),
    m2i = sum(!is.na(res2i$error)), m3i = sum(!is.na(res3i$error)))

  # Residual scanner effects; after ComBat these should be small, near zero
  summary(res1$scannerb)

  # Models 1 and 3 both estimate Fetal vs Preterm. Model 1 weights subgroups by
  # scan count; Model 3 uses their unweighted average. They should correlate
  # highly if subgroups are roughly balanced.
  cor(res1$scanTypeb, res3$scanTypeb, use = "complete.obs")

  # Same comparison for the slope differences in Models 1i and 3i
  cor(res1i$xScanTypeb, res3i$xScanTypeb, use = "complete.obs")

  # Pairwise contrasts check: m2 and m2i are still the LAST connection's models,
  # so compare against emmeans on that connection. Estimates and SEs should
  # match exactly (p-values may differ slightly if emmeans picks a different df).
  # Contrasts in emmeans are labeled like "VPT - EPT".
  pairs(emmeans(m2, ~ BirthCat, data = roiAnalysisFrame), adjust = "none")
  pair2[n_conn, ]

  pairs(emtrends(m2i, ~ BirthCat, var = "ScanPMA_c", data = roiAnalysisFrame), adjust = "none")
  pair2i_slope[n_conn, ]


  # ---- Significant connections ------------------------------------------------
  # Main-effect models
  subset(res1,  scanTypep.Adj  < .05)   # Fetal vs Preterm
  subset(res2,  birthCatp.Adj  < .05)   # any difference across 4 groups
  subset(res3,  scanTypep.Adj  < .05)   # Fetal vs preterm average
  subset(res3,  subp.Adj       < .05)   # preterm subgroups differ

  # Interaction models
  subset(res1i, xScanTypep.Adj < .05)   # Fetal vs Preterm slopes differ
  subset(res2i, xBirthCatp.Adj < .05)   # any slope difference across 4 groups
  subset(res3i, xScanTypep.Adj < .05)   # preterm-average slope differs from Fetal
  subset(res3i, xSubp.Adj      < .05)   # preterm subgroup slopes differ

  # Pairwise subgroup differences, reported only where the omnibus test is
  # significant (protected approach). Row order matches res2 / res2i.
  for (nm in names(pair_defs)) {
    padj <- paste0(nm, ".p.Adj")
    print(paste("Level difference (Model 2):", nm))
    print(pair2[which(res2$birthCatp.Adj < .05 & pair2[[padj]] < .05), ])
    print(paste("Slope difference (Model 2i):", nm))
    print(pair2i_slope[which(res2i$xBirthCatp.Adj < .05 & pair2i_slope[[padj]] < .05), ])
  }
}
