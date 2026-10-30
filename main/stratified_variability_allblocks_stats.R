## T2TVAR BRAIN-STATE-RESOLVED VARIABILITY STATISTICS
# By Timo van Hattem (timo.van-hattem@uni-tuebingen.de)
# Brain Networks & Plasticity Lab, University of Tuebingen
#
# Description: Run statistics on brain-state-resolved variability in motor-evoked
# potentials (MEPs) using different EEG features and variability metrics.
# Analysis settings allow different quantification methods to be selected.

################
# Required packages
library(tidyverse)
library(lme4)        # for mixed models
library(lmerTest)    # adds p-values for lmer
library(emmeans)     # post-hoc tests
library(effectsize)  # effect size calculations
library(ez)
library(afex)
library(dplyr)
library(tidyr)
library(stringr)

# User input
response        <- "MEP"
var_method      <- "IQR"
prestim_feature <- "Phase"
site        <- "Tuebingen"
splits  <- "4phasebins"
bands_of_interest <- 'alpha'
dv <- "VARrel" #MEAN OR VARrel or VAR
if (dv == "MEAN") {
  var_method <- "SD"
}  

save_output     <- FALSE
out_dir <- "Y:/VANHATTEM/2025-08 T2TVAR/Data/REFTEP_M1/extracted_ep_info/"
out_file <- paste0(
  "T_VAR_allblocks_M1_Hjorth_",
  prestim_feature, "_",
  response, "_",
  var_method, "_",
  site, "_",
  splits, "_",
  dv,
  "_stats_log.txt"
)
out_path <- file.path(out_dir, out_file)
if (save_output) {
  sink(out_path, split = TRUE)
}

# Construct filename
file_path <- paste0(
  "Y:/VANHATTEM/2025-08 T2TVAR/Data/REFTEP_M1/extracted_ep_info/",
  "T_VAR_allblocks_M1_Hjorth_", prestim_feature, "_", response, "_", var_method, "_", site, "_", splits, "_", bands_of_interest, ".csv"
)

# Load data
T_all <- read_csv(file_path)

# Reshape to long format
if (prestim_feature %in% c('Power', 'Connectivity')) {
  
  # Convert to long format for Power/Connectivity
  T_long <- T_all %>%
    pivot_longer(
      cols = c(starts_with("VARrel"), starts_with("MEAN_"), 
               VAR_high, VAR_low, VAR_rand),
      names_to = c(".value", "state"),
      names_pattern = "(.+)_(high|low|rand)"
    ) %>%
    filter(state != "all") %>%  # Remove the "all" rows if you don't want them
    select(subject, site, type, band, VAR_method, state, VARrel, VAR, MEAN) %>%
    mutate(
      subject = factor(subject),
      band    = factor(band, levels = c("theta", "alpha", "beta", "gamma")),
      state   = factor(state, levels = c("low", "high", "rand"))
    )
  
} else if (prestim_feature %in% c('Phase', 'Phase_alternative_preprocessing_1', 'Phase_alternative_preprocessing_2')) {
  
  if (splits == "4phasebins") {
    phasebin_cols  <- paste0("phasebins_", 1:4)
    phasebin_levels <- c("Trough", "Rising", "Peak", "Falling")
  } else if (splits == "8phasebins") {
    phasebin_cols  <- paste0("phasebins_", 1:8)
    phasebin_levels <- c("Trough", "Early rising", "Rising", "Late rising",
                         "Peak", "Early falling", "Falling", "Late falling")
  }
  phasebin_pattern <- paste(c(phasebin_cols, "rand"), collapse = "|")
  recode_map <- setNames(phasebin_levels, phasebin_cols)
  recode_map["rand"] <- "Random"
  
  # Convert to long format for Phase
  T_long <- T_all %>%
    pivot_longer(
      cols = c(starts_with("VARrel"), starts_with("MEAN_"),
               all_of(paste0("VAR_", phasebin_cols)), VAR_rand),
      names_to = c(".value", "phase"),
      names_pattern = paste0("(.+)_(", phasebin_pattern, ")")
    ) %>%
    filter(phase != "all") %>%
    select(subject, site, type, band, VAR_method, phase, VARrel, VAR, MEAN) %>%
    mutate(
      subject = factor(subject),
      band    = factor(band, levels = c("theta", "alpha", "beta", "gamma")),
      phase   = recode_factor(phase, !!!recode_map),
      phase   = factor(phase, levels = c(phasebin_levels, "Random"))
    )
  
} else if (prestim_feature == 'PowerConnectivity') {
  
  # Convert to long format for Power × Connectivity
  T_long <- T_all %>%
    pivot_longer(
      cols = c(starts_with("VARrel"), starts_with("MEAN_"), 
               VAR_combos_1, VAR_combos_2, VAR_combos_3, VAR_combos_4, VAR_rand),
      names_to = c(".value", "combos"),
      names_pattern = "(.+)_(combos_1|combos_2|combos_3|combos_4|rand)"
    ) %>%
    select(subject, site, type, band, VAR_method, combos, VARrel, VAR, MEAN) %>%
    mutate(
      subject = factor(subject),
      band    = factor(band, levels = c("theta", "alpha", "beta", "gamma")),
      combos  = recode_factor(combos,
                              "combos_1" = "low Power - low FC",
                              "combos_2" = "low Power - high FC",
                              "combos_3" = "high Power - low FC",
                              "combos_4" = "high Power - high FC",
                              "rand" = "Random")
    ) %>%
    # Separate Power and Connectivity into separate columns
    mutate(
      Power = case_when(
        str_detect(combos, "low Power") ~ "Low",
        str_detect(combos, "high Power") ~ "High",
        combos == "Random" ~ "Random"
      ),
      Connectivity = case_when(
        str_detect(combos, "low FC") ~ "Low",
        str_detect(combos, "high FC") ~ "High",
        combos == "Random" ~ "Random"
      ),
      Power = factor(Power, levels = c("Low", "High", "Random")),
      Connectivity = factor(Connectivity, levels = c("Low", "High", "Random"))
    )
  T_long <- T_long %>%
      filter(combos != "Random")
  
} else if (prestim_feature == 'PowerPhase') {
  
  # Convert to long format for Power × Phase
  T_long <- T_all %>%
    mutate(subject = factor(subject)) %>%
    # Ensure the combos columns are numeric
    mutate(across(matches("^(VARrel|VAR|MEAN)_combos_\\d+$"), as.numeric)) %>%
    pivot_longer(
      cols = c(matches("^VARrel_combos_\\d+$"), matches("^VAR_combos_\\d+$"), matches("^MEAN_combos_\\d+$")),
      names_to = c(".value", "combo_n"),
      names_pattern = "(.+)_combos_(\\d+)"
    ) %>%
    select(subject, site, type, band, VAR_method, combo_n, VARrel, VAR, MEAN) %>%
    mutate(
      combo_n = as.integer(combo_n),
      Power = if_else(combo_n <= 4, "Low", "High"),
      Phase = case_when(
        combo_n %in% c(1, 5) ~ "Trough",
        combo_n %in% c(2, 6) ~ "Rising",
        combo_n %in% c(3, 7) ~ "Peak",
        combo_n %in% c(4, 8) ~ "Falling",
        TRUE ~ NA_character_
      ),
      band  = factor(band, levels = c("theta", "alpha", "beta", "gamma")),
      Power = factor(Power, levels = c("Low", "High")),
      Phase = factor(Phase, levels = c("Trough", "Rising", "Peak", "Falling"))
    ) %>%
    filter(!is.na(VARrel), !is.na(Phase))

} else if (prestim_feature == 'Power_betweensessions') {
  
  T_long <- T_all %>%
    # Create base subject ID and session number
    mutate(
      session = if_else(str_detect(subject, "_rep$"), 2L, 1L),
      subject_id = str_remove(subject, "_rep$")  # base subject without _rep
    ) %>%
    pivot_longer(
      cols = c(starts_with("VARrel"), starts_with("MEAN_"), 
               VAR_high, VAR_low, VAR_rand),
      names_to = c(".value", "state"),
      names_pattern = "(.+)_(high|low|rand)"
    ) %>%
    filter(state != "all") %>%  # Remove the "all" rows if you don't want them
    select(subject_id, site, type, band, VAR_method, session, state, VARrel, VAR, MEAN) %>%
    mutate(
      session = factor(session, levels = c(1, 2)),
      state = factor(state, levels = c("low", "high", "rand"))
    )
  
}

glimpse(T_long)

#############################################################################
# Phase

# T_long <- T_long %>% filter(!phase %in% c("Rising", "Falling"))

T_long$subject <- factor(T_long$subject)
T_long$phase   <- factor(T_long$phase)

aov_afex <- aov_ez(
  id = "subject",
  dv = dv,
  within = "phase",
  data = T_long,
  type = 3
)

aov_afex

# Pairwise paired t-tests
pairwise.t.test(
  x = T_long[[dv]],
  g = T_long$phase,
  paired = TRUE,
  p.adjust.method = "fdr"   # or "bonferroni", "fdr", "none"
)
emmeans_res <- emmeans(aov_afex, pairwise ~ phase, adjust = "fdr")
emmeans_res
summary(emmeans_res$contrasts)

if (save_output) {
  sink()
}
 
# #############################################################################
# Power / Connectivity

T_long$subject <- factor(T_long$subject)
T_long$state   <- factor(T_long$state)

aov_afex <- aov_ez(
  id = "subject",
  dv = dv,
  within = "state",
  data = T_long,
  type = 3
)

aov_afex

# Pairwise paired t-tests
pairwise.t.test(
  x = T_long[[dv]],
  g = T_long$state,
  paired = TRUE,
  p.adjust.method = "fdr"   # or "bonferroni", "fdr", "none"
)

emmeans_res <- emmeans(aov_afex, pairwise ~ state, adjust = "fdr")
emmeans_res
summary(emmeans_res$contrasts)

if (save_output) {
  sink()
}

# ################################################################
# # PowerConnectivity
# 
anova_combos <- ezANOVA(
  data = T_long,
  dv = .(VARrel),
  wid = .(subject),
  within = .(Power, Connectivity),
  type = 3,
  detailed = TRUE
)
print(anova_combos)

aov_pp <- aov(
  VARrel ~ Power * Connectivity + Error(subject/(Power*Connectivity)),
  data = T_long
)

emm <- emmeans(aov_pp, ~ Power * Connectivity)

# All 8-condition pairwise comparisons (6? actually 8 levels -> 28 comparisons)
pairs(emm, adjust = "fdr")

pairs(emmeans(aov_pp, ~ Connectivity | Power), adjust = "bonferroni")

pairs(emmeans(aov_pp, ~ Power | Connectivity), adjust = "bonferroni")

joint_tests(emm)  # gives tests of simple effects in a structured way

if (save_output) {
  sink()
}

# 
# #############################################################################
# ## PowerPhase

#!!! MANUALLY CHANGE TO MEAN OR VARREL
anova <- ezANOVA(
  data = T_long,
  dv = .(MEAN),
  wid = .(subject),
  within = .(Power, Phase),
  type = 3,
  detailed = TRUE
)

anova

aov_pp <- aov(
  MEAN ~ Power * Phase + Error(subject/(Power*Phase)),
  data = T_long
)

emm <- emmeans(aov_pp, ~ Power * Phase)

# All 8-condition pairwise comparisons (6? actually 8 levels -> 28 comparisons)
pairs(emm, adjust = "fdr")

pairs(emmeans(aov_pp, ~ Phase | Power), adjust = "fdr")

pairs(emmeans(aov_pp, ~ Power | Phase), adjust = "fdr")

joint_tests(emm)  # gives tests of simple effects in a structured way

if (save_output) {
  sink()
}


#############################################################################
# Phase - SIGNIFICANT/INSIGNIFICANT PHASE EFFECT

out_file <- paste0(
  "T_VAR_allblocks_M1_Hjorth_",
  prestim_feature, "_",
  response, "_",
  var_method, "_",
  site, "_",
  splits, "_",
  dv,
  "_sig_insig_stats_log.txt"
)
out_path <- file.path(out_dir, out_file)
if (save_output) {
  sink(out_path, split = TRUE)
}

# Define subject groups
insig_phase_subjects <- c("REFTEP_018", "REFTEP_019", "REFTEP_020", "REFTEP_021",
                          "REFTEP_023", "REFTEP_024", "REFTEP_025", "REFTEP_027",
                          "REFTEP_028", "REFTEP_029", "REFTEP_031", "REFTEP_034",
                          "REFTEP_035", "REFTEP_036", "REFTEP_038", "REFTEP_039",
                          "REFTEP_040", "REFTEP_041", "REFTEP_042", "REFTEP_043",
                          "REFTEP_045", "REFTEP_046", "REFTEP_048", "REFTEP_049",
                          "REFTEP_050")

sig_phase_subjects <- c("REFTEP_022", "REFTEP_026", "REFTEP_032", "REFTEP_044")

# Run ANOVA for each group
for (group in c("sig", "insig")) {
  
  subj_list <- if (group == "sig") sig_phase_subjects else insig_phase_subjects
  
  T_sub <- T_long %>%
    filter(subject %in% subj_list) %>%
    mutate(subject = factor(subject),
           phase   = factor(phase))
  
  cat("\n====", toupper(group), "PHASE SUBJECTS (n =", length(subj_list), ") ====\n")
  
  aov_afex <- aov_ez(
    id     = "subject",
    dv     = dv,
    within = "phase",
    data   = T_sub,
    type   = 3
  )
  
  print(aov_afex)
  
  # --- Post-hoc: emmeans ---
  cat("\n-- Post-hoc (emmeans, FDR) --\n")
  emmeans_res <- emmeans(aov_afex, pairwise ~ phase, adjust = "fdr")
  print(summary(emmeans_res$contrasts))
  
  # --- Post-hoc: pairwise t-tests ---
  cat("\n-- Pairwise t-tests (FDR) --\n")
  print(pairwise.t.test(
    x = T_sub[[dv]],
    g = T_sub$phase,
    paired = TRUE,
    p.adjust.method = "fdr"
  ))
}

# COMBINE IN MIXED ANOVA
T_long <- T_long %>%
  mutate(phase_group = case_when(
    subject %in% sig_phase_subjects   ~ "sig",
    subject %in% insig_phase_subjects ~ "insig",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(phase_group)) %>%
  mutate(subject    = factor(subject),
         phase      = factor(phase),
         phase_group = factor(phase_group, levels = c("sig", "insig")))

# Mixed ANOVA: phase (within) x phase_group (between)
aov_mixed <- aov_ez(
  id      = "subject",
  dv      = dv,
  within  = "phase",
  between = "phase_group",
  data    = T_long,
  type    = 3
)
print(aov_mixed)
emmeans_res <- emmeans(aov_mixed, pairwise ~ phase * phase_group, adjust = "fdr")
print(summary(emmeans_res$contrasts))


if (save_output) {
  sink()
}


#############################################################################
# Power between sessions

m1 <- aov_ez(
  id      = "subject_id",
  dv      = "VARrel",
  data    = T_long,
  within  = c("state", "session"),  # both are within-subject factors
  type    = 3
)

m1
summary(m1)

if (save_output) {
  sink()
}

#############################################################################
# Power all frequency bands

T_long <- T_long %>% filter(state != 'rand')
T_long$subject <- factor(T_long$subject)
T_long$state   <- factor(T_long$state)

aov_afex <- aov_ez(
  id = "subject",
  dv = dv,
  within = c("state", "band"),
  data = T_long,
  type = 3
)

aov_afex

bands <- levels(T_long$band)

results <- lapply(bands, function(b) {
  dat <- T_long %>% filter(band == b)
  high <- dat %>% filter(state == "high") %>% arrange(subject) %>% pull(VARrel)
  low  <- dat %>% filter(state == "low")  %>% arrange(subject) %>% pull(VARrel)
  t.test(high, low, paired = TRUE)
})

names(results) <- bands
results

p_vals <- sapply(results, function(x) x$p.value)
p.adjust(p_vals, method = "fdr")

if (save_output) {
  sink()
}
