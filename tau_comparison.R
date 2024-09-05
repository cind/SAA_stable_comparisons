# brain maps cross sectional analyses
library(dplyr)
library(ggseg)
library(stringr)
library(ggplot2)
library(patchwork)

####################################
# Getting Relevant RID's
####################################
#loading in stable RID's
stables <- read.csv("~/data/subject_saa_group_assignment.csv") %>%
  dplyr::filter(!is.na(estaget0_ABETA42)) %>%
  dplyr::select(-X) %>%
  dplyr::mutate(stable = "stable")

stable_rids <- unique(stables$RID)
table(stables$group)

amyloid_onset_ages <- read.csv("~/data/saa_estimated_csf_ages.csv") %>%
  dplyr::filter(!(estaget0_ABETA42 == "NaN")) %>%
  dplyr::select(RID, estaget0_ABETA42) %>%
  dplyr::distinct()

amyloid_onset_ages <- amyloid_onset_ages %>%
  dplyr::filter(!(RID == "6816" & estaget0_ABETA42 > 91.2)) #this RID has 2 amyloid onset ages in the file - both 91.something

#since the stable RID's are few, getting data from RID's who are not stable and limiting data from those individuals
amprion_neg <- read.csv("~/data/AMPRION_ASYN_SAA_07Jul2024.csv") %>%
  dplyr::filter(Result == "Not_Detected",
                !(RID %in% stable_rids)) %>%
  dplyr::mutate(last_saa_negative_date = EXAMDATE,
                group = "SAA- Stable",
                stable = "not stable") %>%
  dplyr::select(RID, group, stable, last_saa_negative_date) %>%
  dplyr::left_join(amyloid_onset_ages)

amprion_pos <- read.csv("~/data/AMPRION_ASYN_SAA_07Jul2024.csv") %>%
  dplyr::filter(Result == "Detected-1",
                !(RID %in% stable_rids)) %>%
  dplyr::mutate(first_saa_positive_date = EXAMDATE,
                group = "SAA+ Stable",
                stable = "not stable") %>%
  dplyr::select(RID, group, stable, first_saa_positive_date) %>%
  dplyr::left_join(amyloid_onset_ages)

stables <- merge(stables, amprion_neg, all = TRUE)
stables <- merge(stables, amprion_pos, all = TRUE)
stable_rids <- unique(stables$RID)

####################################
# Labelled Steps:
####################################
# Step 1: Adding in tau positivity (filling everything after 1st tau positive - positive, filling everything before last tau negative - negative)
####################################
# Step 2: Getting tau right before EXAMDATE (SAA-) and right after EXAMDATE (SAA+)
####################################
# Step 3: Adding nearest diagnoses
####################################
# Step 4: Getting demographic information
####################################
# Step 5: Getting p-values
####################################
# Step 6: Organizing dictionary and variable names
####################################
# Step 7: Determining which variables to plot by significance
####################################
# Step 8: Plotting
####################################
# Step 9: p_value adjustment
####################################

#############################################################################################################
#############################################################################################################
#############################################################################################################
# getting diagnosis ready, APOE, education, and birth year into a dataset
#############################################################################################################
#############################################################################################################
#############################################################################################################
##adding diagnoses as CN, MCI, or AD
adni_diagnoses <-
  readr::read_delim("~/data/DXSUM_PDXCONV_07Jun2024.csv") %>% # change to your file name and location here 
  dplyr::select(RID,EXAMDATE, DIAGNOSIS, EXAMDATE, VISCODE2) %>%
  dplyr::mutate(DX = case_when(
    (DIAGNOSIS == 1) ~ "CU",
    (DIAGNOSIS == 2) ~ "MCI",
    (DIAGNOSIS == 3) ~ "Dementia")
  ) %>%
  dplyr::select(RID, EXAMDATE, DX, VISCODE2) %>%
  dplyr::filter(!is.na(DX))
adni_diagnoses <- adni_diagnoses %>%
  dplyr::rename(DX.DATE = EXAMDATE) %>%
  dplyr::filter(RID %in% stable_rids)

#getting APOE type
apoeres <- read.csv("~/data/APOERES.csv") %>% 
  dplyr::filter(RID %in% stable_rids) %>%
  dplyr::select(RID, APGEN1, APGEN2) %>%
  dplyr::mutate(apoe = paste(paste("E", APGEN1, sep = ""), paste("E", APGEN2, sep = ""), sep = "/"),
                RID = as.character(RID)) %>%
  dplyr::distinct()

#getting information to calculate age
dem <- read.csv("~/data/demographics.csv") %>%
  dplyr::mutate(day = 1,
                birthdate = as.character(paste(day, PTDOB, sep = "/")))

dem <- dem %>%
  dplyr::mutate(birthdate = as.Date(birthdate, format = "%d/%m/%Y")) %>%
  dplyr::select(RID, PTGENDER, PTEDUCAT, birthdate) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(birthdate))

dem <- merge(dem, apoeres, all = TRUE) %>%
  dplyr::mutate(apoe = case_when(apoe == "E2/E3" | apoe == "E2/E2" ~ "E2",
                                 apoe == "E3/E3" ~ "E3",
                                 apoe == "E3/E4" | apoe == "E4/E4" ~ "E4"),
                RID = as.character(RID))

dem <- dem %>%
  dplyr::filter(RID %in% stable_rids,
                !(RID == "672" & is.na(PTEDUCAT)))

#getting amyloid positivity status
amyloid_pet <- readr::read_delim("~/data/UCBERKELEY_AMY_6MM_31Jul2024.csv") [,1:16]
amyloid_pet <- amyloid_pet %>% dplyr::rename(suvr_summary=SUMMARY_SUVR,
                                             Centiloid=CENTILOIDS,
                                             AmyloidPosPET=AMYLOID_STATUS,
                                             EXAMDATE_pet = SCANDATE) %>%
  dplyr::mutate(EXAMDATE_pet = as.Date(EXAMDATE_pet),
                RID = as.character(RID)) %>%
  dplyr::select(RID, EXAMDATE_pet, AmyloidPosPET, Centiloid, suvr_summary)

amyloid_pet <- amyloid_pet %>%
  dplyr::filter(RID %in% stable_rids)

rm(apoeres)

#############################################################################################################
#############################################################################################################
#############################################################################################################
#                                           Regional tau
#############################################################################################################
#############################################################################################################
#############################################################################################################

#regional pet values for tau
tau_roi_cs <- read.csv("~/data/UCBERKELEY_TAU_6MM_17Jun2024.csv") %>% 
  dplyr::filter(qc_flag == 2 | qc_flag == 1,
                RID %in% stable_rids) %>%
  dplyr::mutate(SCANDATE = as.Date(SCANDATE))

#adding in the stable cases
tau_roi_cs <- merge(tau_roi_cs, stables, by = "RID")

####################################
# Step 1: Adding in tau positivity (filling everything after 1st tau positive - positive, filling everything before last tau negative - negative)
####################################

#adding in tau pet
#getting tau negative cases with last tau negative dates
last_a_negative_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 0) %>%
  dplyr::mutate(last_a_neg_date_pet = as.Date(max(EXAMDATE_pet)),
                AmyNeg_Centiloid = Centiloid) %>%
  dplyr::select(RID, Centiloid, EXAMDATE_pet, last_a_neg_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
last_a_negative_date_pet <- last_a_negative_date_pet %>%
  dplyr::filter(EXAMDATE_pet == last_a_neg_date_pet) %>%
  dplyr::mutate(AmyNeg_Centiloid = Centiloid) %>%
  dplyr::select(-Centiloid, -EXAMDATE_pet, - AmyNeg_Centiloid)

# getting tau positive cases with first tau positive dates
first_a_positive_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 1) %>%
  dplyr::mutate(first_a_pos_date_pet = as.Date(min(EXAMDATE_pet))) %>%
  dplyr::select(RID, Centiloid, EXAMDATE_pet, first_a_pos_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
first_a_positive_date_pet <- first_a_positive_date_pet %>%
  dplyr::filter(EXAMDATE_pet == first_a_pos_date_pet) %>%
  dplyr::mutate(AmyPos_Centiloid = Centiloid) %>%
  dplyr::select(-Centiloid, -EXAMDATE_pet, - AmyPos_Centiloid)

#merging tau info back into dataset so that I can use the first and last tau status dates
tau_roi_cs <- merge(tau_roi_cs, last_a_negative_date_pet, by = "RID", all.x = TRUE)
tau_roi_cs <- merge(tau_roi_cs, first_a_positive_date_pet, by = "RID", all.x = TRUE) %>%
  dplyr::mutate(last_saa_negative_date = as.Date(last_saa_negative_date),
                first_saa_positive_date = as.Date(first_saa_positive_date))

tau_roi_cs <- tau_roi_cs %>%
  dplyr::mutate(AmyloidPosPET = case_when(SCANDATE <= last_a_neg_date_pet ~ 0,
                                          SCANDATE >= first_a_pos_date_pet ~ 1)) %>%
  dplyr::select(RID, SCANDATE, last_a_neg_date_pet, first_a_pos_date_pet, AmyloidPosPET, estaget0_ABETA42, group, first_saa_positive_date, last_saa_negative_date, ends_with("SUVR"))

rm(first_a_positive_date_pet, last_a_negative_date_pet)

####################################
# Step 2: getting SCAN right before EXAMDATE (SAA-) and right after EXAMDATE (SAA+)
####################################

tau_roi_cs_saa_neg <- tau_roi_cs %>%
  dplyr::filter(group == "SAA- Stable")
tau_roi_cs_saa_neg <- tau_roi_cs_saa_neg %>%
  dplyr::mutate(RID = as.character(RID))
tau_roi_cs_saa_pos <- tau_roi_cs %>%
  dplyr::filter(group == "SAA+ Stable")
tau_roi_cs_saa_pos <- tau_roi_cs_saa_pos %>%
  dplyr::mutate(RID = as.character(RID))

#now pulling the EXAMDATE right before the last SAA negative date
tau_roi_cs_saa_neg <- tau_roi_cs_saa_neg %>%
  dplyr::left_join(dem) %>%
  dplyr::mutate(age = lubridate::time_length(difftime(as.Date(SCANDATE), as.Date(birthdate)), "years"))
tau_roi_cs_saa_neg <- tau_roi_cs_saa_neg %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(SCANDATE), last_saa_negative_date), "years")) %>%
  dplyr::select(RID, SCANDATE, last_a_neg_date_pet, first_a_pos_date_pet, age, PTGENDER, PTEDUCAT, apoe, time_diff, AmyloidPosPET, estaget0_ABETA42, group, first_saa_positive_date, last_saa_negative_date, ends_with("SUVR"))
tau_roi_cs_saa_neg <- tau_roi_cs_saa_neg %>%
  dplyr::filter(time_diff <= 0)
tau_roi_cs_saa_neg <- tau_roi_cs_saa_neg %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_diff == max(time_diff)) %>%
  dplyr::ungroup()
tau_roi_cs_saa_neg <- tau_roi_cs_saa_neg %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::arrange(SCANDATE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

#now pulling the EXAMDATE right after the first SAA positive date
tau_roi_cs_saa_pos <- tau_roi_cs_saa_pos %>%
  dplyr::left_join(dem) %>%
  dplyr::mutate(age = lubridate::time_length(difftime(as.Date(SCANDATE), as.Date(birthdate)), "years"))
tau_roi_cs_saa_pos <- tau_roi_cs_saa_pos %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(SCANDATE), first_saa_positive_date), "years")) %>%
  dplyr::select(RID, SCANDATE, last_a_neg_date_pet, first_a_pos_date_pet, age, PTGENDER, PTEDUCAT, apoe, time_diff, AmyloidPosPET, estaget0_ABETA42, group, first_saa_positive_date, last_saa_negative_date, ends_with("SUVR"))
tau_roi_cs_saa_pos <- tau_roi_cs_saa_pos %>%
  dplyr::filter(time_diff >= 0)
tau_roi_cs_saa_pos <- tau_roi_cs_saa_pos %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
tau_roi_cs_saa_pos <- tau_roi_cs_saa_pos %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::arrange(SCANDATE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

#putting the two datasets with tau status and correct visit dates together
tau_roi_cs <- rbind(tau_roi_cs_saa_neg, tau_roi_cs_saa_pos) #604 RID's

####################################
# Step 3: adding nearest diagnoses
####################################
#adding diagnosis in
tau_roi_cs <- tau_roi_cs %>%
  dplyr::left_join(adni_diagnoses %>% 
                     dplyr::mutate(RID = as.character(RID))) %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(SCANDATE), DX.DATE), "years")) %>%
  dplyr::distinct() %>%
  dplyr::select(RID, SCANDATE, last_a_neg_date_pet, first_a_pos_date_pet, age, PTGENDER, PTEDUCAT, apoe, time_diff, AmyloidPosPET, estaget0_ABETA42, group, first_saa_positive_date, last_saa_negative_date, DX.DATE, DX, ends_with("SUVR"))
tau_roi_cs <- tau_roi_cs %>%
  dplyr::filter(!is.na(time_diff)) %>%
  dplyr::mutate(time_diff = abs(time_diff))
tau_roi_cs <- as.data.frame(tau_roi_cs)
tau_roi_cs <- tau_roi_cs %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
tau_roi_cs <- tau_roi_cs %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::arrange(SCANDATE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

####################################
# Step 4: getting demographic information
####################################
#first getting rid of rows that have no information on features
tau_roi_cs <- tau_roi_cs[which(rowMeans(!is.na(tau_roi_cs)) > 0.8), ]

#getting distributions of diagnoses, gender, apoe, tau status
table(tau_roi_cs$group, tau_roi_cs$DX)
table(tau_roi_cs$group, tau_roi_cs$PTGENDER)
table(tau_roi_cs$group, tau_roi_cs$apoe)
table(tau_roi_cs$group, tau_roi_cs$AmyloidPosPET)
table(tau_roi_cs$group)

#getting averages of age and education
mean(tau_roi_cs_saa_pos$age)
sd(tau_roi_cs_saa_pos$age)
mean(tau_roi_cs_saa_neg$age)
sd(tau_roi_cs_saa_neg$age)
mean(tau_roi_cs_saa_pos$PTEDUCAT)
sd(tau_roi_cs_saa_pos$PTEDUCAT)
mean(tau_roi_cs_saa_neg$PTEDUCAT)
sd(tau_roi_cs_saa_neg$PTEDUCAT)

rm(tau_roi_cs_saa_neg, tau_roi_cs_saa_pos, amyloid_pet)

####################################
# Step 5: getting p-values
####################################
tau_roi_cs <- tau_roi_cs %>%
  dplyr::select(-VENTRICLE_5TH_SUVR, - INFERIORCEREBELLUM_SUVR, -NON_WM_HYPOINTENSITIES_SUVR)
features_tau <- tau_roi_cs %>%
  dplyr::select(ends_with("SUVR"))
featurenames_tau <- names(features_tau)

aseg_data <- as.data.frame(aseg)
aseg_data$label
aparc_data <- as.data.frame(dk)
aparc_data$label

#stratifying data
#pairwise comparisons of each region and saving to a list then data frame
tau_roi_cs_cn <- tau_roi_cs %>%
  dplyr::filter(DX == "CU")
tau_roi_cs_mci <- tau_roi_cs %>%
  dplyr::filter(DX == "MCI")
tau_roi_cs_ad <- tau_roi_cs %>%
  dplyr::filter(DX == "Dementia")

#full data
tau_p_values <- vector()
tau_roi_names <- vector()
tau_effect_sizes <- vector()
tau_n_all <- vector()

for (f in 1:length(featurenames_tau)){
  feature <- featurenames_tau[f]
  print(paste("tau data regression: ", feature))
  lm_feature <- lm(scale(tau_roi_cs[, 16+f]) ~ group + scale(age) + scale(PTGENDER) + PTEDUCAT + apoe + DX + AmyloidPosPET, data = tau_roi_cs, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  tau_effect_sizes <- append(tau_effect_sizes, estimate)
  tau_p_values <- append(tau_p_values, group_p_value)
  tau_roi_names <- append(tau_roi_names, feature)
  tau_n_all <- append(tau_n_all, length(lm_feature$residuals))
}

p_values_data <- data.frame(tau_p_values)
tau_roi_names <- gsub("CTX_", "", tau_roi_names) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data$ROI <- tau_roi_names #creating column for the roi names for each p value
p_values_data$Effect_Size <- tau_effect_sizes
p_values_data$n <- tau_n_all

#CN
tau_p_values_cn <- vector()
tau_roi_names_cn <- vector()
tau_effect_sizes_cn <- vector()
tau_n_cn <- vector()

for (f in 1:length(featurenames_tau)){
  feature <- featurenames_tau[f]
  print(paste("tau data regression: ", feature))
  lm_feature <- lm(scale(tau_roi_cs_cn[, 16+f]) ~ group + scale(age) + scale(PTGENDER) + PTEDUCAT + apoe + AmyloidPosPET, data = tau_roi_cs_cn, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  tau_effect_sizes_cn <- append(tau_effect_sizes_cn, estimate)
  tau_p_values_cn <- append(tau_p_values_cn, group_p_value)
  tau_roi_names_cn <- append(tau_roi_names_cn, feature)
  tau_n_cn <- append(tau_n_cn, length(lm_feature$residuals))
}

p_values_data_cn <- data.frame(tau_p_values_cn)
tau_roi_names_cn <- gsub("CTX_", "", tau_roi_names_cn) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data_cn$ROI <- tau_roi_names_cn #creating column for the roi names for each p value
p_values_data_cn$Effect_Size <- tau_effect_sizes_cn
p_values_data_cn$n <- tau_n_cn

#MCI
tau_p_values_mci <- vector()
tau_roi_names_mci <- vector()
tau_effect_sizes_mci <- vector()
tau_n_mci <- vector()

for (f in 1:length(featurenames_tau)){
  feature <- featurenames_tau[f]
  print(paste("tau data regression: ", feature))
  lm_feature <- lm(scale(tau_roi_cs_mci[, 16+f]) ~ group + scale(age) + scale(PTGENDER) + PTEDUCAT + apoe + AmyloidPosPET, data = tau_roi_cs_mci, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  tau_effect_sizes_mci <- append(tau_effect_sizes_mci, estimate)
  tau_p_values_mci <- append(tau_p_values_mci, group_p_value)
  tau_roi_names_mci <- append(tau_roi_names_mci, feature)
  tau_n_mci <- append(tau_n_mci, length(lm_feature$residuals))
}

p_values_data_mci <- data.frame(tau_p_values_mci)
tau_roi_names_mci <- gsub("CTX_", "", tau_roi_names_mci) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data_mci$ROI <- tau_roi_names_mci #creating column for the roi names for each p value
p_values_data_mci$Effect_Size <- tau_effect_sizes_mci
p_values_data_mci$n <- tau_n_mci

#AD
tau_p_values_ad <- vector()
tau_roi_names_ad <- vector()
tau_effect_sizes_ad <- vector()
tau_n_ad <- vector()

for (f in 1:length(featurenames_tau)){
  feature <- featurenames_tau[f]
  print(paste("tau data regression: ", feature))
  lm_feature <- lm(scale(tau_roi_cs_ad[, 16+f]) ~ group + scale(age) + scale(PTGENDER) + PTEDUCAT + apoe + AmyloidPosPET, data = tau_roi_cs_ad, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  tau_effect_sizes_ad <- append(tau_effect_sizes_ad, estimate)
  tau_p_values_ad <- append(tau_p_values_ad, group_p_value)
  tau_roi_names_ad <- append(tau_roi_names_ad, feature)
  tau_n_ad <- append(tau_n_ad, length(lm_feature$residuals))
}

p_values_data_ad <- data.frame(tau_p_values_ad)
tau_roi_names_ad <- gsub("CTX_", "", tau_roi_names_ad) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data_ad$ROI <- tau_roi_names_ad #creating column for the roi names for each p value
p_values_data_ad$Effect_Size <- tau_effect_sizes_ad
p_values_data_ad$n <- tau_n_ad

####################################
# Step 6: organizing dictionary and variable names
####################################

#getting the dictionary so that the ROI names are listed
dict <- as.data.frame(featurenames_tau)
dict <- dict %>%
  dplyr::rename(FLDNAME = featurenames_tau)

# getting the ROI name as an additional column
p_values_data_test <- merge(p_values_data, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE) %>%
  dplyr::mutate(label_aseg = case_when(ROI == "LEFT_THALAMUS_PROPER_SUVR" ~ "Left-Thalamus-Proper",
                                       ROI == "RIGHT_THALAMUS_PROPER_SUVR" ~ "Right-Thalamus-Proper",
                                       ROI == "LEFT_LATERAL_VENTRICLE_SUVR" ~ "Left-Lateral-Ventricle",
                                       ROI == "RIGHT_LATERAL_VENTRICLE_SUVR" ~ "Right-Lateral-Ventricle",
                                       ROI == "LEFT_HIPPOCAMPUS_SUVR" ~ "Left-Hippocampus",
                                       ROI == "RIGHT_HIPPOCAMPUS_SUVR" ~ "Right-Hippocampus",
                                       ROI == "LEFT_PUTAMEN_SUVR" ~ "Left-Putamen",
                                       ROI == "RIGHT_PUTAMEN_SUVR" ~ "Right-Putamen", 
                                       ROI == "LEFT_AMYGDALA_SUVR" ~ "Left-Amygdala",
                                       ROI == "RIGHT_AMYGDALA_SUVR" ~ "Right-Amygdala",
                                       ROI == "LEFT_PALLIDUM_SUVR" ~ "Left-Pallidum",
                                       ROI == "RIGHT_PALLIDUM_SUVR" ~ "Right-Pallidum",
                                       ROI == "LEFT_CAUDATE_SUVR" ~ "Left-Caudate",
                                       ROI == "RIGHT_CAUDATE_SUVR" ~ "Right-Caudate",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "LH_CAUDALMIDDLEFRONTAL_SUVR" ~ "lh_caudalmiddlefrontal",
                                       ROI == "LH_FUSIFORM_SUVR" ~ "lh_fusiform",
                                       ROI == "LH_INFERIORPARIETAL_SUVR" ~ "lh_inferiorparietal",
                                       ROI == "LH_INFERIORTEMPORAL_SUVR" ~ "lh_inferiortemporal",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_MIDDLETEMPORAL_SUVR" ~ "lh_middletemporal",
                                       ROI == "LH_PARSOPERCULARIS_SUVR" ~ "lh_parsopercularis",
                                       ROI == "LH_PARSORBITALIS_SUVR" ~ "lh_parsorbitalis",
                                       ROI == "LH_PARSTRIANGULARIS_SUVR" ~ "lh_parstriangularis",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "LH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "lh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "LH_SUPERIORPARIETAL_SUVR" ~ "lh_superiorparietal",
                                       ROI == "LH_SUPERIORTEMPORAL_SUVR" ~ "lh_superiortemporal",
                                       ROI == "LH_SUPRAMARGINAL_SUVR" ~ "lh_supramarginal",
                                       ROI == "LH_TEMPORALPOLE_SUVR" ~ "lh_temporalpole",
                                       ROI == "LH_TRANSVERSETEMPORAL_SUVR" ~ "lh_transversetemporal",
                                       ROI == "LH_INSULA_SUVR" ~ "lh_insula",
                                       ROI == "LH_CAUDALANTERIORCINGULATE_SUVR" ~ "lh_caudalanteriorcingulate",
                                       ROI == "LH_CUNEUS_SUVR" ~ "lh_cuneus",
                                       ROI == "LH_ENTORHINAL_SUVR" ~ "lh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "LH_ISTHMUSCINGULATE_SUVR" ~ "lh_isthmuscingulate",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_LINGUAL_SUVR" ~ "lh_lingual",
                                       ROI == "LH_MEDIALORBITOFRONTAL_SUVR" ~ "lh_medialorbitofrontal",
                                       ROI == "LH_PARAHIPPOCAMPAL_SUVR" ~ "lh_parahippocampal",
                                       ROI == "LH_PARACENTRAL_SUVR" ~ "lh_paracentral",
                                       ROI == "LH_PERICALCARINE_SUVR" ~ "lh_pericalcarine",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_POSTERIORCINGULATE_SUVR" ~ "lh_posteriorcingulate",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "LH_PRECUNEUS_SUVR" ~ "lh_precuneus",
                                       ROI == "LH_ROSTRALANTERIORCINGULATE_SUVR" ~ "lh_rostralanteriorcingulate",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORPARIETAL_SUVR" ~ "rh_superiorparietal",
                                       ROI == "LH_FRONTALPOLE_SUVR" ~ "lh_frontalpole",
                                       ROI == "RH_TEMPORALPOLE_SUVR" ~ "rh_temporalpole",
                                       ROI == "RH_CAUDALANTERIORCINGULATE_SUVR" ~ "rh_caudalanteriorcingulate",
                                       ROI == "RH_CUNEUS_SUVR" ~ "rh_cuneus",
                                       ROI == "RH_ENTORHINAL_SUVR" ~ "rh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "RH_ISTHMUSCINGULATE_SUVR" ~ "rh_isthmuscingulate",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_LINGUAL_SUVR" ~ "rh_lingual",
                                       ROI == "RH_MEDIALORBITOFRONTAL_SUVR" ~ "rh_medialorbitofrontal",
                                       ROI == "RH_PARAHIPPOCAMPAL_SUVR" ~ "rh_parahippocampal",
                                       ROI == "RH_PARACENTRAL_SUVR" ~ "rh_paracentral",
                                       ROI == "RH_PERICALCARINE_SUVR" ~ "rh_pericalcarine",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_POSTERIORCINGULATE_SUVR" ~ "rh_posteriorcingulate",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "RH_PRECUNEUS_SUVR" ~ "rh_precuneus",
                                       ROI == "RH_ROSTRALANTERIORCINGULATE_SUVR" ~ "rh_rostralanteriorcingulate",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_FRONTALPOLE_SUVR" ~ "rh_frontalpole",
                                       ROI == "LH_TEMPORALPOLE_SUVR"  ~ "Lh_temporalpole",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "RH_CAUDALMIDDLEFRONTAL_SUVR" ~ "rh_caudalmiddlefrontal",
                                       ROI == "RH_BANKSSTS_SUVR" ~ "rh_bankssts",
                                       ROI == "RH_INFERIORPARIETAL_SUVR" ~ "rh_inferiorparietal",
                                       ROI == "RH_INFERIORTEMPORAL_SUVR" ~ "rh_inferiortemporal",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_MIDDLETEMPORAL_SUVR" ~ "rh_middletemporal",
                                       ROI == "RH_PARSOPERCULARIS_SUVR" ~ "rh_parsopercularis",
                                       ROI == "RH_PARSORBITALIS_SUVR" ~ "rh_parsorbitalis",
                                       ROI == "RH_PARSTRIANGULARIS_SUVR" ~ "rh_parstriangularis",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "RH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "rh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORTEMPORAL_SUVR" ~ "rh_superiortemporal",
                                       ROI == "RH_SUPRAMARGINAL_SUVR" ~ "rh_supramarginal",
                                       ROI == "RH_TEMPORALPOLE_SUVR"  ~ "rh_temporalpole",
                                       ROI == "RH_TRANSVERSETEMPORAL_SUVR" ~ "rh_transversetemporal",
                                       ROI == "RH_INSULA_SUVR" ~ "rh_insula"))

p_tau_aseg <- merge(p_values_data_test, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_tau_aparc <- merge(p_values_data_test, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

# CN
p_values_data_test_cn <- merge(p_values_data_cn, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE) %>%
  dplyr::mutate(label_aseg = case_when(ROI == "LEFT_THALAMUS_PROPER_SUVR" ~ "Left-Thalamus-Proper",
                                       ROI == "RIGHT_THALAMUS_PROPER_SUVR" ~ "Right-Thalamus-Proper",
                                       ROI == "LEFT_LATERAL_VENTRICLE_SUVR" ~ "Left-Lateral-Ventricle",
                                       ROI == "RIGHT_LATERAL_VENTRICLE_SUVR" ~ "Right-Lateral-Ventricle",
                                       ROI == "LEFT_HIPPOCAMPUS_SUVR" ~ "Left-Hippocampus",
                                       ROI == "RIGHT_HIPPOCAMPUS_SUVR" ~ "Right-Hippocampus",
                                       ROI == "LEFT_PUTAMEN_SUVR" ~ "Left-Putamen",
                                       ROI == "RIGHT_PUTAMEN_SUVR" ~ "Right-Putamen", 
                                       ROI == "LEFT_AMYGDALA_SUVR" ~ "Left-Amygdala",
                                       ROI == "RIGHT_AMYGDALA_SUVR" ~ "Right-Amygdala",
                                       ROI == "LEFT_PALLIDUM_SUVR" ~ "Left-Pallidum",
                                       ROI == "RIGHT_PALLIDUM_SUVR" ~ "Right-Pallidum",
                                       ROI == "LEFT_CAUDATE_SUVR" ~ "Left-Caudate",
                                       ROI == "RIGHT_CAUDATE_SUVR" ~ "Right-Caudate",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "LH_CAUDALMIDDLEFRONTAL_SUVR" ~ "lh_caudalmiddlefrontal",
                                       ROI == "LH_FUSIFORM_SUVR" ~ "lh_fusiform",
                                       ROI == "LH_INFERIORPARIETAL_SUVR" ~ "lh_inferiorparietal",
                                       ROI == "LH_INFERIORTEMPORAL_SUVR" ~ "lh_inferiortemporal",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_MIDDLETEMPORAL_SUVR" ~ "lh_middletemporal",
                                       ROI == "LH_PARSOPERCULARIS_SUVR" ~ "lh_parsopercularis",
                                       ROI == "LH_PARSORBITALIS_SUVR" ~ "lh_parsorbitalis",
                                       ROI == "LH_PARSTRIANGULARIS_SUVR" ~ "lh_parstriangularis",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "LH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "lh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "LH_SUPERIORPARIETAL_SUVR" ~ "lh_superiorparietal",
                                       ROI == "LH_SUPERIORTEMPORAL_SUVR" ~ "lh_superiortemporal",
                                       ROI == "LH_SUPRAMARGINAL_SUVR" ~ "lh_supramarginal",
                                       ROI == "LH_TEMPORALPOLE_SUVR" ~ "lh_temporalpole",
                                       ROI == "LH_TRANSVERSETEMPORAL_SUVR" ~ "lh_transversetemporal",
                                       ROI == "LH_INSULA_SUVR" ~ "lh_insula",
                                       ROI == "LH_CAUDALANTERIORCINGULATE_SUVR" ~ "lh_caudalanteriorcingulate",
                                       ROI == "LH_CUNEUS_SUVR" ~ "lh_cuneus",
                                       ROI == "LH_ENTORHINAL_SUVR" ~ "lh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "LH_ISTHMUSCINGULATE_SUVR" ~ "lh_isthmuscingulate",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_LINGUAL_SUVR" ~ "lh_lingual",
                                       ROI == "LH_MEDIALORBITOFRONTAL_SUVR" ~ "lh_medialorbitofrontal",
                                       ROI == "LH_PARAHIPPOCAMPAL_SUVR" ~ "lh_parahippocampal",
                                       ROI == "LH_PARACENTRAL_SUVR" ~ "lh_paracentral",
                                       ROI == "LH_PERICALCARINE_SUVR" ~ "lh_pericalcarine",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_POSTERIORCINGULATE_SUVR" ~ "lh_posteriorcingulate",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "LH_PRECUNEUS_SUVR" ~ "lh_precuneus",
                                       ROI == "LH_ROSTRALANTERIORCINGULATE_SUVR" ~ "lh_rostralanteriorcingulate",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORPARIETAL_SUVR" ~ "rh_superiorparietal",
                                       ROI == "LH_FRONTALPOLE_SUVR" ~ "lh_frontalpole",
                                       ROI == "RH_TEMPORALPOLE_SUVR" ~ "rh_temporalpole",
                                       ROI == "RH_CAUDALANTERIORCINGULATE_SUVR" ~ "rh_caudalanteriorcingulate",
                                       ROI == "RH_CUNEUS_SUVR" ~ "rh_cuneus",
                                       ROI == "RH_ENTORHINAL_SUVR" ~ "rh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "RH_ISTHMUSCINGULATE_SUVR" ~ "rh_isthmuscingulate",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_LINGUAL_SUVR" ~ "rh_lingual",
                                       ROI == "RH_MEDIALORBITOFRONTAL_SUVR" ~ "rh_medialorbitofrontal",
                                       ROI == "RH_PARAHIPPOCAMPAL_SUVR" ~ "rh_parahippocampal",
                                       ROI == "RH_PARACENTRAL_SUVR" ~ "rh_paracentral",
                                       ROI == "RH_PERICALCARINE_SUVR" ~ "rh_pericalcarine",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_POSTERIORCINGULATE_SUVR" ~ "rh_posteriorcingulate",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "RH_PRECUNEUS_SUVR" ~ "rh_precuneus",
                                       ROI == "RH_ROSTRALANTERIORCINGULATE_SUVR" ~ "rh_rostralanteriorcingulate",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_FRONTALPOLE_SUVR" ~ "rh_frontalpole",
                                       ROI == "LH_TEMPORALPOLE_SUVR"  ~ "Lh_temporalpole",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "RH_CAUDALMIDDLEFRONTAL_SUVR" ~ "rh_caudalmiddlefrontal",
                                       ROI == "RH_BANKSSTS_SUVR" ~ "rh_bankssts",
                                       ROI == "RH_INFERIORPARIETAL_SUVR" ~ "rh_inferiorparietal",
                                       ROI == "RH_INFERIORTEMPORAL_SUVR" ~ "rh_inferiortemporal",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_MIDDLETEMPORAL_SUVR" ~ "rh_middletemporal",
                                       ROI == "RH_PARSOPERCULARIS_SUVR" ~ "rh_parsopercularis",
                                       ROI == "RH_PARSORBITALIS_SUVR" ~ "rh_parsorbitalis",
                                       ROI == "RH_PARSTRIANGULARIS_SUVR" ~ "rh_parstriangularis",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "RH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "rh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORTEMPORAL_SUVR" ~ "rh_superiortemporal",
                                       ROI == "RH_SUPRAMARGINAL_SUVR" ~ "rh_supramarginal",
                                       ROI == "RH_TEMPORALPOLE_SUVR"  ~ "rh_temporalpole",
                                       ROI == "RH_TRANSVERSETEMPORAL_SUVR" ~ "rh_transversetemporal",
                                       ROI == "RH_INSULA_SUVR" ~ "rh_insula"))

p_tau_aseg_cn <- merge(p_values_data_test_cn, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_tau_aparc_cn <- merge(p_values_data_test_cn, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

# MCI
p_values_data_test_mci <- merge(p_values_data_mci, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE) %>%
  dplyr::mutate(label_aseg = case_when(ROI == "LEFT_THALAMUS_PROPER_SUVR" ~ "Left-Thalamus-Proper",
                                       ROI == "RIGHT_THALAMUS_PROPER_SUVR" ~ "Right-Thalamus-Proper",
                                       ROI == "LEFT_LATERAL_VENTRICLE_SUVR" ~ "Left-Lateral-Ventricle",
                                       ROI == "RIGHT_LATERAL_VENTRICLE_SUVR" ~ "Right-Lateral-Ventricle",
                                       ROI == "LEFT_HIPPOCAMPUS_SUVR" ~ "Left-Hippocampus",
                                       ROI == "RIGHT_HIPPOCAMPUS_SUVR" ~ "Right-Hippocampus",
                                       ROI == "LEFT_PUTAMEN_SUVR" ~ "Left-Putamen",
                                       ROI == "RIGHT_PUTAMEN_SUVR" ~ "Right-Putamen", 
                                       ROI == "LEFT_AMYGDALA_SUVR" ~ "Left-Amygdala",
                                       ROI == "RIGHT_AMYGDALA_SUVR" ~ "Right-Amygdala",
                                       ROI == "LEFT_PALLIDUM_SUVR" ~ "Left-Pallidum",
                                       ROI == "RIGHT_PALLIDUM_SUVR" ~ "Right-Pallidum",
                                       ROI == "LEFT_CAUDATE_SUVR" ~ "Left-Caudate",
                                       ROI == "RIGHT_CAUDATE_SUVR" ~ "Right-Caudate",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "LH_CAUDALMIDDLEFRONTAL_SUVR" ~ "lh_caudalmiddlefrontal",
                                       ROI == "LH_FUSIFORM_SUVR" ~ "lh_fusiform",
                                       ROI == "LH_INFERIORPARIETAL_SUVR" ~ "lh_inferiorparietal",
                                       ROI == "LH_INFERIORTEMPORAL_SUVR" ~ "lh_inferiortemporal",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_MIDDLETEMPORAL_SUVR" ~ "lh_middletemporal",
                                       ROI == "LH_PARSOPERCULARIS_SUVR" ~ "lh_parsopercularis",
                                       ROI == "LH_PARSORBITALIS_SUVR" ~ "lh_parsorbitalis",
                                       ROI == "LH_PARSTRIANGULARIS_SUVR" ~ "lh_parstriangularis",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "LH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "lh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "LH_SUPERIORPARIETAL_SUVR" ~ "lh_superiorparietal",
                                       ROI == "LH_SUPERIORTEMPORAL_SUVR" ~ "lh_superiortemporal",
                                       ROI == "LH_SUPRAMARGINAL_SUVR" ~ "lh_supramarginal",
                                       ROI == "LH_TEMPORALPOLE_SUVR" ~ "lh_temporalpole",
                                       ROI == "LH_TRANSVERSETEMPORAL_SUVR" ~ "lh_transversetemporal",
                                       ROI == "LH_INSULA_SUVR" ~ "lh_insula",
                                       ROI == "LH_CAUDALANTERIORCINGULATE_SUVR" ~ "lh_caudalanteriorcingulate",
                                       ROI == "LH_CUNEUS_SUVR" ~ "lh_cuneus",
                                       ROI == "LH_ENTORHINAL_SUVR" ~ "lh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "LH_ISTHMUSCINGULATE_SUVR" ~ "lh_isthmuscingulate",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_LINGUAL_SUVR" ~ "lh_lingual",
                                       ROI == "LH_MEDIALORBITOFRONTAL_SUVR" ~ "lh_medialorbitofrontal",
                                       ROI == "LH_PARAHIPPOCAMPAL_SUVR" ~ "lh_parahippocampal",
                                       ROI == "LH_PARACENTRAL_SUVR" ~ "lh_paracentral",
                                       ROI == "LH_PERICALCARINE_SUVR" ~ "lh_pericalcarine",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_POSTERIORCINGULATE_SUVR" ~ "lh_posteriorcingulate",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "LH_PRECUNEUS_SUVR" ~ "lh_precuneus",
                                       ROI == "LH_ROSTRALANTERIORCINGULATE_SUVR" ~ "lh_rostralanteriorcingulate",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORPARIETAL_SUVR" ~ "rh_superiorparietal",
                                       ROI == "LH_FRONTALPOLE_SUVR" ~ "lh_frontalpole",
                                       ROI == "RH_TEMPORALPOLE_SUVR" ~ "rh_temporalpole",
                                       ROI == "RH_CAUDALANTERIORCINGULATE_SUVR" ~ "rh_caudalanteriorcingulate",
                                       ROI == "RH_CUNEUS_SUVR" ~ "rh_cuneus",
                                       ROI == "RH_ENTORHINAL_SUVR" ~ "rh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "RH_ISTHMUSCINGULATE_SUVR" ~ "rh_isthmuscingulate",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_LINGUAL_SUVR" ~ "rh_lingual",
                                       ROI == "RH_MEDIALORBITOFRONTAL_SUVR" ~ "rh_medialorbitofrontal",
                                       ROI == "RH_PARAHIPPOCAMPAL_SUVR" ~ "rh_parahippocampal",
                                       ROI == "RH_PARACENTRAL_SUVR" ~ "rh_paracentral",
                                       ROI == "RH_PERICALCARINE_SUVR" ~ "rh_pericalcarine",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_POSTERIORCINGULATE_SUVR" ~ "rh_posteriorcingulate",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "RH_PRECUNEUS_SUVR" ~ "rh_precuneus",
                                       ROI == "RH_ROSTRALANTERIORCINGULATE_SUVR" ~ "rh_rostralanteriorcingulate",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_FRONTALPOLE_SUVR" ~ "rh_frontalpole",
                                       ROI == "LH_TEMPORALPOLE_SUVR"  ~ "Lh_temporalpole",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "RH_CAUDALMIDDLEFRONTAL_SUVR" ~ "rh_caudalmiddlefrontal",
                                       ROI == "RH_BANKSSTS_SUVR" ~ "rh_bankssts",
                                       ROI == "RH_INFERIORPARIETAL_SUVR" ~ "rh_inferiorparietal",
                                       ROI == "RH_INFERIORTEMPORAL_SUVR" ~ "rh_inferiortemporal",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_MIDDLETEMPORAL_SUVR" ~ "rh_middletemporal",
                                       ROI == "RH_PARSOPERCULARIS_SUVR" ~ "rh_parsopercularis",
                                       ROI == "RH_PARSORBITALIS_SUVR" ~ "rh_parsorbitalis",
                                       ROI == "RH_PARSTRIANGULARIS_SUVR" ~ "rh_parstriangularis",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "RH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "rh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORTEMPORAL_SUVR" ~ "rh_superiortemporal",
                                       ROI == "RH_SUPRAMARGINAL_SUVR" ~ "rh_supramarginal",
                                       ROI == "RH_TEMPORALPOLE_SUVR"  ~ "rh_temporalpole",
                                       ROI == "RH_TRANSVERSETEMPORAL_SUVR" ~ "rh_transversetemporal",
                                       ROI == "RH_INSULA_SUVR" ~ "rh_insula"))

p_tau_aseg_mci <- merge(p_values_data_test_mci, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_tau_aparc_mci <- merge(p_values_data_test_mci, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

# AD
p_values_data_test_ad <- merge(p_values_data_ad, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE) %>%
  dplyr::mutate(label_aseg = case_when(ROI == "LEFT_THALAMUS_PROPER_SUVR" ~ "Left-Thalamus-Proper",
                                       ROI == "RIGHT_THALAMUS_PROPER_SUVR" ~ "Right-Thalamus-Proper",
                                       ROI == "LEFT_LATERAL_VENTRICLE_SUVR" ~ "Left-Lateral-Ventricle",
                                       ROI == "RIGHT_LATERAL_VENTRICLE_SUVR" ~ "Right-Lateral-Ventricle",
                                       ROI == "LEFT_HIPPOCAMPUS_SUVR" ~ "Left-Hippocampus",
                                       ROI == "RIGHT_HIPPOCAMPUS_SUVR" ~ "Right-Hippocampus",
                                       ROI == "LEFT_PUTAMEN_SUVR" ~ "Left-Putamen",
                                       ROI == "RIGHT_PUTAMEN_SUVR" ~ "Right-Putamen", 
                                       ROI == "LEFT_AMYGDALA_SUVR" ~ "Left-Amygdala",
                                       ROI == "RIGHT_AMYGDALA_SUVR" ~ "Right-Amygdala",
                                       ROI == "LEFT_PALLIDUM_SUVR" ~ "Left-Pallidum",
                                       ROI == "RIGHT_PALLIDUM_SUVR" ~ "Right-Pallidum",
                                       ROI == "LEFT_CAUDATE_SUVR" ~ "Left-Caudate",
                                       ROI == "RIGHT_CAUDATE_SUVR" ~ "Right-Caudate",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "LH_CAUDALMIDDLEFRONTAL_SUVR" ~ "lh_caudalmiddlefrontal",
                                       ROI == "LH_FUSIFORM_SUVR" ~ "lh_fusiform",
                                       ROI == "LH_INFERIORPARIETAL_SUVR" ~ "lh_inferiorparietal",
                                       ROI == "LH_INFERIORTEMPORAL_SUVR" ~ "lh_inferiortemporal",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_MIDDLETEMPORAL_SUVR" ~ "lh_middletemporal",
                                       ROI == "LH_PARSOPERCULARIS_SUVR" ~ "lh_parsopercularis",
                                       ROI == "LH_PARSORBITALIS_SUVR" ~ "lh_parsorbitalis",
                                       ROI == "LH_PARSTRIANGULARIS_SUVR" ~ "lh_parstriangularis",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "LH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "lh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "LH_SUPERIORPARIETAL_SUVR" ~ "lh_superiorparietal",
                                       ROI == "LH_SUPERIORTEMPORAL_SUVR" ~ "lh_superiortemporal",
                                       ROI == "LH_SUPRAMARGINAL_SUVR" ~ "lh_supramarginal",
                                       ROI == "LH_TEMPORALPOLE_SUVR" ~ "lh_temporalpole",
                                       ROI == "LH_TRANSVERSETEMPORAL_SUVR" ~ "lh_transversetemporal",
                                       ROI == "LH_INSULA_SUVR" ~ "lh_insula",
                                       ROI == "LH_CAUDALANTERIORCINGULATE_SUVR" ~ "lh_caudalanteriorcingulate",
                                       ROI == "LH_CUNEUS_SUVR" ~ "lh_cuneus",
                                       ROI == "LH_ENTORHINAL_SUVR" ~ "lh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "LH_ISTHMUSCINGULATE_SUVR" ~ "lh_isthmuscingulate",
                                       ROI == "LH_LATERALOCCIPITAL_SUVR" ~ "lh_lateraloccipital",
                                       ROI == "LH_LATERALORBITOFRONTAL_SUVR" ~ "lh_lateralorbitofrontal",
                                       ROI == "LH_LINGUAL_SUVR" ~ "lh_lingual",
                                       ROI == "LH_MEDIALORBITOFRONTAL_SUVR" ~ "lh_medialorbitofrontal",
                                       ROI == "LH_PARAHIPPOCAMPAL_SUVR" ~ "lh_parahippocampal",
                                       ROI == "LH_PARACENTRAL_SUVR" ~ "lh_paracentral",
                                       ROI == "LH_PERICALCARINE_SUVR" ~ "lh_pericalcarine",
                                       ROI == "LH_POSTCENTRAL_SUVR" ~ "lh_postcentral",
                                       ROI == "LH_POSTERIORCINGULATE_SUVR" ~ "lh_posteriorcingulate",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "LH_PRECUNEUS_SUVR" ~ "lh_precuneus",
                                       ROI == "LH_ROSTRALANTERIORCINGULATE_SUVR" ~ "lh_rostralanteriorcingulate",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORPARIETAL_SUVR" ~ "rh_superiorparietal",
                                       ROI == "LH_FRONTALPOLE_SUVR" ~ "lh_frontalpole",
                                       ROI == "RH_TEMPORALPOLE_SUVR" ~ "rh_temporalpole",
                                       ROI == "RH_CAUDALANTERIORCINGULATE_SUVR" ~ "rh_caudalanteriorcingulate",
                                       ROI == "RH_CUNEUS_SUVR" ~ "rh_cuneus",
                                       ROI == "RH_ENTORHINAL_SUVR" ~ "rh_entorhinal",
                                       ROI == "RH_FUSIFORM_SUVR" ~ "rh_fusiform",
                                       ROI == "RH_ISTHMUSCINGULATE_SUVR" ~ "rh_isthmuscingulate",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_LINGUAL_SUVR" ~ "rh_lingual",
                                       ROI == "RH_MEDIALORBITOFRONTAL_SUVR" ~ "rh_medialorbitofrontal",
                                       ROI == "RH_PARAHIPPOCAMPAL_SUVR" ~ "rh_parahippocampal",
                                       ROI == "RH_PARACENTRAL_SUVR" ~ "rh_paracentral",
                                       ROI == "RH_PERICALCARINE_SUVR" ~ "rh_pericalcarine",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_POSTERIORCINGULATE_SUVR" ~ "rh_posteriorcingulate",
                                       ROI == "LH_PRECENTRAL_SUVR" ~ "lh_precentral",
                                       ROI == "RH_PRECUNEUS_SUVR" ~ "rh_precuneus",
                                       ROI == "RH_ROSTRALANTERIORCINGULATE_SUVR" ~ "rh_rostralanteriorcingulate",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_FRONTALPOLE_SUVR" ~ "rh_frontalpole",
                                       ROI == "LH_TEMPORALPOLE_SUVR"  ~ "Lh_temporalpole",
                                       ROI == "LH_BANKSSTS_SUVR" ~ "lh_bankssts",
                                       ROI == "RH_CAUDALMIDDLEFRONTAL_SUVR" ~ "rh_caudalmiddlefrontal",
                                       ROI == "RH_BANKSSTS_SUVR" ~ "rh_bankssts",
                                       ROI == "RH_INFERIORPARIETAL_SUVR" ~ "rh_inferiorparietal",
                                       ROI == "RH_INFERIORTEMPORAL_SUVR" ~ "rh_inferiortemporal",
                                       ROI == "RH_LATERALOCCIPITAL_SUVR" ~ "rh_lateraloccipital",
                                       ROI == "RH_LATERALORBITOFRONTAL_SUVR" ~ "rh_lateralorbitofrontal",
                                       ROI == "RH_MIDDLETEMPORAL_SUVR" ~ "rh_middletemporal",
                                       ROI == "RH_PARSOPERCULARIS_SUVR" ~ "rh_parsopercularis",
                                       ROI == "RH_PARSORBITALIS_SUVR" ~ "rh_parsorbitalis",
                                       ROI == "RH_PARSTRIANGULARIS_SUVR" ~ "rh_parstriangularis",
                                       ROI == "RH_POSTCENTRAL_SUVR" ~ "rh_postcentral",
                                       ROI == "RH_PRECENTRAL_SUVR" ~ "rh_precentral",
                                       ROI == "RH_ROSTRALMIDDLEFRONTAL_SUVR" ~ "rh_rostralmiddlefrontal",
                                       ROI == "LH_SUPERIORFRONTAL_SUVR" ~ "lh_superiorfrontal",
                                       ROI == "RH_SUPERIORFRONTAL_SUVR" ~ "rh_superiorfrontal",
                                       ROI == "RH_SUPERIORTEMPORAL_SUVR" ~ "rh_superiortemporal",
                                       ROI == "RH_SUPRAMARGINAL_SUVR" ~ "rh_supramarginal",
                                       ROI == "RH_TEMPORALPOLE_SUVR"  ~ "rh_temporalpole",
                                       ROI == "RH_TRANSVERSETEMPORAL_SUVR" ~ "rh_transversetemporal",
                                       ROI == "RH_INSULA_SUVR" ~ "rh_insula"))

p_tau_aseg_ad <- merge(p_values_data_test_ad, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_tau_aparc_ad <- merge(p_values_data_test_ad, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

####################################
# Step 7: Determining which variables to plot by significance
####################################
#ASEG rois
p_tau_aseg <- p_tau_aseg %>%
  dplyr::mutate(p_sig = case_when(tau_p_values >= .05 ~ "not sig",
                                  tau_p_values < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_tau_aseg_cn <- p_tau_aseg_cn %>%
  dplyr::mutate(p_sig = case_when(tau_p_values_cn >= .05 ~ "not sig",
                                  tau_p_values_cn < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_tau_aseg_mci <- p_tau_aseg_mci %>%
  dplyr::mutate(p_sig = case_when(tau_p_values_mci >= .05 ~ "not sig",
                                  tau_p_values_mci < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_tau_aseg_ad <- p_tau_aseg_ad %>%
  dplyr::mutate(p_sig = case_when(tau_p_values_ad >= .05 ~ "not sig",
                                  tau_p_values_ad < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))

#APARC rois
p_tau_aparc <- p_tau_aparc %>%
  dplyr::mutate(p_sig = case_when(tau_p_values >= .05 ~ "not sig",
                                  tau_p_values < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_tau_aparc_cn <- p_tau_aparc_cn %>%
  dplyr::mutate(p_sig = case_when(tau_p_values_cn >= .05 ~ "not sig",
                                  tau_p_values_cn < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_tau_aparc_mci <- p_tau_aparc_mci %>%
  dplyr::mutate(p_sig = case_when(tau_p_values_mci >= .05 ~ "not sig",
                                  tau_p_values_mci < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_tau_aparc_ad <- p_tau_aparc_ad %>%
  dplyr::mutate(p_sig = case_when(tau_p_values_ad >= .05 ~ "not sig",
                                  tau_p_values_ad < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))

####################################
# Step 8: plotting
####################################
#aseg
b <- c(-.7, -.35, 0, .35, .7) # b <- c(-1.25, -.625,  0, .625, 1.25)

aseg_full <- ggplot(p_tau_aseg %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) +
  scale_fill_gradientn(limits = c(-.7,.7), #  scale_fill_gradientn(limits = c(-1.3,1.3),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("All Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aseg_cn <- ggplot(p_tau_aseg_cn %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) + 
  scale_fill_gradientn(limits = c(-.7,.7),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("CU Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aseg_mci <- ggplot(p_tau_aseg_mci %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size))  + 
  scale_fill_gradientn(limits = c(-.7,.7),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("MCI Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aseg_ad <- ggplot(p_tau_aseg_ad %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) + 
  scale_fill_gradientn(limits = c(-.7,.7),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("AD Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()

aseg_full + aseg_cn + aseg_mci + aseg_ad + 
  plot_layout(ncol = 1, guides = "collect")

#aparc
b2 <- c(-.5, -.25, 0, .25, .5) # b2 <- c(-1, -.5, 0, .5, 1)


aparc_full <- ggplot(p_tau_aparc %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.5, .5), #  scale_fill_gradientn(limits = c(-1.03,1.03),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("All Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aparc_cn <- ggplot(p_tau_aparc_cn %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.5, .5),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("CU Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aparc_mci <- ggplot(p_tau_aparc_mci %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size))  + # significant rois include: CC_Anterior, Left-Amygdala, Left-Hippocampus, Right-Hippocampus
  scale_fill_gradientn(limits = c(-.5, .5),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("MCI Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aparc_ad <- ggplot(p_tau_aparc_ad %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #no significant rois
  scale_fill_gradientn(limits = c(-.5, .5),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("AD Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()

aparc_full + aparc_cn + aparc_mci + aparc_ad + 
  plot_layout(ncol = 1, guides = "collect")

####################################
# Step 9: p_value adjustment
####################################

p_tau_all <- rbind(p_tau_aseg %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values, Effect_Size, p_sig), p_tau_aparc %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)
p_tau_all_cn <- rbind(p_tau_aseg_cn %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values_cn, Effect_Size, p_sig), p_tau_aparc_cn %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values_cn, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)
p_tau_all_mci <- rbind(p_tau_aseg_mci %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values_mci, Effect_Size, p_sig), p_tau_aparc_mci %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values_mci, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)
p_tau_all_ad <- rbind(p_tau_aseg_ad %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values_ad, Effect_Size, p_sig), p_tau_aparc_ad %>% dplyr::select(ROI, atlas, hemi, side, region, tau_p_values_ad, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)

p_tau_all$adjusted_p <- p.adjust(p_tau_all$tau_p_values, method = "holm")
p_tau_all_cn$adjusted_p <- p.adjust(p_tau_all_cn$tau_p_values, method = "holm")
p_tau_all_mci$adjusted_p <- p.adjust(p_tau_all_mci$tau_p_values, method = "holm")
p_tau_all_ad$adjusted_p <- p.adjust(p_tau_all_ad$tau_p_values, method = "holm")
