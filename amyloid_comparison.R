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
# Step 1: Adding in Amyloid positivity (filling everything after 1st amyloid positive - positive, filling everything before last amyloid negative - negative)
####################################
# Step 2: getting MRI right before EXAMDATE (SAA-) and right after EXAMDATE (SAA+)
####################################
# Step 3: adding nearest diagnoses
####################################
# Step 4: getting demographic information
####################################
# Step 5: getting p-values
####################################
# Step 6: ASEG - organizing dictionary and variable names
####################################
# Step 7: ASEG - determining which variables to plot by significance
####################################
# Step 8: APARC - organizing dictionary and variable names
####################################
# Step 9: APARC - determining which variables to plot by significance
####################################
# Step 10: Plotting
####################################
# Step 11: p-value corrections
####################################


####################################
# Labelled Steps:
####################################
# Step 1: Adding in Amyloid positivity (filling everything after 1st amyloid positive - positive, filling everything before last amyloid negative - negative)
####################################
# Step 2: getting amyloid right before EXAMDATE (SAA-) and right after EXAMDATE (SAA+)
####################################
# Step 3: adding nearest diagnoses
####################################
# Step 4: getting demographic information
####################################
# Step 5: getting p-values
####################################
# Step 6: organizing dictionary and variable names
####################################
# Step 7: determining which variables to plot by significance
####################################
# Step 8: plotting
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

rm(apoeres)












####################################
# Labelled Steps:
####################################
# Step 1: Adding in Amyloid positivity (filling everything after 1st amyloid positive - positive, filling everything before last amyloid negative - negative)
####################################
# Step 2: getting amyloid right before EXAMDATE (SAA-) and right after EXAMDATE (SAA+)
####################################
# Step 3: adding nearest diagnoses
####################################
# Step 4: getting demographic information
####################################
# Step 5: getting p-values
####################################
# Step 6: organizing dictionary and variable names
####################################
# Step 7: determining which variables to plot by significance
####################################
# Step 8: plotting
####################################
# Step 9: p_value adjustment
####################################



#############################################################################################################
#############################################################################################################
#############################################################################################################
#                                           Regional Amyloid
#############################################################################################################
#############################################################################################################
#############################################################################################################

#regional pet values for amyloid
amyloid_roi_cs <- read.csv("C:\\Work Folder\\paper data longitudinal phases\\UCBERKELEY_AMY_6MM_31Jul2024.csv") %>% 
  dplyr::filter(qc_flag == 2 | qc_flag == 1,
                RID %in% stable_rids) %>%
  dplyr::mutate(SCANDATE = as.Date(SCANDATE))

#adding in the stable cases
amyloid_roi_cs <- merge(amyloid_roi_cs, stables, by = "RID")

####################################
# Step 2: getting SCAN right before EXAMDATE (SAA-) and right after EXAMDATE (SAA+)
####################################

amyloid_roi_cs_saa_neg <- amyloid_roi_cs %>%
  dplyr::filter(group == "SAA- Stable")
amyloid_roi_cs_saa_neg <- amyloid_roi_cs_saa_neg %>%
  dplyr::mutate(RID = as.character(RID))
amyloid_roi_cs_saa_pos <- amyloid_roi_cs %>%
  dplyr::filter(group == "SAA+ Stable")
amyloid_roi_cs_saa_pos <- amyloid_roi_cs_saa_pos %>%
  dplyr::mutate(RID = as.character(RID))

#now pulling the EXAMDATE right before the last SAA negative date
amyloid_roi_cs_saa_neg <- amyloid_roi_cs_saa_neg %>%
  dplyr::left_join(dem) %>%
  dplyr::mutate(age = lubridate::time_length(difftime(as.Date(SCANDATE), as.Date(birthdate)), "years"))
amyloid_roi_cs_saa_neg <- amyloid_roi_cs_saa_neg %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(SCANDATE), last_saa_negative_date), "years")) %>%
  dplyr::select(RID, SCANDATE, age, PTGENDER, PTEDUCAT, apoe, time_diff, estaget0_ABETA42, group, first_saa_positive_date, last_saa_negative_date, ends_with("SUVR"))
amyloid_roi_cs_saa_neg <- amyloid_roi_cs_saa_neg %>%
  dplyr::filter(time_diff <= 0)
amyloid_roi_cs_saa_neg <- amyloid_roi_cs_saa_neg %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_diff == max(time_diff)) %>%
  dplyr::ungroup()
amyloid_roi_cs_saa_neg <- amyloid_roi_cs_saa_neg %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::arrange(SCANDATE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

#now pulling the EXAMDATE right after the first SAA positive date
amyloid_roi_cs_saa_pos <- amyloid_roi_cs_saa_pos %>%
  dplyr::left_join(dem) %>%
  dplyr::mutate(age = lubridate::time_length(difftime(as.Date(SCANDATE), as.Date(birthdate)), "years"))
amyloid_roi_cs_saa_pos <- amyloid_roi_cs_saa_pos %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(SCANDATE), first_saa_positive_date), "years")) %>%
  dplyr::select(RID, SCANDATE, age, PTGENDER, PTEDUCAT, apoe, time_diff, estaget0_ABETA42, group, first_saa_positive_date, last_saa_negative_date, ends_with("SUVR"))
amyloid_roi_cs_saa_pos <- amyloid_roi_cs_saa_pos %>%
  dplyr::filter(time_diff >= 0)
amyloid_roi_cs_saa_pos <- amyloid_roi_cs_saa_pos %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
amyloid_roi_cs_saa_pos <- amyloid_roi_cs_saa_pos %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::arrange(SCANDATE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

#putting the two datasets with amyloid status and correct visit dates together
amyloid_roi_cs <- rbind(amyloid_roi_cs_saa_neg, amyloid_roi_cs_saa_pos) #604 RID's

####################################
# Step 3: adding nearest diagnoses
####################################
#adding diagnosis in
amyloid_roi_cs <- amyloid_roi_cs %>%
  dplyr::left_join(adni_diagnoses %>% 
                     dplyr::mutate(RID = as.character(RID))) %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(SCANDATE), DX.DATE), "years")) %>%
  dplyr::distinct() %>%
  dplyr::select(RID, SCANDATE, age, PTGENDER, PTEDUCAT, apoe, time_diff, estaget0_ABETA42, group, first_saa_positive_date, last_saa_negative_date, DX.DATE, DX, ends_with("SUVR"))
amyloid_roi_cs <- amyloid_roi_cs %>%
  dplyr::filter(!is.na(time_diff)) %>%
  dplyr::mutate(time_diff = abs(time_diff))
amyloid_roi_cs <- as.data.frame(amyloid_roi_cs)
amyloid_roi_cs <- amyloid_roi_cs %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
amyloid_roi_cs <- amyloid_roi_cs %>%
  dplyr::group_by(RID, SCANDATE) %>%
  dplyr::arrange(SCANDATE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

####################################
# Step 4: getting demographic information
####################################
#first getting rid of rows that have no information on features
amyloid_roi_cs <- amyloid_roi_cs[which(rowMeans(!is.na(amyloid_roi_cs)) > 0.8), ]

# amyloid_roi_cs$RID[duplicated(amyloid_roi_cs$RID)]

#getting distributions of diagnoses, gender, apoe, amyloid status
table(amyloid_roi_cs$group, amyloid_roi_cs$DX)
table(amyloid_roi_cs$group, amyloid_roi_cs$PTGENDER)
table(amyloid_roi_cs$group, amyloid_roi_cs$apoe)
table(amyloid_roi_cs$group)

#getting averages of age and education
mean(amyloid_roi_cs_saa_pos$age)
sd(amyloid_roi_cs_saa_pos$age)
mean(amyloid_roi_cs_saa_neg$age)
sd(amyloid_roi_cs_saa_neg$age)
mean(amyloid_roi_cs_saa_pos$PTEDUCAT)
sd(amyloid_roi_cs_saa_pos$PTEDUCAT)
mean(amyloid_roi_cs_saa_neg$PTEDUCAT)
sd(amyloid_roi_cs_saa_neg$PTEDUCAT)




rm(amyloid_roi_cs_saa_neg, amyloid_roi_cs_saa_pos)


####################################
# Step 5: getting p-values
####################################
amyloid_roi_cs <- amyloid_roi_cs %>%
  dplyr::select(-VENTRICLE_5TH_SUVR, -WHOLECEREBELLUM_SUVR)
features_amyloid <- amyloid_roi_cs %>%
  dplyr::select(ends_with("SUVR"))
featurenames_amyloid <- names(features_amyloid)

aseg_data <- as.data.frame(aseg)
aseg_data$label
aparc_data <- as.data.frame(dk)
aparc_data$label

#stratifying data
#pairwise comparisons of each region and saving to a list then data frame
amyloid_roi_cs_cn <- amyloid_roi_cs %>%
  dplyr::filter(DX == "CU")
amyloid_roi_cs_mci <- amyloid_roi_cs %>%
  dplyr::filter(DX == "MCI")
amyloid_roi_cs_ad <- amyloid_roi_cs %>%
  dplyr::filter(DX == "Dementia")

#full data
amyloid_p_values <- vector()
amyloid_roi_names <- vector()
amyloid_effect_sizes <- vector()
amyloid_n_all <- vector()

for (f in 1:length(featurenames_amyloid)){
  feature <- featurenames_amyloid[f]
  print(paste("amyloid data regression: ", feature))
  lm_feature <- lm(scale(amyloid_roi_cs[, 13+f]) ~ group + scale(age) + as.factor(PTGENDER) + scale(PTEDUCAT) + apoe + DX, data = amyloid_roi_cs, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  amyloid_effect_sizes <- append(amyloid_effect_sizes, estimate)
  amyloid_p_values <- append(amyloid_p_values, group_p_value)
  amyloid_roi_names <- append(amyloid_roi_names, feature)
  amyloid_n_all <- append(amyloid_n_all, length(lm_feature$residuals))
}

p_values_data <- data.frame(amyloid_p_values)
amyloid_roi_names <- gsub("CTX_", "", amyloid_roi_names) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data$ROI <- amyloid_roi_names #creating column for the roi names for each p value
p_values_data$Effect_Size <- amyloid_effect_sizes
p_values_data$n <- amyloid_n_all

#CN
amyloid_p_values_cn <- vector()
amyloid_roi_names_cn <- vector()
amyloid_effect_sizes_cn <- vector()
amyloid_n_cn <- vector()

for (f in 1:length(featurenames_amyloid)){
  feature <- featurenames_amyloid[f]
  print(paste("amyloid data regression: ", feature))
  lm_feature <- lm(scale(amyloid_roi_cs_cn[, 13+f]) ~ group + scale(age) + as.factor(PTGENDER) + scale(PTEDUCAT) + apoe, data = amyloid_roi_cs_cn, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  amyloid_effect_sizes_cn <- append(amyloid_effect_sizes_cn, estimate)
  amyloid_p_values_cn <- append(amyloid_p_values_cn, group_p_value)
  amyloid_roi_names_cn <- append(amyloid_roi_names_cn, feature)
  amyloid_n_cn <- append(amyloid_n_cn, length(lm_feature$residuals))
}

p_values_data_cn <- data.frame(amyloid_p_values_cn)
amyloid_roi_names_cn <- gsub("CTX_", "", amyloid_roi_names_cn) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data_cn$ROI <- amyloid_roi_names_cn #creating column for the roi names for each p value
p_values_data_cn$Effect_Size <- amyloid_effect_sizes_cn
p_values_data_cn$n <- amyloid_n_cn

#MCI
amyloid_p_values_mci <- vector()
amyloid_roi_names_mci <- vector()
amyloid_effect_sizes_mci <- vector()
amyloid_n_mci <- vector()

for (f in 1:length(featurenames_amyloid)){
  feature <- featurenames_amyloid[f]
  print(paste("amyloid data regression: ", feature))
  lm_feature <- lm(scale(amyloid_roi_cs_mci[, 13+f]) ~ group + scale(age) + as.factor(PTGENDER) + scale(PTEDUCAT) + apoe, data = amyloid_roi_cs_mci, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  amyloid_effect_sizes_mci <- append(amyloid_effect_sizes_mci, estimate)
  amyloid_p_values_mci <- append(amyloid_p_values_mci, group_p_value)
  amyloid_roi_names_mci <- append(amyloid_roi_names_mci, feature)
  amyloid_n_mci <- append(amyloid_n_mci, length(lm_feature$residuals))
}

p_values_data_mci <- data.frame(amyloid_p_values_mci)
amyloid_roi_names_mci <- gsub("CTX_", "", amyloid_roi_names_mci) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data_mci$ROI <- amyloid_roi_names_mci #creating column for the roi names for each p value
p_values_data_mci$Effect_Size <- amyloid_effect_sizes_mci
p_values_data_mci$n <- amyloid_n_mci

#AD
amyloid_p_values_ad <- vector()
amyloid_roi_names_ad <- vector()
amyloid_effect_sizes_ad <- vector()
amyloid_n_ad <- vector()

for (f in 1:length(featurenames_amyloid)){
  feature <- featurenames_amyloid[f]
  print(paste("amyloid data regression: ", feature))
  lm_feature <- lm(scale(amyloid_roi_cs_ad[, 13+f]) ~ group + scale(age) + as.factor(PTGENDER) + scale(PTEDUCAT) + apoe, data = amyloid_roi_cs_ad, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  amyloid_effect_sizes_ad <- append(amyloid_effect_sizes_ad, estimate)
  amyloid_p_values_ad <- append(amyloid_p_values_ad, group_p_value)
  amyloid_roi_names_ad <- append(amyloid_roi_names_ad, feature)
  amyloid_n_ad <- append(amyloid_n_ad, length(lm_feature$residuals))
}

p_values_data_ad <- data.frame(amyloid_p_values_ad)
amyloid_roi_names_ad <- gsub("CTX_", "", amyloid_roi_names_ad) #getting rid of the CTX part so that it lines up with the aseg dictionaries
p_values_data_ad$ROI <- amyloid_roi_names_ad #creating column for the roi names for each p value
p_values_data_ad$Effect_Size <- amyloid_effect_sizes_ad
p_values_data_ad$n <- amyloid_n_ad

####################################
# Step 6: organizing dictionary and variable names
####################################

#getting the dictionary so that the ROI names are listed
dict <- as.data.frame(featurenames_amyloid)
dict <- dict %>%
  dplyr::rename(FLDNAME = featurenames_amyloid)

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
                                  # ROI == "LEFT_VENTRALDC_SUVR" ~ "Left-VentralDC",
                                  # ROI == "RIGHT_VENTRALDC_SUVR" ~ "Right-VentralDC",
                                  # ROI == "VENTRICLE_3RD_SUVR" ~ "x3rd-ventricle",
                                  # ROI == "VENTRICLE_4TH_SUVR" ~ "x4th-ventricle",
                                  # ROI == "BRAINSTEM_SUVR" ~ "Brain-Stem",
                                  # ROI == "CC_ANTERIOR_SUVR" ~ "CC_Anterior",
                                  # ROI == "CC_CENTRAL_SUVR" ~ "CC_Central",
                                  # ROI == "CC_MID_ANTERIOR_SUVR" ~ "CC_Mid_Anterior",
                                  # ROI == "CC_MID_POSTERIOR_SUVR" ~ "CC_Mid_Posterior",
                                  # ROI == "CC_POSTERIOR_SUVR" ~ "CC_Posterior",
                                  # ROI == "LEFT_CEREBELLUM_CORTEX_SUVR" ~ "Left-Cerebellum-Cortex",
                                  # ROI == "RIGHT_CEREBELLUM_CORTEX_SUVR" ~ "Right-Cerebellum-Cortex",
                                  # ROI == "LEFT_CEREBELLUM_WHITE_MATTER_SUVR" ~ "Left-Cerebellum-White-Matter",
                                  # ROI == "RIGHT_CEREBELLUM_WHITE_MATTER_SUVR" ~ "Right-Cerebellum-White-Matter",
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
                                   # ~ "lh_corpuscallosum",
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
                                  # ~ "rh_corpuscallosum",
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

p_amyloid_aseg <- merge(p_values_data_test, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_amyloid_aparc <- merge(p_values_data_test, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
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

p_amyloid_aseg_cn <- merge(p_values_data_test_cn, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_amyloid_aparc_cn <- merge(p_values_data_test_cn, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
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

p_amyloid_aseg_mci <- merge(p_values_data_test_mci, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_amyloid_aparc_mci <- merge(p_values_data_test_mci, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
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

p_amyloid_aseg_ad <- merge(p_values_data_test_ad, aseg_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

p_amyloid_aparc_ad <- merge(p_values_data_test_ad, aparc_data, by.x = "label_aseg", by.y = "label", all.y = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

####################################
# Step 7: Determining which variables to plot by significance
####################################
#ASEG rois
p_amyloid_aseg <- p_amyloid_aseg %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values >= .05 ~ "not sig",
                                  amyloid_p_values < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_amyloid_aseg_cn <- p_amyloid_aseg_cn %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values_cn >= .05 ~ "not sig",
                                  amyloid_p_values_cn < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_amyloid_aseg_mci <- p_amyloid_aseg_mci %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values_mci >= .05 ~ "not sig",
                                  amyloid_p_values_mci < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_amyloid_aseg_ad <- p_amyloid_aseg_ad %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values_ad >= .05 ~ "not sig",
                                  amyloid_p_values_ad < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))

#APARC rois
p_amyloid_aparc <- p_amyloid_aparc %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values >= .05 ~ "not sig",
                                  amyloid_p_values < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_amyloid_aparc_cn <- p_amyloid_aparc_cn %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values_cn >= .05 ~ "not sig",
                                  amyloid_p_values_cn < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_amyloid_aparc_mci <- p_amyloid_aparc_mci %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values_mci >= .05 ~ "not sig",
                                  amyloid_p_values_mci < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_amyloid_aparc_ad <- p_amyloid_aparc_ad %>%
  dplyr::mutate(p_sig = case_when(amyloid_p_values_ad >= .05 ~ "not sig",
                                  amyloid_p_values_ad < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))


####################################
# Step 8: plotting
####################################
#aseg
b <- c(-.3, 0, .3) # b <- c(-.4, 0, .4)


aseg_full <- ggplot(p_amyloid_aseg %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) +
  scale_fill_gradientn(limits = c(-.34,.34),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("All Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aseg_cn <- ggplot(p_amyloid_aseg_cn %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) + 
  scale_fill_gradientn(limits = c(-.34,.34),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("CU Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aseg_mci <- ggplot(p_amyloid_aseg_mci %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size))  + 
  scale_fill_gradientn(limits = c(-.34,.34),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("MCI Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aseg_ad <- ggplot(p_amyloid_aseg_ad %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) + 
  scale_fill_gradientn(limits = c(-.34,.34),
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
b2 <- c(-.3, 0, .3) # b2 <- c(-.8, -.4, 0, .4, .8)

aparc_full <- ggplot(p_amyloid_aparc %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.3,.3), #   scale_fill_gradientn(limits = c(-.81,.81),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("All Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aparc_cn <- ggplot(p_amyloid_aparc_cn %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.3,.3),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("CU Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aparc_mci <- ggplot(p_amyloid_aparc_mci %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size))  + # significant rois include: CC_Anterior, Left-Amygdala, Left-Hippocampus, Right-Hippocampus
  scale_fill_gradientn(limits = c(-.3,.3),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("MCI Data") + 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) + 
  theme_void()
aparc_ad <- ggplot(p_amyloid_aparc_ad %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #no significant rois
  scale_fill_gradientn(limits = c(-.3,.3),
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

p_amyloid_all <- rbind(p_amyloid_aseg %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values, Effect_Size, p_sig), p_amyloid_aparc %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)
p_amyloid_all_cn <- rbind(p_amyloid_aseg_cn %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values_cn, Effect_Size, p_sig), p_amyloid_aparc_cn %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values_cn, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)
p_amyloid_all_mci <- rbind(p_amyloid_aseg_mci %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values_mci, Effect_Size, p_sig), p_amyloid_aparc_mci %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values_mci, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)
p_amyloid_all_ad <- rbind(p_amyloid_aseg_ad %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values_ad, Effect_Size, p_sig), p_amyloid_aparc_ad %>% dplyr::select(ROI, atlas, hemi, side, region, amyloid_p_values_ad, Effect_Size, p_sig)) %>%
  dplyr::distinct(ROI, .keep_all = TRUE)

p_amyloid_all$adjusted_p <- p.adjust(p_amyloid_all$amyloid_p_values, method = "holm")
p_amyloid_all_cn$adjusted_p <- p.adjust(p_amyloid_all_cn$amyloid_p_values, method = "holm")
p_amyloid_all_mci$adjusted_p <- p.adjust(p_amyloid_all_mci$amyloid_p_values, method = "holm")
p_amyloid_all_ad$adjusted_p <- p.adjust(p_amyloid_all_ad$amyloid_p_values, method = "holm")
