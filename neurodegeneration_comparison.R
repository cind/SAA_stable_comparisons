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
#                                                    MRI
#############################################################################################################
#############################################################################################################
#############################################################################################################

#switch so first dates close to saa dates
#make everything after  amyloid positive - positive
#make everythin before amyloid negative - negative
#saa negative - take MRI image just before last saa- observation (ignore 2 year window)
#for saa positive - take MRI image just after first saa+ observation (ignore 2 year window)

#this file is the output of the MRI_CSV_processing steps plys apoe-, education-, and gender-corrected
mri_roi_cs <- read.csv("~/data/age_apoe_education_gender_SV_CV_corrected_mri.csv") %>%
  dplyr::filter(RID %in% stable_rids) %>%
  dplyr::mutate(EXAMDATE = as.Date(EXAMDATE))

#adding in the stable cases
mri_roi_cs <- merge(mri_roi_cs, stables, by = "RID") %>%
  dplyr::distinct()

####################################
# Step 1: Adding in Amyloid positivity (filling everything after 1st amyloid positive - positive, filling everything before last amyloid negative - negative)
####################################
#adding in amyloid pet
#getting amyloid negative cases with last amyloid negative dates
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

# getting amyloid positive cases with first amyloid positive dates
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

#merging amyloid info back into dataset so that I can use the first and last amyloid status dates
mri_roi_cs <- merge(mri_roi_cs, last_a_negative_date_pet, by = "RID", all.x = TRUE)
mri_roi_cs <- merge(mri_roi_cs, first_a_positive_date_pet, by = "RID", all.x = TRUE) %>%
  dplyr::mutate(last_saa_negative_date = as.Date(last_saa_negative_date),
                first_saa_positive_date = as.Date(first_saa_positive_date))

mri_roi_cs <- mri_roi_cs %>%
  dplyr::mutate(AmyloidPosPET = case_when(EXAMDATE <= last_a_neg_date_pet ~ 0,
                                          EXAMDATE >= first_a_pos_date_pet ~ 1)) %>%
  dplyr::select(RID, EXAMDATE, last_a_neg_date_pet, first_a_pos_date_pet, AmyloidPosPET, age, estaget0_ABETA42, STUDY, PTEDUCAT, PTGENDER, apoe, group, first_saa_positive_date, last_saa_negative_date, ends_with("_adjusted"), ends_with("TA"),ends_with("SV"), ends_with("CV"), ends_with("SA"))

rm(first_a_positive_date_pet, last_a_negative_date_pet)

####################################
# Step 2: getting MRI right before EXAMDATE (SAA-) and right after EXAMDATE (SAA+)
####################################
mri_roi_cs_saa_neg <- mri_roi_cs %>%
  dplyr::filter(group == "SAA- Stable")
mri_roi_cs_saa_pos <- mri_roi_cs %>%
  dplyr::filter(group == "SAA+ Stable")

#now pulling the EXAMDATE right before the last SAA negative date
mri_roi_cs_saa_neg <- mri_roi_cs_saa_neg %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(EXAMDATE), last_saa_negative_date), "years")) %>%
  dplyr::select(RID, EXAMDATE, last_a_neg_date_pet, first_a_pos_date_pet, time_diff, AmyloidPosPET, age, estaget0_ABETA42, STUDY, PTEDUCAT, PTGENDER, apoe, group, first_saa_positive_date, last_saa_negative_date, ends_with("_adjusted"), ends_with("SV"), ends_with("TA"), ends_with("CV"), ends_with("SA"))
mri_roi_cs_saa_neg <- mri_roi_cs_saa_neg %>%
  dplyr::filter(time_diff <= 0)
mri_roi_cs_saa_neg <- mri_roi_cs_saa_neg %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_diff == max(time_diff)) %>%
  dplyr::ungroup()

#now pulling the EXAMDATE right after the first SAA positive date
mri_roi_cs_saa_pos <- mri_roi_cs_saa_pos %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(EXAMDATE), first_saa_positive_date), "years")) %>%
  dplyr::select(RID, EXAMDATE, last_a_neg_date_pet, first_a_pos_date_pet, time_diff, AmyloidPosPET, age, estaget0_ABETA42, STUDY, PTEDUCAT, PTGENDER, apoe, group, first_saa_positive_date, last_saa_negative_date, ends_with("_adjusted"), ends_with("SV"), ends_with("TA"), ends_with("CV"), ends_with("SA"))
mri_roi_cs_saa_pos <- mri_roi_cs_saa_pos %>%
  dplyr::filter(time_diff >= 0)
mri_roi_cs_saa_pos <- mri_roi_cs_saa_pos %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()

#putting the two datasets with amyloid status and correct visit dates together
mri_roi_cs <- rbind(mri_roi_cs_saa_neg, mri_roi_cs_saa_pos) #604 RID's

####################################
# Step 3: adding nearest diagnoses
####################################
#adding diagnosis in
mri_roi_cs <- mri_roi_cs %>%
  dplyr::left_join(adni_diagnoses) %>%
  dplyr::mutate(time_diff = lubridate::time_length(difftime(as.Date(EXAMDATE), DX.DATE), "years")) %>%
  dplyr::distinct() %>%
  dplyr::select(RID, EXAMDATE, last_a_neg_date_pet, first_a_pos_date_pet, time_diff, AmyloidPosPET, age, estaget0_ABETA42, STUDY, PTEDUCAT, PTGENDER, apoe, group, first_saa_positive_date, last_saa_negative_date, DX.DATE, DX, ends_with("_adjusted"), ends_with("SV"), ends_with("TA"), ends_with("CV"), ends_with("SA"))
mri_roi_cs <- mri_roi_cs %>%
  dplyr::filter(!is.na(time_diff)) %>%
  dplyr::mutate(time_diff = abs(time_diff))
mri_roi_cs <- as.data.frame(mri_roi_cs)
mri_roi_cs <- mri_roi_cs %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
mri_roi_cs <- mri_roi_cs %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::arrange(EXAMDATE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

####################################
# Step 4: getting demographic information
####################################
#first getting rid of rows that have no information on features
mri_roi_cs <- mri_roi_cs[which(rowMeans(!is.na(mri_roi_cs)) > 0.8), ]

#getting distributions of diagnoses, gender, apoe, amyloid status
table(mri_roi_cs$group, mri_roi_cs$DX)
table(mri_roi_cs$group, mri_roi_cs$PTGENDER)
table(mri_roi_cs$group, mri_roi_cs$apoe)
table(mri_roi_cs$group, mri_roi_cs$AmyloidPosPET)
table(mri_roi_cs$group)

#getting averages of age and education
mean(mri_roi_cs_saa_pos$age)
sd(mri_roi_cs_saa_pos$age)
mean(mri_roi_cs_saa_neg$age)
sd(mri_roi_cs_saa_neg$age)
mean(mri_roi_cs_saa_pos$PTEDUCAT)
sd(mri_roi_cs_saa_pos$PTEDUCAT)
mean(mri_roi_cs_saa_neg$PTEDUCAT)
sd(mri_roi_cs_saa_neg$PTEDUCAT)

rm(mri_roi_cs_saa_neg, mri_roi_cs_saa_pos, amyloid_pet)

####################################
# Step 5: getting p-values
####################################

features_mri <- mri_roi_cs %>%
  dplyr::select(ends_with("_adjusted") | ends_with("SV") | ends_with("TA") | ends_with("CV") | ends_with("SA"))
featurenames_mri <- names(features_mri)

#pairwise comparisons of each region and saving to a list then data frame
mri_roi_cs_cn <- mri_roi_cs %>%
  dplyr::filter(DX == "CU")
mri_roi_cs_mci <- mri_roi_cs %>%
  dplyr::filter(DX == "MCI")
mri_roi_cs_ad <- mri_roi_cs %>%
  dplyr::filter(DX == "Dementia")

mri_p_values_cn <- vector()
mri_roi_names_cn <- vector()
mri_effect_sizes_cn <- vector()
n_cn <- vector()

for (f in 1:length(featurenames_mri)){
  feature <- featurenames_mri[f]
  print(paste("MRI data regression: ", feature))
  lm_feature <- lm(scale(mri_roi_cs_cn[, 17+f]) ~ group + AmyloidPosPET, data = mri_roi_cs_cn, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  mri_effect_sizes_cn <- append(mri_effect_sizes_cn, estimate)
  mri_p_values_cn <- append(mri_p_values_cn, group_p_value)
  mri_roi_names_cn <- append(mri_roi_names_cn, feature)
  n_cn <- append(n_cn, length(lm_feature$residuals))
}

p_values_data_cn <- data.frame(mri_p_values_cn)
mri_roi_names_cn <- gsub("_adjusted", "", mri_roi_names_cn) #getting rid of the adjusted part so that it lines up with the dict
p_values_data_cn$ROI <- mri_roi_names_cn #creating column for the roi names for each p value
p_values_data_cn$Effect_Size <- mri_effect_sizes_cn
p_values_data_cn$n <- n_cn

#MCI
mri_p_values_mci <- vector()
mri_roi_names_mci <- vector()
mri_effect_sizes_mci <- vector()
n_mci <- vector()

for (f in 1:length(featurenames_mri)){
  feature <- featurenames_mri[f]
  print(paste("MRI data regression: ", feature))
  lm_feature <- lm(scale(mri_roi_cs_mci[, 17+f]) ~ group + AmyloidPosPET, data = mri_roi_cs_mci, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  mri_effect_sizes_mci <- append(mri_effect_sizes_mci, estimate)
  mri_p_values_mci <- append(mri_p_values_mci, group_p_value)
  mri_roi_names_mci <- append(mri_roi_names_mci, feature)
  n_mci <- append(n_mci, length(lm_feature$residuals))
}

p_values_data_mci <- data.frame(mri_p_values_mci)
mri_roi_names_mci <- gsub("_adjusted", "", mri_roi_names_mci) #getting rid of the adjusted part so that it lines up with the dict
p_values_data_mci$ROI <- mri_roi_names_mci #creating column for the roi names for each p value
p_values_data_mci$Effect_Size <- mri_effect_sizes_mci
p_values_data_mci$n <- n_mci

#Dementia
mri_p_values_ad <- vector()
mri_roi_names_ad <- vector()
mri_effect_sizes_ad <- vector()
n_ad <- vector()

for (f in 1:length(featurenames_mri)){
  feature <- featurenames_mri[f]
  print(paste("MRI data regression: ", feature))
  lm_feature <- lm(scale(mri_roi_cs_ad[, 17+f]) ~ group + AmyloidPosPET, data = mri_roi_cs_ad, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  mri_effect_sizes_ad <- append(mri_effect_sizes_ad, estimate)
  mri_p_values_ad <- append(mri_p_values_ad, group_p_value)
  mri_roi_names_ad <- append(mri_roi_names_ad, feature)
  n_ad <- append(n_ad, length(lm_feature$residuals))
}

p_values_data_ad <- data.frame(mri_p_values_ad)
mri_roi_names_ad <- gsub("_adjusted", "", mri_roi_names_ad) #getting rid of the adjusted part so that it lines up with the dict
p_values_data_ad$ROI <- mri_roi_names_ad #creating column for the roi names for each p value
p_values_data_ad$Effect_Size <- mri_effect_sizes_ad
p_values_data_ad$n <- n_ad

#full data
mri_p_values <- vector()
mri_roi_names <- vector()
mri_effect_sizes <- vector()
n_all <- vector()

for (f in 1:length(featurenames_mri)){
  feature <- featurenames_mri[f]
  print(paste("MRI data regression: ", feature))
  lm_feature <- lm(scale(mri_roi_cs[, 17+f]) ~ group + DX + AmyloidPosPET, data = mri_roi_cs, na.action=na.exclude)
  summary <- summary(lm_feature)
  coefficients <- as.data.frame(summary$coefficients)
  p_values <- coefficients$`Pr(>|t|)`
  group_p_value <- p_values[2]
  estimate_values <- coefficients$Estimate
  estimate <- estimate_values[2]
  print(coefficients)
  mri_effect_sizes <- append(mri_effect_sizes, estimate)
  mri_p_values <- append(mri_p_values, group_p_value)
  mri_roi_names <- append(mri_roi_names, feature)
  n_all <- append(n_all, length(lm_feature$residuals))
}

p_values_data <- data.frame(mri_p_values)
mri_roi_names <- gsub("_adjusted", "", mri_roi_names) #getting rid of the adjusted part so that it lines up with the dict
p_values_data$ROI <- mri_roi_names #creating column for the roi names for each p value
p_values_data$Effect_Size <- mri_effect_sizes
p_values_data$n <- n_all

####################################
# Step 6: ASEG - organizing dictionary and variable names
####################################

#getting the dictionary so that the ROI names are listed
dict <- read.csv("C:\\Work Folder\\scanner information\\UCSF_scan_dict.csv") %>%
  dplyr::select(FLDNAME, TEXT) %>%
  dplyr::filter(!(FLDNAME == "RID"))

# getting the ROI name as an additional column
p_values_data <- merge(p_values_data, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE)
p_values_data_cn <- merge(p_values_data_cn, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE) 
p_values_data_mci <- merge(p_values_data_mci, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE) 
p_values_data_ad <- merge(p_values_data_ad, dict, by.x = "ROI", by.y = "FLDNAME", all.x = TRUE)

# replacing the column names with the name of the region itself
colnames(mri_roi_cs) <- gsub('_adjusted', '', colnames(mri_roi_cs))
colnames(mri_roi_cs_cn) <- gsub('_adjusted', '', colnames(mri_roi_cs_cn))
colnames(mri_roi_cs_mci) <- gsub('_adjusted', '', colnames(mri_roi_cs_mci))
colnames(mri_roi_cs_ad) <- gsub('_adjusted', '', colnames(mri_roi_cs_ad))

#filtering dictionary to only include variables that are in our dataset
dict_filtered <- dict %>%
  dplyr::filter(FLDNAME %in% names(mri_roi_cs))

test_data <- mri_roi_cs %>%
  rename_with(.cols = dict_filtered$FLDNAME, .fn = function(x) dict_filtered$TEXT[dict_filtered$FLDNAME %in% x])
test_data_cn <- mri_roi_cs_cn %>%
  rename_with(.cols = dict_filtered$FLDNAME, .fn = function(x) dict_filtered$TEXT[dict_filtered$FLDNAME %in% x])
test_data_mci <- mri_roi_cs_mci %>%
  rename_with(.cols = dict_filtered$FLDNAME, .fn = function(x) dict_filtered$TEXT[dict_filtered$FLDNAME %in% x])
test_data_ad <- mri_roi_cs_ad %>%
  rename_with(.cols = dict_filtered$FLDNAME, .fn = function(x) dict_filtered$TEXT[dict_filtered$FLDNAME %in% x])

#pulling aseg variables
volume_aseg <- test_data %>%
  dplyr::select(contains("aseg.stat"))
volume_aseg_cn <- test_data_cn %>%
  dplyr::select(contains("aseg.stat"))
volume_aseg_mci <- test_data_mci %>%
  dplyr::select(contains("aseg.stat"))
volume_aseg_ad <- test_data_ad %>%
  dplyr::select(contains("aseg.stat"))

#now looking at volume (aseg)
aseg_data <- as.data.frame(ggseg::aseg)

p_volume_aseg <- p_values_data %>%
  dplyr::filter(grepl("aseg.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1)) %>%
  dplyr::mutate(label1= case_when(label1 == "Brainstem" ~ "BrainStem",
                                  TRUE ~ label1))
p_volume_aseg_cn <- p_values_data_cn %>%
  dplyr::filter(grepl("aseg.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1)) %>%
  dplyr::mutate(label1= case_when(label1 == "Brainstem" ~ "BrainStem",
                                  TRUE ~ label1))
p_volume_aseg_mci <- p_values_data_mci %>%
  dplyr::filter(grepl("aseg.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1)) %>%
  dplyr::mutate(label1= case_when(label1 == "Brainstem" ~ "BrainStem",
                                  TRUE ~ label1))
p_volume_aseg_ad <- p_values_data_ad %>%
  dplyr::filter(grepl("aseg.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1)) %>%
  dplyr::mutate(label1= case_when(label1 == "Brainstem" ~ "BrainStem",
                                  TRUE ~ label1))

p_volume_aseg <- p_volume_aseg %>%
  dplyr::mutate(label_aseg = case_when(label1 == "RightPallidum" ~ "Right-Pallidum",
                                       label1 == "LeftPallidum" ~ "Left-Pallidum",
                                       label1 == "RightPutamen" ~ "Right-Putamen",
                                       label1 == "LeftPutamen" ~ "Left-Putamen",
                                       label1 == "RightThalamus" ~ "Right-Thalamus-Proper",
                                       label1 == "LeftThalamus" ~ "Left-Thalamus-Proper",
                                       label1 == "RightVentralDC" ~ "Right-VentralDC",
                                       label1 == "LeftVentralDC" ~ "Left-VentralDC",
                                       label1 == "RightLateralVentricle" ~ "Right-Lateral-Ventricle",
                                       label1 == "LeftLateralVentricle" ~ "Left-Lateral-Ventricle",
                                       label1 == "ThirdVentricle" ~ "x3rd-ventricle",
                                       label1 == "FourthVentricle" ~ "x4th-ventricle",
                                       label1 == "LeftAmygdala" ~ "Left-Amygdala",
                                       label1 == "RightAmygdala" ~  "Right-Amygdala",
                                       label1 == "LeftCaudate" ~ "Left-Caudate",
                                       label1 == "RightCaudate" ~ "Right-Caudate",
                                       label1 == "RightHippocampus" ~ "Right-Hippocampus",
                                       label1 == "LeftHippocampus" ~ "Left-Hippocampus",
                                       label1 == "CorpusCallosumMidPosterior" ~ "CC_Mid_Posterior",
                                       label1 == "CorpusCallosumMidAnterior" ~ "CC_Mid_Anterior",
                                       label1 == "CorpusCallosumAnterior" ~ "CC_Anterior",
                                       label1 == "CorpusCallosumPosterior" ~ "CC_Posterior",
                                       label1 == "CorpusCallosumCentral" ~ "CC_Central",
                                       label1 == "BrainStem" ~ "Brain-Stem",
                                       label1 == "RightCerebellumCortex" ~ "Right-Cerebellum-Cortex",
                                       label1 == "LeftCerebellumCortex" ~ "Left-Cerebellum-Cortex",
                                       label1 == "RightCerebellumWM" ~ "Right-Cerebellum-White-Matter",
                                       label1 == "LeftCerebellumWM" ~ "Left-Cerebellum-White-Matter",
                                       TRUE ~ label1))
p_volume_aseg_cn <- p_volume_aseg_cn %>%
  dplyr::mutate(label_aseg = case_when(label1 == "RightPallidum" ~ "Right-Pallidum",
                                       label1 == "LeftPallidum" ~ "Left-Pallidum",
                                       label1 == "RightPutamen" ~ "Right-Putamen",
                                       label1 == "LeftPutamen" ~ "Left-Putamen",
                                       label1 == "RightThalamus" ~ "Right-Thalamus-Proper",
                                       label1 == "LeftThalamus" ~ "Left-Thalamus-Proper",
                                       label1 == "RightVentralDC" ~ "Right-VentralDC",
                                       label1 == "LeftVentralDC" ~ "Left-VentralDC",
                                       label1 == "RightLateralVentricle" ~ "Right-Lateral-Ventricle",
                                       label1 == "LeftLateralVentricle" ~ "Left-Lateral-Ventricle",
                                       label1 == "ThirdVentricle" ~ "x3rd-ventricle",
                                       label1 == "FourthVentricle" ~ "x4th-ventricle",
                                       label1 == "LeftAmygdala" ~ "Left-Amygdala",
                                       label1 == "RightAmygdala" ~  "Right-Amygdala",
                                       label1 == "LeftCaudate" ~ "Left-Caudate",
                                       label1 == "RightCaudate" ~ "Right-Caudate",
                                       label1 == "RightHippocampus" ~ "Right-Hippocampus",
                                       label1 == "LeftHippocampus" ~ "Left-Hippocampus",
                                       label1 == "CorpusCallosumMidPosterior" ~ "CC_Mid_Posterior",
                                       label1 == "CorpusCallosumMidAnterior" ~ "CC_Mid_Anterior",
                                       label1 == "CorpusCallosumAnterior" ~ "CC_Anterior",
                                       label1 == "CorpusCallosumPosterior" ~ "CC_Posterior",
                                       label1 == "CorpusCallosumCentral" ~ "CC_Central",
                                       label1 == "BrainStem" ~ "Brain-Stem",
                                       label1 == "RightCerebellumCortex" ~ "Right-Cerebellum-Cortex",
                                       label1 == "LeftCerebellumCortex" ~ "Left-Cerebellum-Cortex",
                                       label1 == "RightCerebellumWM" ~ "Right-Cerebellum-White-Matter",
                                       label1 == "LeftCerebellumWM" ~ "Left-Cerebellum-White-Matter",
                                       TRUE ~ label1))
p_volume_aseg_mci <- p_volume_aseg_mci %>%
  dplyr::mutate(label_aseg = case_when(label1 == "RightPallidum" ~ "Right-Pallidum",
                                       label1 == "LeftPallidum" ~ "Left-Pallidum",
                                       label1 == "RightPutamen" ~ "Right-Putamen",
                                       label1 == "LeftPutamen" ~ "Left-Putamen",
                                       label1 == "RightThalamus" ~ "Right-Thalamus-Proper",
                                       label1 == "LeftThalamus" ~ "Left-Thalamus-Proper",
                                       label1 == "RightVentralDC" ~ "Right-VentralDC",
                                       label1 == "LeftVentralDC" ~ "Left-VentralDC",
                                       label1 == "RightLateralVentricle" ~ "Right-Lateral-Ventricle",
                                       label1 == "LeftLateralVentricle" ~ "Left-Lateral-Ventricle",
                                       label1 == "ThirdVentricle" ~ "x3rd-ventricle",
                                       label1 == "FourthVentricle" ~ "x4th-ventricle",
                                       label1 == "LeftAmygdala" ~ "Left-Amygdala",
                                       label1 == "RightAmygdala" ~  "Right-Amygdala",
                                       label1 == "LeftCaudate" ~ "Left-Caudate",
                                       label1 == "RightCaudate" ~ "Right-Caudate",
                                       label1 == "RightHippocampus" ~ "Right-Hippocampus",
                                       label1 == "LeftHippocampus" ~ "Left-Hippocampus",
                                       label1 == "CorpusCallosumMidPosterior" ~ "CC_Mid_Posterior",
                                       label1 == "CorpusCallosumMidAnterior" ~ "CC_Mid_Anterior",
                                       label1 == "CorpusCallosumAnterior" ~ "CC_Anterior",
                                       label1 == "CorpusCallosumPosterior" ~ "CC_Posterior",
                                       label1 == "CorpusCallosumCentral" ~ "CC_Central",
                                       label1 == "BrainStem" ~ "Brain-Stem",
                                       label1 == "RightCerebellumCortex" ~ "Right-Cerebellum-Cortex",
                                       label1 == "LeftCerebellumCortex" ~ "Left-Cerebellum-Cortex",
                                       label1 == "RightCerebellumWM" ~ "Right-Cerebellum-White-Matter",
                                       label1 == "LeftCerebellumWM" ~ "Left-Cerebellum-White-Matter",
                                       TRUE ~ label1))
p_volume_aseg_ad <- p_volume_aseg_ad %>%
  dplyr::mutate(label_aseg = case_when(label1 == "RightPallidum" ~ "Right-Pallidum",
                                       label1 == "LeftPallidum" ~ "Left-Pallidum",
                                       label1 == "RightPutamen" ~ "Right-Putamen",
                                       label1 == "LeftPutamen" ~ "Left-Putamen",
                                       label1 == "RightThalamus" ~ "Right-Thalamus-Proper",
                                       label1 == "LeftThalamus" ~ "Left-Thalamus-Proper",
                                       label1 == "RightVentralDC" ~ "Right-VentralDC",
                                       label1 == "LeftVentralDC" ~ "Left-VentralDC",
                                       label1 == "RightLateralVentricle" ~ "Right-Lateral-Ventricle",
                                       label1 == "LeftLateralVentricle" ~ "Left-Lateral-Ventricle",
                                       label1 == "ThirdVentricle" ~ "x3rd-ventricle",
                                       label1 == "FourthVentricle" ~ "x4th-ventricle",
                                       label1 == "LeftAmygdala" ~ "Left-Amygdala",
                                       label1 == "RightAmygdala" ~  "Right-Amygdala",
                                       label1 == "LeftCaudate" ~ "Left-Caudate",
                                       label1 == "RightCaudate" ~ "Right-Caudate",
                                       label1 == "RightHippocampus" ~ "Right-Hippocampus",
                                       label1 == "LeftHippocampus" ~ "Left-Hippocampus",
                                       label1 == "CorpusCallosumMidPosterior" ~ "CC_Mid_Posterior",
                                       label1 == "CorpusCallosumMidAnterior" ~ "CC_Mid_Anterior",
                                       label1 == "CorpusCallosumAnterior" ~ "CC_Anterior",
                                       label1 == "CorpusCallosumPosterior" ~ "CC_Posterior",
                                       label1 == "CorpusCallosumCentral" ~ "CC_Central",
                                       label1 == "BrainStem" ~ "Brain-Stem",
                                       label1 == "RightCerebellumCortex" ~ "Right-Cerebellum-Cortex",
                                       label1 == "LeftCerebellumCortex" ~ "Left-Cerebellum-Cortex",
                                       label1 == "RightCerebellumWM" ~ "Right-Cerebellum-White-Matter",
                                       label1 == "LeftCerebellumWM" ~ "Left-Cerebellum-White-Matter",
                                       TRUE ~ label1))

p_volume_aseg <- merge(p_volume_aseg, aseg_data, by.x = "label_aseg", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)
p_volume_aseg_cn <- merge(p_volume_aseg_cn, aseg_data, by.x = "label_aseg", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)
p_volume_aseg_mci <- merge(p_volume_aseg_mci, aseg_data, by.x = "label_aseg", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)
p_volume_aseg_ad <- merge(p_volume_aseg_ad, aseg_data, by.x = "label_aseg", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) %>%
  dplyr::rename(label = label_aseg)

####################################
# Step 7: ASEG - determining which variables to plot by significance
####################################

p_volume_aseg <- p_volume_aseg %>%
  dplyr::mutate(p_sig = case_when(mri_p_values >= .05 ~ "not sig",
                                  mri_p_values < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_volume_aseg_cn <- p_volume_aseg_cn %>%
  dplyr::mutate(p_sig = case_when(mri_p_values_cn >= .05 ~ "not sig",
                                  mri_p_values_cn < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_volume_aseg_mci <- p_volume_aseg_mci %>%
  dplyr::mutate(p_sig = case_when(mri_p_values_mci >= .05 ~ "not sig",
                                  mri_p_values_mci < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_volume_aseg_ad <- p_volume_aseg_ad %>%
  dplyr::mutate(p_sig = case_when(mri_p_values_ad >= .05 ~ "not sig",
                                  mri_p_values_ad < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))

####################################
# Step 8: APARC - organizing dictionary and variable names
####################################

aparc_data <- as.data.frame(dk)

p_volume_aparc <- p_values_data %>%
  dplyr::filter(grepl("aparc.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1))
p_volume_aparc_cn <- p_values_data_cn %>%
  dplyr::filter(grepl("aparc.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1))
p_volume_aparc_mci <- p_values_data_mci %>%
  dplyr::filter(grepl("aparc.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1))
p_volume_aparc_ad <- p_values_data_ad %>%
  dplyr::filter(grepl("aparc.stat", TEXT)) %>%
  dplyr::mutate(label1 = word(TEXT, -1))

volume_aparc <- test_data %>%
  dplyr::select(contains("aparc.stat"))

p_volume_aparc$label1 <- tolower(p_volume_aparc$label1)
p_volume_aparc_cn$label1 <- tolower(p_volume_aparc_cn$label1)
p_volume_aparc_mci$label1 <- tolower(p_volume_aparc_mci$label1)
p_volume_aparc_ad$label1 <- tolower(p_volume_aparc_ad$label1)

p_volume_aparc$label1 <- gsub("right", "rh_", p_volume_aparc$label1)
p_volume_aparc_cn$label1 <- gsub("right", "rh_", p_volume_aparc_cn$label1)
p_volume_aparc_mci$label1 <- gsub("right", "rh_", p_volume_aparc_mci$label1)
p_volume_aparc_ad$label1 <- gsub("right", "rh_", p_volume_aparc_ad$label1)

p_volume_aparc$label1 <- gsub("left", "lh_", p_volume_aparc$label1)
p_volume_aparc_cn$label1 <- gsub("left", "lh_", p_volume_aparc_cn$label1)
p_volume_aparc_mci$label1 <- gsub("left", "lh_", p_volume_aparc_mci$label1)
p_volume_aparc_ad$label1 <- gsub("left", "lh_", p_volume_aparc_ad$label1)


p_volume_aparc <- merge(p_volume_aparc, aparc_data, by.x = "label1", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) 
p_volume_aparc_cn <- merge(p_volume_aparc_cn, aparc_data, by.x = "label1", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) 
p_volume_aparc_mci <- merge(p_volume_aparc_mci, aparc_data, by.x = "label1", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) 
p_volume_aparc_ad <- merge(p_volume_aparc_ad, aparc_data, by.x = "label1", by.y = "label", all.x = TRUE) %>%
  dplyr::filter(!(is.na(region))) 

####################################
# Step 9: APARC - determining which variables to plot by significance
####################################

p_volume_aparc <- p_volume_aparc %>%
  dplyr::mutate(p_sig = case_when(mri_p_values >= .05 ~ "not sig",
                                  mri_p_values < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_volume_aparc_cn <- p_volume_aparc_cn %>%
  dplyr::mutate(p_sig = case_when(mri_p_values_cn >= .05 ~ "not sig",
                                  mri_p_values_cn < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_volume_aparc_mci <- p_volume_aparc_mci %>%
  dplyr::mutate(p_sig = case_when(mri_p_values_mci >= .05 ~ "not sig",
                                  mri_p_values_mci < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))
p_volume_aparc_ad <- p_volume_aparc_ad %>%
  dplyr::mutate(p_sig = case_when(mri_p_values_ad >= .05 ~ "not sig",
                                  mri_p_values_ad < .05 ~ "sig")) %>%
  dplyr::filter(!(is.na(p_sig)))

####################################
# Step 10: Plotting
####################################
#ASEG
b <- c(-.45, -.30, -.15, 0, .15, .30, .45)

aseg_full <- ggplot(p_volume_aseg %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.51,.51),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("All Data")
aseg_cn <- ggplot(p_volume_aseg_cn %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.51,.51),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("CU Data")
aseg_mci <- ggplot(p_volume_aseg_mci %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size))  + # significant rois include: CC_Anterior, Left-Amygdala, Left-Hippocampus, Right-Hippocampus
  scale_fill_gradientn(limits = c(-.51,.51),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("MCI Data")
aseg_ad <- ggplot(p_volume_aseg_ad %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=aseg, mapping=aes(fill=Effect_Size)) + #no significant rois
  scale_fill_gradientn(limits = c(-.51,.51),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b, labels=format(b)) + 
  ggtitle("AD Data")

aseg_full + aseg_cn + aseg_mci + aseg_ad + 
  plot_layout(ncol = 1, guides = "collect")

#APARC
b2 <- c(-.5, -.25, 0, .25, .5) # b2 <- c(-1, -.5, 0, .5, 1)

aparc_full <- ggplot(p_volume_aparc %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.6,.6), #   scale_fill_gradientn(limits = c(-1.1,1.1),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("All Data")
aparc_cn <- ggplot(p_volume_aparc_cn %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #significant rois include: Left-Lateral-Ventricle
  scale_fill_gradientn(limits = c(-.6,.6),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) +
  ggtitle("CU Data")
aparc_mci <- ggplot(p_volume_aparc_mci %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size))  + # significant rois include: CC_Anterior, Left-Amygdala, Left-Hippocampus, Right-Hippocampus
  scale_fill_gradientn(limits = c(-.6,.6),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2)) + 
  ggtitle("MCI Data")
aparc_ad <- ggplot(p_volume_aparc_ad %>% dplyr::filter(p_sig == "sig")) +
  geom_brain( atlas=dk, mapping=aes(fill=Effect_Size)) + #no significant rois
  scale_fill_gradientn(limits = c(-.6,.6),
                       colours=c("red3", "goldenrod2", "darkgreen", "darkblue"),
                       breaks=b2, labels=format(b2))  + 
  ggtitle("AD Data")

aparc_full + aparc_cn + aparc_mci + aparc_ad + 
  plot_layout(ncol = 1, guides = "collect")

####################################
# Step 12: p-value adjustment
####################################

p_volumes_all <- rbind(p_volume_aseg %>% dplyr::select(atlas, hemi, side, region, mri_p_values, Effect_Size, p_sig), p_volume_aparc %>% dplyr::select(atlas, hemi, side, region, mri_p_values, Effect_Size, p_sig))
p_volumes_all_cn <- rbind(p_volume_aseg_cn %>% dplyr::select(atlas, hemi, side, region, mri_p_values_cn, Effect_Size, p_sig), p_volume_aparc_cn %>% dplyr::select(atlas, hemi, side, region, mri_p_values_cn, Effect_Size, p_sig))
p_volumes_all_mci <- rbind(p_volume_aseg_mci %>% dplyr::select(atlas, hemi, side, region, mri_p_values_mci, Effect_Size, p_sig), p_volume_aparc_mci %>% dplyr::select(atlas, hemi, side, region, mri_p_values_mci, Effect_Size, p_sig))
p_volumes_all_ad <- rbind(p_volume_aseg_ad %>% dplyr::select(atlas, hemi, side, region, mri_p_values_ad, Effect_Size, p_sig), p_volume_aparc_ad %>% dplyr::select(atlas, hemi, side, region, mri_p_values_ad, Effect_Size, p_sig))

p_volumes_all$adjusted_p <- p.adjust(p_volumes_all$mri_p_values, method = "holm")
p_volumes_all_cn$adjusted_p <- p.adjust(p_volumes_all_cn$mri_p_values, method = "holm")
p_volumes_all_mci$adjusted_p <- p.adjust(p_volumes_all_mci$mri_p_values, method = "holm")
p_volumes_all_ad$adjusted_p <- p.adjust(p_volumes_all_ad$mri_p_values, method = "holm")
