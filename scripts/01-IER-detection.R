######################################################
##                                                  ##
##        Unravelling the Capacity Action Gap       ##
##                                                  ##
##                Annika Schubert                   ##
##                                                  ##
##                01-IER-detection                  ##
##                                                  ##
######################################################

# @brief        Data preparation: restrict sample, recode variables 
# @license      GNU General Public License v3.0  
# @authors      Annika Schubert & Felix Bauer

# Code for:
# Schubert, A.; von Streit, A. & Garschagen, M. (2024):
# Unravelling the capacity-action gap in flood risk adaptation. 
# In: EGUsphere [preprint], https://doi.org/X. 


# ###################################################### #
#                                                        #
#### Part 1: Insufficient Effort Responding Detection ####
#                                                        #
# ###################################################### #

### IER can influence data quality, correlations, factor structures and 
# statistical power and may lead to
# - reduce measurement reliability
# - inflated correlation

### 3 measures: 
# IRV: to detect straightlining in item batteries (Dunn et al. 2018)
# RIR: to check individual consistency within item batteries (Brühlmann et al. 2020)
# response time: to detect rushing, fatigue, inattentiveness (Greszki et al. 2015)

library(here)               # for file referencing
library(tidyverse)          # for data manipulation, exploration and visualization
library(car)                # for recoding
library(haven)              # for working with (SPSS) labbelled dataset
library(remotes)            # to install packages from repositories (e.g. careless)
library(careless)           # to detect such careless/insufﬁcient effort responses
                            # remotes::install_github('ryentes/careless')


# 1.1 Load & recode data---------------------------

load(here::here("./data/raw/KARE_town.RData"))

# no. of answered survey items
KARE_IERDetection <- mutate(KARE_IERDetection, ansq = rowSums(!is.na(KARE_IERDetection[,6:288])))
table(KARE_IERDetection$ansq)

# create dataset with measures & flags
KARE_DQ <- data.frame(IDS = KARE_IERDetection$IDS,
                      ansq = KARE_IERDetection$ansq)

# 1.2 Post-hoc measures to detect IER ---------------------------

# 1.2.1 Intra-Individual Response Variability (IRV) ---------------------------

# = standard deviation of responses across a set of consecutive item responses
# to detect straightlining in item batteries (Dunn et al. 2018)

# no recoding of reversed items
# cut-off for flagging depends on no of items, answer patterns, ...

### calculation IRV for every battery

# 1) R16 NEP Scale
KARE_IERDetection$IRV_16<- irv(KARE_IERDetection[, c("R16_1", "R16_2", "R16_3")], na.rm = T)
KARE_IERDetection$IRV16_flag <- ifelse(KARE_IERDetection$IRV_16  < 1, 1, 0) 
table(KARE_IERDetection$IRV16_flag) # 193 flags

# 2) R40 Damage Consequences Loop 1
KARE_IERDetection$IRV_40 <- irv(KARE_IERDetection[, c("R40_1", "R40_2", "R40_3", "R40_4", "R40_5", "R40_6", "R40_7", "R40_8", "R40_9", "R40_10")], na.rm = T)
KARE_IERDetection$IRV40_flag <- ifelse(KARE_IERDetection$IRV_40 < 0.4216370, 1, 0)
table(KARE_IERDetection$IRV40_flag) # 114 flags

# 3) R63 Responsibility State, Household, ...
KARE_IERDetection$IRV_63 <- irv(KARE_IERDetection[,  c("R63_1", "R63_2", "R63_3", "R63_4", "R63_5", "R63_6", "R63_7")], na.rm = T)
KARE_IERDetection$IRV63_flag <- ifelse(KARE_IERDetection$IRV_63 < 0.8, 1, 0)
table(KARE_IERDetection$IRV63_flag) # 112 flags

# 4) R64 mandatory measures for houses
KARE_IERDetection$IRV_64 <- irv(KARE_IERDetection[, c("R64_1", "R64_2", "R64_3", "R64_4")], na.rm = T)
KARE_IERDetection$IRV64_flag <- ifelse(KARE_IERDetection$IRV_64 < 0.51, 1, 0)
table(KARE_IERDetection$IRV64_flag)  # 359 flags

# 5) R80 Effective Protection from different groups (municipality, fire brigade, ...)
KARE_IERDetection$IRV_80 <- irv(KARE_IERDetection[, c("R80a_1", "R80a_2", "R80a_3", "R80b_1", "R80b_2")], na.rm = T)
KARE_IERDetection$IRV80_flag <- ifelse(KARE_IERDetection$IRV_80 < 0.51, 1, 0)
table(KARE_IERDetection$IRV80_flag)  # 442 flags

# 6) R81 Responsibility - own vs municipal administration
KARE_IERDetection$IRV_81 <- irv(KARE_IERDetection[, c("R81_1", "R81_2")], na.rm = T)
KARE_IERDetection$IRV81_flag <- ifelse(KARE_IERDetection$IRV_81 < 0.01, 1, 0)
table(KARE_IERDetection$IRV81_flag) # 380 flags

# 7) R88 Aid
KARE_IERDetection$IRV_88 <- irv(KARE_IERDetection[, c("R88_1", "R88_2", "R88_3", "R88_4")], na.rm = T)
KARE_IERDetection$IRV88_flag <- ifelse(KARE_IERDetection$IRV_88 < 0.8, 1, 0)
table(KARE_IERDetection$IRV88_flag) # 193 flags

# 8) R89 Firms
KARE_IERDetection$IRV_89 <- irv(KARE_IERDetection[, c("R89_1", "R89_2", "R89_3")], na.rm = T)
KARE_IERDetection$IRV89_flag <- ifelse(KARE_IERDetection$IRV_89 < 1, 1, 0)
table(KARE_IERDetection$IRV89_flag) # 281 flags

# 9) R90 Social Network
KARE_IERDetection$IRV_90 <- irv(KARE_IERDetection[, c("R90_1", "R90_3", "R90_5")], na.rm = T)
KARE_IERDetection$IRV90_flag <- ifelse(KARE_IERDetection$IRV_90 < 0.01, 1, 0)
table(KARE_IERDetection$IRV90_flag)  # 307 flags


### Flagging due to IRV

# no. of answered batteries
KARE_IERDetection$batq <- apply(KARE_IERDetection[,c("IRV16_flag", "IRV40_flag", "IRV63_flag", 
                                                     "IRV64_flag", "IRV80_flag", "IRV81_flag", 
                                                     "IRV88_flag", "IRV89_flag", "IRV90_flag")], 
                                1, function(x) sum(!is.na(x)))
table(KARE_IERDetection$batq)

# times flagged due to straightlining
KARE_IERDetection$IRV_sumflags <- rowSums(KARE_IERDetection[,c("IRV16_flag", "IRV40_flag", "IRV63_flag", 
                                                               "IRV64_flag", "IRV80_flag", "IRV81_flag", 
                                                               "IRV88_flag", "IRV89_flag", "IRV90_flag")],
                                          na.rm=TRUE)
table(KARE_IERDetection$IRV_sumflags)

# percentage of flagged batteries
KARE_IERDetection$IRV_flagsperbattery <- KARE_IERDetection$IRV_sumflags / KARE_IERDetection$batq
table(KARE_IERDetection$IRV_flagsperbattery, useNA = "always")
hist(KARE_IERDetection$IRV_flagsperbattery)
plot(density(KARE_IERDetection$IRV_flagsperbattery, na.rm = T))

# flag if straightlining in more than 40% answered batteries or if no battery answered
KARE_IERDetection$IRV_flag <- ifelse(KARE_IERDetection$IRV_flagsperbattery > 0.4 | 
                                       is.na(KARE_IERDetection$IRV_flagsperbattery), 1, 0)
table(KARE_IERDetection$IRV_flag, useNA = "always") # 230 flags

KARE_DQ <-left_join(KARE_DQ, 
                    KARE_IERDetection[,c("IDS", "IRV_flagsperbattery", "IRV_flag")], 
                    by = c("IDS"))


# 1.2.2 Even-odd inconsistency (EO) / Resampled Individual Reliability (RIR) ---------------------------

# to check individual consistency within item batteries

# check internal consistency of batteries with Cronbachs alpha,
# keep scales with alpha >=  0.65
psych::alpha(subset(KARE_IERDetection, select = c(R12_1, R12_2, R13_1, R13_2)), check.keys = TRUE) # 0.84
psych::alpha(subset(KARE_IERDetection, select = c(R13a, R14, R15)), check.keys = TRUE) # 0.65
psych::alpha(subset(KARE_IERDetection, select = c(R64_1, R64_2, R64_3, R64_4)), check.keys = TRUE)  # 0.79
psych::alpha(subset(KARE_IERDetection, select = c(R80a_1, R80a_2, R80a_3, R80b_1, R80b_2)), check.keys = TRUE) # 0.74
psych::alpha(subset(KARE_IERDetection, select = c(R90_1, R90_3, R90_5)), check.keys = TRUE) # 0.7 

# use only batteries which measure one construct
Consvars <- c("R12_1", "R12_2", "R13_1", "R13_2",
              "R13a", "R14", "R15", 
              "R64_1", "R64_2", "R64_3", "R64_4",
              "R80a_1", "R80a_2", "R80a_3", "R80b_1", "R80b_2", 
              "R90_1", "R90_3", "R90_5")

KARE_Cons <- KARE_IERDetection[names(KARE_IERDetection) %in% Consvars]
KARE_Cons[] <- lapply(KARE_Cons, function(x) { attributes(x) <- NULL; x }) # remove attributes
KARE_Cons <- labelled::remove_attributes(KARE_Cons, "attr")     # remove attribute

# change vars names (important for RIR)
KARE_Cons <- rename(KARE_Cons,
                    R12_3 = R13_1,
                    R12_4 = R13_2,
                    R13b = R14,
                    R13c = R15)
 
# 1.2.2.1 Even-odd inconsistency (EO) ---------------------------

# = correlation between scores on odd and even halves across subscales

KARE_IERDetection$EO <- evenodd(KARE_Cons, factors = c(4,3,4,5,3))

hist(KARE_IERDetection$EO, 
     main = "Histogram of the Even-Odd consistency index", 
     xlab = "Even-Odd Index")

# EO flag: > 0 (Curran 2016)
KARE_IERDetection$EO_flag <- ifelse(KARE_IERDetection$EO > 0 | 
                                      is.na(KARE_IERDetection$EO), 1, 0)
table(KARE_IERDetection$EO_flag) # 341 flags


# 1.2.2.2 Resampled Individual Reliability (RIR) ---------------------------

# RIR code from Brühlmann et al. (2020)

# changes: added na action to also calculate values in case of NA(s)
n.iter = 100
boot_mat <- matrix(nrow = nrow(KARE_Cons), ncol = n.iter)
for (i in 1:n.iter) {
  set.seed(i)                                      # added: set.seed to make output reproducible
  res_indiv_reliability <- c()
  for (row in 1:nrow(KARE_Cons)) {
    df_person <- KARE_Cons[row,]
    vector1 <- c()
    vector2 <- c()
    for (SCALE in c("R12", "R13", "R64", "R80", "R90")) {
      df_tmp <- df_person[,grep(SCALE, names(df_person), value = TRUE)]
      half1 <- sample(names(df_tmp),ncol(df_tmp)/2)
      half2 <- !(names(df_tmp) %in% half1)
      h1_i <- df_tmp[,half1]
      h2_i <- df_tmp[,half2]
      vector1 <- c(vector1, ifelse(length(half1) > 1, rowMeans(h1_i, na.rm = T), h1_i))            # added: rowmeans na.rm = T; if only one var, no mean but value
      vector2 <- c(vector2, ifelse(sum(half2, na.rm=TRUE) > 1, rowMeans(h2_i, na.rm = T), h2_i))   # added: na.rm = T; if only one var, no mean but value
    }
    res_indiv_reliability <- c(res_indiv_reliability, cor(unlist(vector1), vector2,  # added: bug fixing in case vector1 is saved as a list
                                                          use = "pairwise.complete.obs")) # added:use = "pairwise.complete.obs" 
  }
  boot_mat[,i] <- res_indiv_reliability
}
save(boot_mat, file = here::here("./data/raw/RIR_mat.RData"))

load(here::here("./data/raw/RIR_mat.RData"))

# calculate mean for each respondent
KARE_IERDetection$RIR <- rowMeans(boot_mat, na.rm = T)
table(KARE_IERDetection$RIR)

# flag respondents with NA values & < 0 (Curran 2016, Brühlmann et al. 2020)
KARE_IERDetection$RIR_flag<- ifelse(KARE_IERDetection$RIR < 0 | is.na(KARE_IERDetection$RIR), 1, 0) 
table(KARE_IERDetection$RIR_flag)
prop.table(table(KARE_IERDetection$RIR_flag))
# 200 flags

# strong negative correlation between EO & RIR (-0.87)
cor(KARE_IERDetection$RIR, KARE_IERDetection$EO, use = "pairwise.complete.obs") 

# add RIR & RIR_flag to the KARE_DQ dataset
KARE_DQ <-left_join(KARE_DQ, KARE_IERDetection[,c("IDS", "RIR", "RIR_flag")], by = c("IDS"))


# 1.2.3 Response Time ---------------------------

# to detect rushing

# average response time (in seconds) per item
KARE_IERDetection$avg_time <- KARE_IERDetection$Befragungsdauer_Sek / KARE_IERDetection$ansq
# respondents who answered not a single question have an Inf (dividing by 0)
KARE_IERDetection$avg_time[KARE_IERDetection$avg_time == Inf] <- NA

hist(KARE_IERDetection$avg_time, xlim = c(0,60), breaks = 200, main = 'Average time per question, showing only values from 0-60')
quantile(KARE_IERDetection$avg_time, na.rm = T, probs = seq(.1, .9, by = .1))
boxplot(KARE_IERDetection$avg_time, main = 'Average time per question')
boxplot(KARE_IERDetection$avg_time, outline = FALSE, main = 'Average time per question, without upper outliers')


### combined approach: statistical + cognitive approach 
# Greszki et al. (2015):
# --> response times are likely to be influenced by individual differences, for example 
# due to cognitive abilities and cognitive aging 
# --> utilize group-specific rather than overall, medians to avoid falsely flagging 
# highly educated young respondents as speeders, who answer particularly quickly 
# due to cognitive abilities and Internet experience
# --> calculated group-specific page medians

# age differences in response time
plot(KARE_IERDetection$R98[KARE_IERDetection$avg_time <100], KARE_IERDetection$avg_time[KARE_IERDetection$avg_time <100])

# age groups:
table(KARE_IERDetection$R98)                  
KARE_IERDetection <- mutate(KARE_IERDetection,                    
                            age_cat = case_when(
                              R98 %in% 19:44 ~ "19-44",
                              R98 %in% 45:69 ~ "45-69",
                              R98 %in% 70:93 ~ "70-93"))
table(KARE_IERDetection$age_cat)

# slight differences between educational levels
boxplot(KARE_IERDetection$avg_time[KARE_IERDetection$avg_time <100] ~ KARE_IERDetection$R104[KARE_IERDetection$avg_time <100])

# education groups
table(KARE_IERDetection$R104)
KARE_IERDetection$edu_cat <-  recode_factor(as.numeric(KARE_IERDetection$R104), 
                                            `1` = "Low",
                                            `2` = "Low",
                                            `3` = "Intermediate/high",
                                            `4` = "Intermediate/high",
                                            `5` = "Intermediate/high",
                                            `6` = "Low")
table(KARE_IERDetection$edu_cat) 


# 3*3*2 = 18 subgroups
subgroups <- table(KARE_IERDetection$Modus, KARE_IERDetection$age_cat, KARE_IERDetection$edu_cat, useNA = "always")
table(subgroups)
# some subgroups are rather small (Low edu + low age)

# 1) calculate median per subgroup
KARE_IERDetection$speed_med <- NA
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CATI" & 
                                        KARE_IERDetection$age_cat == "19-44" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" &  KARE_IERDetection$age_cat == "19-44" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CATI" & 
                                        KARE_IERDetection$age_cat == "45-69" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" &  KARE_IERDetection$age_cat == "45-69" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CATI" & 
                                        KARE_IERDetection$age_cat == "70-93" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" &  KARE_IERDetection$age_cat == "70-93" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "19-44" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" &  KARE_IERDetection$age_cat == "19-44" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "45-69" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" &  KARE_IERDetection$age_cat == "45-69" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "70-93" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" &  KARE_IERDetection$age_cat == "70-93" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "Unsicher" & 
                                        KARE_IERDetection$age_cat == "19-44" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" &  KARE_IERDetection$age_cat == "19-44" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "Unsicher" & 
                                        KARE_IERDetection$age_cat == "45-69" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" &  KARE_IERDetection$age_cat == "45-69" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "Unsicher" & 
                                        KARE_IERDetection$age_cat == "70-93" &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" &  KARE_IERDetection$age_cat == "70-93" &KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)

KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CATI" & 
                                        KARE_IERDetection$age_cat == "19-44" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" &  KARE_IERDetection$age_cat == "19-44" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CATI" & 
                                        KARE_IERDetection$age_cat == "45-69" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" &  KARE_IERDetection$age_cat == "45-69" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CATI" & 
                                        KARE_IERDetection$age_cat == "70-93" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" &  KARE_IERDetection$age_cat == "70-93" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "19-44" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" &  KARE_IERDetection$age_cat == "19-44" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "45-69" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" &  KARE_IERDetection$age_cat == "45-69" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "70-93" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" &  KARE_IERDetection$age_cat == "70-93" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "Unsicher" & 
                                        KARE_IERDetection$age_cat == "19-44" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" &  KARE_IERDetection$age_cat == "19-44" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "Unsicher" & 
                                        KARE_IERDetection$age_cat == "45-69" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" &  KARE_IERDetection$age_cat == "45-69" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(KARE_IERDetection$Modus == "Unsicher" & 
                                        KARE_IERDetection$age_cat == "70-93" &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" &  KARE_IERDetection$age_cat == "70-93" &KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
# in case of NAs: specify subgroup as far as possible
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CATI" & 
                                        is.na(KARE_IERDetection$age_cat)&
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" & KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CATI" & 
                                        KARE_IERDetection$age_cat == "70-93"&
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI" & KARE_IERDetection$age_cat == "70-93"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CAWI" & 
                                        is.na(KARE_IERDetection$age_cat) &
                                        KARE_IERDetection$edu_cat == "Low",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" & KARE_IERDetection$edu_cat == "Low"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CAWI" & 
                                        is.na(KARE_IERDetection$age_cat)&
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" & KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "19-44"&
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" & KARE_IERDetection$age_cat == "19-44"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "45-69"&
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" & KARE_IERDetection$age_cat == "45-69"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CAWI" & 
                                        KARE_IERDetection$age_cat == "70-93"&
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI" & KARE_IERDetection$age_cat == "70-93"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "Unsicher" & 
                                        is.na(KARE_IERDetection$age_cat) &
                                        KARE_IERDetection$edu_cat == "Intermediate/high",
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" & KARE_IERDetection$edu_cat == "Intermediate/high"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "Unsicher" & 
                                        KARE_IERDetection$age_cat == "45-69" &
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher" & KARE_IERDetection$age_cat == "45-69"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CATI" & 
                                        is.na(KARE_IERDetection$age_cat)&
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CATI"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "CAWI" & 
                                        is.na(KARE_IERDetection$age_cat)&
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "CAWI"], na.rm = T),
                                      KARE_IERDetection$speed_med)
KARE_IERDetection$speed_med <- ifelse(is.na(KARE_IERDetection$speed_med) &
                                        KARE_IERDetection$Modus == "Unsicher" & 
                                        is.na(KARE_IERDetection$age_cat)&
                                        is.na(KARE_IERDetection$edu_cat),
                                      median(KARE_IERDetection$avg_time[KARE_IERDetection$Modus == "Unsicher"], na.rm = T),
                                      KARE_IERDetection$speed_med)
table(KARE_IERDetection$speed_med, useNA= "always")


# 2) compare individual response times with group-specific median
KARE_IERDetection$speedindex <- KARE_IERDetection$avg_time / KARE_IERDetection$speed_med
table(KARE_IERDetection$speedindex, useNA= "always")
plot(density(KARE_IERDetection$speedindex, na.rm = T))
boxplot(KARE_IERDetection$speedindex[KARE_IERDetection$speedindex<10])

# cut off:
table(KARE_IERDetection$speedindex[KARE_IERDetection$speedindex < 0.5])   # fast
table(KARE_IERDetection$speedindex[KARE_IERDetection$speedindex > 2])     # slow

# The most exclusive measures: 50% faster than the respective group median
KARE_IERDetection$time_flag <- ifelse(KARE_IERDetection$speedindex < 0.5 |    # 50% faster than the respective group median
                                        KARE_IERDetection$speedindex > 2 |    # but also flag very slow respondents
                                        is.na(KARE_IERDetection$speedindex),  # and respondents with NA value
                                      1, 0) 
table(KARE_IERDetection$time_flag) # 222 flags

KARE_DQ <-left_join(KARE_DQ, KARE_IERDetection[,c("IDS", "avg_time", "time_flag")],
                    by = c("IDS"))



# 1.3 Check Associations of Measures ---------------------------

# 1) Response Time & IRV
DescTools::Assocs(table(KARE_IERDetection$time_flag, KARE_IERDetection$IRV_flag))
# Cramer's V: 0.2197 
chisq.test(table(KARE_IERDetection$time_flag, KARE_IERDetection$IRV_flag), correct = T) 
# Contingency Coeff. 0.2146***

# 2) Response Time & RIR
DescTools::Assocs(table(KARE_IERDetection$time_flag, KARE_IERDetection$RIR_flag))
# Cramer's V: 0.2311
chisq.test(table(KARE_IERDetection$time_flag, KARE_IERDetection$RIR_flag)) 
# Contingency Coeff. 0.2252***

# 3) IRV & RIR
DescTools::Assocs(table(KARE_IERDetection$IRV_flag, KARE_IERDetection$RIR_flag))
# Cramer's V: 0.4866
chisq.test(table(KARE_IERDetection$IRV_flag, KARE_IERDetection$RIR_flag)) 
# Contingency Coeff. 0.4376***


# 1.4 Flagging of IER Responents ---------------------------

# add up time, RIR and IRV flags
KARE_DQ$flagsum <- rowSums(KARE_DQ[ ,c("time_flag", "RIR_flag", "IRV_flag")], na.rm=TRUE)
table(KARE_DQ$flag)

KARE_DQ$IER_flag <- ifelse(KARE_DQ$flagsum > 1, 1, 0)
table(KARE_DQ$IER_flag)  # 131 respondents are flagged due to IER

# save KARE_DQ dataset
save(KARE_DQ, file = here::here("./data/tidy/KARE_DQ.RData"))

