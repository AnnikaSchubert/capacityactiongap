######################################################
##                                                  ##
##        Unravelling the Capacity Action Gap       ##
##                                                  ##
##                Annika Schubert                   ##
##                                                  ##
##            02-sample-and-recoding                ##
##                                                  ##
######################################################

# @brief        Data preparation: restrict sample, recode variables 
# @license      GNU General Public License v3.0  
# @authors      Annika Schubert

# Code for:
# Schubert, A.; von Streit, A. & Garschagen, M. (2024):
# Unravelling the capacity-action gap in flood risk adaptation. 
# In: EGUsphere [preprint], https://doi.org/X. 


# ############################################## #
#                                                #
#### Part 2: Data Preparation                 ####
#                                                #
# ############################################## #

library(here)                # for file referencing
library(tidyverse)           # for data manipulation, exploration and visualization
library(forcats)             # for recoding categorical (factor) variables
library(ggplot2)             # plots


load(here::here("./data/raw/KARE_town.RData"))


# 3.1 Sample ------------

# Sample Restriction
# - drop respondents, who did not finish the questionnaire (3.1.1)
# - drop careless respondents (3.1.2)
# - drop respondents who are neither tenant nor homeowner (3.1.3)


# 3.1.1 Finished Interviews ------------

# keep only finished interviews --> respondents must have answered at least
# parts of the sociodemographic questions at the end of the questionnaire

socio_miss <- c("R97", "R98", "R99", "R100", "R101", "R102", "R104", 
                "R105", "R106", "R109", "R96") 

# calculate number of NAs for each observation across sociodemographics
KARE <- mutate(KARE, miss_socio = rowSums(across(all_of(socio_miss), ~ is.na(.))))
table(KARE$miss_socio)

# sample: drop respondents with more than seven missing on sociodemographic vars
KARE <- subset(KARE, KARE$miss_socio < 7)
# n = 1,614



# 3.1.2 Data Quality Checks ---------------------------

# load dataset with data quality indicators (see 01-IER-detection)
load(here::here("./data/tidy/KARE_DQ.RData"))

# join it with dataset
KARE <-left_join(KARE, KARE_DQ, by = c("IDS"))

# drop respondents who got flagged more than once for careless responding
KARE <- subset(KARE, KARE$flag < 2)
# n = 1,604


# 3.1.3 Drops Respondents due to small cell sizes ------------

# drop "other" type of housing
# due to unclear relationship with other variables (e.g. responsibility)
# and small number of cases (n = 19)
table(KARE$R7)
KARE <- subset(KARE, R7 != "Sonstiges")
# n = 1,585

# drop diverse gender 
# due to the small number of cases (n = 1)
table(KARE$R97, useNA = "always")
KARE <- subset(KARE, R97 != "Divers" | is.na(R97))
# n = 1,584

# remove unused factor levels
KARE <- droplevels(KARE)


# 3.2 Recoding ---------------------------

# recode & transform all variables before imputing missing data
# (von Hippel 2009: transform-then-impute)


# 3.2.1 Y ---------------------------

# 10 adaptation measures
admeas <- c("R61_1", "R61_2", "R61_3", "R61_4", "R61_5", "R61_6", "R61_7", 
            "R61_8", "R61_9", "R61_10")

# 1) Number of implemented measures

# create index, treat NA as 0
score_counter <- function(row) {
  sum(row == "Implementiert", na.rm = T)    # sum up if 2 (= implemented) 
}

# count no. of implemented measures for each row
KARE$numadmeas <- apply(KARE[admeas], 1, score_counter)

# calculate number of NAs per respondent
KARE$numadmeasNA <- rowSums(is.na(KARE[admeas]))
# set respondents with 10 NAs to NA
KARE$numadmeas <-ifelse(KARE$numadmeasNA == 10, NA, KARE$numadmeas)
table(KARE$numadmeas, useNA = "always")


# 2) Yes / No
KARE$admeas <- ifelse(KARE$numadmeas > 0, 1, 0)
table(KARE$admeas, useNA = "always")



# 3.2.2 Independent Variables ---------------------------

# 3.2.2.1 Reversed polarity ---------------------------

### reverse polarity of some variables (for easier interpretation)
KARE$R18 <- forcats::fct_rev(KARE$R18)
KARE$R63_1 <- forcats::fct_rev(KARE$R63_1)
KARE$R63_2 <- forcats::fct_rev(KARE$R63_2)
KARE$R63_4 <- forcats::fct_rev(KARE$R63_4)
KARE$R63_6 <- forcats::fct_rev(KARE$R63_6)
KARE$R90_1 <- forcats::fct_rev(KARE$R90_1)
KARE$R90_5 <- forcats::fct_rev(KARE$R90_5)

# 3.2.2.2 Quasi-metric variables ---------------------------

### 6-level likert scales are assumed to be (quasi-)metric
# assumption: equidistance between categories
# considered as fulfilled as scale was only endpoint labeled 
KARE$R9 <- as.numeric(KARE$R9)
KARE$R63_1 <- as.numeric(KARE$R63_1)
KARE$R63_2 <- as.numeric(KARE$R63_2)
KARE$R63_3 <- as.numeric(KARE$R63_3)
KARE$R63_4 <- as.numeric(KARE$R63_4)
KARE$R63_6 <- as.numeric(KARE$R63_6)
KARE$R90_1 <- as.numeric(KARE$R90_1)
KARE$R90_5 <- as.numeric(KARE$R90_5)


### Midpoint estimator for binned income data
# closed bins: assigns cases to the midpoints of their bins
# bottom bin: treat bottom bin as though its lower bound is zero
# top bin: assigns cases to the harmonic mean of a Pareto distribution

# 1) visualise distribution
income_data <- data.frame (income = c(rep("single", 11), rep("multiple", 14)),
                           bin_min = c(0, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 7500,
                                       0, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 6500, 7500, 8500, 9500),
                           bin_max = c(1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 7500, NA,
                                       1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 6500, 7500, 8500, 9500, NA),
                           n_group = c(15, 48, 71, 44, 42, 16, 25, 10, 10, 4, 7,
                                       8, 12, 36, 89, 112, 143, 130, 105, 179, 119, 55, 34, 20, 50))

ggplot(subset(income_data, income == "single")) +   # single-person households
  geom_rect(aes(xmin = bin_min, xmax = bin_max,
                ymin = 0, ymax = n_group))

windows()
ggplot(subset(income_data, income == "multiple")) +   # multi-person households 
  geom_rect(aes(xmin = bin_min, xmax = bin_max,
                ymin = 0, ymax = n_group))


# 2) calculate a harmonic robust Pareto midpoint estimator (RPME) 
# for the upper boundary based on von Hippel et al. (2015: 219f.)

# single households
s_alpha <- (log(4+7) - log(7)) / (log(7500) - log(5500)) # estimation of alpha parameter for Pareto distribution
s_hb <- 7500 * (1+1/s_alpha)   # harmonic mean
KARE$inc <- NA
KARE$inc <- ifelse(KARE$R109 == 1, 500, KARE$inc)
KARE$inc <- ifelse(KARE$R109 == 2, 1250, KARE$inc)
KARE$inc <- ifelse(KARE$R109 == 3, 1750, KARE$inc)
KARE$inc <- ifelse(KARE$R109 == 4, 2250, KARE$inc) 
KARE$inc <- ifelse(KARE$R109 == 5, 2750, KARE$inc) 
KARE$inc <- ifelse(KARE$R109 == 6, 3250, KARE$inc)
KARE$inc <- ifelse(KARE$R109 == 7, 3750, KARE$inc) 
KARE$inc <- ifelse(KARE$R109 == 8, 4250, KARE$inc) 
KARE$inc <- ifelse(KARE$R109 == 9, 5000, KARE$inc) 
KARE$inc <- ifelse(KARE$R109 == 10, 6500, KARE$inc)
KARE$inc <- ifelse(KARE$R109 == 11, round(s_hb, 0), KARE$inc) 

# multiperson households
m_alpha <- (log(20+50) - log(50)) / (log(9500) - log(8500)) # estimation of alpha parameter for Pareto distribution
m_hb <- 9500 * (1+1/m_alpha) # harmonic mean
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 1, 500, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 2, 1250, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 3, 1750, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 4, 2250, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 5, 2750, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 6, 3250, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 7, 3750, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 8, 4250, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 9, 5000, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 10, 6000, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 11, 7000, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 12, 8000, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 13, 9000, KARE$inc)
KARE$inc <- ifelse(is.na(KARE$inc) & KARE$R108 == 14, round(m_hb, 0), KARE$inc)

KARE$R109 <- KARE$inc
table(KARE$R109, useNA = "always")

# rescale variable to estimate effect of 1000 € more
KARE$R109 <- KARE$R109/1000

### for income robustness check in regression models:
# 1) take into account marginal utility of income
KARE$logR109 <- log(KARE$R109 * 1000)
plot(density(KARE$logR109, na.rm = T))
# 2) income groups
# equivalised disposable income: income / equalised adults
KARE$hhsize_eq <- ifelse(KARE$R106 == 1, 1, NA)
KARE$hhsize_eq <- ifelse(KARE$R106 != 1, 1 + (KARE$R106 - 1)*0.5, KARE$hhsize_eq)
KARE$R109_eq <- (KARE$R109 * 1000) / KARE$hhsize_eq
plot(density(KARE$R109_eq, na.rm = T))
table(KARE$R109_eq, useNA = "always")
# define income groups
quantile(KARE$R109_eq, prob=c(.15,.9), na.rm = T)
KARE <- mutate(KARE,                    
               rich = case_when(
                 R109_eq < 1300 ~ "low", # 15%
                 R109_eq >= 1300 & R109_eq < 4000 ~ "middle",
                 R109_eq >= 4000 ~ "high")) # 10%
KARE$rich <- as.ordered(KARE$rich)
KARE$rich <- factor(KARE$rich, levels=c('poor', 'middle', 'rich'))
table(KARE$rich, useNA = "always")


# 3.2.2.3 Regrouping  due to small cell sizes ---------------------------

### regroup variable as some categories are too small to draw conclusions 

# education
table(KARE$R104)
KARE$R104 <-  recode_factor(as.numeric(KARE$R104), 
                            `1` = "no / lower secondary",
                            `2` = "no / lower secondary",
                            `3` = "intermediate secondary",
                            `4` = "upper secondary",
                            `5` = "upper secondary",
                            `6` = "other degree")
table(KARE$R99, KARE$R104, useNA = "always")
# most respondents with "other" degree are Germans -> set to NA and impute
levels(KARE$R104)[levels(KARE$R104)=='other degree'] <- NA


# future residence
table(KARE$R92)
KARE$R92 <-  recode_factor(as.numeric(KARE$R92), 
                           `1` = "unsure/short-term",
                           `2` = "unsure/short-term",
                           `3` = "medium-term",
                           `4` = "medium-term",
                           `5` = "long-term")
table(KARE$R92, useNA = "always")

# trust in effective protection against flooding from local government
# create dummy to increase cell sizes (esp. for tenants)
table(KARE$R80a_1, KARE$R7, useNA = "always") 
KARE$R80a_1 <-  recode_factor(as.numeric(KARE$R80a_1), 
                              `1` = "Rather yes/yes",
                              `2` = "Rather yes/yes",
                              `3` = "Rather no/no",
                              `4` = "Rather no/no")
table(KARE$R80a_1, useNA = "always")



# 3.2.3 Control Variables ---------------------------

### regroup variable as some categories are too small to draw conclusions 
# housing type
table(KARE$R69)
table(as.numeric(KARE$R69))
KARE$R69 <-  recode_factor(as.numeric(KARE$R69), 
                           `1` = "Single family home",
                           `2` = "Duplexe/terraced house",
                           `3` = "Duplexe/terraced house",
                           `4` = "Apartment building",
                           `5` = "Apartment building",
                           `6` = "Apartment building")

### high correlation between R4 & R5, code into one variable to reduce
#   multicollinearity in imputation models
# past experience
KARE$R4a5 <- "none"
KARE$R4a5 <- ifelse(KARE$R4 == 1, "experience", KARE$R4a5)
KARE$R4a5 <- ifelse(KARE$R5 == 1, "damage", KARE$R4a5)
KARE$R4a5 <- as_factor(KARE$R4a5)
table(KARE$R4a5)


### recode NAs for town into a separate category
# assumption: people from same towns are more similar (not independent)
table(KARE$town, useNA = "always")
# recode 58 respondents who did not answer in a separate category named "noresp"
KARE$town <- as.character(KARE$town)
KARE$town <-ifelse(is.na(KARE$town), "noresp", KARE$town) 
KARE$town <- as_factor(KARE$town)


# 3.2.3 Control Variables ---------------------------

### regroup housing type categories are too small to draw conclusions 
# housing type
table(KARE$R69)
table(as.numeric(KARE$R69))
KARE$R69 <-  recode_factor(as.numeric(KARE$R69), 
                           `1` = "Single family home",
                           `2` = "Duplexe/terraced house",
                           `3` = "Duplexe/terraced house",
                           `4` = "Apartment building",
                           `5` = "Apartment building",
                           `6` = "Apartment building")

# recode year of construction into age of building
KARE$R70 <- 2022 - as.numeric(KARE$R70)

# recode age based on year of birth
KARE$R98 <- 2022 - KARE$R98


# past experience
# high correlation between R4 (experience) & R5 (damage), code into one 
# variable to reduce multicollinearity in imputation models
KARE$R4a5 <- "none"
KARE$R4a5 <- ifelse(KARE$R4 == 1, "experience", KARE$R4a5)
KARE$R4a5 <- ifelse(KARE$R5 == 1, "damage", KARE$R4a5)
KARE$R4a5 <- as_factor(KARE$R4a5)
table(KARE$R4a5)


### recode 58 NAs for town into a separate category ("noresp")
# assumption: people from same towns are more similar (obs not independent)
table(KARE$town, useNA = "always")
KARE$town <- as.character(KARE$town)
KARE$town <-ifelse(is.na(KARE$town), "noresp", KARE$town) 
KARE$town <- as_factor(KARE$town)


save(KARE, file = here::here("./data/tidy/KARE_02_rec.RData"))
