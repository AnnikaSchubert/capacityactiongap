######################################################
##                                                  ##
##        Unravelling the Capacity Action Gap       ##
##                                                  ##
##                Annika Schubert                   ##
##                                                  ##
##            03-descriptive-statistics             ##
##                                                  ##
######################################################

# @brief        Descriptive statistics for variables of interest 
# @license      GNU General Public License v3.0  
# @authors      Annika Schubert

# Code for:
# Schubert, A.; von Streit, A. & Garschagen, M. (2024):
# Unravelling the capacity-action gap in flood risk adaptation. 
# In: EGUsphere [preprint], https://doi.org/X. 


# ############################################## #
#                                                #
#### Part 3: Descriptive Statistics           ####
#                                                #
# ############################################## #

library(here)                # for file referencing
library(ggplot2)             # plots
library(dplyr)               # data transformation, filtering and summarising
library(missMethods)         # compute median of ordered factors
library(ComplexHeatmap)      # create upset plot
library(readxl)              # load excel files


load(here::here("./data/tidy/KARE_02_rec.RData"))

# exclude respondents with NAs on Y
KARE_stats <- subset(KARE, !is.na(KARE$numadmeas))


# 3.1 Compare sociodemographics of survey with census data (Appendix Tab. A1) ---------------------------

# gender
table(KARE_stats$R97) / length(KARE_stats$R97[!is.na(KARE_stats$R97)])
# German nationality
table(KARE_stats$R99) / length(KARE_stats$R99[!is.na(KARE_stats$R99)])

# education
table(KARE_stats$R104) / length(KARE_stats$R99[!is.na(KARE_stats$R104)])

# age
table(KARE_stats$R98)
KARE_stats <- mutate(KARE_stats,                    
                     age_cat = case_when(
                       R98 %in% 19:24 ~ "19-24",
                       R98 %in% 25:44 ~ "25-44",
                       R98 %in% 45:64 ~ "45-64",
                       R98 %in% 65:93 ~ "65-93"))
table(KARE_stats$age_cat) / length(KARE_stats$age_cat[!is.na(KARE_stats$age_cat)])

# income
table(KARE_stats$R109)
KARE_stats <- mutate(KARE_stats,                    
                     inc_cat = case_when(
                       R109 < 1500 ~ "low", 
                       R109 >= 1500 & R109 < 4000 ~ "medium",
                       R109 >= 4000 ~ "high")) 
table(KARE_stats$inc_cat) / length(KARE_stats$inc_cat[!is.na(KARE_stats$inc_cat)])


# 3.2 Overview: adaptive capacity indicators (Table 2) ---------------------------

### full sample
KARE_stats %>%
  select(R109, R71, R98, R92, R106, R104, R91, R9, R63_3, R7, R75a,
         R63_6, R63_1, R63_2, R90_1, R90_5,R80a_1, R18, R63_4,
         R97, R99, R4a5, R70, R69, Modus) %>%
  summary()

KARE_stats %>%
  select(R109, R71, R98, R106, R91, R9, R63_3, R7, R75a,
         R63_6, R63_1, R63_2, R90_1, R90_5, R80a_1, R18, R63_4,
         R97, R99, R4a5, R70, R69, Modus) %>%
  psych::describe()

median(as.ordered(KARE_stats$R92), na.rm = T)
median(as.ordered(KARE_stats$R104), na.rm = T)
median(as.ordered(KARE_stats$R18), na.rm = T)

### for owner 
KARE_stats %>%
  select(R109, R71, R98, R92, R106, R104, R91, R9, R63_3, R7, R75b,
         R63_6, R63_1, R63_2, R90_1, R90_5,R80a_1, R18, R63_4,
         R97, R99, R4a5, R70, R69, Modus) %>%
  filter(KARE_stats$R7 == "Eigentum") %>%
  summary()

KARE_stats %>%
  select(R109, R71, R98, R106, R91, R9, R63_3, R75b,
         R63_6, R63_1, R63_2, R90_1, R90_5, R80a_1, R18, R63_4,
         R97, R99, R4a5, R70, R69, Modus) %>%
  filter(KARE_stats$R7 == "Eigentum") %>%
  psych::describe()

median(as.ordered(KARE_stats$R104[KARE_stats$R7 == "Eigentum"]), na.rm = T)
median(as.ordered(KARE_stats$R92[KARE_stats$R7 == "Eigentum"]), na.rm = T)
median(as.ordered(KARE_stats$R18[KARE_stats$R7 == "Eigentum"]), na.rm = T)


### for tenants 
KARE_stats %>%
  select(R109, R71, R98, R92, R106, R104, R91, R9, R63_3, R7, R75a,
         R63_6, R63_1, R63_2, R90_1, R90_5,R80a_1, R18, R63_4,
         R97, R99, R4a5, R70, R69, Modus) %>%
  filter(KARE_stats$R7 == "Miete") %>%
  summary()

KARE_stats %>%
  select(R109, R71, R98, R106, R91, R9, R63_3, R75a,
         R63_6, R63_1, R63_2, R90_1, R90_5, R80a_1, R18, R63_4,
         R97, R99, R4a5, R70, R69, Modus) %>%
  filter(KARE_stats$R7 == "Miete") %>%
  psych::describe()

median(as.ordered(KARE_stats$R104[KARE_stats$R7 == "Miete"]), na.rm = T)
median(as.ordered(KARE_stats$R92[KARE_stats$R7 == "Miete"]), na.rm = T)
median(as.ordered(KARE_stats$R18[KARE_stats$R7 == "Miete"]), na.rm = T)


# 3.3 Overview: adaptation actions (Fig. 4) ---------------------------

# 3.3.1 Summary Statistics ---------------------------

# recode adaptive measures to dummies for mean, sd
KARE_stats$R61_1_d <- ifelse(KARE_stats$R61_1 == "Implementiert", 1, 0)
KARE_stats$R61_2_d <- ifelse(KARE_stats$R61_2 == "Implementiert", 1, 0)
KARE_stats$R61_3_d <- ifelse(KARE_stats$R61_3 == "Implementiert", 1, 0)
KARE_stats$R61_4_d <- ifelse(KARE_stats$R61_4 == "Implementiert", 1, 0)
KARE_stats$R61_5_d <- ifelse(KARE_stats$R61_5 == "Implementiert", 1, 0)
KARE_stats$R61_6_d <- ifelse(KARE_stats$R61_6 == "Implementiert", 1, 0)
KARE_stats$R61_7_d <- ifelse(KARE_stats$R61_7 == "Implementiert", 1, 0)
KARE_stats$R61_8_d <- ifelse(KARE_stats$R61_8 == "Implementiert", 1, 0)
KARE_stats$R61_9_d <- ifelse(KARE_stats$R61_9 == "Implementiert", 1, 0)
KARE_stats$R61_10_d <- ifelse(KARE_stats$R61_10 == "Implementiert", 1, 0)

### full sample
KARE_stats %>%
  select(numadmeas, admeas, R61_1_d, R61_2_d, R61_3_d, R61_4_d, 
         R61_5_d, R61_6_d, R61_7_d, R61_8_d, R61_9_d,R61_10_d) %>%
  psych::describe()

KARE_stats %>%
  select(numadmeas, admeas, R61_1_d, R61_2_d, R61_3_d, R61_4_d, 
         R61_5_d, R61_6_d, R61_7_d, R61_8_d, R61_9_d,R61_10_d) %>%
  summary()


### for owner 
KARE_stats %>%
  select(numadmeas, admeas, R61_1_d, R61_2_d, R61_3_d, R61_4_d, 
         R61_5_d, R61_6_d, R61_7_d, R61_8_d, R61_9_d, R61_10_d) %>%
  filter(KARE_stats$R7 == "Eigentum") %>%
  psych::describe()

KARE_stats %>%
  select(numadmeas, admeas, R61_1_d, R61_2_d, R61_3_d, R61_4_d, 
         R61_5_d, R61_6_d, R61_7_d, R61_8_d, R61_9_d,R61_10_d) %>%
  filter(KARE_stats$R7 == "Eigentum") %>%
  summary()


### and tenants
KARE_stats %>%
  select(numadmeas, admeas, R61_2_d, R61_4_d, R61_7_d, R61_10_d) %>%
  filter(KARE_stats$R7 == "Miete") %>%
  psych::describe()

KARE_stats %>%
  select(numadmeas, admeas, R61_1_d, R61_2_d, R61_3_d, R61_4_d, 
         R61_5_d, R61_6_d, R61_7_d, R61_8_d, R61_9_d,R61_10_d) %>%
  filter(KARE_stats$R7 == "Miete") %>%
  summary()

## Two-Proportions Z-Test - Owners vs. tenants
table(KARE_stats$R7, KARE_stats$admeas)
prop.test(x = c(1049, 200), n = c(1157, 414),
          alternative = "two.sided",
          correct = T)


# distribution: adaptive actions
ggplot(data = KARE_stats, aes(x = admeas))+
  geom_bar(fill = "#1C71A6") 

ggplot(data = KARE_stats, aes(x = numadmeas))+
  geom_bar(fill = "#1C71A6")
quantile(KARE_stats$numadmeas, c(0.25, 0.5, 0.75), type = 1)
# skew: 0.92 --> slightly right skewed

# owner
ggplot(subset(KARE_stats, R7 == "Eigentum"), aes(x = numadmeas)) +
  geom_bar(fill = "#1C71A6") +
  labs(title = "Distribution: No. of implemeted measures by homeowners")
quantile(KARE_stats$numadmeas[KARE_stats$R7 == "Eigentum"], c(0.25, 0.5, 0.75), type = 1)
mean(KARE_stats$numadmeas[KARE_stats$R7 == "Eigentum"])
table(KARE_stats$numadmeas[KARE_stats$R7 == "Eigentum"]) / 1157

# tenants
ggplot(subset(KARE_stats, R7 == "Miete"), aes(x = numadmeas)) +
  geom_bar(fill = "#1C71A6") +
  labs(title = "Distribution: No. of implemeted measures by tenants")
quantile(KARE_stats$numadmeas[KARE_stats$R7 == "Miete"], c(0.25, 0.5, 0.75), type = 1)
mean(KARE_stats$numadmeas[KARE_stats$R7 == "Miete"])
table(KARE_stats$numadmeas[KARE_stats$R7 == "Miete"]) / 417


# 3.3.2 Combinations of measures - Upset plot (Fig. 4) ---------------------------

# create dataframe with private adaptation measures
measures <- c("R61_1_d", "R61_2_d", "R61_3_d", "R61_4_d", "R61_5_d", "R61_6_d", 
              "R61_7_d", "R61_8_d", "R61_9_d", "R61_10_d", "R7")
KARE_measures <- KARE_stats[names(KARE_stats) %in% measures]

# set NAs to 0
KARE_measures[is.na(KARE_measures)] <- 0

# subset data
KARE_owner <- subset(KARE_measures, R7 == "Eigentum")
KARE_owner <- KARE_owner[,-1]
KARE_tenants <- subset(KARE_measures, R7 == "Miete")
KARE_tenants <- KARE_tenants[,-c(1:2, 4, 6:7, 9:10)] # drop mbuilding measures


# 3.3.2.1 Property owners Upset plot (Fig. 4a) ---------------------------

owner_m <- make_comb_mat(KARE_owner)
rownames(owner_m) <- c("Backflow preventer", "Mobile water barriers", "Flood-proof windows/doors",
                       "Information seeking", "Flood-proof heating system", "Water-resistant plaster", "Insurance",
                       "Water-resistant flooring", "Structural adjuments", "Moved furniture in higher floors")
owner_m <- owner_m[comb_size(owner_m) >= 11] # display only combinations with n >= 11 (1 % of tenants)

# change absolute to relative frequency
ss_owner <- set_size(owner_m)
cs_owner <-  comb_size(owner_m)
rel_comb_size <- (cs_owner/nrow(KARE_owner)) * 100
rel_comb_size <- round(rel_comb_size, 0) # in %
rel_set_size <- (ss_owner/nrow(KARE_owner)) * 100
rel_set_size <- round(rel_set_size, 0) # in %

sum(cs_owner) / nrow(KARE_owner)
# other combinations: 55 %

lab_cap_own <- expression(paste("Source: own calculation, data from KARE household survey 2022 (", n[Owners], " = 1,157)"))

# save plot as pdf
Cairo::Cairo(
  38.5, 17,
  units = "in",
  file = "D:/LRZ Sync+Share/Transfer_Hiwi (Anne von Streit)/Luna Ciesla/Grafik/Grafik 2/Plots_final/measures_owner.pdf",
  bg = "transparent",
  type = "pdf",
  dpi = 300
)

# plot
UpSet(owner_m, # display only intersections with n > 10
      bg_col = "transparent",  # set background behind set combinations transparent
      bg_pt_col = "transparent", # change unmarked dots to transparent
      row_title = "Implemented measures", # add label for measures
      comb_order = order(comb_size(owner_m), decreasing=F), # sort by decreasing frequency 
      top_annotation = HeatmapAnnotation("Percentage of each combination type" = anno_barplot(rel_comb_size, 
                                                                                              border = FALSE,  # no border line
                                                                                              height = unit(1.6, "cm"),
                                                                                              add_numbers = TRUE, # add numbers to bars
                                                                                              gp = gpar(col = "coral3", fill = "coral3"),  # change color of bars
                                                                                              annotation_name_side = "right",   # label at top
                                                                                              axis = F,                   # remove axis
                                                                                              bar_width = 0.05)), 
      right_annotation = rowAnnotation("Percentage of households which\nimplemented the measure" = anno_barplot(rel_set_size,
                                                                                                                border = FALSE,  # no border line
                                                                                                                add_numbers = TRUE,
                                                                                                                gp = gpar(col = "cadetblue", fill = "cadetblue"), # change color of bars
                                                                                                                axis = F,              # remove axis
                                                                                                                bar_width = 0.05,# change width of bars
                                                                                                                width = unit(6, "cm")), 
                                       annotation_name_side = "top",   # label at top
                                       annotation_name_rot = 0)       # rotate label
)

invisible(dev.off())

# 3.3.2.2 Tenants Upset plot (Fig. 4a) ---------------------------

# create combination matrix
mt <- make_comb_mat(KARE_tenants)
rownames(mt) <- c("Mobile water barriers", "Information seeking", "Insurance",
                  "Moved furniture in higher floors")

# display only combinations with n >= 5
tenants_m <- mt[comb_size(mt) >= 1]


# change absolute to relative frequency
ss_tenants <- set_size(tenants_m)
cs_tenants <-  comb_size(tenants_m)
rel_comb_size_ten <- (cs_tenants/nrow(KARE_tenants)) * 100
rel_comb_size_ten <- round(rel_comb_size_ten, 0) # in %
rel_set_size_ten <- (ss_tenants/nrow(KARE_tenants)) * 100
rel_set_size_ten <- round(rel_set_size_ten, 0) # in %

sum(cs_tenants) / nrow(KARE_tenants)
# other combinations: 0 % --> all covered

lab_cap_ten <- expression(paste("Source: own calculation, data from KARE household survey 2022 (", n[Tenants], " = 414)"))

# save plot as pdf
Cairo::Cairo(
  38.5, 13,
  units = "in",
  file = "D:/LRZ Sync+Share/Transfer_Hiwi (Anne von Streit)/Luna Ciesla/Grafik/Grafik 2/Plots_final/measures_tenants.svg",
  bg = "transparent",
  type = "svg", # svg
  dpi = 300
)

UpSet(tenants_m, # display only intersections with n > 10
      bg_col = "transparent",  # set background behind set combinations transparent
      bg_pt_col = "transparent", # change unmarked dots to transparent
      row_title = "Implemented measures", # add label for measures
      comb_order = order(comb_size(tenants_m), decreasing=F), # sort by decreasing frequency 
      top_annotation = HeatmapAnnotation("Percentage of each combination type" = anno_barplot(rel_comb_size_ten, 
                                                                                              border = FALSE,  # no border line
                                                                                              height = unit(4, "cm"),
                                                                                              add_numbers = TRUE, # add numbers to bars
                                                                                              gp = gpar(col = "coral3", fill = "coral3"),  # change color of bars
                                                                                              annotation_name_side = "right",   # label at top
                                                                                              axis = F,                   # remove axis
                                                                                              bar_width = 0.05)), 
      right_annotation = rowAnnotation("Percentage of households which\nimplemented the measure" = anno_barplot(rel_set_size_ten,
                                                                                                                border = FALSE,  # no border line
                                                                                                                add_numbers = TRUE,
                                                                                                                gp = gpar(col = "cadetblue", fill = "cadetblue"), # change color of bars
                                                                                                                axis = F,              # remove axis
                                                                                                                bar_width = 0.05,# change width of bars
                                                                                                                width = unit(2.5, "cm")), 
                                       annotation_name_side = "top",   # label at top
                                       annotation_name_rot = 0)       # rotate label
)

invisible(dev.off())


# 3.4 Correlation / Measures of Association ---------------------------

### type of correlation depends on level of measurement:
# dummy - dummy --> Phi
# dummy - ordinal --> Cramer's V
# dummy - nominal --> Cramer's V
# metric - metric --> Pearson cor 
# metric - ordinal --> Spearman's rank
# metric - dummy --> Point biserial
# metric - nominal --> Cramer's V

KARE_stats_own <- subset(KARE_stats, KARE_stats$R7 == "Eigentum")
KARE_stats_ten <- subset(KARE_stats, KARE_stats$R7 == "Miete")

# 3.4.1 Adaptation yes/no ---------------------------

# binary - binary: Phi coefficient (== Pearson Corr)
psych::phi(table(KARE_stats$admeas, KARE_stats$R7))
psych::phi(table(KARE_stats$admeas, KARE_stats$R80a_1))
cor(as.numeric(KARE_stats$admeas), as.numeric(KARE_stats$R80a_1),  use = "complete.obs")
cor(as.numeric(KARE_stats$admeas), as.numeric(KARE_stats$R75b),  use = "complete.obs") # owner

# binary - ordinal: Cramer's V
DescTools::CramerV(table(KARE_stats$R104, KARE_stats$admeas))
DescTools::CramerV(table(KARE_stats$R92, KARE_stats$admeas))
DescTools::CramerV(table(KARE_stats$R18, KARE_stats$admeas))
DescTools::CramerV(table(KARE_stats$R4a5, KARE_stats$admeas))
DescTools::CramerV(table(KARE_stats$R75a, KARE_stats$admeas))

# metric - binary: Point-Biserial 
KARE_stats$admeas_m <- as.numeric(KARE_stats$admeas)
metricvars <- c("admeas_m", "R109", "R71", "R91", "R90_1", "R90_5", "R9",
                "R63_6", "R63_1", "R63_3", "R63_2", "R63_4")
round(cor(KARE_stats[metricvars], use = "pairwise.complete.obs"), 4)


# 3.4.2 No. of implemented measures ---------------------------

# nominal - metric: Cramer's V
DescTools::CramerV(table(KARE_stats$R75a, KARE_stats$numadmeas))

# ordinal - metric: Spearman's Rank Correlation Coefficient 
cor(as.numeric(KARE_stats$R104), KARE_stats$numadmeas, 
    method="spearman", use = "complete.obs")
cor(as.numeric(KARE_stats$R92), KARE_stats$numadmeas, 
    method="spearman", use = "complete.obs")
cor(as.numeric(KARE_stats$R18), KARE_stats$numadmeas, 
    method="spearman", use = "complete.obs")
cor(as.numeric(KARE_stats$R4a5), KARE_stats$numadmeas, 
    method="spearman", use = "complete.obs")

# metric - metric
metricvars <- c("numadmeas", "R109", "R71", "R91", "R90_1", "R90_5", "R9",
                "R63_6", "R63_1", "R63_3", "R63_2", "R63_4")
round(cor(KARE_stats[metricvars], use = "pairwise.complete.obs"), 4)

# metric - binary: Point-Biserial 
cor(as.numeric(KARE_stats$R7), KARE_stats$numadmeas, use = "complete.obs")
cor(as.numeric(KARE_stats$R80a_1), KARE_stats$numadmeas, use = "complete.obs")
cor(as.numeric(KARE_stats$R75b), KARE_stats$numadmeas, use = "complete.obs") # owner



# 3.4.3 Correlation Heatmap (Fig. 5) ---------------------------

# load excel with correlations
cormat <- readxl::read_excel(here::here("./data/tidy/Cor_Heatmap.xlsx"),
                             col_names = FALSE, col_types = "numeric", sheet = 2)

# reshape dataframe into matrix
cormat <- data.matrix(cormat)

# variable names as rownames
rownames(cormat) <- c("Competing concerns", "Motivation","Self-efficacy", 
                      "Public protection is sufficient", "Trust in authorities",
                      "Expectation in authorities", "Main responsibility", 
                      "Previous experience", "Risk perception", "Future risk perception",
                      "Social cohesion", "Social network", "Place attachment", 
                      "Duration of residence", "Property ownership", "Living area",
                      "Income","Education")

# reshape back into data frame
melted_cor_mat <- reshape2::melt(cormat)

# create label with subscript for caption
lab_cap <- expression(paste("Source: own calculation, data from KARE household survey 2022 (n = 1,571; ", 
                            n[Owners], " = 1,157; ", n[Tenants], " = 414)"))

# plot
windows()
ggplot(data = melted_cor_mat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#A63446", high = "#0C6291", mid = "#FBFEF9", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  labs(x = NULL, y = NULL, title = "Correlation heatmap", subtitle = "", 
       caption = lab_cap) + 
  theme_minimal()+ 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.subtitle = element_text(size=50),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  geom_text(aes(Var2, Var1, label = round(value,2)), 
            color = "black", size = 4) +
  coord_fixed(ratio = 0.2) 

# save plot
ggsave(filename = here::here("./output/cor_heatmap.svg"), dpi = 300)
ggsave(filename = here::here("./output/cor_heatmap.svg"))
