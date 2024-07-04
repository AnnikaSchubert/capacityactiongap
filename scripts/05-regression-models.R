######################################################
##                                                  ##
##        Unravelling the Capacity Action Gap       ##
##                                                  ##
##                Annika Schubert                   ##
##                                                  ##
##             05-regression-models                 ##
##                                                  ##
######################################################

# @brief        Run logistic and Poisson regression models on imputed datasets
# @license      GNU General Public License v3.0  
# @authors      Annika Schubert

# Code for:
# Schubert, A.; von Streit, A. & Garschagen, M. (2024):
# Unravelling the capacity-action gap in flood risk adaptation. 
# In: EGUsphere [preprint], https://doi.org/X. 


# ############################################## #
#                                                #
#### Part 5: Regression Analysis              ####
#                                                #
# ############################################## #

library(here)                # for file referencing
library(haven)               # to work with the SPSS labbelled dataset
library(tidyverse)           # for data manipulation
library(mice)                # for data preparation
library(sandwich)            # cluster-robust standard error
library(marginaleffects)     # for AME calculation
library(DHARMa)              # residual diagnostics for generalized linear (mixed) models 
library(stargazer)           # to produce tables
library(ggplot2)             # to create forest plots
library(cowplot)             # to arrange plots 


# 5.1 Data Preparation ---------------------------

# load data and rename dataset
load(here::here("./data/tidy/KARE_imp_inck.RData"))
KARE_MI <- KARE_imp

### MID (multiple imputation, then deletion): 
# use Y to impute X, but drop obs with missing Y from the analysis (von Hippel 2007)

# create a long format dataset  
KARE_MI_long <- complete(KARE_MI, action = "long", include = TRUE)

# drop obs with NA on Y
KARE_MI_long <- subset(KARE_MI_long, !(IDS %in% KARE_MI_long$IDS[is.na(KARE_MI_long$numadmeas)]))

# recode R75 from imputed R75a (tenants - polyreg) and R75b (owners - logreg)
table(KARE_MI_long$R7, KARE_MI_long$R75a, useNA = "always")
table(KARE_MI_long$R7, KARE_MI_long$R75b)
KARE_MI_long$R75 <- ifelse(KARE_MI_long$R7 == "Miete", KARE_MI_long$R75a, 
                           ifelse(KARE_MI_long$R75b == "Öffentliche Stellen (d.h. Bund, Länder oder Gemeinden)", 3, 1))
KARE_MI_long$R75 <- as_factor(KARE_MI_long$R75)
table(KARE_MI_long$R7, KARE_MI_long$R75, useNA = "always")

KARE_MI <- mice::as.mids(KARE_MI_long)
# n = 1,571 obs


### turn ordered factors into unordered to get pairwise comparisons 
# instead of polynomial contrasts
KARE_MI$data$R18 <- factor(KARE_MI$data$R18,
                           ordered = F)
KARE_MI$data$R80a_1 <- factor(KARE_MI$data$R80a_1,
                              ordered = F)
KARE_MI$data$R92 <- factor(KARE_MI$data$R92,
                           ordered = F)
KARE_MI$data$R104 <- factor(KARE_MI$data$R104,
                            ordered = F)
KARE_MI$data$Modus <- factor(KARE_MI$data$Modus,
                             ordered = F)

### change reference levels
KARE_MI$data$Modus <- relevel(KARE_MI$data$Modus, ref = "2")
KARE_MI$data$R92 <- relevel(KARE_MI$data$R92, ref = "unsure/short-term")
KARE_MI$data$R18 <- relevel(KARE_MI$data$R18, ref = "Überhaupt nicht wahrscheinlich")
KARE_MI$data$R4a5 <- relevel(KARE_MI$data$R4a5, ref = "none")
KARE_MI$data$R7 <- relevel(KARE_MI$data$R7, ref = "Miete")
KARE_MI$data$R80a_1 <- relevel(KARE_MI$data$R80a_1, ref = "Rather no/no")


### change display of results
options(scipen=999)
options(digits = 4)


# 7.2 Log Reg: ACs --> Adaptation action yes/no ---------------------------

# binary outcome --> binary logistic regression

# 7.2.1 Logistic Reg with cluster-robust SEs - full sample ---------------------------

# fit logistic regression to each imputed data set 
glm_multiimp <- with(KARE_MI,
                     glm(admeas ~ R104 + R109 + R71 +    # generic capacity
                           R7 + R91 + R92 +
                           R90_1 + R90_5 +
                           R9 + R18 + R4a5 + R75 +       # flood-specific capacity
                           R63_6 + R80a_1 + R63_1 +
                           R63_3 + R63_2 + R63_4 +
                           R97 + R98 + R99 + R106 +      # CVs: gender, age, mig, hhsize
                           R69 + R70 +                   # CVs: house
                           Modus,                        # Survey design
                         family = binomial("logit")))

# calculate AME
marginal <- avg_slopes(glm_multiimp, vcov = ~town)

# save results
marginal$order <- c(33, 32, 1, 2, 28, 3, 12, 13, 14, 15, 16, 21, 
                    23, 22, 24, 19, 30, 29, 5, 31, 4, 17, 18, 20, 11,
                    9, 10, 6, 8, 7, 25, 26, 27)
marginal <- marginal[order(marginal$order, decreasing = F),]
b_log_AME <- c(NA, marginal$estimate)
se_log_AME <- c(NA, marginal$std.error)
p_log_AME <- c(NA, marginal$p.value)


### Goodness of fit measures
#R^2
Nagelkerke_log_all <- numeric(30)
for (i in 1:30) {
  Nagelkerke_log_all_res <- c(DescTools::PseudoR2(glm_multiimp$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_log_all[i] <- Nagelkerke_log_all_res
}
median(Nagelkerke_log_all)

# BIC
BIC_log_all <- numeric(30)
for (i in 1:30) {
  BIC_log_all_res <- c(BIC(glm_multiimp$analyses[[i]]))
  BIC_log_all[i] <- BIC_log_all_res
}
median(BIC_log_all)


# check assumptions:
# 1) Linearity 
plot(glm_multiimp$analyses[[7]], 1)
# not useful to plot the raw residuals, examine binned residual plots
arm::binnedplot(fitted(glm_multiimp$analyses[[7]]), 
                residuals(glm_multiimp$analyses[[7]], type = "response"), 
                nclass = NULL, 
                xlab = "Expected Values", 
                ylab = "Average residual", 
                main = "Binned residual plot", 
                cex.pts = 0.8, 
                col.pts = 1, 
                col.int = "gray")
KARE_data6 <- complete(KARE_MI, action=6)
arm::binnedplot(KARE_data6$R109, 
                residuals(glm_multiimp$analyses[[6]], type = "response"), 
                nclass = NULL, 
                xlab = "Expected Values", 
                ylab = "Average residual", 
                main = "Binned residual plot", 
                cex.pts = 0.8, 
                col.pts = 1, 
                col.int = "gray")


# 2) Independence 
# --> Respondents from the same town might be more similar to one another on the 
#     outcome measure, on average, than they are with respondents across towns

# plot residuals vs spatial variable
glm_data6 <- glm(data = KARE_data6,
                 admeas ~ R104 + R109 + R71 +    # generic capacity
                   R7 + R91 + R92 +
                   R90_1 + R90_5 +
                   R9 + R18 + R4a5 + R75 +       # flood-specific capacity
                   R63_6 + R80a_1 + R63_1 +
                   R63_3 + R63_2 + R63_4 +
                   R97 + R98 + R99 + R106 +      # CVs: gender, age, mig, hhsize
                   R69 + R70 +                   # CVs: house
                   Modus,                        # Survey design
                 family = binomial("logit"))                      

# by town
ggplot(data = data.frame(x = KARE_data6$town, y = glm_data6$residuals), aes(x = x, y = y)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(-20, 20))
# some patterns --> use cluster-robust SEs
# maybe also due to small cluster sizes (e.g. 2, 5, 10 respondents)


# 3) No multicollinearity 
car::vif(glm_multiimp$analyses[[10]]) 
# GVIF (1/(2*Df) < 2 for all --> no multicollinearity


# 4) Exogeneity 
# --> all relevant CVs are included in the model


# 7.2.2 Logistic Reg with cluster-robust SEs - Owners ---------------------------

# fit logistic regression to each imputed data set 
glm_multiimp_own <- with(KARE_MI,
                         glm(admeas ~ R104 + R109 + R71 +   # generic capacity
                               R91 + R92 +
                               R90_1 + R90_5 +
                               R9 + R18 + R4a5 + R75 +      # flood-specific capacity
                               R63_6 + R80a_1 + R63_1 +
                               R63_3 + R63_2 + R63_4 +
                               R97 + R98 + R99 + R106 +     # CVs: gender, age, mig, hhsize
                               R69 + R70 +                  # CVs: house
                               Modus,                       # Survey design
                             family = binomial("logit"),
                             subset = (R7 == "Eigentum")))

# calculate AME
marginal_own <- avg_slopes(glm_multiimp_own, vcov = ~town)

# save results
marginal_own$order <- c(33, 32, 1, 2, 28, 3, 12, 13, 14, 15, 16, 21, 
                        23, 22, 24, 19, 30, 29, 31, 4, 18, 20, 11,
                        9, 10, 6, 8, 7, 25, 26, 27)
marginal_own <- marginal_own[order(marginal_own$order, decreasing = F),]
b_log_AME_own <- c(NA, marginal_own$estimate)
se_log_AME_own <- c(NA, marginal_own$std.error)
p_log_AME_own <- c(NA, marginal_own$p.value)


### Goodness of fit measures
#R^2
Nagelkerke_log_own <- numeric(30)
for (i in 1:30) {
  Nagelkerke_log_own_res <- c(DescTools::PseudoR2(glm_multiimp_own$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_log_own[i] <- Nagelkerke_log_own_res
}
median(Nagelkerke_log_own)

# BIC
BIC_log_own <- numeric(30)
for (i in 1:30) {
  BIC_log_own_res <- c(BIC(glm_multiimp_own$analyses[[i]]))
  BIC_log_own[i] <- BIC_log_own_res
}
median(BIC_log_own)


# 7.2.3 Logistic Reg with cluster-robust SEs - Tenants ---------------------------

# fit logistic regression to each imputed data set 
glm_multiimp_ten <- with(KARE_MI,
                         glm(admeas ~ R104 + R109 + R71 +        # generic capacity
                               R91 + R92 +
                               R90_1 + R90_5 +
                               R9 + R18 + R4a5 + R75 +    # flood-specific capacity
                               R63_6 + R80a_1 + R63_1 +
                               R63_3 + R63_2 + R63_4 +
                               R97 + R98 + R99 + R106 +   # CVs: gender, age, mig, hhsize
                               R69 + R70 +                # CVs: house
                               Modus,                     # Survey design
                             family = binomial("logit"),
                             subset = (R7 == "Miete")))

# calculate AME
marginal_ten <- avg_slopes(glm_multiimp_ten, vcov = ~town)

# save results
marginal_ten$order <- c(33, 32, 1, 2, 28, 3, 12, 13, 14, 15, 16, 21, 
                        23, 22, 24, 19, 30, 29, 31, 4, 17, 18, 20, 11,
                        9, 10, 6, 8, 7, 25, 26, 27)
marginal_ten <- marginal_ten[order(marginal_ten$order, decreasing = F),]
b_log_AME_ten <- c(NA, marginal_ten$estimate)
se_log_AME_ten <- c(NA, marginal_ten$std.error)
p_log_AME_ten <- c(NA, marginal_ten$p.value)


### Goodness of fit measures
#R^2
Nagelkerke_log_ten <- numeric(30)
for (i in 1:30) {
  Nagelkerke_log_ten_res <- c(DescTools::PseudoR2(glm_multiimp_ten$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_log_ten[i] <- Nagelkerke_log_ten_res
}
median(Nagelkerke_log_ten)

# BIC
BIC_log_ten <- numeric(30)
for (i in 1:30) {
  BIC_log_ten_res <- c(BIC(glm_multiimp_ten$analyses[[i]]))
  BIC_log_ten[i] <- BIC_log_ten_res
}
median(BIC_log_ten)


# 7.2.4 Create table (App. A2) ---------------------------
 
stargazer(log_model, log_model_own, log_model_ten,                # use complete case models as model objects as list objects are not supported
          coef = list(b_log_AME, b_log_AME_own, b_log_AME_ten),   # but display coefficients etc. from imputed models 
          se = list(se_log_AME, se_log_AME_own, se_log_AME_ten),
          p = list(p_log_AME, p_log_AME_own, p_log_AME_ten),
          digits = 4,
          report=('vc*s'),
          title="Logistic Regressions",
          type = "html",
          out = here::here("./output/LogReg_inck.doc"))


# 7.2.5 Create forest plot (Fig.6) ---------------------------

# label variables
marginal_own <- marginal_own %>%
  mutate(label = c("Education - intermediate secondary", "Education - upper secondary", "Income", 
                   "Living area", "Duration of residence","Place attachment - medium-term", "Place attachment - long-term",
                   "Social network", "Social cohesion", "Future risk perception",
                   "Risk perception - not likely", "Risk perception - likely", "Risk perception - very likely",
                   "Previous experience - damage", "Previous experience - experience",  
                   "Main responsibility state", "Expecation in authorities", "Trust in authorities - yes", 
                   "Public protection is sufficient", "Self-efficacy", "Motivation", 
                   "Competing concerns",
                   NA, NA, NA,  NA, NA, NA, NA, NA, NA ), 
         .after = term)

marginal_ten <- marginal_ten %>%
  mutate(label = c("Education - intermediate secondary", "Education - upper secondary", "Income", 
                   "Living area", "Duration of residence","Place attachment - medium-term", "Place attachment - long-term",
                   "Social network", "Social cohesion", "Future risk perception",
                   "Risk perception - not likely", "Risk perception - likely", "Risk perception - very likely",
                   "Previous experience - damage", "Previous experience - experience", "Main responsibility landlord",  
                   "Main responsibility state", "Expecation in authorities", "Trust in authorities - yes", 
                   "Public protection is sufficient", "Self-efficacy", "Motivation", 
                   "Competing concerns",
                   NA, NA, NA,  NA, NA, NA, NA, NA, NA ), 
         .after = term)

# drop control variables
marginal_own2 <- marginal_own %>% drop_na(label) 
marginal_ten2 <- marginal_ten %>% drop_na(label)

# merge results owners & tenants into one dataframe
marginal_own_ten <- merge(marginal_ten2, marginal_own2, by = "label", all = T)
marginal_own_ten <- marginal_own_ten[order(marginal_own_ten$estimate.y, decreasing = TRUE),]
# for plotting: add value for missing landlord which is not within plotted boundaries
marginal_own_ten$estimate.y[is.na(marginal_own_ten$estimate.y)] <- -2


# change colors 
colors_log <- c("#588BAE", "#588BAE", "#588BAE", "#588BAE","#588BAE", 
                "black", "#588BAE", "black", "black", "black",
                "#588BAE", "#588BAE", "black", "black", "#588BAE",
                "black", "black", "black",  "#588BAE","#588BAE",
                "#588BAE", "#588BAE", "#588BAE")
colors_log_label <- c("#588BAE", "#588BAE", "#588BAE", "#588BAE", "#588BAE",
                      "black", "black", "black", "#588BAE", "black", 
                      "black", "#588BAE", "#588BAE", "black", "black", 
                      "black", "#588BAE", "black", "#588BAE", "#588BAE", 
                      "#588BAE", "#588BAE","#588BAE")
# add source label
lab_cap <- expression(paste("Source: own calculation, data from KARE household survey 2022 (",
                            n[Owners], " = 1,157; ", n[Tenants], " = 414)"))

# plot owners
log_own <- ggplot(marginal_own_ten, aes(y = reorder(label, estimate.y), x = estimate.y, xmin = conf.low.y, xmax = conf.high.y)) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "#758DA3", size=0.5) +
  geom_pointrange(color = colors_log) +
  labs(x ="AME", y = "",
       title = "Owners", caption = " ") + 
  xlim(-0.3, 0.6) + 
  theme_minimal() +
  theme(axis.text.y = element_text(colour = colors_log_label, size = 14),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

# plot tenants
log_ten <- ggplot(marginal_own_ten, aes(y = reorder(label, estimate.y), x = estimate.x, xmin = conf.low.x, xmax = conf.high.x)) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "#758DA3", size=0.5) +
  geom_pointrange(color = colors_log) + 
  labs(x ="AME", y = "",
       title = "Tenants", caption = lab_cap) + 
  xlim(-0.3, 0.6) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

# create title
title_log <- ggdraw() + 
  draw_label(
    "Effect sizes of adaptive capacity indicators explaing household adaptation (yes/no)",
    fontface = 'bold', size = 16,x = 0, hjust = 0)  +
  theme(plot.margin = margin(0, 0, 0, 7))


plot_log <- plot_grid(log_own, log_ten, rel_widths = c(1.85,1.12))
windows()
plot_grid(title_log, plot_log,ncol = 1, rel_heights = c(0.1, 1))

ggsave(filename = here::here("./output/logplot.pdf"), dpi = 300)
ggsave(filename = here::here("./output/logplot.svg"))



# 7.3 Poisson Reg: ACs --> Adaptation Action No. ---------------------------

# count data --> poisson regression

# linear regression was tested, suffered from heteroscedasticity
# results were similar with heteroskedasticity-robust standard errors 


# 7.3.1 Poisson Reg. with cluster-robust SEs - full sample ---------------------------

# Fit poisson regression to each imputed data set 
pois_multiimp <- with(KARE_MI, 
                      glm(numadmeas ~ R104 + R109 + R71 +        # generic capacity
                            R7 + R91 + R92 +
                            R90_1 + R90_5 +
                            R9 + R18 + R4a5 + R75 +    # flood-specific capacity
                            R63_6 + R80a_1 + R63_1 +
                            R63_3 + R63_2 + R63_4 +
                            R97 + R98 + R99 + R106 +   # CVs: gender, age, mig, hhsize
                            R69 + R70 +                # CVs: house
                            Modus,                     # Survey design
                          family = 'poisson'))


# calculate cluster-robust AME
marginal_pois <- avg_slopes(pois_multiimp, vcov = ~town)

# save results
marginal_pois$order <- c(33, 32, 1, 2, 28, 3, 12, 13, 14, 15, 16, 21, 
                         23, 22, 24, 19, 30, 29, 5, 31, 4, 17, 18, 20, 11,
                         9, 10, 6, 8, 7, 25, 26, 27)
marginal_pois <- marginal_pois[order(marginal_pois$order, decreasing = F),]
b_pois_AME <- c(NA, marginal_pois$estimate)
se_pois_AME <- c(NA, marginal_pois$std.error)
p_pois_AME <- c(NA, marginal_pois$p.value)


### Test Assumptions with residual diagnostics
sim_fmp <- simulateResiduals(pois_multiimp$analyses[[5]], nSim = 1000) 

# under- and overdispersion
testDispersion(sim_fmp)
# ratio close to 1 --> Poisson model fits well to the data, no over-/underdispersion
# not sign.

# Zero inflation
testZeroInflation(sim_fmp)
# not sign.

# Heteroscedasticity
testQuantiles(sim_fmp)
# not sign.

# KS test for correct distribution of residuals
testUniformity(sim_fmp)
# not sign.

# test if outliers are a concern, bootstrap for integer-valued distributions 
testOutliers(sim_fmp, type = "bootstrap")
# not sign.

### Goodness of fit measures
# R^2
Nagelkerke_pois_all <- numeric(30)
for (i in 1:30) {
  Nagelkerke_pois_all_res <- c(DescTools::PseudoR2(pois_multiimp$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_pois_all[i] <- Nagelkerke_pois_all_res
}
median(Nagelkerke_pois_all)

# BIC
BIC_pois_all <- numeric(30)
for (i in 1:30) {
  BIC_pois_all_res <- c(BIC(pois_multiimp$analyses[[i]]))
  BIC_pois_all[i] <- BIC_pois_all_res
}
median(BIC_pois_all)


# 7.3.2 Poisson Reg. with cluster-robust SEs - Owners ---------------------------

# Fit poisson regression to each imputed data set 
pois_multiimp_own <- with(KARE_MI, 
                          glm(numadmeas ~ R104 + R109 + R71 +        # generic capacity
                                R91 + R92 +
                                R90_1 + R90_5 +
                                R9 + R18 + R4a5 + R75 +    # flood-specific capacity
                                R63_6 + R80a_1 + R63_1 +
                                R63_3 + R63_2 + R63_4 +
                                R97 + R98 + R99 + R106 +   # CVs: gender, age, mig, hhsize
                                R69 + R70 +                # CVs: house
                                Modus,                     # Survey design
                              family = 'poisson',
                              subset = (R7 == "Eigentum")))

# calculate cluster-robust AME
marginal_pois_own <- avg_slopes(pois_multiimp_own, vcov = ~town)

# save results
marginal_pois_own$order <- c(33, 32, 1, 2, 28, 3, 12, 13, 14, 15, 16, 21, 
                             23, 22, 24, 19, 30, 29, 31, 4, 18, 20, 11,
                             9, 10, 6, 8, 7, 25, 26, 27)
marginal_pois_own <- marginal_pois_own[order(marginal_pois_own$order, decreasing = F),]
b_pois_AME_own <- c(NA, marginal_pois_own$estimate)
se_pois_AME_own <- c(NA, marginal_pois_own$std.error)
p_pois_AME_own <- c(NA, marginal_pois_own$p.value)


### Test Assumptions with residual diagnostics
sim_fmp_own <- simulateResiduals(pois_multiimp_own$analyses[[16]], nSim = 1000) 

# under- and overdispersion
testDispersion(sim_fmp_own)
# ratio close to 1 --> Poisson model fits well to the data, no over-/underdispersion
# not sign.

# Zero inflation
testZeroInflation(sim_fmp_own)
# not sign.

# Heteroscedasticity
testQuantiles(sim_fmp_own)
# not sign.

# KS test for correct distribution of residuals
testUniformity(sim_fmp_own)
# not sign.

# test if outliers are a concern, bootstrap for integer-valued distributions 
testOutliers(sim_fmp_own, type = "bootstrap")
# sign., but not a problem (no underdispersion, to few outliers)


### Goodness of fit measures
# R^2
Nagelkerke_pois_own <- numeric(30)
for (i in 1:30) {
  Nagelkerke_pois_own_res <- c(DescTools::PseudoR2(pois_multiimp_own$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_pois_own[i] <- Nagelkerke_pois_own_res
}
median(Nagelkerke_pois_own)

# BIC
BIC_pois_own <- numeric(30)
for (i in 1:30) {
  BIC_pois_own_res <- c(BIC(pois_multiimp_own$analyses[[i]]))
  BIC_pois_own[i] <- BIC_pois_own_res
}
median(BIC_pois_own)


# 7.3.3 Poisson Reg. with cluster-robust SEs - Tenants ---------------------------

# Fit poisson regression to each imputed data set 
pois_multiimp_ten <- with(KARE_MI, 
                          glm(numadmeas ~ R104 + R109 + R71 +        # generic capacity
                                R91 + R92 +
                                R90_1 + R90_5 +
                                R9 + R18 + R4a5 + R75 +    # flood-specific capacity
                                R63_6 + R80a_1 + R63_1 +
                                R63_3 + R63_2 + R63_4 +
                                R97 + R98 + R99 + R106 +   # CVs: gender, age, mig, hhsize
                                R69 + R70 +                # CVs: house
                                Modus,                     # Survey design
                              family = 'poisson',
                              subset = (R7 == "Miete")))

# calculate cluster-robust AME
marginal_pois_ten <- avg_slopes(pois_multiimp_ten, vcov = ~town)

# save results
marginal_pois_ten$order <- c(33, 32, 1, 2, 28, 3, 12, 13, 14, 15, 16, 21, 
                             23, 22, 24, 19, 30, 29, 31, 4, 17, 18, 20, 11,
                             9, 10, 6, 8, 7, 25, 26, 27)
marginal_pois_ten <- marginal_pois_ten[order(marginal_pois_ten$order, decreasing = F),]
b_pois_AME_ten <- c(NA, marginal_pois_ten$estimate)
se_pois_AME_ten <- c(NA, marginal_pois_ten$std.error)
p_pois_AME_ten <- c(NA, marginal_pois_ten$p.value)


### Test Assumptions with residual diagnostics
sim_fmp_ten <- simulateResiduals(pois_multiimp_ten$analyses[[16]], nSim = 1000) 

# under- and overdispersion
testDispersion(sim_fmp_ten)
# ratio close to 1 --> Poisson model fits well to the data, no over-/underdispersion
# not sign.

# Zero inflation
testZeroInflation(sim_fmp_ten)
# not sign.

# Heteroscedasticity
testQuantiles(sim_fmp_ten)
# not sign.

# KS test for correct distribution of residuals
testUniformity(sim_fmp_ten)
# not sign.

# test if outliers are a concern, bootstrap for integer-valued distributions 
testOutliers(sim_fmp_ten, type = "bootstrap")
# not sign.


### Goodness of fit measures
# R^2
Nagelkerke_pois_ten <- numeric(30)
for (i in 1:30) {
  Nagelkerke_pois_ten_res <- c(DescTools::PseudoR2(pois_multiimp_ten$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_pois_ten[i] <- Nagelkerke_pois_ten_res
}
median(Nagelkerke_pois_ten)

# BIC
BIC_pois_ten <- numeric(30)
for (i in 1:30) {
  BIC_pois_ten_res <- c(BIC(pois_multiimp_ten$analyses[[i]]))
  BIC_pois_ten[i] <- BIC_pois_ten_res
}
median(BIC_pois_ten)

# 7.3.4 Create table (App. A2) ---------------------------

stargazer(pois_model, pois_model_own, pois_model_ten,               # use complete case models as model objects as list objects are not supported
          coef = list(b_pois_AME, b_pois_AME_own, b_pois_AME_ten),  # but display coefficients etc. from imputed models 
          se = list(se_pois_AME, se_pois_AME_own, se_pois_AME_ten),
          p = list(p_pois_AME, p_pois_AME_own, p_pois_AME_ten),
          digits = 4,
          report=('vc*s'),
          title="Poisson Regressions",
          type = "html",
          out = here::here("./output/PoisReg_inck.doc"))


# 7.3.5 Create forest plot (Fig.7) ---------------------------

# label variables
marginal_pois_own <- marginal_pois_own %>%
  mutate(label = c("Education - intermediate secondary", "Education - upper secondary", "Income", 
                   "Living area", "Duration of residence","Place attachment - medium-term", "Place attachment - long-term",
                   "Social network", "Social cohesion", "Future risk perception",
                   "Risk perception - not likely", "Risk perception - likely", "Risk perception - very likely",
                   "Previous experience - damage", "Previous experience - experience",  
                   "Main responsibility state", "Expecation in authorities", "Trust in authorities - yes", 
                   "Public protection is sufficient", "Self-efficacy", "Motivation", 
                   "Competing concerns",
                   NA, NA, NA,  NA, NA, NA, NA, NA, NA ), 
         .after = term)

marginal_pois_ten <- marginal_pois_ten %>%
  mutate(label = c("Education - intermediate secondary", "Education - upper secondary", "Income", 
                   "Living area", "Duration of residence","Place attachment - medium-term", "Place attachment - long-term",
                   "Social network", "Social cohesion", "Future risk perception",
                   "Risk perception - not likely", "Risk perception - likely", "Risk perception - very likely",
                   "Previous experience - damage", "Previous experience - experience", "Main responsibility landlord",  
                   "Main responsibility state", "Expecation in authorities", "Trust in authorities - yes", 
                   "Public protection is sufficient", "Self-efficacy", "Motivation", 
                   "Competing concerns",
                   NA, NA, NA,  NA, NA, NA, NA, NA, NA ), 
         .after = term)

#drop control variables
marginal_pois_own2 <- marginal_pois_own %>% drop_na(label) 
marginal_pois_ten2 <- marginal_pois_ten %>% drop_na(label)

# merge results owners & tenants into one dataframe
marginal_pois_own_ten <- merge(marginal_pois_ten2, marginal_pois_own2, by = "label", all = T)
marginal_pois_own_ten <- marginal_pois_own_ten[order(marginal_pois_own_ten$estimate.y, decreasing = TRUE),]
# for plotting: add value for missing landlord which is not within plotted boundaries
marginal_pois_own_ten$estimate.y[is.na(marginal_pois_own_ten$estimate.y)] <- -2


# change colors 
colors_pois <- c("#588BAE", "#588BAE", "#588BAE", "black", "#588BAE",
                 "black", "#588BAE", "black", "#588BAE", "black", 
                 "#588BAE", "#588BAE", "black", "black", "#588BAE", 
                 "black", "black", "#588BAE", "black", "#588BAE", 
                 "#588BAE", "#588BAE", "#588BAE")
colors_pois_label <- c("#588BAE", "#588BAE", "#588BAE", "#588BAE", "black", 
                       "#588BAE", "black", "black", "#588BAE", "black", 
                       "black", "#588BAE", "#588BAE", "black", "#588BAE", 
                       "black", "#588BAE", "black", "#588BAE", "black",
                       "#588BAE", "#588BAE", "#588BAE")

lab_cap <- expression(paste("Source: own calculation, data from KARE household survey 2022 (",
                            n[Owners], " = 1,157; ", n[Tenants], " = 414)"))


pois_own <- ggplot(marginal_pois_own_ten, aes(y = reorder(label, estimate.y), x = estimate.y, xmin = conf.low.y, xmax = conf.high.y)) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "#758DA3", size=0.5) +
  geom_pointrange(color = colors_pois) + 
  labs(x ="AME", y = "",
       title = "Owners", caption = " ") + 
  xlim(-0.6, 1.3) + 
  theme_minimal() +
  theme(axis.text.y = element_text(colour = colors_pois_label, size = 14),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))


pois_ten <- ggplot(marginal_pois_own_ten, aes(y = reorder(label, estimate.y), x = estimate.x, xmin = conf.low.x, xmax = conf.high.x)) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "#758DA3", size=0.5) +
  geom_pointrange(color = colors_pois) + 
  labs(x ="AME", y = "",
       title = "Tenants", caption = lab_cap) + 
  xlim(-0.6, 1.3) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

title_pois <- ggdraw() + 
  draw_label("Effect sizes of adaptive capacity indicators explaing the number of implemented pluvial flood adaptation measures",
    fontface = 'bold', size = 16, x = 0, hjust = 0)  +
  theme( plot.margin = margin(0, 0, 0, 7))


plot_pois <- plot_grid(pois_own, pois_ten, rel_widths = c(1.85,1.12))
windows()
plot_grid(title_pois, plot_pois,  ncol = 1, rel_heights = c(0.1, 1))

ggsave(filename = here::here("./output/poisplot.pdf"), dpi = 300)
ggsave(filename = here::here("./output/poisplot.svg"))



# 7.4 Robustness check: Income / Income Groups ---------------------------

# load data
load(here::here("./data/tidy/KARE_imp_richlog.RData"))
KARE_MI <- KARE_imp

### MID (multiple imputation, then deletion): 
# use Y to impute X, but drop obs with missing Y from the analysis (von Hippel 2007)

# create a long format dataset  
KARE_MI_long <- complete(KARE_MI, action = "long", include = TRUE)

# drop obs with NA on Y
KARE_MI_long <- subset(KARE_MI_long, !(IDS %in% KARE_MI_long$IDS[is.na(KARE_MI_long$numadmeas)]))

# recode R75 from imputed R75a (tenants - polyreg) and R75b (owners - logreg)
table(KARE_MI_long$R7, KARE_MI_long$R75a, useNA = "always")
table(KARE_MI_long$R7, KARE_MI_long$R75b)
KARE_MI_long$R75 <- ifelse(KARE_MI_long$R7 == "Miete", KARE_MI_long$R75a, 
                           ifelse(KARE_MI_long$R75b == "Öffentliche Stellen (d.h. Bund, Länder oder Gemeinden)", 3, 1))
KARE_MI_long$R75 <- as_factor(KARE_MI_long$R75)
table(KARE_MI_long$R7, KARE_MI_long$R75, useNA = "always")

KARE_MI <- mice::as.mids(KARE_MI_long)
# n = 1,571 obs


### turn ordered factors into unordered to get pairwise comparisons 
# instead of polynomial contrasts
KARE_MI$data$R18 <- factor(KARE_MI$data$R18,
                           ordered = F)
KARE_MI$data$R80a_1 <- factor(KARE_MI$data$R80a_1,
                              ordered = F)
KARE_MI$data$R92 <- factor(KARE_MI$data$R92,
                           ordered = F)
KARE_MI$data$R104 <- factor(KARE_MI$data$R104,
                            ordered = F)
KARE_MI$data$Modus <- factor(KARE_MI$data$Modus,
                             ordered = F)

### change reference levels
KARE_MI$data$Modus <- relevel(KARE_MI$data$Modus, ref = "2")
KARE_MI$data$R92 <- relevel(KARE_MI$data$R92, ref = "unsure/short-term")
KARE_MI$data$R18 <- relevel(KARE_MI$data$R18, ref = "Überhaupt nicht wahrscheinlich")
KARE_MI$data$R4a5 <- relevel(KARE_MI$data$R4a5, ref = "none")
KARE_MI$data$R7 <- relevel(KARE_MI$data$R7, ref = "Miete")
KARE_MI$data$R80a_1 <- relevel(KARE_MI$data$R80a_1, ref = "Rather no/no")
KARE_MI$data$rich <- relevel(KARE_MI$data$rich, ref = "middle")


# 7.4.1 Log Reg: ACs --> Adaptation action yes/no ---------------------------

# binary outcome --> binary logistic regression

# 7.4.1.1 Logistic Reg with cluster-robust SEs - full sample ---------------------------

# fit logistic regression to each imputed data set 
glm_multiimp <- with(KARE_MI,
                     glm(admeas ~ R104 + logR109 + rich + R71 +        # generic capacity
                           R7 + R91 + R92 +
                           R90_1 + R90_5 +
                           R9 + R18 + R4a5 + R75 +                     # flood-specific capacity
                           R63_6 + R80a_1 + R63_1 +
                           R63_3 + R63_2 + R63_4 +
                           R97 + R98 + R99 + R106 +                    # CVs: gender, age, mig, hhsize
                           R69 + R70 +                                 # CVs: house
                           Modus,                                      # Survey design
                         family = binomial("logit")))

# calculate AME
marginal <- avg_slopes(glm_multiimp, vcov = ~town)

# save results
marginal$order <- c(34, 35, 1, 2, 30, 14, 15, 16, 17, 18, 23, 
                    25, 24, 26, 21, 32, 31, 7, 33, 6, 19, 20, 22, 13, 
                    11, 12, 8, 10, 9, 27, 28, 29, 3, 4, 5)  
marginal <- marginal[order(marginal$order, decreasing = F),]
b_log_rich_AME <- c(NA, marginal$estimate)
se_log_rich_AME <- c(NA, marginal$std.error)
p_log_rich_AME <- c(NA, marginal$p.value)


### Goodness of fit measures
#R^2
Nagelkerke_log_all <- numeric(30)
for (i in 1:30) {
  Nagelkerke_log_all_res <- c(DescTools::PseudoR2(glm_multiimp$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_log_all[i] <- Nagelkerke_log_all_res
}
median(Nagelkerke_log_all)

# BIC
BIC_log_all <- numeric(30)
for (i in 1:30) {
  BIC_log_all_res <- c(BIC(glm_multiimp$analyses[[i]]))
  BIC_log_all[i] <- BIC_log_all_res
}
median(BIC_log_all)


# check assumptions:
# 1) Linearity 
plot(glm_multiimp$analyses[[7]], 1)
# not useful to plot the raw residuals, examine binned residual plots
arm::binnedplot(fitted(glm_multiimp$analyses[[7]]), 
                residuals(glm_multiimp$analyses[[7]], type = "response"), 
                nclass = NULL, 
                xlab = "Expected Values", 
                ylab = "Average residual", 
                main = "Binned residual plot", 
                cex.pts = 0.8, 
                col.pts = 1, 
                col.int = "gray")
KARE_data6 <- complete(KARE_MI, action=6)
arm::binnedplot(KARE_data6$logR109, 
                residuals(glm_multiimp$analyses[[6]], type = "response"), 
                nclass = NULL, 
                xlab = "Expected Values", 
                ylab = "Average residual", 
                main = "Binned residual plot", 
                cex.pts = 0.8, 
                col.pts = 1, 
                col.int = "gray")


# 2) Independence 
# --> Respondents from the same town might be more similar to one another on the 
#     outcome measure, on average, than they are with respondents across towns

# plot residuals vs spatial variable
glm_data6 <- glm(data = KARE_data6,
                 admeas ~ R104 + logR109 + rich + R71 +        # generic capacity
                   R7 + R91 + R92 +
                   R90_1 + R90_5 +
                   R9 + R18 + R4a5 + R75 +                    # flood-specific capacity
                   R63_6 + R80a_1 + R63_1 +
                   R63_3 + R63_2 + R63_4 +
                   R97 + R98 + R99 + R106 +                   # CVs: gender, age, mig, hhsize
                   R69 + R70 +                                # CVs: house
                   Modus,                                     # Survey design
                 family = binomial("logit"))                                  

# by town
ggplot(data = data.frame(x = KARE_data6$town, y = glm_data6$residuals), aes(x = x, y = y)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(-20, 20))
# some patterns --> use cluster-robust SEs
# maybe also due to small cluster sizes (e.g. 2, 5, 10 respondents)


# 3) No multicollinearity 
car::vif(glm_multiimp$analyses[[10]]) 
# GVIF to the power of 1/2df makes the value of the GVIF comparable across 
# different number of parameters --> 
# GVIF (1/(2*Df) < 2 for all --> no multicollinearity


# 4) Exogeneity 
# --> all relevant CVs are included in the model


# 7.4.1.2 Logistic Reg with cluster-robust SEs - Owners ---------------------------

# fit logistic regression to each imputed data set 
glm_multiimp_own <- with(KARE_MI,
                         glm(admeas ~ R104 + logR109 + rich + R71 +   # generic capacity
                               R91 + R92 +
                               R90_1 + R90_5 +
                               R9 + R18 + R4a5 + R75 +      # flood-specific capacity
                               R63_6 + R80a_1 + R63_1 +
                               R63_3 + R63_2 + R63_4 +
                               R97 + R98 + R99 + R106 +     # CVs: gender, age, mig, hhsize
                               R69 + R70 +                  # CVs: house
                               Modus,                       # Survey design
                             family = binomial("logit"),
                             subset = (R7 == "Eigentum")))

# calculate AME
marginal_own <- avg_slopes(glm_multiimp_own, vcov = ~town)

# save results
marginal_own$order <- c(34, 35, 1, 2, 30, 14, 15, 16, 17, 18, 23, 
                        25, 24, 26, 21, 32, 31, 33, 6, 20, 22, 13, 
                        11, 12, 8, 10, 9, 27, 28, 29, 3, 4, 5) 
marginal_own <- marginal_own[order(marginal_own$order, decreasing = F),]
b_log_rich_AME_own <- c(NA, marginal_own$estimate)
se_log_rich_AME_own <- c(NA, marginal_own$std.error)
p_log_rich_AME_own <- c(NA, marginal_own$p.value)


### Goodness of fit measures
#R^2
Nagelkerke_log_own <- numeric(30)
for (i in 1:30) {
  Nagelkerke_log_own_res <- c(DescTools::PseudoR2(glm_multiimp_own$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_log_own[i] <- Nagelkerke_log_own_res
}
median(Nagelkerke_log_own)

# BIC
BIC_log_own <- numeric(30)
for (i in 1:30) {
  BIC_log_own_res <- c(BIC(glm_multiimp_own$analyses[[i]]))
  BIC_log_own[i] <- BIC_log_own_res
}
median(BIC_log_own)


# 7.4.1.3 Logistic Reg with cluster-robust SEs - Tenants ---------------------------

# fit logistic regression to each imputed data set 
glm_multiimp_ten <- with(KARE_MI,
                         glm(admeas ~ R104 + logR109 + rich + R71 +       # generic capacity
                               R91 + R92 +
                               R90_1 + R90_5 +
                               R9 + R18 + R4a5 + R75 +                    # flood-specific capacity
                               R63_6 + R80a_1 + R63_1 +
                               R63_3 + R63_2 + R63_4 +
                               R97 + R98 + R99 + R106 +                   # CVs: gender, age, mig, hhsize
                               R69 + R70 +                                # CVs: house
                               Modus,                                     # Survey design
                             family = binomial("logit"),
                             subset = (R7 == "Miete")))

# calculate AME
marginal_ten <- avg_slopes(glm_multiimp_ten, vcov = ~town)

# save results
marginal_ten$order <- c(34, 35, 1, 2, 30, 14, 15, 16, 17, 18, 23, 
                        25, 24, 26, 21, 32, 31, 33, 6, 19, 20, 22, 13, 
                        11, 12, 8, 10, 9, 27, 28, 29, 3, 4, 5) 
marginal_ten <- marginal_ten[order(marginal_ten$order, decreasing = F),]
b_log_rich_AME_ten <- c(NA, marginal_ten$estimate)
se_log_rich_AME_ten <- c(NA, marginal_ten$std.error)
p_log_rich_AME_ten <- c(NA, marginal_ten$p.value)


### Goodness of fit measures
#R^2
Nagelkerke_log_ten <- numeric(30)
for (i in 1:30) {
  Nagelkerke_log_ten_res <- c(DescTools::PseudoR2(glm_multiimp_ten$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_log_ten[i] <- Nagelkerke_log_ten_res
}
median(Nagelkerke_log_ten)

# BIC
BIC_log_ten <- numeric(30)
for (i in 1:30) {
  BIC_log_ten_res <- c(BIC(glm_multiimp_ten$analyses[[i]]))
  BIC_log_ten[i] <- BIC_log_ten_res
}
median(BIC_log_ten)


# 7.4.1.4 Create table (App. A3) ---------------------------

stargazer(log_rich_model, log_rich_model_own, log_rich_model_ten,               # use complete case models as model objects as list objects are not supported
          coef = list(b_log_rich_AME, b_log_rich_AME_own, b_log_rich_AME_ten),  # but display coefficients etc. from imputed models 
          se = list(se_log_rich_AME, se_log_rich_AME_own, se_log_rich_AME_ten),
          p = list(p_log_rich_AME, p_log_rich_AME_own, p_log_rich_AME_ten),
          digits = 4,
          report=('vc*s'),
          title="Logistic Regressions",
          type = "html",
          out = here::here("./output/LogReg_inclog_rich.doc"))


# 7.4.2 Poisson Reg: ACs --> Adaptation Action No. ---------------------------

# count data --> poisson regression

# 7.4.2.1 Poisson Reg. with cluster-robust SEs - full sample ---------------------------

# Fit poisson regression to each imputed data set 
pois_multiimp <- with(KARE_MI, 
                      glm(numadmeas ~ R104 + logR109 + rich + R71 +        # generic capacity
                            R7 + R91 + R92 +
                            R90_1 + R90_5 +
                            R9 + R18 + R4a5 + R75 +    # flood-specific capacity
                            R63_6 + R80a_1 + R63_1 +
                            R63_3 + R63_2 + R63_4 +
                            R97 + R98 + R99 + R106 +   # CVs: gender, age, mig, hhsize
                            R69 + R70 +                # CVs: house
                            Modus,                     # Survey design
                          family = 'poisson'))


# calculate cluster-robust AME
marginal_pois <- avg_slopes(pois_multiimp, vcov = ~town)

# save results
marginal_pois$order <- c(34, 35, 1, 2, 30,14, 15, 16, 17, 18, 23, 
                         25, 24, 26, 21, 32, 31, 7, 33, 6, 19, 20, 22, 13, 
                         11, 12, 8, 10, 9, 27, 28, 29, 3, 4, 5)
marginal_pois <- marginal_pois[order(marginal_pois$order, decreasing = F),]
b_pois_rich_AME <- c(NA, marginal_pois$estimate)
se_pois_rich_AME <- c(NA, marginal_pois$std.error)
p_pois_rich_AME <- c(NA, marginal_pois$p.value)


### Test Assumptions with residual diagnostics
sim_fmp <- simulateResiduals(pois_multiimp$analyses[[5]], nSim = 1000) 

# under- and overdispersion
testDispersion(sim_fmp)
# ratio close to 1 --> Poisson model fits well to the data, no over-/underdispersion
# not sign.

# Zero inflation
testZeroInflation(sim_fmp)
# not sign.

# Heteroscedasticity
testQuantiles(sim_fmp)
# not sign.

# KS test for correct distribution of residuals
testUniformity(sim_fmp)
# not sign.

# test if outliers are a concern, bootstrap for integer-valued distributions 
testOutliers(sim_fmp, type = "bootstrap")
# not sign.

### Goodness of fit measures
# R^2
Nagelkerke_pois_all <- numeric(30)
for (i in 1:30) {
  Nagelkerke_pois_all_res <- c(DescTools::PseudoR2(pois_multiimp$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_pois_all[i] <- Nagelkerke_pois_all_res
}
median(Nagelkerke_pois_all)

# BIC
BIC_pois_all <- numeric(30)
for (i in 1:30) {
  BIC_pois_all_res <- c(BIC(pois_multiimp$analyses[[i]]))
  BIC_pois_all[i] <- BIC_pois_all_res
}
median(BIC_pois_all)


# 7.4.2.2 Poisson Reg. with cluster-robust SEs - Owners ---------------------------

# Fit poisson regression to each imputed data set 
pois_multiimp_own <- with(KARE_MI, 
                          glm(numadmeas ~ R104 + logR109 + rich + R71 +        # generic capacity
                                R91 + R92 +
                                R90_1 + R90_5 +
                                R9 + R18 + R4a5 + R75 +    # flood-specific capacity
                                R63_6 + R80a_1 + R63_1 +
                                R63_3 + R63_2 + R63_4 +
                                R97 + R98 + R99 + R106 +   # CVs: gender, age, mig, hhsize
                                R69 + R70 +                # CVs: house
                                Modus,                     # Survey design
                              family = 'poisson',
                              subset = (R7 == "Eigentum")))

# calculate cluster-robust AME
marginal_pois_own <- avg_slopes(pois_multiimp_own, vcov = ~town)

# save results
marginal_pois_own$order <- c(34, 35, 1, 2, 30, 14, 15, 16, 17, 18, 23, 
                             25, 24, 26, 21, 32, 31, 33, 6, 20, 22, 13, 
                             11, 12, 8, 10, 9, 27, 28, 29, 3, 4, 5)   
marginal_pois_own <- marginal_pois_own[order(marginal_pois_own$order, decreasing = F),]
b_pois_rich_AME_own <- c(NA, marginal_pois_own$estimate)
se_pois_rich_AME_own <- c(NA, marginal_pois_own$std.error)
p_pois_rich_AME_own <- c(NA, marginal_pois_own$p.value)


### Test Assumptions with residual diagnostics
sim_fmp_own <- simulateResiduals(pois_multiimp_own$analyses[[16]], nSim = 1000) 

# under- and overdispersion
testDispersion(sim_fmp_own)
# ratio close to 1 --> Poisson model fits well to the data, no over-/underdispersion
# not sign.

# Zero inflation
testZeroInflation(sim_fmp_own)
# not sign.

# Heteroscedasticity
testQuantiles(sim_fmp_own)
# not sign.

# KS test for correct distribution of residuals
testUniformity(sim_fmp_own)
# not sign.

# test if outliers are a concern, bootstrap for integer-valued distributions 
testOutliers(sim_fmp_own, type = "bootstrap")
# sign., but not a problem (no underdispersion, to few outliers)


### Goodness of fit measures
# R^2
Nagelkerke_pois_own <- numeric(30)
for (i in 1:30) {
  Nagelkerke_pois_own_res <- c(DescTools::PseudoR2(pois_multiimp_own$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_pois_own[i] <- Nagelkerke_pois_own_res
}
median(Nagelkerke_pois_own)

# BIC
BIC_pois_own <- numeric(30)
for (i in 1:30) {
  BIC_pois_own_res <- c(BIC(pois_multiimp_own$analyses[[i]]))
  BIC_pois_own[i] <- BIC_pois_own_res
}
median(BIC_pois_own)


# 7.4.2.3 Poisson Reg. with cluster-robust SEs - Tenants ---------------------------

# Fit poisson regression to each imputed data set 
pois_multiimp_ten <- with(KARE_MI, 
                          glm(numadmeas ~ R104 + logR109 + rich + R71 +        # generic capacity
                                R91 + R92 +
                                R90_1 + R90_5 +
                                R9 + R18 + R4a5 + R75 +    # flood-specific capacity
                                R63_6 + R80a_1 + R63_1 +
                                R63_3 + R63_2 + R63_4 +
                                R97 + R98 + R99 + R106 +   # CVs: gender, age, mig, hhsize
                                R69 + R70 +                # CVs: house
                                Modus,                     # Survey design
                              family = 'poisson',
                              subset = (R7 == "Miete")))

# calculate cluster-robust AME
marginal_pois_ten <- avg_slopes(pois_multiimp_ten, vcov = ~town)

# save results
marginal_pois_ten$order <- c(34, 35, 1, 2, 30, 14, 15, 16, 17, 18, 23, 
                             25, 24, 26, 21, 32, 31, 33, 6, 19, 20, 22, 13, 
                             11, 12, 8, 10, 9, 27, 28, 29, 3, 4, 5)  
marginal_pois_ten <- marginal_pois_ten[order(marginal_pois_ten$order, decreasing = F),]
b_pois_rich_AME_ten <- c(NA, marginal_pois_ten$estimate)
se_pois_rich_AME_ten <- c(NA, marginal_pois_ten$std.error)
p_pois_rich_AME_ten <- c(NA, marginal_pois_ten$p.value)


### Test Assumptions with residual diagnostics
sim_fmp_ten <- simulateResiduals(pois_multiimp_ten$analyses[[16]], nSim = 1000) 

# under- and overdispersion
testDispersion(sim_fmp_ten)
# ratio close to 1 --> Poisson model fits well to the data, no over-/underdispersion
# not sign.

# Zero inflation
testZeroInflation(sim_fmp_ten)
# not sign.

# Heteroscedasticity
testQuantiles(sim_fmp_ten)
# not sign.

# KS test for correct distribution of residuals
testUniformity(sim_fmp_ten)
# not sign.

# test if outliers are a concern, bootstrap for integer-valued distributions 
testOutliers(sim_fmp_ten, type = "bootstrap")
# not sign.


### Goodness of fit measures
# R^2
Nagelkerke_pois_ten <- numeric(30)
for (i in 1:30) {
  Nagelkerke_pois_ten_res <- c(DescTools::PseudoR2(pois_multiimp_ten$analyses[[i]], which = "Nagelkerke"))
  Nagelkerke_pois_ten[i] <- Nagelkerke_pois_ten_res
}
median(Nagelkerke_pois_ten)

# BIC
BIC_pois_ten <- numeric(30)
for (i in 1:30) {
  BIC_pois_ten_res <- c(BIC(pois_multiimp_ten$analyses[[i]]))
  BIC_pois_ten[i] <- BIC_pois_ten_res
}
median(BIC_pois_ten)


# 7.4.2.4 Create table (App. A3) ---------------------------

stargazer(pois_rich_model, pois_rich_model_own, pois_rich_model_ten,              # use complete case models as model objects as list objects are not supported
          coef = list(b_pois_rich_AME, b_pois_rich_AME_own, b_pois_rich_AME_ten), # but display coefficients etc. from imputed models
          se = list(se_pois_rich_AME, se_pois_rich_AME_own, se_pois_rich_AME_ten),
          p = list(p_pois_rich_AME, p_pois_rich_AME_own, p_pois_rich_AME_ten),
          digits = 4,
          report=('vc*s'),
          title="Poisson Regressions",
          type = "html",
          out = here::here("./output/PoisReg_inclog_rich.doc"))


