######################################################
##                                                  ##
##        Unravelling the Capacity Action Gap       ##
##                                                  ##
##                Annika Schubert                   ##
##                                                  ##
##           04-imputation-missing-data             ##
##                                                  ##
######################################################

# @brief        Exploration of missing data patterns and imputation
# @license      GNU General Public License v3.0  
# @authors      Annika Schubert

# Code for:
# Schubert, A.; von Streit, A. & Garschagen, M. (2024):
# Unravelling the capacity-action gap in flood risk adaptation. 
# In: EGUsphere [preprint], https://doi.org/X. 


# ####################################################### #
#                                                         #
#### Part 4: Exploration & Imputation of Missing Data  ####
#                                                         #
# ####################################################### #

library(here)                # for file referencing
library(naniar)              # summaries for missing data 
library(VIM)                 # visualisation of missing values
library(ggplot2)             # plots
library(RColorBrewer)        # color palettes for plots
library(mice)                # imputation of missing values
library(remotes)             # to install packages from repositories (e.g. rNuggets)
library(rNuggets)            # calculation of the (extended) predmat (remotes::install_github("LukasWallrich/rNuggets"))
library(howManyImputations)  # to determine no. of multiple imputations (m)
library(reshape2)            # for reshaping data


# load dataset
load(here::here("./data/tidy/KARE_02_rec.RData"))


# 4.1 Explore Missing Data Patterns & Mechanisms ---------------------------

# 4.1.1 Overview Missings per Observation / Variable ---------------------------

# NAs per observation
obs_mis <- miss_case_summary(KARE)
plot(density(obs_mis$pct_miss, na.rm = T))

# % of NAs per var
var_mis <- miss_var_summary(KARE)
plot(density(var_mis$pct_miss, na.rm = T))
# high no. of missing data  
# - sensitive variables: political orientation (R96), income (R109)
# - filtered questions (R73, R31, R84, R62, etc.)
# - factual questions: natural hazard insurance coverage (R83, R83a), year of house construction (R70)


# 4.1.2 Missing Data Patterns & Mechanisms  ---------------------------

### 1) Frequency & occurring missingness patterns
windows()
aggr(KARE)

windows()
KARE %>% 
  select(R4, R5, R7, R9, R18, R61_1, R61_2, R61_3, R61_4, R61_5, R61_6, R61_7, 
         R61_8, R61_9, R61_10, R63_1, R63_3, R63_4, R63_6, R69, R70, R71, R75a,
         R80a_1, R90_1, R90_5, R91, R92, R97, R98, R99, R104, R106, R109, admeas,
         numadmeas) %>% 
  aggr()


### 2) Matrix plots: relationships between variable values and the propensity to be missing
# for relevant variables with high amount of missing data (> 5%) 

# income (R109)
KARE %>% 
  select(R97, R98, R99, R104, R105, R106, R109) %>% 
  matrixplot(sortby = "R109")
# no indication for relationships between variables and propensity to be missing

# competing concerns (R63_4)
KARE %>% 
  select(R63_4, R97, R98, R99, R101, R104, R105, R106, R109) %>% 
  matrixplot(sortby = "R101")
# no indication for relationships between variables and propensity to be missing

# expectation in authorities (R63_6)
KARE %>% 
  select(R63_6, R4, R5, R9, R12_1, R12_2, R65, R63_7, R97, R98, R99, R104, R105, R109) %>% 
  matrixplot(sortby = "R63_6")
# overall, no indication for relationships between variables and propensity to be missing
# however, monotone pattern for R63_6 & R63_7; high propensity to be missing when the other one is missing


### 3) Margin plot: scatterplot matrix with information about missing values in the plot margins 

windows()
KARE %>% 
  select(R90_3, R109) %>% 
  marginplot(delimiter = "imp")
KARE %>% 
  select(R98, R109) %>% 
  marginplot(delimiter = "imp")
# No indication for relationships between (numerical) variables and propensity for 
# income to be missing 


### 4) Mosaic plot

windows()
KARE %>% 
  select(R104, R97) %>% 
  mosaicMiss()
# No indication for relationships between variables and propensity to be missing 
# (checked for several variables)


# 4.2 Identify Variables for Imputation  ---------------------------

# based on
# Influx: connection of vars missingness indicator with observed data on other vars
# Outflux: connection of a vars observed data with missing data on other vars

fluxplot(KARE)
options(max.print=1000000)
flux(KARE)

# remove variables with low outflux / with high fractions of missing data.
# --> most variables from the two damage loops (R21 - R60), cellar usage (R72), 
#     and other filtered questions     

# Select Variables for Imputation
KARE <- KARE[names(KARE) %in% c("IDS", "Quelle", "Modus", "flag", "R1", "R2", "R4a5", "R6", "R7", 
                                "R8_1", "R8_2", "R8_3", "R8_4", "R8_5", "R8_6", "R9", 
                                "R10_1", "R10_2", "R10_3", "R10a", "R11", "R13_1",  "R13a", "R14",
                                "R15", "R16_1", "R16_2", "R16_3", "R17", "R18", "R19", "R20", 
                                "R26", "R27_2", "R27_3", "R30", "R36", "R37", "R39", 
                                "R40_1", "R40_2", "R40_3", "R40_4", "R40_5", "R40_6", "R40_7", "R40_8", "R40_9", "R40_10",
                                "R61_1", "R61_2", "R61_3", "R61_4", "R61_5", "R61_6", "R61_7", "R61_8", 
                                "R61_9", "R61_10", "admeas", "numadmeas",
                                "R63_1", "R63_2", "R63_3", "R63_4", "R63_5", "R63_6", "R63_7", 
                                "R69", "R70", "R71", "R75a", "R75b", "R76_3", "R77_1",  
                                "R77_2", "R77_3", "R77_4", "R78", "R80a_1", "R80a_2", "R80a_3", "R80b_1", "R80b_2", 
                                "R81_1", "R81_2", "R83", "R83a", "R85_1", "R85_2", "R85_3", 
                                "R85_4", "R85_5", "R86", "R87", "R88_1", "R88_2", "R88_3", "R88_4", 
                                "R90_1", "R90_3", "R90_5", "R91", "R92", "R95", "R97", "R98",  
                                "R99", "R100", "R104", "R105", "R106", "R109", "R96", "town")]
#                               "hhsize_eq", "R109_eq", "rich")]  # for robustness check

# 4.3 Predictor Matrix  ---------------------------

# 4.3.1 Calculate Predictor Matrix ---------------------------

# for full sample
pred_mat <- rNuggets::quickpred_ext(    # ext: tests each category of unordered factors
  KARE,                                             
  mincor = 0.1,  # use 0.8 to check for potential collinearity
  method = "spearman")

# for owner subset (needed for owner-only measures)
KARE_owner <- subset(KARE, KARE$R7 == "Eigentum")
KARE_owner <- droplevels(KARE_owner)
pred_mat_owner <- quickpred_ext(
   KARE_owner,
   mincor = 0.1,
   method = "spearman")
# save pred_mat as excel


# 4.3.2 Refine Predictor Matrix ---------------------------

# R61: Variables predict R61, but R61 not other variables
pred_mat[, c("R61_1", "R61_2", "R61_3", "R61_4", "R61_5", "R61_6", "R61_7", "R61_8", "R61_9", "R61_10")] <- 0

# use owner cors for prediction of owner-only measures (to reduce multicollinarity & no. of vars)
pred_mat[c("R61_1", "R61_3", "R61_5", "R61_6", "R61_8", "R61_9" ), ] <- 0            # owners 
pred_mat["R61_1", c("R2", "R4a5", "R17", "R30", "R40_6", "R40_8", "R40_9", "R63_2", "R83a")] <- 1           # owners
pred_mat["R61_3", c("Modus", "R2", "R4a5", "R18", "R26", "R27_2", "R27_3", "R30", "R36", "R37", "R40_4", "R40_5", 
                    "R40_7", "R40_10", "R63_2", "R70", "R83", "R85_3", "R85_5", "R90_1")] <- 1 # owners
pred_mat["R61_5", c("R27_3", "R37", "R40_1", "R40_3", "R40_9", "R40_10", "R63_2", "R85_5")] <- 1          # owners
pred_mat["R61_6", c("R30", "R63_2", "R63_3")] <- 1            # owners
pred_mat["R61_8", c("Quelle", "Modus", "R17", "R26", "R27_2", "R30", "R36", "R37", "R40_4", "R40_5", "R40_9", 
                    "R63_2", "R63_3", "R85_1", "R85_5", "R99")] <- 1  # owners
pred_mat["R61_9", c("town", "R4a5", "R8_4", "R9", "R10_1", "R10_2", "R17", "R18", "R20", "R26", "R27_2", 
                    "R27_3", "R30", "R36", "R37", "R40_3", "R40_4", "R40_5", "R40_6", "R40_7", "R40_8", 
                    "R40_9", "R40_10", "R63_1", "R63_2", "R63_3", "R63_4", "R63_6", "R69", "R71", 
                    "R83a", "R85_1", "R85_3", "R85_5", "R98", "R99", "R105" )] <- 1  # owners

# R61 predict each other (cor >= 0.1,identify correlation for owner-only
# measures based on cors from owner dataset)
pred_mat["R61_1", c("R61_2", "R61_3", "R61_4", "R61_5", "R61_8", "R61_10")] <- 1           # owners
pred_mat["R61_2", c("R61_1", "R61_3", "R61_4", "R61_5", "R61_6", "R61_7", "R61_8", "R61_9", "R61_10")] <- 1     
pred_mat["R61_3", c("R61_1", "R61_2", "R61_4", "R61_5", "R61_6", "R61_7", "R61_8", "R61_9", "R61_10")] <- 1 # owners 
pred_mat["R61_4", c("R61_1", "R61_2", "R61_3", "R61_5", "R61_6", "R61_7", "R61_8", "R61_9", "R61_10")] <- 1
pred_mat["R61_5", c("R61_1", "R61_3", "R61_4", "R61_6",  "R61_8", "R61_10")] <- 1          # owners 
pred_mat["R61_6", c("R61_3", "R61_4", "R61_5", "R61_7", "R61_8", "R61_9")] <- 1            # owners
pred_mat["R61_7", c("R61_1", "R61_2", "R61_3", "R61_4", "R61_5", "R61_6", "R61_8", "R61_9")] <- 1         
pred_mat["R61_8", c("R61_3", "R61_4", "R61_5", "R61_6", "R61_7", "R61_9", "R61_10")] <- 1  # owners
pred_mat["R61_9", c("R61_2", "R61_3", "R61_4", "R61_6", "R61_8", "R61_10")] <- 1  # owners 
pred_mat["R61_10", c("R61_1", "R61_2", "R61_3", "R61_4", "R61_5",  "R61_8", "R61_9")] <- 1  

# R61 predicts numadmeas, numadmeas predicts admeas
pred_mat["numadmeas", ] <- 0 
pred_mat["admeas", ] <- 0 
pred_mat["numadmeas", c("R61_1", "R61_2", "R61_3", "R61_4", "R61_5", "R61_6", "R61_7", "R61_8", "R61_9", "R61_10")] <- 1
pred_mat["admeas", c("numadmeas")] <- 1

# R106: refine prediction to avoid multicollinearity (dropped: R8, R10, R75, R77, R81, R83)
pred_mat["R106", ] <- 0 
pred_mat["R106", c("R1", "R7", "R63_3", "R69", "R71", "R85_3", "R90_1", "R90_3", "R98", "R100", "R105", "R109", "numadmeas", "admeas")] <- 1

# R75: only R75a is used for predicting other vars, not R75b (to avoid multicollinearity)
pred_mat[, c("R75b")] <- 0
pred_mat["R75b", c("R75a")] <- 0   # R75a does not predict R75b
pred_mat["R75b", c("R7","R76_3")] <- 0   # no ownership-related variables for prediction

# numadmeas & admeas predict all adaptive capacity indicators (= regression model)
pred_mat[c("R4a5", "R7", "R9", "R18", "R63_1", "R63_2", "R63_3", "R63_4", "R63_6", 
           "R69", "R70", "R71", "R75a", "R75b", "R80a_1", "R90_1", "R90_5", "R91", "R92", 
           "R97", "R98", "R99", "R104", "R106", "R109", "Modus"), c("numadmeas", "admeas") ] <- 1
           # "R109log", "rich"  for robustness check

### robustness check: income & rich 
#pred_mat["hhsize_eq", ] <- 0 
#pred_mat[, c("hhsize_eq")] <- 0
#pred_mat["R109_eq", ] <- 0 
#pred_mat[, c("R109_eq")] <- 0
#pred_mat["rich", ] <- 0 
#pred_mat["hhsize_eq", c("R106")] <- 1
#pred_mat["R109_eq", c("R106", "R109")] <- 1
#pred_mat["rich", c("R109_eq")] <- 1

# check no of predictors per variable
rowSums(pred_mat)
mean(rowSums(pred_mat))  # 22 is okay

# times var is used as predictor
colSums(pred_mat)


# 4.4 MICE Imputation  ---------------------------

# 4.4.1 Imputation Set up ---------------------------

KARE_imp <- mice(KARE,                         
                 m = 10,                       
                 predictorMatrix = pred_mat,
                 seed = 248,
                 maxit = 0,
                 defaultMethod = c("rf", "logreg", "polyreg", "polr"),
                 nnet.MaxNWts = 2000)


### Imputation method

meth <- KARE_imp$method
# pmm for count data
meth["R106"] <- "pmm"
# count no. of implemented measures (calculate numadmeas based on imputed R61_X)
# for tenants, consider only the four implementable measures
meth["numadmeas"] <- "~I(ifelse(R7 == 'Eigentum',
                        (ifelse(R61_1 == 'Implementiert', 1, 0)+ifelse(R61_2 == 'Implementiert', 1, 0)+ifelse(R61_3 == 'Implementiert', 1, 0)+ifelse(R61_4 == 'Implementiert', 1, 0)+ifelse(R61_5 == 'Implementiert', 1, 0)+ifelse(R61_6 == 'Implementiert', 1, 0)+ifelse(R61_7 == 'Implementiert', 1, 0)+ifelse(R61_8 == 'Implementiert', 1, 0)+ifelse(R61_9 == 'Implementiert', 1, 0)+ifelse(R61_10 == 'Implementiert', 1, 0)),
                        (ifelse(R61_2 == 'Implementiert', 1, 0)+ifelse(R61_4 == 'Implementiert', 1, 0)+ifelse(R61_7 == 'Implementiert', 1, 0)+ifelse(R61_10 == 'Implementiert', 1, 0))))" 
# Adaption yes if at least one implemented measure
meth["admeas"] <- "~ifelse(numadmeas > 0, 1, 0)" 
# for robustness check:
# household size equivalent
#meth["hhsize_eq"] <- "~ifelse(R106 == 1, 1, 1 + (R106 - 1)*0.5)" 
# income equivalent
#meth["R109_eq"] <- "~I(R109/hhsize_eq)" 
# rich indicator
#meth["rich"] <- "~ifelse(R109_eq < 1300, 'poor', 
#                  ifelse(R109_eq >= 4000, 'rich',
#                  'middle'))"


### Post processing ###

post <- KARE_imp$post
# duration of housing can not be greater than age (R91 =< R98)
post["R91"] <- "imp[[j]][, i] <- pmin(imp[[j]][, i], data[(!r[, j]) & where[, j], 'R98'])"
# age has to be greater than duration of residence (R98 >= R91)
post["R98"] <- "imp[[j]][, i] <- pmax(imp[[j]][, i], data[(!r[, j]) & where[, j], 'R91'])"


### Visiting Sequence 

# adjust visiting sequence due to passive imputation
visSeq <- KARE_imp$vis 
# change visiting sequence so that numadmeas & admeas are calculated after R61
# and R98 is updated after imputation of R91
vis_update <- miceadds::visitSequence.determine(impMethod=meth, vis = visSeq, data=KARE, maxit = 15)


### Determine No. of Imputations
# based on von Hippel (2020)

# Log reg (run regressions beforehand)
howManyImputations::how_many_imputations(glm_rob_own_multiimp, cv = 0.05, alpha = 0.05) # 20
howManyImputations::how_many_imputations(glm_rob_ten_multiimp, cv = 0.05, alpha = 0.05) # 25

# Poisson reg (run regressions beforehand)
howManyImputations::how_many_imputations(pois_rob_own_multiimp, cv = 0.05, alpha = 0.05) # 9
howManyImputations::how_many_imputations(pois_rob_ten_multiimp, cv = 0.05, alpha = 0.05) # 30


# 4.4.2 Imputation with MICE  ---------------------------

KARE_imp <- mice(KARE,                         
                 m = 30,                   
                 predictorMatrix = pred_mat,
                 seed = 248,
                 maxit = 15,   
                 method = meth,
                 post = post,
                 vis = visSeq,
                 nnet.MaxNWts = 2000)    

# save imputation
save(KARE_imp, file = here::here("./data/tidy/KARE_imp_inck.RData"))  
# save(KARE_imp, file = here::here("./data/tidy/KARE_imp_richlog.RData"))  # for robustness check

load(here::here("./data/tidy/KARE_imp_inck.RData"))
# load(here::here("./data/tidy/KARE_imp_richlog.RData")) # for robustness check


# 4.5 Evaluation of Imputations  ---------------------------

# 1) Logged Events
KARE_imp$loggedEvents
head(KARE_imp$loggedEvents, 80)


# 2) Convergence: Traceplots
windows()
plot(KARE_imp, layout = c(15,12))
# no trends detectable


# 3) Plausibility checks (values)
# display summary statistics and NAs
my_summary <- function(v){
  if(!any(is.na(v))){
    res <- c(summary(v),"NA's"=0)
  } else{
    res <- summary(v)
  }
  return(res)
}

longDF <- complete(KARE_imp, "long")
summary(longDF)
lapply(KARE_imp$imp, function(x)   # get the summary of the "raw" imputed values
  my_summary(unlist(x))
)

# check: age not larger than duration of residence
longDF$ageresi_controll <- longDF$R98 - longDF$R91
table(longDF$ageresi_controll)

# check: tenants not more than 4 implemented measures
table(longDF$R7, longDF$numadmeas)



# 4) Univariate distribution checks

## density plot (only metric vars)
windows()
densityplot(KARE_imp)


## bar charts (all vars)
# use propplot() function from Nicole Erler (GitHub: NErler)
propplot <- function(x, formula, facet = "wrap", ...) {
  library(ggplot2)
  
  cd <- data.frame(mice::complete(x, "long", include = TRUE))
  cd$.imp <- factor(cd$.imp)
  
  r <- as.data.frame(is.na(x$data))
  
  impcat <- x$meth != "" & sapply(x$data, is.factor)
  vnames <- names(impcat)[impcat]
  
  if (missing(formula)) {
    formula <- as.formula(paste(paste(vnames, collapse = "+",
                                      sep = ""), "~1", sep = ""))
  }
  
  tmsx <- terms(formula[-3], data = x$data)
  xnames <- attr(tmsx, "term.labels")
  xnames <- xnames[xnames %in% vnames]
  
  if (paste(formula[3]) != "1") {
    wvars <- gsub("[[:space:]]*\\|[[:print:]]*", "", paste(formula)[3])
    # wvars <- all.vars(as.formula(paste("~", wvars)))
    wvars <- attr(terms(as.formula(paste("~", wvars))), "term.labels")
    if (grepl("\\|", formula[3])) {
      svars <- gsub("[[:print:]]*\\|[[:space:]]*", "", paste(formula)[3])
      svars <- all.vars(as.formula(paste("~", svars)))
    } else {
      svars <- ".imp"
    }
  } else {
    wvars <- NULL
    svars <- ".imp"
  }
  
  for (i in seq_along(xnames)) {
    xvar <- xnames[i]
    select <- cd$.imp != 0 & !r[, xvar]
    cd[select, xvar] <- NA
  }
  
  
  for (i in which(!wvars %in% names(cd))) {
    cd[, wvars[i]] <- with(cd, eval(parse(text = wvars[i])))
  }
  
  meltDF <- reshape2::melt(cd[, c(wvars, svars, xnames)], id.vars = c(wvars, svars))
  meltDF <- meltDF[!is.na(meltDF$value), ]
  
  
  wvars <- if (!is.null(wvars)) paste0("`", wvars, "`")
  
  a <- plyr::ddply(meltDF, c(wvars, svars, "variable", "value"), plyr::summarize,
                   count = length(value))
  b <- plyr::ddply(meltDF, c(wvars, svars, "variable"), plyr::summarize,
                   tot = length(value))
  mdf <- merge(a,b)
  mdf$prop <- mdf$count / mdf$tot
  
  plotDF <- merge(unique(meltDF), mdf)
  plotDF$value <- factor(plotDF$value,
                         levels = unique(unlist(lapply(x$data[, xnames], levels))),
                         ordered = T)
  
  p <- ggplot(plotDF, aes(x = value, fill = get(svars), y = prop)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme(legend.position = "bottom", ...) +
    ylab("proportion") +
    scale_fill_manual(name = "",
                      values = c("black",
                                 colorRampPalette(
                                   RColorBrewer::brewer.pal(9, "Blues"))(x$m + 3)[1:x$m + 3])) +
    guides(fill = guide_legend(nrow = 1))
  
  if (facet == "wrap")
    if (length(xnames) > 1) {
      print(p + facet_wrap(c("variable", wvars), scales = "free"))
    } else {
      if (is.null(wvars)) {
        print(p)
      } else {
        print(p + facet_wrap(wvars, scales = "free"))
      }
    }
  
  if (facet == "grid")
    if (!is.null(wvars)) {
      print(p + facet_grid(paste(paste(wvars, collapse = "+"), "~ variable"),
                           scales = "free"))
    }
}

windows()
par(mfrow = c(6, 4))
propplot(KARE_imp)    

