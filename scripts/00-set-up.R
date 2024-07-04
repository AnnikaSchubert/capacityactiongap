######################################################
##                                                  ##
##        Unravelling the Capacity Action Gap       ##
##                                                  ##
##                Annika Schubert                   ##
##                                                  ##
##                 00-set-up                        ##
##                                                  ##
######################################################

# @brief        Cleans the data, imputes missing data and analyses the data 
# @license      GNU General Public License v3.0  
# @authors      Annika Schubert

# Code for:
# Schubert, A.; von Streit, A. & Garschagen, M. (2024):
# Unravelling the capacity-action gap in flood risk adaptation. 
# In: EGUsphere [preprint], https://doi.org/X. 

# Scripts:
# 00-set-up
# 01-IER-detection
# 02-sample-and-recoding
# 03-descriptive-statistics
# 04-imputation-missing-data
# 05-regression-models



# ############################################## #
#                                                #
#### Part 0: Setup                            ####
#                                                #
# ############################################## #

# 1.1 UTF-8 Encoding   ---------------------------

# reload script for a correct display of German umlauts (ä, ö, ü):
# reopen with UTF-8 encoding


# 1.2 Libraries  ---------------------------

# 00-set-up
library(here)                # for file referencing
library(renv)                # create reproducible environment

# 01-IER-detection
library(here)                # for file referencing
library(tidyverse)           # for data manipulation, exploration and visualization
library(car)                 # for recoding
library(haven)               # for working with (SPSS) labbelled dataset
library(remotes)             # to install packages from repositories (e.g. careless)
library(careless)            # to detect such careless/insufﬁcient effort responses (remotes::install_github('ryentes/careless'))

# 02-sample-and-recoding
library(here)                # for file referencing
library(tidyverse)           # for data manipulation, exploration and visualization
library(ggplot2)             # plots


# 03-descriptive-statistics
library(here)                # for file referencing
library(ggplot2)             # plots
library(dplyr)               # data transformation, filtering and summarising
library(missMethods)         # compute median of ordered factors
library(ComplexHeatmap)      # create upset plot
library(readxl)              # load excel files

# 04-imputation-missing-data
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

# 05-regression-models
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


# 1.2 Set Path ---------------------------

# set path
here::set_here()

# check path
here::here()


# 1.3 Create Reproducible Environment ---------------------------

# initialize renv
renv::init()

# to create snapshot of used packages
renv::snapshot()

