# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                  LIBRAIRIES                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(tidyverse)
library(here)
library(readr)
library(forcats)
theme_set(theme_classic2())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                  CONSTANTS                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cols = c("#262D37",
         "#3D5A5D",
         "#64A395",
         "#E1B782",
         "#E67C37",
         "#8F2F0E")


folderFig = here('OutputFigures')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                 LOADING DATA                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Temporal dynamics
df = read.csv(here('Data/temporalDynamicsData.csv'))

# Adding slenderness
df = df %>%
  mutate(
    Date = factor(Date),
    Pop = factor(Pop),
    Slenderness = Total_height/Plant_Biomass)


# Within population analysis
df.pop5 = read.csv(here('Data/withinPopData.csv'))

# Adding standardized male and female fecundities
df.pop5 = df.pop5 %>% 
  mutate(nbFruit_sd = nbFruit/subsample_dry_mass,
         nbMfl_sd = nbMfl/subsample_dry_mass)





