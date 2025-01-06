
#  This script is directly derived from the template for linear modelling 
#  developped by Matteo Santon, Fraenzi Korner-Nievergelt, Nico Michiels, Nils Anthes
#  Version date: 27 March 2023

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                         TEMPORAL DYNAMICS ANALYSIS                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Setup -----------------------------------------------------------------------

rm(list = ls()) 

# Activate the support 'functions file' for this routine:
source(here("LinearModellingWorkflow_supportFunctions.R")) 
# This script also contains librairies that must be available

# Data importation and first data treatment:
source(here("DataImportation.R")) 

## Definition of variables -----------------------------------------------------

# Response variable
var_resp <- "Weight_mfl_net"

# Fixed predictors
# Factors (none here)
var_fac <- NA                           
# Continuous
var_num <- c("Gen", "Weight_seeds_Top_net", "Slenderness", "Total_height")

# Random predictors
var_rand <- c("Pop", "Gen_Pop", "Date")

## Temporal and spatial data structure -----------------------------------------

# Time variable (generations)
var_time <- "Gen"

# Group structure
var_time_groups <- "Pop"

# Space variable (none here)
var_space <- NA


## Missing values --------------------------------------------------------------

df.pr <- remove_NAs(data = df, variables = c(var_num, var_fac, var_rand, var_resp))
df.NAs <- anti_join(df,  df.pr)

# Keep only complete cases in the dataset for further analysis: 
df <- df.pr


## Extreme values --------------------------------------------------------------
# Graphical inspection of extreme values.

# Dotplots for all numerical predictors plus the response variable.
dotplot_num(data = df, variables = c(var_num, var_resp))

# 1 slenderness value is extreme:
df %>% filter(Slenderness > 60)
# Removing it out of prudence
df = df %>% filter(Slenderness<60)

# Barplots for all factor predictors and random terms.
barplot_fac(data = df, variables = c(var_fac, var_rand))

## Predictor collinearity ------------------------------------------------------
# Graphical and numerical inspection variables for predictor collinearity.

# Pairwise scatterplots for all numeric predictors.
coll_num(data = df, predictors = var_num) 


## Variance Inflation Factors --------------------------------------------------

# VIF are calculated for all fixed predictors potentially included in the model.
corvif(data = df, variables = c(var_num, var_fac))
# Derived from Zuur, Ieno & Elphick (2010).


## Predictor-response relationships --------------------------------------------

# Plots of the response variable against each of the predictor variables.
relat_single(data = df,
             response = var_resp,
             predictors = c(var_num, var_fac))


## Response distribution -------------------------------------------------------

distr_response(data = df, response = var_resp) 
# A Gamma model seems to be a good fit


## Model formulation -----------------------------------------------------------

# Standardisation of numeric predictors:
df <- cbind(df, z_transform(data = df,
                            predictors = var_num))


# Model: 
mod <- glmmTMB(Weight_mfl_net ~ 
                 Weight_seeds_Top_net_z + Total_height_z + Gen + Slenderness_z + 
                 (1|Pop) + (1|Date) + (1|Gen_Pop),
               data = df,
               ziformula=~1,
               family = ziGamma(link = "log"))


## Model assessment ------------------------------------------------------------

# Distribution of residuals
residual_plots(data = df,
               modelTMB = mod,
               response = var_resp)

# Residuals against possible two-way interactions among predictors
residual_plots_interactions(data = df, 
                            modelTMB = mod, 
                            predictors = c(var_num, var_fac))


# Compare the variance of the observed data with the variance distribution in model-simulated dataframes:
dispersion_simulation(data = df, 
                      modelTMB = mod,
                      response = var_resp,
                      n.sim = 500)

# The function compares the observed number of zeros to the zero-distribution in model-simulated datasets. 
zero_simulation(data = df, 
                modelTMB = mod,
                response = var_resp,
                n.sim = 500)

# The function compares the observed raw data distribution with that observed in model-simulated dataframes.
ppcheck_fun(data = df,
            modelTMB = mod,
            response = var_resp,
            n.sim = 500)

# We check for temporal correlation patterns in model residuals (= autocorrelation) 
# using standardised semivariograms.
autocor_check(data = df, 
              modelTMB = mod,
              variable = var_time,
              n.sim = 500)


## Model results ---------------------------------------------------------------

# Extract coefficient estimates and 95% compatibility intervals:
comp_int(modelTMB = mod,
         ci_range = 0.95,
         effects = "all",
         component = "all")

# Extract marginal and conditional R-squared value.
r2(mod, tolerance = 1e-10)



## Final plotting --------------------------------------------------------------


library(smatr)

modSma = sma(Weight_mfl_net ~ Weight_seeds_Top_net, data = df %>% mutate(Weight_mfl_net = Weight_mfl_net * 1000,
                                                                         Weight_seeds_Top_net = Weight_seeds_Top_net * 1000))

df.plot = df %>% 
  mutate(
    mflMass = mean(Weight_mfl_net, na.rm = T) * 1000,
    seedMass = mean(Weight_seeds_Top_net, na.rm = T) * 1000
  )


df.byPop = df %>% 
  group_by(Gen, Pop) %>% 
  summarise(
    mflMass = mean(Weight_mfl_net, na.rm = T) * 1000,
    seedMass = mean(Weight_seeds_Top_net, na.rm = T) * 1000,
    mflMass_se = sd(Weight_mfl_net, na.rm = T)/sqrt(n()) * 1000,
    seedMass_se = sd(Weight_seeds_Top_net, na.rm = T)/sqrt(n()) * 1000,
    sexAlloc = mean(relMRE, na.rm = T),
    nbInds = n()
  ) %>% 
  mutate(
    Pop = plyr::revalue(Pop, c("5" = "Population 5", 
                               "3" = "Population 3", 
                               "1" = "Population 1"))
  )


df.plot = 
  left_join(df.plot, df.byPop, by = c("Gen", "Pop"), suffix = c(".ind", ".pop"))


fg1 = function(df, noYlab = F, noXtitle = F) {
  g = df %>% ggplot() +
    geom_abline(slope = coef(modSma)[2], intercept = coef(modSma)[1], 
                show.legend = T, linetype = "dashed", linewidth = 1, color = "grey50") + 
    geom_path(aes(x = seedMass, y = mflMass), linewidth = 1.1, color = "grey70", lineend = "round") +
    geom_point(aes(x = seedMass, y = mflMass, col = sexAlloc), size = 6) +
    geom_text(aes(x = seedMass, y = mflMass, label = paste0("G", Gen)), col = "white", fontface = "bold", size = 3) + 
    scale_color_gradientn(colours = colorRampPalette(cols[c(3,3,4,5,6,6)])(6), limits = c(0, 1)) +
    scale_x_continuous(expand = c(0.01, 0), limits = c(0, 150), breaks = c(0, 0.05, 0.1, 0.15)*1000) +
    scale_y_continuous(expand = c(0.01, 0), limits = c(0, 200), breaks = c(0, 0.05, 0.1, 0.15, 0.20)*1000) + 
    theme(legend.position="none",
          plot.title = element_text(size=14, face="bold", hjust = 0.5)) +
    labs(title = df$Pop[1]) +
    ylab("Biomass of male flowers (mg)") +
    xlab("Biomass of seeds (mg)")
  if (noYlab){
    g = g + theme(axis.title.y = element_blank(),
                  axis.text.y = element_blank())
  }
  if (noXtitle){
    g = g + theme(axis.title.x = element_blank())
  }
  return(g)
}

g1.1 = fg1(df.byPop %>% filter(Pop == "Population 1"), noXtitle = T)
g1.2 = fg1(df.byPop %>% filter(Pop == "Population 3"), noYlab = T)
g1.3 = fg1(df.byPop %>% filter(Pop == "Population 5"), noYlab = T, noXtitle = T)

g1.4 = df %>% 
  group_by(Gen) %>% 
  mutate(Gen = as.integer(as.character(Gen))) %>% 
  summarise(
    mflMass = mean(Weight_mfl_net, na.rm = T) * 1000,
    seedMass = mean(Weight_seeds_Top_net, na.rm = T) * 1000,
    mflMass_se = sd(Weight_mfl_net, na.rm = T)/sqrt(n()) * 1000,
    seedMass_se = sd(Weight_seeds_Top_net, na.rm = T)/sqrt(n()) * 1000,
    sexAlloc = mean(relMRE, na.rm = T)
  ) %>% 
  mutate(Pop = "All populations") %>% 
  ggplot() +
  geom_abline(slope = coef(modSma)[2], intercept = coef(modSma)[1], 
              show.legend = T, linetype = "dashed", lineend = "round", linewidth = 1, color = "grey50") + 
  geom_path(aes(x = seedMass, y = mflMass), linewidth = 1.1, color = "grey70") +
  geom_point(aes(x = seedMass, y = mflMass, col = sexAlloc), size = 6) +
  geom_text(aes(x = seedMass, y = mflMass, label = paste0("G", Gen)), col = "white", fontface = "bold", size = 3) + 
  scale_color_gradientn(colours = colorRampPalette(cols[c(3,3,4,5,6,6)])(6), limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(40, 105), breaks = c(50, 75, 100)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(15, 115), breaks = c(25, 50, 75, 100)) + 
  theme(plot.title = element_text(size=14, face="bold", hjust = 0.5), 
        legend.title = element_text(margin = margin(b = 15))) +
  labs(title = "All populations", color = "Relative male\n allocation") +
  ylab("Biomass of male flowers (mg)") +
  xlab("Biomass of seeds (mg)") +
  guides(color = guide_colourbar(theme = theme(
    legend.ticks = element_blank()
  )))

g1 = cowplot::plot_grid(
  cowplot::plot_grid(
    g1.1,
    g1.2,
    g1.3,
    align = "h",
    nrow = 1, 
    rel_widths = c(1.17, 1, 1),
    labels = LETTERS[1:3],
    label_size = 15),
  g1.4 + theme(plot.margin = margin(1, 1, 1, 3.5, "cm")),
  ncol = 1, 
  rel_heights = c(1, 1.3), 
  labels = c("", LETTERS[4]),
  label_size = 15, hjust = c(0, -9), vjust = c(0, 3.54)
  )

ggsave(paste0(folderFig, "temporalDynamicsFigure.jpg"), 
       plot = g1,
       width = 7.5,
       height = 10)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                         WITHIN POPULATION ANALYSIS                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Compute relative male allocation through calibration and analysing the trade-off:

mod.pop5 = sma(nbMfl_sd ~ nbFruit_sd, data = df.pop5)
summary(mod.pop5)

df.pop5 = df.pop5 %>% 
  mutate(relMRE = nbMfl_sd/(nbMfl_sd + nbFruit_sd*-coef(mod.pop5)[2]))


## Final plotting --------------------------------------------------------------

g2 = df.pop5 %>% 
  ggplot(aes(x = nbFruit/subsample_dry_mass, y = nbMfl/subsample_dry_mass, col = relMRE)) +
  geom_abline(slope = coef(mod.pop5)[2], intercept = coef(mod.pop5)[1],
              linewidth = 1, linetype = "dashed", col = "grey30") +
  geom_point() +
  scale_color_gradientn(colours = cols[3:6], name = "Rellative male\nallocation") +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.85, 0.75),
        legend.text.position = "left",
        text = element_text(size=15),
        legend.title = element_text(margin = margin(b = 15), hjust = 1)) +
  ylab("Number of male flower") +
  xlab("Number of fruits")


ggsave(paste0(folderFig, "withinPopFigure.jpg"), 
       plot = g2,
       width = 7.5,
       height = 5.75)





