##########
# This script is used to reproduce the results from the manuscript entitled
# "Convergent evolution of fern nectaries facilitated independent recruitment of
# ant-bodyguard from flowering plants" authored by JS Suissa, FW Li, and C Moreau,
# Nature Communications 2024
##########

# Load libraries
remotes::install_github("joelnitta/taxastand")
install.packages("taxon-tools")
install.packages("ftolr")
library(tidyverse)
library(taxastand)
library(devtools)
library(phytools)
library(ape)
library(reshape2)
library(ftolr)
library(geiger)
library(Taxonstand)
library(ggsci)
library(scales)
library(taxize)
library(ggtree)
library(BAMMtools)
library(hisse)
library(diversitree)
library(taxize)
library(taxizedb)
library(brew)
library(corHMM)
library(MEDUSA)

# Set working directory (replace with your desired directory)
setwd("path/to/working/directory")

# Build the fern phylogeny or upload from GitHub repository. Note that the FTOL is continuously updated and will be different from the one used in this manuscript which was downloaded August 2022.
#### Use the 5k fern phylogeny from the Fern Tree Of Life Project ####

#Generate an ancestral character state reconstruction for the fern genus tree
fern_genus_tree <- read.tree("/fern_5k_genus_tree.tre")

# Bring in the nectary data
genus_nectary.dat <- read.csv("data/fern_5k_genus_nectary_data")

# Filter tips in the tree to just match tips in the dataframe
rownames(genus_nectary.dat) <- genus_nectary.dat[, 6]

# select nectaries
fern_genus_nectar_reduced.dat <- genus_nectary.dat %>%
  select(Nectary)

# Check species in tree
chk <- name.check(fern_genus_tree, fern_genus_nectar_reduced.dat)

# Now prune the tree based on the data available
tree.pruned <- drop.tip(fern_genus_tree, chk$tree_not_data)
fern.Data.pruned <-
  fern_genus_nectar_reduced.dat[!(rownames(fern_genus_nectar_reduced.dat) %in% chk$data_not_tree), , drop = FALSE]

# Format the data
nectary.mode <-
  setNames(fern.Data.pruned[, 1], rownames(fern.Data.pruned))
nectary.mode <- as.factor(nectary.mode)

# Run ancestral character estimation, specify the root node as non-nectary bearing
nectary_genus_simmap <- make.simmap(
  tree = tree.pruned,
  model = "ARD",
  x = nectary.mode,
  Q = "empirical",
  nsim = 1000,
  burnin = 500,
  pi = c(1, 0)
)

# Summarize data
pd_genus <- describe.simmap(nectary_genus_simmap, plot = FALSE)

# Make a density map
dMapGenus <- densityMap(nectary_genus_simmap, plot = FALSE)

# Change the colors
n <- length(dMapGenus$cols)
dMapGenus$cols[1:n] <-
  colorRampPalette(c("gray", "#00A08799"), space = "Lab")(n)

# Plot the density map
plot(
  dMapGenus,
  lwd = 0.8,
  fsize = c(0.15, 0.7),
  ftype = "i",
  colors = cols,
  outline = FALSE
)

setEnv = TRUE
nodelabels(
  node = as.numeric(row.names(pd_genus$ace)),
  pie = pd_genus$ace,
  piecol = c("gray", "#00A08799"),
  cex = 0.15,
  frame = "none"
)


#Anc recon for fern species.
fern_5k_tree <- read.tree("/fern_5k_tree.tre")

# Bring in species nectary data
species_nectar.dat <-
  read.csv("/Fern_Nectary_Database_Nov19_2022.csv")

# Filter tips in the tree to just match tips in the dataframe
rownames(species_nectar.dat) <- species_nectar.dat$Species

# Drop all tips that are not No or Yes
ss <- species_nectar.dat %>%
  filter(Nectary == "Yes" | Nectary == "No") %>%
  select(Nectary)

# Check species in tree
chk <- name.check(fern_5k_tree, ss)

# Now prune the tree based on the data available
tree.pruned <- drop.tip(fern_5k_tree, chk$tree_not_data)
fern.Data.pruned <-
  ss[!(rownames(ss) %in% chk$data_not_tree), , drop = FALSE]

nectary.mode <-
  setNames(fern.Data.pruned[, 1], rownames(fern.Data.pruned))
nectary.mode <- as.factor(nectary.mode)


# Test models
fit_nectary_ER <- fitMk(tree.pruned,
                        nectary.mode,
                        model = "ER",
                        pi = c(1, 0))

fit_nectary_ARD <- fitMk(tree.pruned,
                         nectary.mode,
                         model = "ARD",
                         pi = c(1, 0))

aic <- c(AIC(fit_nectary_ARD),
         AIC(fit_nectary_ER))

model.fit <- data.frame(
  model = c("fit_nectary_ARD", "fit_nectary_ER"),
  logL = c(logLik(fit_nectary_ARD), logLik(fit_nectary_ER)),
  AIC = aic,
  delta.AIC = aic - min(aic)
)

# Arrange the models according to their delta AIC
model.fit <- model.fit %>%
  arrange(delta.AIC)

mod_ard <- fit_nectary_ARD$rates

# Fit the data that is just for all species that have nectaries and only 1 species for each genus that lacks it
nectary_simmap <- make.simmap(
  tree = tree.pruned,
  model = "ARD",
  x = nectary.mode,
  Q = "empirical",
  nsim = 1000,
  burnin = 500,
  pi = c(1, 0)
)

pd <- describe.simmap(nectary_simmap, plot = FALSE)


# LTT plots

# Make a density map for ferns
dMap <- densityMap(nectary_simmap, plot = FALSE)


# Substitute the colors in dmap
n <- length(dMap$cols)
dMap$cols[1:n] <-
  colorRampPalette(c("gray", "black"), space = "Lab")(n)
dMap$cols[1:n] <-
  colorRampPalette(c("gray", "#26cf0e99"), space = "Lab")(n)


# Sample 100 nectary simmaps
simap_sample <- sample(nectary_simmap, 100, replace = FALSE)
fern_nectaries.ltts <- ltt(simap_sample)
ltt_fern_nec_phytools <-
  ltt.plot.coords(tree.pruned_yes, type = "s")
ltt_fern_nec_phytools[, 1] <- abs(ltt_fern_nec_phytools[, 1])

# Set colors
cols <- setNames(c("gray", "#26cf0e99"), levels(nectary.mode))






#Code for generating Ant figures
# Bring in the ant tree
anttree <-
  read.tree("/Nelson_anttree_ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

# Read in the csv and drop all tips that don't have plant associations
antPlant.dat <- read.csv("/Ant_plant_data.csv")

# Filter out species that only have plant interactions
antPlant1.dat <- antPlant.dat %>%
  filter(
    DietModanyplant == 1 |
      ForagingModanycanopy == 1 |
      NestingModanycanopy == 1 | DietModanyomni == 1,
    .preserve = TRUE
  ) %>%
  dplyr::select(Taxon,
                Genus,
                Family,
                Subfamily,
                Tribe,
                Order,
                NestingModanycanopy)

# Make a reduced ant tree with only species that are plant-associated
# Filter tips in the tree to just match tips in the dataframe
rownames(antPlant1.dat) <- antPlant1.dat[, 1]

# Just select nectary column
antPlant2.dat <- antPlant1.dat %>%
  dplyr::select(NestingModanycanopy)

# Check species in tree
chk_ant <- name.check(anttree, antPlant2.dat)

# Now prune the tree based on the data available
anttreePruned <- drop.tip(anttree, chk_ant$tree_not_data)
ltt_ants_phytools <- ltt.plot.coords(anttreePruned, type = "s")
ltt_ants_phytools[, 1] <- abs(ltt_ants_phytools[, 1])
ant_ult <- force.ultrametric(anttreePruned)


# Make a new dataset if species interact with plants give them a 1 and if not give them a 0
antPlant2.dat <- antPlant.dat %>%
  mutate(
    plant_interact = ifelse(
      DietModanyplant == 1 |
        ForagingModanycanopy == 1 |
        NestingModanycanopy == 1 |
        DietModanyomni == 1,
      1,
      0
    )
  ) %>%
  select(Taxon, plant_interact)

# Check species in tree
# Filter tips in the tree to just match tips in the dataframe
rownames(antPlant2.dat) <- antPlant2.dat[, 1]

# Just select nectary column
antPlant3.dat <- antPlant2.dat %>%
  dplyr::select(plant_interact)

chk_ant <- name.check(anttree, antPlant3.dat)

ant.mode <- setNames(antPlant3.dat[, 1], rownames(antPlant3.dat))
ant.mode <- as.factor(ant.mode)

# Fit the data that is just for all species that have nectaries and only 1 species for each genus that lacks it
ant_plant_simmap <- make.simmap(
  tree = anttree,
  model = "ARD",
  x = ant.mode,
  Q = "empirical",
  nsim = 1000,
  burnin = 500
)

#Ant Simmaps
plotSimmap(
  ant_plant_simmap[[50]],
  ftype = "i",
  type = "fan",
  fsize = 0.1,
  lwd = 0.5,
  outline = FALSE
)

#Ant portion
ant_simap_sample <- sample(ant_plant_simmap, 100, replace = FALSE)
ant_plant.ltts <- ltt(ant_simap_sample)

# Set colors
cols <- setNames(c("gray", "#ff0e85"), c(0, 1))

# Create graph
plot(
  ant_plant.ltts,
  colors = cols,
  show.total = FALSE,
  xlim = c(0, max(nodeHeights(ant_plant_simmap[[1]]))),
  ylim = c(1, 2000),
  bty = "n",
  cex.axis = 0.8,
  las = 1,
  xlab = "millions of years (from the root)",
  ylab = "lineages",
  axes = FALSE,
  log = "y"
)

axis(1, at = round(seq(0, max(
  nodeHeights(ant_plant_simmap[[1]])
), length.out = 5), 0), cex.axis = 0.8)
axis(2, las = 1, cex.axis = 0.8)
grid()

antpd <- describe.simmap(ant_plant_simmap, plot = FALSE)

#Plot densimap with arc labels
dmap_ants <- densityMap(ant_plant_simmap, plot = FALSE)

plot(
  dmap_ants,
  lwd = 0.4,
  type = "fan",
  fsize = c(0.15, 0.7),
  ftype = "off",
  outline = FALSE
)

setEnv = TRUE
nodelabels(
  node = as.numeric(row.names(antpd$ace)),
  pie = antpd$ace,
  piecol = c("blue", "red"),
  cex = 0.1
)





#Code for generating Angiosperm figures
# Bring in angiosperm nectary data
tree_sp_nectary.dat <-
  read.csv("/Angiosperm_nectaries_with_tree.csv")

# Swap space with _
tree_sp_nectary.dat$Species <-
  gsub(" ", "_", tree_sp_nectary.dat$Species)

ang_nec_YES.dat <- tree_sp_nectary.dat %>%
  filter(EFN == "Yes") %>%
  dplyr::select(Species, EFN) %>%
  distinct(Species, .keep_all = TRUE) %>%
  column_to_rownames(var = "Species")

ang_nec.dat <- tree_sp_nectary.dat %>%
  mutate(EFN = ifelse(is.na(EFN), "No", "Yes")) %>%
  dplyr::select(Species, EFN) %>%
  distinct(Species, .keep_all = TRUE) %>%
  column_to_rownames(var = "Species")

ang_nect.mode <- setNames(ang_nec.dat$EFN, rownames(ang_nec.dat))
ang_nect.mode <- as.factor(ang_nect.mode)

# Bring in angiosperm tree
ang_tree <- read.tree("/Smith_and_Brown_angiosprm_ALLOTB.tre")

# Subset the tree to only have taxa in ang_nec.dat
# Check species in tree
chk_ang <- name.check(ang_tree, ang_nec.dat)

# Now prune the tree based on the data available and make LTT plots
angtreePruned <- drop.tip(ang_tree, chk_ang$tree_not_data)
ult_tree <- force.ultrametric(angtreePruned)
ltt_nectary_angios <- ltt.plot.coords(ult_tree, type = "s")
ltt_nectary_angios[, 1] <- abs(ltt_nectary_angios[, 1])
ltt(ult_tree, log.lineages = FALSE)

# Test the diversification rate of nectary-bearing angiosperms

# Now prune the tree based on the data available
ff <- branching.times(ult_tree)
diversi.time(ff, census = NULL, censoring.codes = c(1, 0))


# Make angiosperm simmaps figure S5
ang_tips <- read.csv("/new_reduced_ang_tipsmay21.csv")
ang_tree <- read.tree("/filtered_90_pruned_ult_ang_treemay1923.tre")

# Remove acrogymnosperm tips
gym.tips <- c(
  "Cycas_multipinnata",
  "Tsuga_heterophylla",
  "Baiera_bidens",
  "Encephalartos_hirsutus",
  "Pachylepis_commersonia",
  "Sequoiadendron_giganteum",
  "Chamaecyparis_lawsoniana",
  "Cupressus_sempervirens",
  "Juniperus_phoenicea",
  "Phyllocladus_trichomanoides",
  "Halocarpus_biformis",
  "Saxegothaea_conspicua"
)

# Drop the tips from the tree
tree_pruned_3 <- drop.tip(tree_pruned_2, gym.tips)

# Now run a series of code that prunes the tree to include only 1 tip per clade with the same trait
tree_pruned_resolved <- multi2di(tree_pruned_3)

# Force ultrametric (this is just a branch length rounding issue)
ult_ang_tree <- force.ultrametric(tree_pruned_resolved)

# Drop corresponding rows from data frame
df_final <- df_sub[df_sub$Species %in% ult_ang_tree$tip.label,]

df_1 <- df_final %>%
  remove_rownames() %>%
  column_to_rownames(var = "Species")

# Check species in tree
chk_ang <- name.check(ult_ang_tree, df_1)

ang_nect.mode <- setNames(df_1[, 1], rownames(df_1))
ang_nect.mode <-
  as.factor(ang_nect.mode) # This is essential now for some reason...

angnectary_species_simmap_ARD <- make.simmap(
  tree = ult_ang_tree,
  model = "ARD",
  x = ang_nect.mode,
  Q = "empirical",
  nsim = 1000,
  burnin = 500,
  pi = c(1, 0)
)


simap_sample_ang <-
  sample(angnectary_species_simmap_ARD, 100, replace = FALSE)

# Angiosperm LTT
ang_nectaries_subsample_pruned.ltts <-
  phytools::ltt.multiSimmap(simap_sample_ang)


# Set colors
cols <- setNames(c("blue", "#ff1f12"), c("No", "Yes"))

# Create graph
plot(
  ang_nectaries_subsample_pruned.ltts,
  colors = cols,
  show.total = FALSE,
  ylim = c(0, 8),
  bty = "n",
  cex.axis = 0.8,
  las = 1,
  xlab = "millions of years (from the root)",
  ylab = "lineages",
  log = "y"
) +
  xlim = c(150, 317) +
  axis(1, cex.axis = 0.8) +
  axis(2, las = 1, cex.axis = 0.8) +
  grid()

# Plot simmap
angpd <-
  describe.simmap(angnectary_species_simmap_ARD, plot = FALSE)

plot(
  dMap_ang_nectaries_Total,
  lwd = 0.4,
  type = "fan",
  fsize = c(0.15, 0.7),
  ftype = "off",
  outline = FALSE
)

setEnv = TRUE
nodelabels(
  node = as.numeric(row.names(angpd$ace)),
  pie = angpd$ace,
  piecol = c("blue", "red"),
  cex = 0.1
)





# Bring in herbivore data to make LTT plots
herbGenTree <- read.tree("/fern_herb_genus_NCBI_syn.nwk")

herbgenTree2 <- drop.tip(
  herbGenTree,
  tip = c(
    "Hyposmocoma_carnivora",
    "Hyposmocoma_kaikuono",
    "Hyposmocoma_kaupo"
  )
)

ltt_herb_genera <- ltt.plot.coords(herbgenTree2, type = "s")
ltt_herb_genera[, 1] <- abs(ltt_herb_genera[, 1])


# Make the LTT plots for all species
show_col(pal_npg("npg", alpha = 0.6)(10))
show_col(pal_npg("nrc", alpha = 0.6)(2))
pal_jama("default", alpha = 0.6)(7)[c(2, 6)]

plot(
  NA,
  xlim = rev(range(times)),
  ylim = c(1, 8),
  xlab = "time",
  ylab = "lineages",
  bty = "n",
  las = 1
)

# Plot the ltt
cols_ferns <- setNames(c("gray", "#26cf0e99"), c("No", "Yes"))
plot(
  fern_nectaries.ltts,
  colors = cols_ferns,
  show.total = FALSE,
  axes = FALSE,
  xlim = c(0, (max(nodeHeights(
    dMap[[1]]
  )))),
  bty = "n",
  cex.axis = 0.8,
  las = 1,
  xlab = "millions of years (from the root)",
  log.lineages = TRUE,
  lwd = 3,
  ylim = c(1, 8)
)

axis(1, at = round(seq(0, max(nodeHeights(
  dMap[[1]]
)), length.out = 5), 0), cex.axis = 0.8)
axis(2, las = 1, cex.axis = 0.8)
plot.window(xlim = rev(c(0, max(nodeHeights(
  dMap[[1]]
)))), ylim = c(1, 8))

# Plot angios
cols_ang <- setNames(c("gray", "#4dbbd599"), c("No", "Yes"))
plot(
  ang_nect_plant.ltts,
  colors = cols_ang,
  show.total = FALSE,
  add = TRUE,
  axes = FALSE,
  xlim = c(0, max(nodeHeights(dMap[[1]]))),
  bty = "n",
  cex.axis = 0.8,
  las = 1,
  log.lineages = TRUE,
  lwd = 3,
  ylim = c(1, 8)
)

axis(1, at = round(seq(0, max(
  nodeHeights(angnectary_species_simmap_ARD[[1]])
), length.out = 5), 0), cex.axis = 0.8)
axis(2, las = 1, cex.axis = 0.8)

# Plot ants
cols_ants <- setNames(c("gray", "#e64b3599"), c(0, 1))
plot(
  ant_plant.ltts,
  colors = cols_ants,
  show.total = FALSE,
  xlim = c(0, max(nodeHeights(dMap[[1]]))),
  add = TRUE,
  axes = FALSE,
  bty = "n",
  cex.axis = 0.8,
  las = 1,
  log.lineages = TRUE,
  lwd = 3,
  ylim = c(1, 8)
)

axis(1, at = round(seq(0, max(
  nodeHeights(ant_plant_simmap[[1]])
), length.out = 5), 0), cex.axis = 0.8)
axis(2, las = 1, cex.axis = 0.8)

# Plot the density map
plot(
  dMap$tree,
  dMap$cols,
  0.5,
  lwd = 0.5,
  ftype = "off",
  add = TRUE,
  mar = par()$mar,
  tips = setNames(1:length(dMap[[1]]$tip.label), dMap[[1]]$tip.label)
)


plot.window(xlim = rev(c(0, max(nodeHeights(
  dMap[[1]]
)))), ylim = c(1, 8))

plot(
  NA,
  xlim = rev(c(0, max(nodeHeights(
    dMap[[1]]
  )))),
  ylim = c(1, 8),
  xlab = "time",
  ylab = "lineages",
  bty = "n",
  las = 1
)

# Plot species herbivores
lines(
  ltt_herb_genera[, 1],
  log(ltt_herb_genera[, 2]),
  type = "s",
  col = "#b09c8599",
  lwd = 3
)

# Plot ant lines
lines(
  (ltt_ants_phytools[, 1]),
  log(ltt_ants_phytools[, 2]),
  type = "s",
  col = pal_npg("nrc", alpha = 0.6)(1),
  lwd = 3
) # plot ants

# Plot flowering plant lines
lines(
  ltt_nectary_angios[, 1],
  log(ltt_nectary_angios[, 2]),
  type = "s",
  col = pal_npg("nrc", alpha = 0.6)(2)[2],
  lwd = 3
) # flowering plants with nectaries

# Plot fern species
lines(
  ltt_fern_nec_phytools[, 1],
  log(ltt_fern_nec_phytools[, 2]),
  type = "s",
  col = "#26cf0e99",
  lwd = 3
)




#Code for generating corHMM plots

# Bring in the tree
fern_5k_tree <- read.tree("/fern_5k_tree.tre")

# Bring in the data
hab.dat <- read.csv("/fern_trait_database.csv")

data <- hab.dat %>%
  filter(Species %in% fern_5k_tree$tip.label, .preserve = TRUE) %>%
  select(Species, Growth_habit_E1_T0_Arb2_climb_3, Nectary) %>%
  filter(Nectary == "Yes" | Nectary == "No") %>%
  filter(!is.na(Growth_habit_E1_T0_Arb2_climb_3)) %>%
  distinct(Species, .keep_all = TRUE)

# Check species in tree
dat_names <- column_to_rownames(data, var = "Species")
chk <- name.check(fern_5k_tree, dat_names)

# This is a crazy
tree.pruned <- drop.tip(fern_5k_tree, chk$tree_not_data)
data$Growth_habit_E1_T0_Arb2_climb_3 <-
  as.character(data$Growth_habit_E1_T0_Arb2_climb_3)

# Test for correlated evolution
a <-
  fitCorrelationTest(tree.pruned, data, simplified_models = FALSE)


plot(a$hidden_Markov_independent_model_fit)
corHMM::plotMKmodel(a$hidden_Markov_independent_model_fit,
                    color = "col.blind",
                    display = "square")

# Try another way
par(mar = rep(2, 4))
rate_cat_mat <- corHMM:::getRateCatMat(2)
indep_model_1 <-
  getStateMat4Dat(data, "ARD", collapse = FALSE, indep = TRUE)$rate.mat
indep_model_2 <-
  getFullMat(list(indep_model_1, indep_model_1), RateClassMat = rate_cat_mat)
hidden_independent_model_fit <-
  corHMM(
    phy = tree.pruned,
    data = data,
    rate.cat = 2,
    rate.mat = indep_model_2,
    upper.bound = 1000
  )

hidden_independent_model_fit_10 <-
  corHMM(
    phy = tree.pruned,
    data = data,
    rate.cat = 2,
    rate.mat = indep_model_2,
    upper.bound = 10000
  )


corHMM::plotMKmodel(hidden_independent_model_fit_10,
                    color = "col.blind",
                    display = "square")

# Rerun the corhmm with growth habit as "canopy" and "understory"
data <- hab.dat %>%
  filter(Species %in% fern_5k_tree$tip.label, .preserve = TRUE) %>%
  filter(Nectary == "Yes" | Nectary == "No", .preserve = TRUE) %>%
  filter(!is.na(Growth_habit_E1_T0_Arb2_climb_3), .preserve = TRUE) %>%
  mutate(GH_mut = ifelse(Growth_habit_E1_T0_Arb2_climb_3 == 0, 0, 1)) %>% # mutate the data so that it's either understory or canopy
  select(Species, GH_mut, Nectary) %>%
  distinct(Species, .keep_all = TRUE)

# Check species in tree
dat_names <- column_to_rownames(data, var = "Species")
chk <- name.check(fern_5k_tree, dat_names)

# This is a crazy
tree.pruned <- drop.tip(fern_5k_tree, chk$tree_not_data)
data$GH_mut <- as.character(data$GH_mut)

# Test for correlated evolution
b <-
  fitCorrelationTest(tree.pruned, data, simplified_models = FALSE)

# Plot it
corHMM::plotMKmodel(b$hidden_Markov_correlated_model_fit,
                    color = "col.blind",
                    display = "square")

# Ask about the total number of changes across the states
# Make simmaps and then coplot them
model <- b$hidden_Markov_correlated_model_fit$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)

# Change model column names
phy <- b$hidden_Markov_correlated_model_fit$phy
states <- b$hidden_Markov_correlated_model_fit$states
tip.states <- b$hidden_Markov_correlated_model_fit$tip.states
rawdat <- b$hidden_Markov_correlated_model_fit$data

# Run get simmap
colnames(model) <- c("1", "2", "3", "4", "5", "6", "7", "8")
rownames(model) <- c("1", "2", "3", "4", "5", "6", "7", "8")
colnames(tip.states) <- c("1", "2", "3", "4", "5", "6", "7", "8")

nec_growth_simmap <- make.simmap(
  tree = b$hidden_Markov_correlated_model_fit$phy,
  x = tip.states,
  Q = model,
  nsim = 100,
  pi = c(1, 0, 0, 0, 1, 0, 0, 0)
)

# Plot it
plotSimmap(nec_growth_simmap, ftype = "off")

# Plot all simmaps on top of each other
cols <- setNames(
  c(
    "darkblue",
    "navajowhite4",
    "green4",
    "goldenrod",
    "darkblue",
    "navajowhite4",
    "green4",
    "goldenrod"
  ),
  c("1", "2", "3", "4", "5", "6", "7", "8")
)

pd <- describe.simmap(nec_growth_simmap)
simmap_counts <- countSimmap(nec_growth_simmap)


# Pull out 100 random trees:
sampled_simmaps <- sample(nec_growth_simmap, 100, replace = FALSE)
plotTree(
  nec_growth_simmap[[50]],
  ftype = "off",
  type = "fan",
  fsize = 0.5,
  offset = 0.5,
  lwd = 0.6
)
par(fg = "transparent", lend = 1)
plotTree(
  nec_growth_simmap[[50]],
  ftype = "off",
  type = "fan",
  fsize = 0.5,
  offset = 0.5,
  lwd = 0.4,
  color = "white",
  add = TRUE
)

for (i in 1:length(sampled_simmaps)) {
  plot(
    sampled_simmaps[[i]],
    ftype = "off",
    type = "fan",
    colors = sapply(cols, make.transparent, alpha = 0.1),
    add = TRUE,
    lwd = 0.4,
    fsize = 0.5,
    offset = 0.5
  )
}

setEnv <-
  TRUE 

# Add node labels
par(fg = "black") # this makes the piecharts have a black ring around them
par(fg = "transparent") # this might help with eliminating the black ring on the piecharts

# Make another plot with the node labels
plot(
  sampled_simmaps[[50]],
  ftype = "off",
  type = "fan",
  colors = cols,
  lwd = 0.4,
  fsize = 0.5,
  offset = 0.5
)
setEnv <- TRUE
nodelabels(
  pie = pd$ace,
  piecol = cols,
  cex = 0.15,
  frame = "none"
)

# Plot tips on tree
library(diversitree)
rawdat <- b$hidden_Markov_correlated_model_fit$data
rawdat1 <- rawdat %>%
  column_to_rownames(var = "Species") %>%
  mutate(GH_bin = as.numeric(GH_mut)) %>%
  mutate(nec_bin = ifelse(Nectary == "No", 0, 1)) %>%
  select(GH_bin, nec_bin)

trait.plot(
  sampled_simmaps[[50]],
  rawdat1,
  cols = list(
    GH_bin = c("pink", "red"),
    nec_bin = c("lightblue", "blue")
  ),
  legend = FALSE,
  cex.lab = 0.000000000000000000000001
)





#MEDUSA analyses

#read fern tree
fern_5k_tree <-
  read.tree("/fern_5k_tree.tre")

#bring in species nectary data
species_nectar.dat <-
  read.csv("/Fern_Nectary_Database_Nov19_2022.csv")

# filter tips in the tree to just match tips in the dataframe
rownames(species_nectar.dat) <- species_nectar.dat$Species

#drop all tips that are not No or Yes
ss <- species_nectar.dat %>%
  filter(Nectary == "Yes" | Nectary ==  "No") %>%
  select(Nectary)

# check species in tree
chk <-
  name.check(fern_5k_tree, ss)# this is a crazy awesome and useful function

# now prune the tre based on the data available
tree.pruned <- drop.tip(fern_5k_tree, chk$tree_not_data)

# now prune the data based on the data available
fern.Data.pruned_labs <-
  species_nectar.dat[!(rownames(species_nectar.dat) %in%
                         chk$data_not_tree), , drop = FALSE]

#filter by the medusa tree
chk_1 <-
  name.check(summ$summary.tree, species_nectar.dat)# this is a crazy awesome and useful function

tree.pruned_1 <- drop.tip(summ$summary.tree, chk_1$tree_not_data)


# now prune the data based on the data available
fern.Data.pruned_labs <-
  species_nectar.dat[!(rownames(species_nectar.dat) %in%
                         chk_1$data_not_tree), , drop =
                       FALSE]


#RUn Medusa

fie <- c(tree.pruned, tree.pruned)
res1 <- MEDUSA(fie,
               mc = TRUE,
               model = "bd",
               criterion = "aicc")
summ <-
  multiMedusaSummary(res1,
                     tree.pruned,
                     plotModelSizes = TRUE,
                     plotTree = FALSE)

#grab all the nodes where there is a positive rate shift
medusa_spp_upshifts <-
  as.data.frame(summ$shift.summary[summ$shift.summary[, 4] > 0, ])

#Grab all the negative shifts
medusa_spp_downshifts <-
  as.data.frame(summ$shift.summary[summ$shift.summary[, 4] < 0, ])

medusa_spp_rate_shifts <- medusa_spp_rate_shifts %>%
  mutate(Node = as.character(shift.node))

#Make total shifts datafram

total_medusa_rateshifts <-
  rbind(medusa_spp_upshifts, medusa_spp_downshifts)


plotMultiMedusa(
  summ,
  type = "fan",
  time = FALSE,
  annotateShift = TRUE,
  annotateRate = "r.mean",
  show.tip.label = FALSE
)


#Run code using a all fern genera and the species richness per genera based on Pteridophyte Phylogeny Group

#Read in the genus tree
phy <-
  read.tree("/fern_5k_genus_tree.tre")

#file specifying species richness for the genera used as tips
richness <-
  read.csv("/fern_genera_species_richness.csv")

#tricking MEDUSA to run two runs to make the summary easier for the big tree
fie2 <- c(phy, phy)
res_gen <-
  MEDUSA(fie2,
         richness,
         mc = TRUE,
         model = "bd",
         criterion = "aicc")
sum_gen <-
  multiMedusaSummary(res_gen, phy, plotModelSizes = TRUE, plotTree = FALSE)

plotMultiMedusa(
  sum_gen,
  annotateShift = TRUE,
  annotateRate = "r.mean",
  tip.cex = 0.2
)

#get upshifts and downshifts at the genus level
medusa_gen_upshifts <-
  as.data.frame(sum_gen$shift.summary[sum_gen$shift.summary[, 4] > 0, ])

#Grab all the negative shifts
medusa_gen_downshifts <-
  as.data.frame(sum_gen$shift.summary[sum_gen$shift.summary[, 4] < 0, ])


#Grab all the descendent nodes from the simmap tree using the getDescendants functioni in phytools following Landis et al., 2017


#decribe fern nectary summary
pd <- describe.simmap(nectary_simmap)


a <- as.data.frame(pd$ace[, 2])

b <- a %>%
  mutate(pp = as.numeric(`pd$ace[, 2]`)) %>%
  select(pp) %>%
  filter(pp >= 0.9)

#Grab node names
nodes_d <- rownames(b)


# Function to get the descendants of a node
#get the descendents

result_df <-
  data.frame(matrix(ncol = 5))  # Create an empty dataframe with 5 columns
colnames(result_df) <-
  c("Node",
    "Descendant1",
    "Descendant2",
    "Descendant3",
    "Descendant4")  # Set column names

for (i in 1:length(nodes_d)) {
  z <- getDescendants(pd$tree[[1]], node = nodes_d[i])
  
  if (length(z) >= 1) {
    z <- rep(NA, 4)  # Create a vector of NAs for descendants
    z[1:length(z)] <-
      getDescendants(pd$tree[[1]], node = nodes_d[i])[1:4]  # Assign first 1-4 descendants
    nod.name <- nodes_d[i]
    row <- c(nod.name, z)
    result_df <- rbind(result_df, row)
  }
}

result_df <- result_df[-1, ]  # Remove the first empty row

#Make a dataframe of only the exact nodes where

exact_nectary_gain_nodes <- result_df$Node

#rename medusa_spp_upshifts#shift.node
medusa_spp_upshifts_corr <- medusa_spp_upshifts %>%
  mutate(Node = as.character(shift.node))

#make a list of the exact shifts and nectary gain nodes
exact_shifts_and_nect_gains <-
  left_join(result_df, medusa_spp_upshifts_corr, by = "Node")

#Only one node has a correlation, this is Cyatheaceae

#Check the single node with a correlation

indexed.tree <- extract.clade(summ$summary.tree, node = 8655)


#plot the phylogeny and the nodes

nectary_nodes = as.numeric(exact_nectary_gain_nodes)
down_nodes = medusa_spp_downshifts$shift.node
up_nodes = medusa_spp_upshifts$shift.node

plot(
  tree.pruned ,
  type = "fan",
  font = 1,
  adj = 0.5,
  no.margin = TRUE,
  edge.width = 0.5,
  show.tip.label = FALSE
)


plotMultiMedusa(
  summ,
  type = "fan",
  time = FALSE,
  annotateShift = TRUE,
  annotateRate = "r.mean",
  show.tip.label = FALSE,
  shift.leg.pos = "bottomleft"
)


setEnv = TRUE

for (i in 1:length(nectary_nodes))
  nodelabels(
    "W",
    nectary_nodes[i],
    frame = "circle" ,
    col = "springgreen",
    bg = "springgreen",
    cex = 0.25
  )
for (i in 1:length(down_nodes))
  nodelabels(
    "d",
    down_nodes[i],
    frame = "circle" ,
    col = "blue",
    bg = "blue",
    cex = 0.3
  )
for (i in 1:length(up_nodes))
  nodelabels(
    "u",
    up_nodes[i],
    frame = "circle" ,
    col = "red",
    bg = "red",
    cex = 0.3
  )



####
#take the subset of the fern tree showing the lineages with nectaries and ask when the rate shifts occur and then take the herbivore tree and ask when the herbivore shifts occur and see if they correlate.
####


#Read in the fern tree and prune to only keep nectary bearing lineages.
fern.tot.tree <-
  read.tree("/fern_5k_tree.tre")

#bring in species nectary data
species_nectar.dat <-
  read.csv("/Fern_Nectary_Database_Nov19_2022.csv")

# filter tips in the tree to just match tips in the dataframe
rownames(species_nectar.dat) <- species_nectar.dat$Species

#drop all tips that are not No or Yes
ss <- species_nectar.dat %>%
  filter(Nectary == "Yes") %>%
  select(Nectary)

# check species in tree
chk <-
  name.check(fern.tot.tree, ss)

# now prune the tre based on the data available
nec.tre <- drop.tip(fern.tot.tree, chk$tree_not_data)

# now prune the data based on the data available
fern.Data.pruned_labs <- ss[!(rownames(ss) %in%
                                chk$data_not_tree), , drop =
                              FALSE]

fie_fern_nec <- c(nec.tre, nec.tre)
res1_fern_nec <-
  MEDUSA(fie_fern_nec,
         mc = TRUE,
         model = "bd",
         criterion = "aicc")
summ_fern_nec <-
  multiMedusaSummary(res1_fern_nec,
                     nec.tre,
                     plotModelSizes = TRUE,
                     plotTree = FALSE)

plotMultiMedusa(
  summ_fern_nec,
  type = "fan",
  time = FALSE,
  annotateShift = TRUE,
  annotateRate = "r.mean",
  show.tip.label = FALSE,
  shift.leg.pos =  "bottomleft"
)


# Get the time of each node in the tree
node_times <- node.depth.edgelength(nec.tre)

# Map each shift node to its corresponding time
shift_times <-
  356.1544 - node_times[summ_fern_nec$shift.summary[, 1]] #max age of tree

# Define the different window sizes and step sizes you want to use
window_sizes <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
step_sizes <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

# Number of bootstrap replicates
n_boot <- 100

# Initialize a list to hold the data frames for each window size
dfs <- list()

# Iterate over the different window sizes
for (j in seq_along(window_sizes)) {
  window_size <- window_sizes[j]
  step_size <- step_sizes[j]
  start_times <- seq(0, max(shift_times), by = step_size)
  
  # Initialize a matrix to hold the counts of rate shifts for each bootstrap replicate
  rate_shift_counts <-
    matrix(nrow = n_boot, ncol = length(start_times))
  
  # For each bootstrap replicate
  for (b in 1:n_boot) {
    # Resample the shift times with replacement
    resampled_shift_times <- sample(shift_times, replace = TRUE)
    
    # For each sliding window, count the number of resampled shift times that fall within the window
    for (i in seq_along(start_times)) {
      start_time <- start_times[i]
      end_time <- start_time + window_size
      rate_shift_counts[b, i] <-
        sum(resampled_shift_times >= start_time &
              resampled_shift_times < end_time)
    }
  }
  
  # Convert the matrix to a data frame in long format
  df <- as.data.frame(rate_shift_counts)
  names(df) <- start_times
  df$bootstrap_replicate <- 1:n_boot
  df$window_size <- window_size
  df_long <-
    gather(df,
           start_time,
           rate_shift_count,-bootstrap_replicate,-window_size)
  
  # Convert the start_time column to numeric
  df_long$start_time <- as.numeric(as.character(df_long$start_time))
  
  # Add the data frame to the list
  dfs[[j]] <- df_long
}

# Combine all the data frames into one
df_combined <- bind_rows(dfs)


#Now we need to summarize the data and get the percent of total shifts but we need to group by window size and count the total number of shifts

# Group by bootstrap replicate and calculate the total number of shifts for each replicate
df_summ <- df_combined %>%
  group_by(bootstrap_replicate, window_size) %>%
  mutate(total_shifts = sum(rate_shift_count)) %>%
  mutate(percentage_of_total_shifts = rate_shift_count / total_shifts * 100)# Calculate the percentage of total shifts for each time window within each bootstrap replicate


# Calculate the mean percentage of total shifts for each time window across the bootstrap replicates
df_summary <- df_summ %>%
  group_by(start_time) %>%
  summarise(
    mean_percentage_of_total_shifts = mean(percentage_of_total_shifts),
    .groups = "drop"
  )


#plot barplots 
ggplot(df_summary,
       aes(x = start_time, y = mean_percentage_of_total_shifts)) +
  geom_bar(
    stat = "identity",
    fill = "steelblue",
    alpha = 0.5,
    width = step_size
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "red",
    size = 1.5
  ) +
  labs(x = "Time", y = "Mean total shifts (%)") +
  theme_minimal() +
  scale_x_reverse()



#Bring in the herbivore tree and run medusa
herbGenTree <- read.tree("/correctedherbGentree.tre")

#RUn Medusa

fie <- c(herbGenTree, herbGenTree)
res1 <- MEDUSA(fie,
               mc = TRUE,
               model = "bd",
               criterion = "aicc")
summ <-
  multiMedusaSummary(res1,
                     herbGenTree,
                     plotModelSizes = TRUE,
                     plotTree = FALSE)


plotMultiMedusa(
  summ,
  type = "fan",
  time = FALSE,
  annotateShift = TRUE,
  annotateRate = "r.mean",
  show.tip.label = FALSE,
  shift.leg.pos =  "bottomleft"
)



# Get the time of each node in the tree
node_times <- node.depth.edgelength(herbGenTree)

upshifts <-
  as.data.frame(summ$shift.summary[summ$shift.summary[, 4] > 0, ])


# Map each shift node to its corresponding time
shift_times <- 380.601 - node_times[upshifts[, 1]]

# Define the different window sizes and step sizes you want to use
window_sizes <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
step_sizes <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

# Number of bootstrap replicates
n_boot <- 100

# Initialize a list to hold the data frames for each window size
dfs <- list()

# Iterate over the different window sizes
for (j in seq_along(window_sizes)) {
  window_size <- window_sizes[j]
  step_size <- step_sizes[j]
  start_times <- seq(0, max(shift_times), by = step_size)
  
  # Initialize a matrix to hold the counts of rate shifts for each bootstrap replicate
  rate_shift_counts <-
    matrix(nrow = n_boot, ncol = length(start_times))
  
  # For each bootstrap replicate
  for (b in 1:n_boot) {
    # Resample the shift times with replacement
    resampled_shift_times <- sample(shift_times, replace = TRUE)
    
    # For each sliding window, count the number of resampled shift times that fall within the window
    for (i in seq_along(start_times)) {
      start_time <- start_times[i]
      end_time <- start_time + window_size
      rate_shift_counts[b, i] <-
        sum(resampled_shift_times >= start_time &
              resampled_shift_times < end_time)
    }
  }
  
  # Convert the matrix to a data frame in long format
  df <- as.data.frame(rate_shift_counts)
  names(df) <- start_times
  df$bootstrap_replicate <- 1:n_boot
  df$window_size <- window_size
  df_long <-
    gather(df,
           start_time,
           rate_shift_count,-bootstrap_replicate,-window_size)
  
  # Convert the start_time column to numeric
  df_long$start_time <- as.numeric(as.character(df_long$start_time))
  
  # Add the data frame to the list
  dfs[[j]] <- df_long
}

# Combine all the data frames into one
df_combined <- bind_rows(dfs)


#Now we need to summarize the data and get the percent of total shifts but we need to group by window size and count the total number of shifts

# Group by bootstrap replicate and calculate the total number of shifts for each replicate
df_summ <- df_combined %>%
  group_by(bootstrap_replicate, window_size) %>%
  mutate(total_shifts = sum(rate_shift_count)) %>%
  mutate(percentage_of_total_shifts = rate_shift_count / total_shifts * 100)# Calculate the percentage of total shifts for each time window within each bootstrap replicate


# Calculate the mean percentage of total shifts for each time window across the bootstrap replicates
df_summary <- df_summ %>%
  group_by(start_time) %>%
  summarise(
    mean_percentage_of_total_shifts = mean(percentage_of_total_shifts),
    .groups = "drop"
  )

# Plot the data in points
ggplot(df_summary,
       aes(x = start_time, y = mean_percentage_of_total_shifts)) +
  geom_point(color = "goldenrod4") +
  labs(x = "Time", y = "Mean total shifts (%)") +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "red",
    size = 1.5
  ) +
  theme_minimal() + scale_x_reverse()


#plot barplots
ggplot(df_summary,
       aes(x = start_time, y = mean_percentage_of_total_shifts)) +
  geom_bar(stat = "identity",
           fill = "steelblue",
           alpha = 0.5) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "red",
    size = 1.5
  ) +
  labs(x = "Time", y = "Mean total shifts (%)") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 250, by = 50), trans = "reverse")


ltt_dfs <- lapply(seq_along(fern_nectaries.ltts), function(i) {
  ltt_data <- fern_nectaries.ltts[[i]]
  time <- max(ltt_data$times) - ltt_data$times
  lineages <- ltt_data$ltt[, 2]
  data.frame(run = i,
             time = time,
             lineages = lineages)
})

# Combine all data frames into one
ltt_df_combined <- do.call(rbind, ltt_dfs)


# Calculate max values for both datasets
max_rate_shifts <-
  max(df_summary$mean_percentage_of_total_shifts, na.rm = TRUE)
max_lineages <- max(ltt_df_combined$lineages, na.rm = TRUE)

# Create a scaling factor
scaling_factor <- max_rate_shifts / max_lineages

# Scale the lineage data
ltt_df_combined$lineages_scaled <-
  ltt_df_combined$lineages * scaling_factor

# Plot the data
ggplot(df_summary,
       aes(x = start_time, y = mean_percentage_of_total_shifts)) +
  geom_bar(stat = "identity",
           fill = "gray",
           alpha = 0.5) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "black",
    size = 1.5
  ) + geom_line(
    data = ltt_df_combined,
    aes(x = time, y = lineages_scaled, group = factor(run)),
    color = "#26cf0e99"
  ) +
  labs(x = "Time", y = "Percentage of Total Shifts") +
  scale_y_continuous(sec.axis = sec_axis(~ . / scaling_factor, name = "Number of Lineages")) +
  scale_x_reverse(limits = c(155, -3)) + theme_classic()







# HiSSE/BiSSE analyses

fern_5k_tree <- read.tree("/fern_5k_tree.tre")


#bring in species nectary data
species_nectar.dat <-
  read.csv("/Fern_Nectary_Database_Nov19_2022.csv")

# filter tips in the tree to just match tips in the dataframe
rownames(species_nectar.dat) <- species_nectar.dat$Species

#drop all tips that are not No or Yes
ss <- species_nectar.dat %>%
  filter(Nectary == "Yes" | Nectary ==  "No") %>%
  select(Nectary)


# check species in tree
chk <-
  name.check(fern_5k_tree, ss)# this is a crazy awesome and useful function

# now prune the tre based on the data available

tree.pruned <- drop.tip(fern_5k_tree, chk$tree_not_data)

tree.pruned <- multi2di(tree.pruned)

#Make ultrametric
tree.pruned_ult <- force.ultrametric(tree.pruned)

fern.Data.pruned <-
  ss[!(rownames(ss) %in% chk$data_not_tree), , drop = FALSE]

df <- fern.Data.pruned %>%
  rownames_to_column(., var = "Species") %>%
  mutate(Nectary = ifelse(Nectary == "Yes", 1, 0))


nectary.mode <- setNames(as.numeric(df$Nectary), df$Species)


total_tree <- length(tree.pruned_ult$tip.label)

# total species in group
# based on estimate
total_sp <- 11000

#global:
samp1 <- c(total_tree / total_sp)

##Do species with nectaries diversify faster than those without?

# convert data to hisse format
taxa2 <- as.matrix(df)

## BiSSE like HiSSE model
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits = 0)
bisse_like_hisse <-
  hisse(
    tree.pruned_ult,
    taxa2,
    hidden.states = FALSE,
    f = c(0.4, 0.8),
    turnover = c(1, 2),
    eps = c(1, 2),
    trans.rate = trans.rates.bisse
  )


## null BiSSE model
null_bisse_like_hisse <-
  hisse(
    tree.pruned_ult,
    taxa2,
    hidden.states = FALSE,
    f = c(0.4, 0.8),
    turnover = c(1, 1),
    eps = c(1, 1),
    trans.rate = trans.rates.bisse
  )


## CID-2 HiSSE model
trans.rates.hisse.1 <-
  TransMatMakerHiSSE(hidden.traits = 1, make.null = TRUE)
null_two_hisse <-
  hisse(
    tree.pruned_ult,
    taxa2,
    f = c(0.4, 0.8),
    hidden.states = TRUE,
    turnover = c(1, 1, 2, 2),
    eps = c(1, 1, 2, 2),
    trans.rate = trans.rates.hisse.1
  )

## CID-4 HiSSE model
trans.rates.hisse.2 <-
  TransMatMakerHiSSE(hidden.traits = 3, make.null = TRUE)
null_four_hisse <-
  hisse(
    tree.pruned_ult,
    taxa2,
    f = c(0.4, 0.8),
    hidden.states = TRUE,
    turnover = c(1, 1, 2, 2, 3, 3, 4, 4),
    eps = c(1, 1, 2, 2, 3, 3, 4, 4),
    trans.rate = trans.rates.hisse.2
  )

## full HiSSE model
trans.rates.hisse.3 <- TransMatMakerHiSSE(hidden.traits = 1)
full_hisse <-
  hisse(
    tree.pruned_ult,
    taxa2,
    f = c(0.4, 0.8),
    hidden.states = TRUE,
    turnover = c(1, 2, 3, 4),
    eps = c(1, 2, 3, 4),
    trans.rate = trans.rates.hisse.3
  )

## HiSSE with one hidden state only
hisse_1hidden <-
  hisse(
    tree.pruned_ult,
    taxa2,
    f = c(0.4, 0.8),
    hidden.states = TRUE,
    turnover = c(1, 2, 1, 2),
    eps = c(1, 2, 1, 2),
    trans.rate = trans.rates.hisse.3
  )

hisse_1hidden2 <-
  hisse(
    tree.pruned_ult,
    taxa2,
    f = c(0.4, 0.8),
    hidden.states = TRUE,
    turnover = c(1, 2, 0, 3),
    eps = c(1, 2, 0, 3),
    trans.rate = trans.rates.hisse.3
  )
