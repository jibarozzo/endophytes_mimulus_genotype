### Distance based dereundancy Model (dbRDA) for the MIM3 project ###
###2024-04-04###
###Bolívar Aponte Rolón###

########################################
### Package installation and loading ###
########################################
if (!requireNamespace(c("vegan", "doParallel"), quietly = TRUE))
    install.packages(c("vegan", "doParallel")
    )

########################################
########### Parallel jobs ##############
########################################

# Parallelization through Rmpi
# install.packages(
#     "Rmpi", 
#     configure.args = c(
#         "--with-Rmpi-include=/usr/include/", # This is where OPENMPI's mpi.h is located
#         "--with-Rmpi-libpath=/usr/lib64/",     # This is where liblam.so is located
#         "--with-Rmpi-type=OPENMPI"               # This says that the type is OPENMPI
#     ))
# library(Rmpi)

# Initialize MPI (Rmpi)
#mpi.spawn.Rslaves(nslaves=8)
#mpi.remote.exec(paste("library(vegan)")) # Loading vegan in the remote slave


# Parallelization through doParallel
library("doParallel")

# use the environment variable SLURM_CPUS_PER_TASK to set the number of cores
registerDoParallel(cores=(Sys.getenv("SLURM_CPUS_PER_TASK")))
getDoParWorkers()

########################################
########### Load libraries #############
########################################
library("vegan")

# Session info
sessionInfo()

################################################################
#### dbRDA modelling for all samples: parentals and hybrids ####
################################################################
rrfy_hell_matrix <- readRDS("clean_data/statistics/rrfy_hell_matrix.rds")

# Traits for dbRDA
##### Hellinger transformed rarefied data
# All samples
# dbrda_hell_matrix <- rrfy_hell_matrix |>
#     column_to_rownames(var = "X") |>
#     select(-c(1, 2:17))

#saveRDS(dbrda_hell_matrix, file = "clean_data/statistics/dbrda_hell_matrix.rds")
dbrda_hell_matrix <- readRDS("clean_data/statistics/dbrda_hell_matrix.rds")


# dbrda_traits <- rrfy_hell_matrix |>
#     select(X, Unique_ID, Site, Habitat, Genotype, logLBI) |>
#     unite(Habitat_Genotype,
#           Habitat,
#           Genotype,
#           sep = "_",
#           remove = FALSE) # To model interactions
#saveRDS(dbrda_traits, file = "clean_data/statistics/dbrda_traits.rds")
dbrda_traits <- readRDS("clean_data/statistics/dbrda_traits.rds")

# Model with intercept only ####
# Using rrfy_hell_matrix as the distance matrix
# m0_hell <- dbrda(
#     dbrda_hell_matrix ~ 1,
#     distance = "bray",
#     dfun = vegdist,
#     data = dbrda_traits,
#     parallel = 20,
#     #Passes parallelization to metaMDS function # Do not use when using Rmpi
#     na.action = na.omit
# ) #Model with intercept only.
# saveRDS(m0_hell, file = "clean_data/statistics/m0_hell.rds")
# 
# m1_hell <- dbrda(
#     dbrda_hell_matrix ~ .,
#     distance = "bray",
#     dfun = vegdist,
#     data = dbrda_traits,
#     parallel = 20,
#     na.action = na.omit
# ) # Model with all explanatory variables.
# saveRDS(m1_hell, file = "clean_data/statistics/m1_hell.rds")

# Model with species, site, leaf traits and elevation. ####
m2_hell <-
    dbrda(
        dbrda_hell_matrix ~ logLBI + Habitat + Genotype,
        distance = "bray",
        dfun = vegdist,
        data = dbrda_traits,
        parallel = 20,
        na.action = na.omit
    )
saveRDS(m2_hell, file = "clean_data/statistics/m2_hell.rds")

#Anovas for m2
#Margins
set.seed(123)
anova_m2_hell_margin <- anova.cca(m2_hell,
                                  by = "margin",
                                  permutations = 999,
                                  parallel = 10)
saveRDS(anova_m2_hell_margin, file = "clean_data/statistics/anova_m2_hell_margin.rds")

#Terms
set.seed(123)
anova_m2_hell_terms <- anova.cca(m2_hell,
                                 by = "terms",
                                 permutations = 999,
                                 parallel = 10)
saveRDS(anova_m2_hell_terms, file = "clean_data/statistics/anova_m2_hell_terms.rds")

# Model with interactions between site/habitat and genotype
m3_hell <-
    dbrda(
        dbrda_hell_matrix ~ logLBI + Habitat * Genotype,
        distance = "bray",
        dfun = vegdist,
        data = dbrda_traits,
        parallel = 20,
        na.action = na.omit
    )
saveRDS(m3_hell, file = "clean_data/statistics/m3_hell.rds")

#Anovas for m3
#Margins
set.seed(123)
anova_m3_hell_margin <- anova.cca(m3_hell,
                                  by = "margin",
                                  permutations = 999,
                                  parallel = 20)
saveRDS(anova_m3_hell_margin, file = "clean_data/statistics/anova_m3_hell_margin.rds")

#Terms
set.seed(123)
anova_m3_hell_term <- anova.cca(m3_hell,
                                by = "terms",
                                permutations = 999,
                                parallel = 20)
saveRDS(anova_m3_hell_term, file = "clean_data/statistics/anova_m3_hell_term.rds")

#################################
#### PERMDISP ####
#################################
#Analysis of multivariate homogeneity of group dispersions.

##################################
#### m2_hell model #####
##################################

# By Genotype
set.seed(123)
beta_dis1 <-
    betadisper(
        vegdist(dbrda_hell_matrix, method = "bray"),
        # Matrix of Hellinger transformed data rarefied
        dbrda_traits$Genotype,
        type = "median",
        sqrt.dist = FALSE
    )
save(beta_dis1, file = "clean_data/statistics/beta_dis1.rda")
load("clean_data/statistics/beta_dis1.rda")

## Permutest
## Margin
set.seed(123)
beta_perm1_margin <-
    permutest(beta_dis1,
              parallel = 8,
              permutations = 999,
              by = "margin")
saveRDS(beta_perm1_margin, file = "clean_data/statistics/beta_perm1_margin.rds")
beta_perm1_margin <- readRDS("clean_data/statistics/beta_perm1_margin.rds")

anova(beta_dis1)
TukeyHSD(beta_dis1)

## Terms
set.seed(123)
beta_perm1_terms <-
    permutest(beta_dis1,
              parallel = 8,
              permutations = 999,
              by = "terms")
saveRDS(beta_perm1_terms, file = "clean_data/statistics/beta_perm1_terms.rds")
beta_perm1_terms <- readRDS("clean_data/statistics/beta_perm1_terms.rds")

# By Site
set.seed(123)
beta_dis2 <- betadisper(
    vegdist(dbrda_hell_matrix, method = "bray"),
    dbrda_traits$Site,
    type = "median",
    sqrt.dist = FALSE
)
saveRDS(beta_dis2, file = "clean_data/statistics/beta_dis2.rds")

## Permutest
## Margin
set.seed(123)
beta_perm2_margin <-
    permutest(beta_dis2,
              parallel = 8,
              permutations = 999,
              by = "margin")
saveRDS(beta_perm2_margin, file = "clean_data/statistics/beta_perm2_margin.rds")
beta_perm2_margin <- readRDS("clean_data/statistics/beta_perm2_margin.rds")

## Terms
set.seed(123)
beta_perm2_terms <-
    permutest(beta_dis2,
              parallel = 8,
              permutations = 999,
              by = "terms")
saveRDS(beta_perm2_terms, file = "clean_data/statistics/beta_perm2_terms.rds")
beta_perm2_terms <- readRDS("clean_data/statistics/beta_perm2_terms.rds")

##################################
#### m3_hell model #####
##################################

# By Genotype
set.seed(123)
beta_dis3 <-
    betadisper(
        vegdist(dbrda_hell_matrix, method = "bray"),
        # Matrix of Hellinger transformed data rarefied
        dbrda_traits$Genotype,
        type = "median",
        sqrt.dist = FALSE
    )
saveRDS(beta_dis3, file = "clean_data/statistics/beta_dis3.rds")

anova(beta_dis3)
TukeyHSD(beta_dis3)

## Permutest
## Margin
set.seed(123)
beta_perm3_margin <-
    permutest(beta_dis3,
              parallel = 10,
              permutations = 999,
              by = "margin")
saveRDS(beta_perm3_margin, file = "clean_data/statistics/beta_perm3_margin.rds")

## Terms
set.seed(123)
beta_perm3_terms <-
    permutest(beta_dis3,
              parallel = 10,
              permutations = 999,
              by = "terms")
saveRDS(beta_perm3_margin, file = "clean_data/statistics/beta_perm3_terms.rds")


# By Habitat
set.seed(123)
beta_dis4 <- betadisper(
    vegdist(dbrda_hell_matrix, method = "bray"),
    dbrda_traits$Habitat,
    type = "median",
    sqrt.dist = FALSE
)
saveRDS(beta_dis4, file = "clean_data/statistics/beta_dis4.rds")

## Permutest
set.seed(123)
beta_perm4_margin <-
    permutest(beta_dis4,
              parallel = 10,
              permutations = 999,
              by = "margin")
saveRDS(beta_perm4_margin, file = "clean_data/statistics/beta_perm4_margin.rds")

## Terms
set.seed(123)
beta_perm4_margin <-
    permutest(beta_dis4,
              parallel = 10,
              permutations = 999,
              by = "terms")
saveRDS(beta_perm4_terms, file = "clean_data/statistics/beta_perm4_terms.rds")

# By Habitat*Genotype
set.seed(123)
beta_dis5 <- betadisper(
    vegdist(dbrda_hell_matrix, method = "bray"),
    dbrda_traits$Habitat_Genotype,
    type = "median",
    sqrt.dist = FALSE
)
saveRDS(beta_dis5, file = "clean_data/statistics/beta_dis5.rds")

## Permutest
## Margin
set.seed(123)
beta_perm5_margin <-
    permutest(beta_dis5,
              parallel = 10,
              permutations = 999,
              by = "margin")
saveRDS(beta_perm5_margin, file = "clean_data/statistics/beta_perm5_margin.rds")
beta_perm5

## Terms
set.seed(123)
beta_perm5_terms <-
    permutest(beta_dis5,
              parallel = 10,
              permutations = 999,
              by = "terms")
saveRDS(beta_perm5_terms, file = "clean_data/statistics/beta_perm5_terms.rds")


#################################################################
#### dbRDA modelling for all samples: hybrids ONLY ####
#################################################################

##### Hellinger transformed rarefied data
# dbrda_hybrids_matrix <- rrfy_hell_matrix |>
#     column_to_rownames(var = "X") |>
#     select(-c(1, 2:17)) |>
#     t() |>
#     as.data.frame() |>
#     select(contains(c("F2WY", "F2YW"))) |>
#     t() |>
#     as.data.frame()
# 
# saveRDS(dbrda_hybrids_matrix, file = "clean_data/statistics/dbrda_hybrids_matrix.rds")
dbrda_hybrids_matrix <- readRDS(file.path(path,"clean_data/statistics/dbrda_hybrids_matrix.rds"))


# dbrda_hybrids_traits <- rrfy_hell_matrix |>
#     select(X, Unique_ID, Site, Habitat, Genotype, logLBI) |>
#     filter(X %in% rownames(dbrda_hybrids_matrix)) |>
#     unite(Habitat_Genotype,
#           Habitat,
#           Genotype,
#           sep = "_",
#           remove = FALSE) # To model interactions
# saveRDS(dbrda_hybrids_traits, file = "clean_data/statistics/dbrda_hybrids_traits.rds")
dbrda_hybrids_traits <- readRDS(file.path(path, "clean_data/statistics/dbrda_hybrids_traits.rds"))


# Model with intercept only ####
# Using rrfy_hell_matrix as the distance matrix
# m0_hybrids_hell <- dbrda(
#     dbrda_hybrids_matrix ~ 1,
#     distance = "bray",
#     dfun = vegdist,
#     data = dbrda_hybrids_traits,
#     parallel = 20, #Passes parallelization to metaMDS function
#     na.action = na.omit
# ) #Model with intercept only.
# saveRDS(m0_hybrids_hell, file = "clean_data/statistics/m0_hybrids_hell.rds")
# 
# m1_hybrids_hell <- dbrda(
#     dbrda_hybrids_matrix ~ .,
#     distance = "bray",
#     dfun = vegdist,
#     data = dbrda_hybrids_traits,
#     parallel = 20,
#     na.action = na.omit
# ) # Model with all explanatory variables.
# saveRDS(m1_hybrids_hell, file = "clean_data/statistics/m1_hybrids_hell.rds")

# Model with species, site, leaf traits and elevation. ####
m2_hybrids_hell <-
    dbrda(
        dbrda_hybrids_matrix ~ logLBI + Habitat + Genotype,
        distance = "bray",
        dfun = vegdist,
        data = dbrda_hybrids_traits,
        parallel = 20,
        na.action = na.omit
    )
save(m2_hybrids_hell, file = "clean_data/statistics/m2_hybrids_hell.rda")

#Anovas for m2
## Margin
set.seed(123)
anova_m2_hybrids_margin <- anova.cca(
    m2_hybrids_hell,
    by = "margin",
    permutations = 999,
    parallel = 10
)
saveRDS(anova_m2_hybrids_hell_margin, file = "clean_data/statistics/anova_m2_hybrids_margin.rds")

## Terms
set.seed(123)
anova_m2_hybrids_terms <- anova.cca(
    m2_hybrids_hell,
    by = "terms",
    permutations = 999,
    parallel = 10
)
saveRDS(anova_m2_hybrids_terms, file = "clean_data/statistics/anova_m2_hybrids_terms.rds")


# m3 hybrids
m3_hybrids_hell <-
    dbrda(
        dbrda_bybrids_matrix ~ logLBI + Habitat * Genotype,
        distance = "bray",
        dfun = vegdist,
        data = dbrda_traits,
        parallel = 20,
        na.action = na.omit
    )
saveRDS(m3_hell, file = "clean_data/statistics/m3_hybrids_hell.rds")

#Anovas for m3_hybrids
#Margins
set.seed(123)
anova_m3_hybrids_margin <- anova.cca(m3_hell,
                                     by = "margin",
                                     permutations = 999,
                                     parallel = 10)
saveRDS(anova_m3_hybrids_margin, file = "clean_data/statistics/anova_m3_hybrids_margin.rds")

#Terms
set.seed(123)
anova_m3_hybrids_term <- anova.cca(m3_hell,
                                   by = "terms",
                                   permutations = 999,
                                   parallel = 10)
saveRDS(anova_m3_hybrids_term, file = "clean_data/statistics/anova_m3_hybrids_term.rds")



### PERMDISP

#Analysis of multivariate homogeneity of group dispersions.

# By Genotype
set.seed(123)
beta_hybrids_dis1 <-
    betadisper(
        vegdist(dbrda_hybrids_matrix, method = "bray"), # Matrix of Hellinger transformed data rarefied
        dbrda_hybrids_traits$Genotype,
        type = "median",
        sqrt.dist = FALSE
    )
saveRDS(beta_hybrids_dis1, file = "clean_data/statistics/beta_hybrids_dis1.rds")

## Permutest
## Margin
set.seed(123)
beta_hybrids_perm1_margin <-
    permutest(beta_hybrids_dis1,
              parallel = 10,
              permutations = 999,
              by = "margin")
saveRDS(beta_hybrids_perm1_margin, file = "clean_data/statistics/beta_hybrids_perm1_margin.rds")

anova(beta_hybrids_dis1_margin)
TukeyHSD(beta_hybrids_dis1_margin)

## Terms
set.seed(123)
beta_hybrids_perm1_terms <-
    permutest(beta_hybrids_dis1,
              parallel = 10,
              permutations = 999,
              by = "terms")
saveRDS(beta_hybrids_perm1_terms, file = "clean_data/statistics/beta_hybrids_perm1_terms.rds")

# By Habitat
set.seed(123)
beta_hybrids_dis2 <- betadisper(
    vegdist(dbrda_hybrids_matrix, method = "bray"),
    dbrda_traits$Habitat,
    type = "median",
    sqrt.dist = FALSE
)
saveRDS(beta_hybrids_dis2, file = "clean_data/statistics/beta_hybrids_dis2.rds")

## Permutest
## Margin
set.seed(123)
beta_hybrids_perm2_margin <-
    permutest(beta_dis2,
              parallel = 10,
              permutations = 999,
              by = "margin")
saveRDS(beta_hybrids_perm2_margin, file = "clean_data/statistics/beta_hybrids_perm2_margin.rds")

## Terms
set.seed(123)
beta_hybrids_perm2_terms <-
    permutest(beta_dis2,
              parallel = 10,
              permutations = 999,
              by = "terms")
saveRDS(beta_hybrids_perm2_terms, file = "clean_data/statistics/beta_hybrids_perm2_terms.rds")

# By Habitat*Genotype
set.seed(123)
beta_hybrids_dis3 <- betadisper(
    vegdist(dbrda_hybrids_matrix, method = "bray"),
    dbrda_traits$Habitat_Genotype,
    type = "median",
    sqrt.dist = FALSE
)
saveRDS(beta_hybrids_dis3, file = "clean_data/statistics/beta_hybrids_dis3.rds")

## Permutest
## Margin
set.seed(123)
beta_hybrids_perm3_margin <-
    permutest(beta_dis3,
              parallel = 10,
              permutations = 999,
              by = "margin")
saveRDS(beta_hybrids_perm3_margin, file = "clean_data/statistics/beta_hybrids_perm3_margin.rds")

## Terms
set.seed(123)
beta_hybrids_perm3_terms <-
    permutest(beta_dis3,
              parallel = 10,
              permutations = 999,
              by = "terms")
saveRDS(beta_hybrids_perm3_terms, file = "clean_data/statistics/beta_hybrids_perm3_terms.rds")

### Finalize MPI (Rmpi) ###
# mpi.close.Rslaves()
# mpi.quit()

# Stop the parallel backend
stopCluster(cores=(Sys.getenv("SLURM_CPUS_PER_TASK")))