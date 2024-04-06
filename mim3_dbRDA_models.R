### Distance based dereundancy Model (dbRDA) for the MIM3 project ###
###2024-04-04###
###Bolívar Aponte Rolón###

# Package installation and loading
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("vegan", "Rmpi")) #, "tidyverse"))


# Load libraries
#library(vegan)
library(Rmpi)

# Initialize MPI
mpi.spawn.Rslaves(nslaves=8)
mpi.remote.exec(paste("library(vegan)")) # Loading vegan in the remote slaves

############################################################################################################
#### dbRDA modelling for all samples: parentals and hybrids ####
############################################################################################################

### Load files ###
rrfy_hell_matrix <- readRDS("clean_data/statistics/rrfy_hell_matrix.rds")

##### Hellinger transformed rarefied data
# Traits for dbRDA
# dbrda_hell_matrix <- rrfy_hell_matrix |>
#     column_to_rownames(var = "X") |>
#     select(-c(1, 534:553))
dbrda_hell_matrix <- readRDS("clean_data/statistics/dbrda_hell_matrix.rds")

# dbrda_traits <- rrfy_hell_matrix |>
#     select(X, Unique_ID, Site, Genotype, logLBI)
dbrda_traits <- readRDS("clean_data/statistics/dbrda_traits.rds")

############################################################################################################


# Model with intercept only ####
# Using rrfy_hell_matrix as the distance matrix
m0_hell <- dbrda(
    dbrda_hell_matrix ~ 1,
    distance = "bray",
    dfun = vegdist,
    data = dbrda_traits,
    #parallel = 8, #Passes parallelization to metaMDS function # Do not use when using Rmpi
    na.action = na.omit
) #Model with intercept only.
save(m0_hell, file = "clean_data/statistics/m0_hell.rda")

m1_hell <- dbrda(
    dbrda_hell_matrix ~ .,
    distance = "bray",
    dfun = vegdist,
    data = dbrda_traits,
    #parallel = 8,
    na.action = na.omit
) # Model with all explanatory variables.
save(m1_hell, file = "clean_data/statistics/m1_hell.rda")


# Model with species, site, leaf traits and elevation. ####
m2_hell <-
    dbrda(
        dbrda_hell_matrix ~ logLBI + Site + Genotype,
        distance = "bray",
        dfun = vegdist,
        data = dbrda_traits,
        #parallel = 8,
        na.action = na.omit
    )
save(m2_hell, file = "clean_data/statistics/m2_hell.rda")

#Anovas for m2
set.seed(123)
anova_m2_hell <- anova.cca(m2_hell,
                           by = "margin",
                           permutations = 999)
                           #parallel = 8)
save(anova_m2_hell, file = "clean_data/statistics/m2_hell.rda")


### PERMDISP

#Analysis of multivariate homogeneity of group dispersions.

# By Genotype
set.seed(123)
beta_dis1 <-
    betadisper(
        vegdist(dbrda_hell_matrix, method = "bray"), # Matrix of Hellinger transformed data rarefied
        dbrda_traits$Genotype,
        type = "median",
        sqrt.dist = FALSE
    )
save(beta_dis1, file = "clean_data/statistics/beta_dis1.rda")

## Permutest
set.seed(123)
beta_perm1 <-
    permutest(beta_dis1,
              #parallel = 8,
              permutations = 999,
              by = "margin")
save(beta_perm1, file = "clean_data/statistics/beta_perm1.rda")


# By Site
set.seed(123)
beta_dis2 <- betadisper(
    vegdist(dbrda_hell_matrix, method = "bray"),
    dbrda_traits$Site,
    type = "median",
    sqrt.dist = FALSE
)
save(beta_dis2, file = "clean_data/statistics/beta_dis2.rda")

## Permutest
set.seed(123)
beta_perm2 <-
    permutest(beta_dis2,
              #parallel = 8,
              permutations = 999,
              by = "margin")
save(beta_perm2, file = "clean_data/statistics/beta_perm2.rda")

############################################################################################################
#### dbRDA modelling for all samples: hybrids ONLY ####
############################################################################################################

##### Hellinger transformed rarefied data
# dbrda_hybrids_matrix <- rrfy_hell_matrix |>
#     column_to_rownames(var = "X") |>
#     select(-c(1, 534:553)) |>
#     t() |>
#     as.data.frame() |>
#     select(contains(c("F2WY", "F2YW"))) |>
#     t() |>
#     as.data.frame()
dbrda_hybrids_matrix <-readRDS("clean_data/statistics/dbrda_hybrids_matrix.rds")


# dbrda_hybrids_traits <- rrfy_hell_matrix |>
#     select(X, Unique_ID, Site, Genotype, logLBI)|>
#     filter(X %in% rownames(dbrda_hybrids_matrix))
dbrda_hybrids_traits<-readRDS("clean_data/statistics/dbrda_hybrids_traits.rds")


# Model with intercept only ####
# Using rrfy_hell_matrix as the distance matrix
m0_hybrids_hell <- dbrda(
    dbrda_hybrids_matrix ~ 1,
    distance = "bray",
    dfun = vegdist,
    data = dbrda_hybrids_traits,
    #parallel = 8, #Passes parallelization to metaMDS function
    na.action = na.omit
) #Model with intercept only.
save(m0_hybrids_hell, file = "clean_data/statistics/m0_hybrids_hell.rda")

m1_hybrids_hell <- dbrda(
    dbrda_hybrids_matrix ~ .,
    distance = "bray",
    dfun = vegdist,
    data = dbrda_hybrids_traits,
    #parallel = 8,
    na.action = na.omit
) # Model with all explanatory variables.
save(m1_hybrids_hell, file = "clean_data/statistics/m1_hybrids_hell.rda")


# Model with species, site, leaf traits and elevation. ####
m2_hybrids_hell <-
    dbrda(
        dbrda_hybrids_matrix ~ logLBI + Site + Genotype,
        distance = "bray",
        dfun = vegdist,
        data = dbrda_hybrids_traits,
        #parallel = 8,
        na.action = na.omit
    )
save(m2_hybrids_hell, file = "clean_data/statistics/m2_hybrids_hell.rda")

#Anovas for m2
set.seed(123)
anova_m2_hybrids_hell <- anova.cca(
    m2_hybrids_hell,
    by = "margin",
    permutations = 999,
    #parallel = 8
)
save(anova_m2_hybrids_hell, file = "clean_data/statistics/anova_m2_hybrids_hell.rda")


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
save(beta_hybrids_dis1, file = "clean_data/statistics/beta_hybrids_dis1.rda")

## Permutest
set.seed(123)
beta_hybrids_perm1 <-
    permutest(beta_hybrids_dis1,
              #parallel = 8,
              permutations = 999,
              by = "margin")
save(beta_hybrids_perm1, file = "clean_data/statistics/beta_hybrids_perm1.rda")


# By Site
set.seed(123)
beta_hybrids_dis2 <- betadisper(
    vegdist(dbrda_hybrids__matrix, method = "bray"),
    dbrda_traits$Site,
    type = "median",
    sqrt.dist = FALSE
)
save(beta_hybrids_dis2, file = "clean_data/statistics/beta_hybrids_dis2.rda")

## Permutest
set.seed(123)
beta_hybrids_perm2 <-
    permutest(beta_dis2,
              #parallel = 8,
              permutations = 999,
              by = "margin")
save(beta_hybrids_perm2, file = "clean_data/statistics/beta_hybrids_perm2.rda")


### Finalize MPI ###
mpi.close.Rslaves()
mpi.quit()
