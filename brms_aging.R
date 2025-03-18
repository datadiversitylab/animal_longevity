library(megatrees)
library(rtrees)
library(readr)
library(utils)
library(ape)
library(dplyr)
library(here)
library(brms)
library(mice)
library(cmdstanr)
library(rethinking)
options(scipen = 999) 

#age data----
agedat_combined1 <- read.csv(here("data", "raw", "agedat_combined1.csv"))

#diet data, Price et al. 2012----
price2012_diet <- read.delim(here("data", "raw", "price_mammaldiet.txt")) #1515 species
sum(price2012_diet$Species.Name...Wilson...Reeder.2005 %in% agedat_combined1$species) #1218
colnames(price2012_diet)[1] <- "species"

#combining pantheria mass data in g, averaged longevity data, diet
mammal_agediet <- merge(agedat_combined1[c(8,24,26)], price2012_diet, by="species", all.x=F, all.y=F)
mammal_agediet$logage <- log(mammal_agediet$mean_long)
mammal_agediet$logmass <- log(mammal_agediet$panth_bodymass_g)
plot(mammal_agediet$logmass~mammal_agediet$logage) 

#change ?s to NAs and diet variable to factors
class(mammal_agediet$Diet)
mammal_agediet$diet.new <- replace(mammal_agediet[,4], mammal_agediet[,4]=="?", NA)
mammal_agediet$diet.new <- as.factor(mammal_agediet$diet.new)
class(mammal_agediet$diet.new)

#does this: https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r
mammal_agediet$mass.st <- rethinking::standardize(mammal_agediet$logmass)
mammal_agediet$age.st <- rethinking::standardize(mammal_agediet$logage)

#tree----
mammaltree <- read.tree(file="/Users/kiranbasava/DDL/AnAge/BayesTraitsV4/agemammals.tre")
diet.comp <- comparative.data(mammaltree[[77]], mammal_agediet, names.col="species")

cormat <- vcv(diet.comp$phy, corr=TRUE) #correlation matrix
#https://discourse.mc-stan.org/t/covariance-matrix-phylogenetic-models/20477/2
diet.comp[["data"]]$species <- rownames(diet.comp[["data"]])

View(diet.comp$data)

#pgls----
#dropping not imputing mass data
#priors are 'weakly regularizing' for logged and standardized vars (Rethinking, McElreath 2019)
m.mammal_agediet1 <- 
  brm(bf(logage ~ logmass + diet.new + (1 | gr(species, cov = cormat))),
      prior = c(prior(normal(0,1), class = Intercept),
                prior(normal(0,0.5), class = b)),
      data = diet.comp$data,
      data2 = list(cormat=cormat),
      iter = 2000, chains = 2,
      backend = "cmdstanr",
      cores = 2)
summary(m.mammal_agediet1)

#running with default priors (I think just flat) almost identical results
summary(brm(bf(logage ~ logmass + diet.new + (1 | gr(species, cov = cormat))),
    data = diet.comp$data,
    data2 = list(cormat=cormat),
    iter = 2000, chains = 2,
    backend = "cmdstanr",
    cores = 2))

#imputing mass data with mi() 
#this takes a while to run
m.mammal_agediet <- 
      brm(bf(logage ~ mi(logmass) + diet.new + (1 | gr(species, cov = cormat))) +
          bf(logmass | mi() ~ 1 + (1 | gr(species, cov = cormat))),
          prior = c(prior(normal(0,1), class = Intercept),
                    prior(normal(0,0.5), class = b)),
             data = diet.comp$data,
             data2 = list(cormat=cormat),
             iter = 2000, chains = 2,
             backend = "cmdstanr",
             cores = 2)
summary(m.mammal_agediet)
