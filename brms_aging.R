library(readr)
library(caper)
library(utils)
library(ape)
library(dplyr)
library(here)
library(brms)
library(mice)
library(cmdstanr)
options(scipen = 999) 

#age data----
agedat_filtered <- read.csv(here("DDL", "AnAge", "animal_longevity", "data", "raw", "agedat_filtered.csv"))

#mammal diet data, Price et al. 2012----
price2012_diet <- read.delim(here("DDL", "AnAge", "animal_longevity", "data", "raw", "price_mammaldiet.txt")) #1515 species
colnames(price2012_diet)[1] <- "species"

#combining pantheria mass data in g, averaged longevity data, diet
mammal_agediet <- merge(agedat_filtered[c(1,9,25)], price2012_diet, by="species", all.x=F, all.y=F)
mammal_agediet$logage <- log(mammal_agediet$mean_long)
mammal_agediet$logmass <- log(mammal_agediet$panth_bodymass_g)
plot(mammal_agediet$logmass~mammal_agediet$logage) 

#change ?s to NAs and diet variable to factors
class(mammal_agediet$Diet)
mammal_agediet$diet.new <- replace(mammal_agediet[,4], mammal_agediet[,4]=="?", NA)
mammal_agediet$diet.new <- as.factor(mammal_agediet$diet.new)
class(mammal_agediet$diet.new)

#tree----
mammaltree <- read.tree(file=here("DDL/AnAge/animal_longevity/data/processed/agemammals.tre"))
diet.comp <- comparative.data(mammaltree[[77]], mammal_agediet, names.col="species")

cormat <- vcv(diet.comp$phy, corr=TRUE) #correlation matrix
#https://discourse.mc-stan.org/t/covariance-matrix-phylogenetic-models/20477/2
diet.comp[["data"]]$species <- rownames(diet.comp[["data"]])

View(diet.comp$data)
mammal_dietcomp <- diet.comp$data

#mammals diet----
m.mammal_agediet <- brms::brm(bf(logage ~ logmass + diet.new + (1 | gr(species, cov = cormat))),
            data = mammal_dietcomp,
            data2 = list(cormat=cormat),
            iter = 4000, chains = 4,
            backend = "cmdstanr",
            cores = 4)
summary(m.mammal_agediet)
saveRDS(m.mammal_agediet, here("DDL", "AnAge", "animal_longevity", "models", "m.mammal_agediet.rds"))

##herbivory vs. others----

mammal_dietcomp$herb <- ifelse(mammal_dietcomp$diet.new == "Herbivore", 1, 0)
mammal_dietcomp$herb <- as.factor(mammal_dietcomp$herb)
class(mammal_dietcomp$herb)

m.mammal_ageherb <- brm(bf(logage ~ logmass + herb + (1 | gr(species, cov = cormat))),
            data = mammal_dietcomp,
            data2 = list(cormat=cormat),
            iter = 4000, chains = 4,
            backend = "cmdstanr",
            cores = 4)
summary(m.mammal_ageherb)
saveRDS(m.mammal_ageherb, here("DDL", "AnAge", "animal_longevity", "models", "m.mammal_ageherb.rds"))



#diel data, Maor et al. 2017----
diel <- agedat_filtered[c(1,4,5,6,9,18,20,22,25,26)]
diel <- diel[!is.na(diel$diel),]
class(diel$diel)
diel$diel <- as.factor(diel$diel)
diel$diel #levels are ARR CRE DIU NOC
table(diel$diel)
table(diel$class) 

diel$logage <- log(diel$mean_long)
saveRDS(diel, file=here("DDL/AnAge/diel.rds"))

diel <- readRDS(file=here("DDL/AnAge/diel.rds"))
write_csv(diel, here("DDL", "AnAge", "animal_longevity", "data", "raw", "diel.csv"))

##diel amphibians
diel_amphi <- diel[diel$class == "Amphibia",]
table(diel_amphi$diel) #only 14 not nocturnal

##diel mammals----
diel_mam <- diel[diel$class == "Mammalia",]
table(diel_mam$diel)
diel_mam$noct <- ifelse(diel_mam$diel == "NOC", 1, 0)
diel_mam$diu <- ifelse(diel_mam$diel == "DIU", 1, 0)

diel_mam$logmass <- log(diel_mam$panth_bodymass_g)

mammaltree <- read.tree(here("DDL", "AnAge", "animal_longevity", "data", "processed", "agemammals.tre"))
sum(mammaltree[[50]]$tip.label %in% diel_mam$species) #155

dieL.comp <- comparative.data(mammaltree[[77]], diel_mam, names.col="species", na.omit=FALSE, vcv=TRUE)

dieL.comp$data$species <- rownames(dieL.comp$data)

View(dieL.comp$data)

saveRDS(dieL.comp, here("DDL", "AnAge", "animal_longevity", "data", "processed", "dielmam.comp.rds"))

### mammal diel data for brms----
dielmam.comp <- readRDS(here("DDL", "AnAge", "animal_longevity", "data", "processed", "dielmam.comp.rds"))
cormat <- vcv(dielmam.comp$phy, corr=TRUE) #correlation matrix
View(dielmam.comp$data)


### mammal nocturnal ----
m.mammal_noct <-
  brm(bf(logage ~ logmass + as.factor(noct) + (1 | gr(species, cov = cormat))),
    data = dielmam.comp$data,
    data2 = list(cormat=cormat),
    iter = 4000, chains = 4,
    backend = "cmdstanr",
    cores = 4)
#no effect
summary(m.mammal_noct)
saveRDS(m.mammal_noct, here("DDL", "AnAge", "animal_longevity", "models", "m.mammal_noct.rds"))

### mammal diurnal-----
m.mammal_diu <-
  brm(bf(logage ~ logmass + as.factor(diu) + (1 | gr(species, cov = cormat))),
      data = dielmam.comp$data,
      data2 = list(cormat=cormat),
      iter = 4000, chains = 4,
      backend = "cmdstanr",
      cores = 4)
#no effect
summary(m.mammal_diu)
saveRDS(m.mammal_diu, here("DDL", "AnAge", "animal_longevity", "models", "m.mammal_diu.rds"))


mean(dielmam.comp$data[dielmam.comp$data$diel == "DIU",]$mean_long)
mean(dielmam.comp$data[dielmam.comp$data$diel == "NOC",]$mean_long)

##reptiles----
diel_reptiles <- diel[diel$class == "Reptilia",]
table(diel_reptiles$diel) 
mean(diel_reptiles[diel_reptiles$diel == "DIU",]$mean_long)
mean(diel_reptiles[diel_reptiles$diel == "NOC",]$mean_long)
diel_reptiles$logmass <- log(diel_reptiles$kup_AdultMass)

rep.tree <- read.nexus(here("DDL", "AnAge", "animal_longevity", "data", "processed", "agesquamates.trees"))

dieL_reptiles <- comparative.data(rep.tree, diel_reptiles, names.col="species", na.omit=FALSE, vcv=TRUE)
dieL_reptiles$phy
dieL_reptiles$data$species <- rownames(dieL_reptiles$data)


dieL_reptiles$data$noct <- ifelse(dieL_reptiles$data$diel == "NOC", 1, 0)
dieL_reptiles$data$diu <- ifelse(dieL_reptiles$data$diel == "DIU", 1, 0)

saveRDS(dieL_reptiles, here("DDL", "AnAge", "animal_longevity", "data", "processed", "dielsquam.comp.rds"))

### reptiles comparative data for brms----
dielsquam.comp <- readRDS(here("DDL", "AnAge", "animal_longevity", "data", "processed", "dielsquam.comp.rds"))
cormat <- vcv(dielsquam.comp$phy, corr=TRUE) #correlation matrix

m.reptile_noct <-
  brm(bf(logage ~ logmass + as.factor(noct) + (1 | gr(species, cov = cormat))),
      data = dielsquam.comp$data,
      data2 = list(cormat=cormat),
      iter = 4000, chains = 4,
      backend = "cmdstanr",
      cores = 4)
summary(m.reptile_noct)
saveRDS(m.reptile_noct, here("DDL", "AnAge", "animal_longevity", "models", "m.reptile_noct.rds"))

m.reptile_diu <-
  brm(bf(logage ~ logmass + as.factor(diu) + (1 | gr(species, cov = cormat))),
      data = dielsquam.comp$data,
      data2 = list(cormat=cormat),
      iter = 4000, chains = 4,
      backend = "cmdstanr",
      cores = 4)
summary(m.reptile_diu)
saveRDS(m.reptile_diu, here("DDL", "AnAge", "animal_longevity", "models", "m.reptile_diu.rds"))


##birds----
diel_birds <- diel[diel$class == "Aves",]
table(diel_birds$diel)

table(diel$class)
mean(diel_birds$mean_long)
mean(diel_birds[diel_birds$diel == "DIU",]$mean_long)
mean(diel_birds[diel_birds$diel == "NOC",]$mean_long)



#metabolic rate----
#aa_metrate_w
#AT_metrate

metrate <- agedat_filtered[c(1,4,9,17,18,20,22,23,25)]
checkcols <- c("aa_metrate_w","AT_metrate")

metrate <- metrate %>% filter(!if_all(all_of(checkcols), is.na))

metrate$bmr_avg <- rowMeans(metrate[c(4,8)], na.rm=TRUE)
metrate$logage <- log(metrate$mean_long)
metrate$mass_avg <- rowMeans(metrate[c(3,5,6,7)], na.rm=TRUE)
metrate$logmass <- log(metrate$mass_avg)
metrate$logmr <- log(metrate$bmr_avg)

saveRDS(metrate, file=here("DDL/AnAge/animal_longevity/data/metrate.rda"))

#probably enough mammals and birds, not enough amphs or reptiles

mammaltree <- read.tree(file=here("DDL/AnAge/animal_longevity/data/processed/agemammals.tre"))
sum(mammaltree[[77]]$tip.label %in% metrate$species)
sum(metrate$species %in% mammaltree[[77]]$tip.label)

metrate_mam.comp <- comparative.data(mammaltree[[77]], metrate, names.col="species", na.omit=FALSE)
cormat <- vcv(metrate_mam.comp$phy, corr=TRUE) #correlation matrix
metrate_mam.comp[["data"]]$species <- rownames(metrate_mam.comp[["data"]])

#mammals
m.mammal_metrate <-
  brm(bf(logage ~ logmr + logmass + (1 | gr(species, cov = cormat))),
      data = metrate_mam.comp$data,
      data2 = list(cormat=cormat),
      iter = 4000, chains = 4,
      backend = "cmdstanr",
      cores = 4)
summary(m.mammal_metrate)
saveRDS(m.mammal_metrate, file=here("DDL", "AnAge", "animal_longevity", "models", "m.mammal_metrate.rds"))
#slightly negative but includes 0

#birds
birdtrees <- read.nexus(here("DDL/AnAge/animal_longevity/data/processed/agebirds.trees"))
sum(metrate$species %in% birdtrees$tip.label) #172
mr_bird.comp <- comparative.data(birdtrees, metrate, names.col="species", na.omit=FALSE)
cormat <- vcv(mr_bird.comp$phy, corr=TRUE) #correlation matrix
mr_bird.comp[["data"]]$species <- rownames(mr_bird.comp[["data"]])
bird_mr_dat <- mr_bird.comp$data 

m.bird_metrate <-
  brm(bf(logage ~ logmr + logmass + (1 | gr(species, cov = cormat))),
      data = bird_mr_dat,
      data2 = list(cormat=cormat),
      iter = 4000, chains = 4,
      backend = "cmdstanr",
      cores = 4)
summary(m.bird_metrate)
saveRDS(m.bird_metrate, file=here("DDL", "AnAge", "animal_longevity", "models", "m.bird_metrate.rds"))

#quite similar to mammals, slightly negative but includes 0

