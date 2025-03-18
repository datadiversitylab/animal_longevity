#options(repos = c(
# rtrees = 'https://daijiang.r-universe.dev',
# CRAN = 'https://cloud.r-project.org'))
#install.packages("rtrees")

library(megatrees)
library(rtrees)
library(readr)
library(utils)
library(ape)
library(caper)
library(dplyr)
library(phytools)
library(treedata.table)
library(geiger)
library(here)

options(scipen = 999) #turn off scientific notation

#age data----
agedat_combined1 <- read.csv(here("data", "raw", "agedat_combined1.csv"))

#mammals----
agedat_mammals <- agedat_combined1[agedat_combined1$class == "Mammalia",]
agemammals <- get_tree(sp_list = agedat_mammals,
                       taxon = "mammal",
                       scenario = "at_basal_node",
                       show_grafted = TRUE)

#sum(test_tree[[1]]$tip.label %in% agedat_mammals$species) #1252 species match

d0 <- agedat_mammals[c(26,24)]
rownames(d0) <- d0$species
d0$logage <- log(d0$mean_long)
d0 <- d0[c(1,3)]

mammalsgeig <- geiger::treedata(agemammals[[1]], d0, sort = TRUE, warnings = TRUE)
mammalsgeig$phy
#View(mammalsgeig$data)

write.nexus(mammalsgeig$phy, file = here("data", "processed", "agemammals.trees"))
write_delim(as.data.frame(mammalsgeig$data), file = here("data", "processed", "mammals_agedat.txt"))

## test signal mammals----
mammals.comp <- comparative.data(phy = agemammals[[50]], 
                                 data = agedat_mammals[c(26,24)], 
                                 names.col = "species")

mammals.agedat <- mammals.comp$data
mammals.agedat$logage <- log(mammals.agedat$mean_long)

#using phylosig
mammalssigL <- phylosig(mammals.comp$phy, 
                        mammals.agedat$logage, 
                        method = "lambda", 
                        test = TRUE)
mammalssigK <- phylosig(mammals.comp$phy, 
                        mammals.agedat$logage, 
                        method = "K", 
                        test = TRUE)

#using fitContinuous
logageM <- as.vector(mammals.agedat$logage)
names(logageM) <- rownames(mammals.agedat)

fC_mammalsL <- fitContinuous(multi2di(mammals.comp$phy), logageM, model = "lambda")
fC_mammalsK <- fitContinuous(multi2di(mammals.comp$phy), logageM, model = "kappa")

# birds----
agedat_birds <- agedat_combined1[agedat_combined1$class == "Aves",]
agebirds = get_tree(sp_list = agedat_birds,
                    taxon = "bird",
                    scenario = "at_basal_node",
                    show_grafted = TRUE)

sum(agebirds[[1]]$tip.label %in% agedat_birds$species) #1266 species match

d <- agedat_birds[c(26,24)] 
rownames(d) <- agedat_birds$species
d$logage <- log(d$mean_long)
d <- d[c(1,3)]

birdgeig <- geiger::treedata(agebirds[[1]], d, sort = TRUE, warnings = TRUE)

write.nexus(birdgeig$phy, file = here("data", "processed", "agebirds.trees"))
write_delim(as.data.frame(birdgeig$data), file = here("data", "processed", "birds_agedat.txt"))

## birds signal----
birds.comp <- comparative.data(agebirds[[77]], agedat_birds[c(26,24)], "species", vcv = TRUE)
birds.agedat <- birds.comp$data
birds.agedat$logage <- log(birds.agedat$mean_long)

birdsigL <- phylosig(birds.comp$phy, birds.agedat$logage, method = "lambda", test = T)
birdsigK <- phylosig(birds.comp$phy, birds.agedat$logage, method = "K", test = T)

logageB <- as.vector(birds.agedat$logage)
names(logageB) <- rownames(birds.agedat)

fC_birdsL <- fitContinuous(multi2di(birds.comp$phy), logageB, model = "lambda")
fC_birdsK <- fitContinuous(multi2di(birds.comp$phy), logageB, model = "kappa")

# fish----
agedat_fish <- agedat_combined1[agedat_combined1$class == "Teleostei" | agedat_combined1$class == "Actinopterygii",]
agefish = get_tree(sp_list = agedat_fish,
                   taxon = "fish",
                   scenario = "at_basal_node",
                   show_grafted = TRUE)

d1 <- agedat_fish[c(26,24)] 
rownames(d1) <- d1$species

d1$logage <- log(d1$mean_long)
d1 <- d1[c(1,3)]
#d1

fishgeig <- geiger::treedata(agefish, d1, sort = TRUE, warnings = TRUE)
write.nexus(fishgeig$phy, file = here("data", "processed", "agefish.trees"))
write_delim(as.data.frame(fishgeig$data), file = here("data", "processed", "fish_agedat.txt"))

## fish signal----
fish.comp <- comparative.data(fishgeig$phy, d1, names.col = "species", vcv = TRUE)
fish.comp$phy
fish.comp$data
fish.agedat <- fish.comp$data

fishsigL <- phylosig(fish.comp$phy, fish.agedat$logage, method = "lambda", test = T)
fishsigK <- phylosig(fish.comp$phy, fish.agedat$logage, method = "K", test = T)

#using fitContinuous
logageF <- as.vector(fish.agedat$logage)
names(logageF) <- rownames(fish.agedat)
fC_fishL <- fitContinuous(multi2di(fish.comp$phy), logageF, model = "lambda")
fC_fishK <- fitContinuous(multi2di(fish.comp$phy), logageF, model = "kappa")

# amphibians----
agedat_amphi <- agedat_combined1[agedat_combined1$class == "Amphibia",]
ageamphi = get_tree(sp_list = agedat_amphi,
                    taxon = "amphibian",
                    scenario = "at_basal_node",
                    show_grafted = TRUE)

d2 <- agedat_amphi[c(26,24)] 
rownames(d2) <- agedat_amphi$species
d2$logage <- log(d2$mean_long)
d2 <- d2[c(1,3)]

amphigeig <- geiger::treedata(ageamphi[[1]], d2, sort = TRUE, warnings = TRUE)
write.nexus(amphigeig$phy, file = here("data", "processed", "ageamphi.trees"))
write_delim(as.data.frame(amphigeig$data), file = here("data", "processed", "amphi_agedat.txt"))

## amphibians signal----
amphi.comp <- comparative.data(amphigeig$phy, d2, names.col = "species", vcv = TRUE)

amphi.agedat <- amphi.comp$data


amphisigL <- phylosig(amphi.comp$phy, amphi.agedat$logage, method = "lambda", test = T)
amphisigL
amphisigK <- phylosig(amphi.comp$phy, amphi.agedat$logage, method = "K", test = T)
amphisigK

logageA <- as.vector(amphi.agedat$logage)
names(logageA) <- rownames(amphi.agedat)

fC_amphiL <- fitContinuous(multi2di(amphi.comp$phy), logageA, model = "lambda")
fC_amphiK <- fitContinuous(multi2di(amphi.comp$phy), logageA, model = "kappa")

## ASR amphibians----
ASR_amphi <- contMap(amphi.comp$phy, logageA, plot = FALSE)
ASR_amphi$cols[1:1001]<-colorRampPalette(c("#feba2c","#d6556d","#2a0593"))(1001)

# Plot the mapped characters with the new colors
plot(ASR_amphi, type = "fan", outline = FALSE, legend = 0.7*max(nodeHeights(mammals.comp$phy)),
     fsize = c(0.4, 0.7))

#squamates----
agedat_squamates <- agedat_combined1[agedat_combined1$class == "Reptilia",]

agesquamates = get_tree(sp_list = agedat_squamates,
                        taxon = "reptile",
                        scenario = "at_basal_node",
                        show_grafted = TRUE)
agesquamates[[1]]

d3 <- agedat_squamates[c(26,24)] 
rownames(d3) <- agedat_squamates$species
d3$logage <- log(d3$mean_long)
d3 <- d3[c(1,3)]

repgeig <- geiger::treedata(agesquamates[[1]], d3, sort = TRUE, warnings = TRUE)
write.nexus(repgeig$phy, file = here("data", "processed", "agesquamates.trees"))
write_delim(as.data.frame(repgeig$data), file = here("data", "processed", "rep_agedat.txt"))

squam.comp <- comparative.data(agesquamates[[1]], d3, names.col = "species", vcv = TRUE)
squam.comp$phy

squam.agedat <- squam.comp$data

squamsigL <- phylosig(squam.comp$phy, squam.agedat$logage, method = "lambda", test = T)
squamsigL
squamsigK <- phylosig(squam.comp$phy, squam.agedat$logage, method = "K", test = T)
squamsigK

logageS <- as.vector(squam.agedat$logage)
names(logageS) <- rownames(squam.agedat)

fC_squamL <- fitContinuous(multi2di(squam.comp$phy), logageS, model = "lambda")
fC_squamL
fC_squamK <- fitContinuous(multi2di(squam.comp$phy), logageS, model = "kappa")
fC_squamK

#bats----
agedat_bats <-  agedat_combined1[agedat_combined1$order == "Chiroptera",]
agebats = get_tree(sp_list = agedat_bats,
                        taxon = "mammal",
                        scenario = "at_basal_node",
                        show_grafted = TRUE)

d4 <- agedat_bats[c(26,24)] 
rownames(d4) <- agedat_bats$species
d4$logage <- log(d4$mean_long)
d4 <- d4[c(1,3)]

bats.comp <- comparative.data(agebats[[50]], d4, names.col = "species", vcv = TRUE)
bats.agedat <- bats.comp$data

batsigL <- phylosig(bats.comp$phy, bats.agedat$logage, method = "lambda", test = T)
batsigL
batsigK <- phylosig(bats.comp$phy, bats.agedat$logage, method = "K", test = T)
batsigK


#molluscs----
moll_tree <- read.tree(here("data", "raw", "Stoger_2013.tre"))
agedat_moll <- agedat_combined1[agedat_combined1$phylum == "Mollusca",]
row.names(agedat_moll) <- agedat_moll$species
matched_moll <- geiger::treedata(moll_tree, agedat_moll, sort = TRUE, warnings = TRUE)

#table----
library(sjPlot)
signals <- 
  as.data.frame(cbind(lambda = c(mammalssigL$lambda, birdsigL$lambda, fishsigL$lambda, amphisigL$lambda, squamsigL$lambda),
                      logL = c(mammalssigL$logL, birdsigL$logL, fishsigL$logL, amphisigL$logL, squamsigL$logL), 
                      p = c(mammalssigL$P, birdsigL$P, fishsigL$P, amphisigL$P, squamsigL$P)), 
                row.names = c("mammals","birds","fish","amphibians","reptiles"))
signals$p <- round(signals$p, digits = 4)

signals$count <- c(1252, 1266, 736, 152, 348)

tab_df(signals, show.rownames = TRUE)

write.csv(signals, here("data", "processed", "phylosignal.csv"))

