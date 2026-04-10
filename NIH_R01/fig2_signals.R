### FIGURE 2----
# Phylogenetic signal in longevity across vertebrate clades, make a bar chart of signal from table

#rerun signals OR skip to read in table of signals----
agemammals <- read.nexus((here("data", "processed", "agemammals.trees")))
mammals.comp <- comparative.data(phy = agemammals, 
                                 data = agedat_mammals[c(26,24)], 
                                 names.col = "species")

mammals.agedat <- mammals.comp$data
mammals.agedat$logage <- log(mammals.agedat$mean_long)

mammalssigL <- phylosig(mammals.comp$phy, 
                        mammals.agedat$logage, 
                        method = "lambda", 
                        test = TRUE)

## birds signal----
agedat_birds <- agedat_combined1[agedat_combined1$class == "Aves",]

agebirds <- read.nexus((here("data", "processed", "agebirds.trees")))

birds.comp <- comparative.data(agebirds, agedat_birds[c(26,24)], "species")
birds.agedat <- birds.comp$data
birds.agedat$logage <- log(birds.agedat$mean_long)

birdsigL <- phylosig(birds.comp$phy, birds.agedat$logage, method = "lambda", test = T)

## fish signal----
agedat_fish <- agedat_combined1[agedat_combined1$class == "Teleostei" | agedat_combined1$class == "Actinopterygii",]
agefish <- read.nexus((here("data", "processed", "agefish.trees")))
fish.comp <- comparative.data(agefish,  agedat_fish[c(26,24)], names.col = "species")
fish.agedat <- fish.comp$data
fish.agedat$logage <- log(fish.agedat$mean_long)
fishsigL <- phylosig(fish.comp$phy, fish.agedat$logage, method = "lambda", test = T)

#amphibians signal----
agedat_amphi <- agedat_combined1[agedat_combined1$class == "Amphibia",] 
ageamphi <- read.nexus((here("data", "processed", "ageamphi.trees")))

amphi.comp <- comparative.data(ageamphi, agedat_amphi[c(26,24)], names.col = "species")

amphi.agedat <- amphi.comp$data
amphi.agedat$logage <- log(amphi.agedat$mean_long)

amphisigL <- phylosig(amphi.comp$phy, amphi.agedat$logage, method = "lambda", test = T)

#squamates signal----
agedat_squamates <- agedat_combined1[agedat_combined1$class == "Reptilia",]
agesquamates <- read.nexus((here("data", "processed", "agesquamates.trees")))

squam.comp <- comparative.data(agesquamates, agedat_squamates[c(26,24)], names.col = "species")
squam.agedat <- squam.comp$data
squam.agedat$logage <- log(squam.agedat$mean_long)
squamsigL <- phylosig(squam.comp$phy, squam.agedat$logage, method = "lambda", test = T)

library(sjPlot)

signals <- 
  as.data.frame(cbind(lambda = c(mammalssigL$lambda, birdsigL$lambda, fishsigL$lambda, amphisigL$lambda, squamsigL$lambda),
                      logL = c(mammalssigL$logL, birdsigL$logL, fishsigL$logL, amphisigL$logL, squamsigL$logL), 
                      p = c(mammalssigL$P, birdsigL$P, fishsigL$P, amphisigL$P, squamsigL$P)))
signals$p <- round(signals$p, digits = 4)
signals$lambda <- round(signals$lambda, digits = 3)
signals$group <- c("mammals", "birds", "fish", "amphibians", "reptiles")


#or read in table of signals----
library(readr)
signals <- read_csv(here("NIH_R01", "phylosignal.csv"))

#bar plot signals----
library(ggplot2)
signalsplot <- ggplot(signals, aes(x=group, y = lambda, fill=group)) + 
  geom_bar(stat="identity")

ggsave("signalsplot.jpeg", signalsplot, width=10, units="in")
