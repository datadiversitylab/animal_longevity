library(here)
library(utils)
library(ape)
library(dplyr)
library(phytools)
library(ggtree)
library(devtools)
library(readr)
library(ggplot2)


options(scipen=999) #turn off scientific notation

#current as of march 2025
agedat_combined1 <- read.csv(file="/Users/kiranbasava/DDL/AnAge/agedat_combined1.csv")
#may 2025 filtering out low and questionable quality records from AnAge
agedat_filtered <- subset(agedat_combined1, agedat_combined1$AA_Data.quality != 'low')
agedat_filtered <- subset(agedat_filtered, agedat_filtered$AA_Data.quality != 'questionable')

here()
write_excel_csv(agedat_filtered, here("DDL", "AnAge", "animal_longevity", "data", "raw", "agedat_filtered.csv"))
#current as of May 2025


#combining process----
#reading in most recent version feb 18 2025
anagefeb25 <- read.delim(file="/Users/kiranbasava/DDL/AnAge/anage_data.txt", sep="\t")
anagefeb25$binomial <- paste(anagefeb25$Genus, anagefeb25$Species, sep=" ")

anagefeb25$species <- sub(" ", "_", anagefeb25$binomial)

#reading in most recent version of pantheria----
#pretty sure it's the same on since nothings been updated since 2008?
pantheria0508 <- read.delim(file="/Users/kiranbasava/DDL/AnAge/PanTHERIA_1-0_WR05_Aug2008.txt", sep="\t")
colnames(pantheria0508)

panth_anage <- merge(pantheria0508[c(1,2,3,5,7,13,23,28)], anagefeb25[c(3,4,5,6,7,9,10,11,21,28,29,32)], by.x="MSW05_Binomial", by.y="binomial", all.x=TRUE, all.y=TRUE)

#rename variables 
colnames(panth_anage)
colnames(panth_anage)[1] <- "binomial"
colnames(panth_anage)[5] <- "panth_bodymass_g"
colnames(panth_anage)[6] <- "panth_metrate_g"
colnames(panth_anage)[7] <- "panth_maxLong_m"
colnames(panth_anage)[8] <- "panth_Mat_d"

colnames(panth_anage)[15] <-  "aa_femMat_d"
colnames(panth_anage)[16] <-  "aa_maleMat_d"
colnames(panth_anage)[17] <-  "aa_maxLong_y"
colnames(panth_anage)[18] <-  "aa_metrate_w"
colnames(panth_anage)[19] <-  "aa_bodymass_g"

panth_anage$aa_maxLong_m <- panth_anage$aa_maxLong_y*12 #aa variable equivalent to the pantheria one for max long

#change -999s in pantheria to NAs
panth_anage[c(5:8)] <- replace(panth_anage[c(5:8)], panth_anage[c(5:8)]==-999.00, NA) 

write_excel_csv(panth_anage, "/Users/kiranbasava/DDL/AnAge/newfeb21/panth_anage.csv") 
panth_anage <- read.csv("/Users/kiranbasava/DDL/AnAge/newfeb21/panth_anage.csv")

#remove observations without at least one maximum age
chekcols <-  c("panth_maxLong_m", "aa_maxLong_m")
panth_anage_comp <- panth_anage %>% filter(!if_all(all_of(chekcols), is.na)) 

write_excel_csv(panth_anage_comp, "/Users/kiranbasava/DDL/AnAge/newfeb21/panth_anage_comp.csv") 
#cebidae in pantheria --> callitrichidae (AnAge/GBIF)

#filled out missing taxonomic info in AnAge columns and deleted original pantheria columns
panthanage_clean <- read.csv("/Users/kiranbasava/DDL/AnAge/panthanage_clean.csv")

panthanage_clean1 <- merge(panthanage_clean, panth_anage_comp[c(1,18,19)], by="binomial", all.x=T, all.y=F)

#longevity data from Kuparinen et al. ----
#https://www.sciencedirect.com/science/article/pii/S1146609X23000292
kup_dat <- read.delim(file="/Users/kiranbasava/DDL/AnAge/kuparinen_data/final_data.txt", sep = " ")
#LMaxLife is in years
kup_dat$binomial <- paste(kup_dat$Genus, kup_dat$Species, sep=" ")
kup_dat$kup_maxLong_m <- kup_dat$MaxLife*12

#so some of this overlaps with AnAge and some is from other places
panthanagek <- merge(panthanage_clean1, kup_dat, by="binomial", all.x=T, all.y=T)
write_excel_csv(panthanagek, "/Users/kiranbasava/DDL/AnAge/newfeb21/panthanagek.csv") 
#want to keep records from kuparinen that are species not already in current
#I think there are only 4 that don't overlap: Andrias davidianus, Crotalus mitchelli, Cryptomys mechowi, and Scomberomorus commerson
#but this is so much fewer than all the species they list in the supplement from FishBase and Amphiobio

data_mar4 <- read.csv("/Users/kiranbasava/DDL/AnAge/newfeb21/panthanagek.csv")

#now import AnimalTraits data----
AT <- read.csv(file="/Users/kiranbasava/DDL/AnAge/animaltraits_2.25.25.csv")
#body mass is in kg
#metabolic rate is in W
#brain size is kg
#mass-specific metabolic rate is W/kg
#since there are multiple measures for some species, average together 
length(unique(AT$species)) 

AT1 <- AT %>%
  group_by(species) %>%
  summarise(AT_bodymass_g = mean(body.mass))
AT1$AT_bodymass_g <- AT1$AT_bodymass_g*1000

AT2 <- AT %>%
  group_by(species) %>%
  summarise(AT_metrate = mean(metabolic.rate))

AT3 <- AT %>%
  group_by(species) %>%
  summarise(AT_brain_g = mean(brain.size))
AT3$AT_brain_g <- AT3$AT_brain_g*1000

data_mar5 <- merge(data_mar4, AT1, by.x="binomial", by.y="species", all.x=T, all.y=F)
data_mar5 <- merge(data_mar5, AT2, by.x="binomial", by.y="species", all.x=T, all.y=F)
data_mar5 <- merge(data_mar5, AT3, by.x="binomial", by.y="species", all.x=T, all.y=F)

#AnAge has A. dux as having a max longevity of 35 years!!! I think they overestimated from 'estimate of over 25 years' 
#to be fair, this data quality is marked as 'questionable'
data_mar5 <- merge(data_mar5, anagefeb25[c(25,32)], by="binomial", all.x=T, all.y=F) #adding original data quality notations from anage
write_excel_csv(data_mar5, "/Users/kiranbasava/DDL/AnAge/data_mar5.csv") 

#manually filling out a little missing taxonomic data
#changing these fish orders https://genomics.senescence.info/species/browser.php?type=2&name=Actinopterygii to be consistent with AnAge taxonomy
#although AnAge site has them classified as Teleostei as well gahhh

#March 5, 2025 data updated-----
agedat_combined <- read.csv(file="/Users/kiranbasava/DDL/AnAge/data_mar5.csv")

#mean longevity measure----
#create mean measure for max longevity across panth_maxLong_m, aa_maxLong_m
agedat_combined$mean_long <- apply(agedat_combined[c(10,15)], 1, function(row) {
  if (is.na(row[1]) & !is.na(row[2])) {
    return(row[2])
  } else if (!is.na(row[1]) & is.na(row[2])) {
    return(row[1])
  } else if (is.na(row[1]) & is.na(row[2])) {
    return(NA)
  } else {
    return(mean(c(row[1], row[2])))
  }
})

write_excel_csv(agedat_combined, "/Users/kiranbasava/DDL/AnAge/agedat_combined.csv") #copied kuparinen measures for 4 NAs 
agedat_combined <- read.csv(file="/Users/kiranbasava/DDL/AnAge/agedat_combined.csv")

diel <- read.csv("/Users/kiranbasava/DDL/AnAge/Supp_File_10_appendix1914.csv")

agedat_combined$binomial %in% diel$Species

agedat_combined1 <- merge(agedat_combined, diel[c(1,2)], by.x="binomial", by.y="Species", all.x=T, all.y=F)

agedat_combined1[228,2] <- "Chordata"

write_excel_csv(agedat_combined1, file="/Users/kiranbasava/DDL/AnAge/agedat_combined1.csv")

agedat_combined1$treesp <- sub(" ", "_", agedat_combined1$binomial)


#diet, for family level----
diet <- read.csv("/Users/kiranbasava/DDL/AnAge/diet.csv")
length(unique(diet$Family)) #1048 families
# mean of max longevity for each family
family_maxlong <- agedat_combined1[c(2:5)]
family_maxlong <- unique(family_maxlong)

family_maxlong1 <- agedat_combined1[c(2:5,24)] %>% group_by(Family) %>% summarise(maxlongmean = mean(mean_long))

family_maxlong <- merge(family_maxlong, family_maxlong1, by="Family")

write_excel_csv(family_maxlong, "/Users/kiranbasava/DDL/AnAge/family_maxlong.csv")

family_maxlong2 <- merge(family_maxlong, diet[c(5,8,9)], by="Family", all.x=T, all.y=F)

diet$Family %in%  agedat_combined1$Family
#ok that's like 14 matches, not great

#using mammal-specific data, Price et al. 2012----
price2012_diet <- read.delim("/Users/kiranbasava/DDL/AnAge/price_mammaldiet.txt") #1515 species
sum(price2012_diet$Species.Name...Wilson...Reeder.2005 %in% agedat_combined1$species) #1218




#taxonomic frequency estimates----
View(table(agedat_combined$Family))
View(table(agedat_combined$Class))


View(table(agedat_combined$Order))
ggplot(agedat_combined, aes(x = reorder(as.factor(Order), -table(Order)[Order]))) +  
  geom_bar() +
  labs(x = "Order", y = "Count in database") +
  theme(axis.text.x = element_text(size=7, angle=75, hjust=1))

ggplot(agedat_combined, aes(x = reorder(as.factor(Class), -table(Class)[Class]))) +  
  geom_bar() +
  labs(x = "Class", y = "Count in database") +
  theme(axis.text.x = element_text(size=7, angle=75, hjust=1))  +
  scale_y_continuous(breaks=seq(0,1400,100))


ggplot(agedat_combined, aes(x = reorder(as.factor(Phylum), -table(Phylum)[Phylum]))) +  
  geom_bar() +
  labs(x = "Phylum", y = "Count in database") +
  theme(axis.text.x = element_text(size=7, angle=75, hjust=1)) 


#means for taxonomic levels
mean(agedat_combined$mean_long[agedat_combined$Phylum == "Chordata"], na.rm = T) #217.7092
mean(agedat_combined$mean_long[agedat_combined$Phylum == "Chordata"], na.rm = T)

#climateNiche https://github.com/RS-eco/climateNiche----
library(rgbif)
library(raster)
library(sf)
library(ggplot2)
library(rebird)
library(remotes)
library(outliers)

packages_github <- c("rasterSp", "climateNiche", "ggmap2", "rISIMIP")

remotes::install_github("RS-eco/rasterSp")
remotes::install_github("RS-eco/climateNiche")
remotes::install_github("RS-eco/ggmap2")

library(rasterSp)
library(climateNiche)
library(ggmap2)
install.packages("geodata")
library(geodata)

data(amphibians_dist, package="rasterSp")
data(ter_birds_dist, package="rasterSp")
data(ter_mammals_dist, package="rasterSp")
data(reptiles_dist, package="rasterSp")

filedir <- "/Users/kiranbasava/DDL/AnAge/filedir/"

bioclim <- getData(name="worldclim", download=TRUE, path=filedir, res=10, var="bio") 

#this didn't work so I downloaded it here https://www.worldclim.org/data/worldclim21.html
bioclim <- worldclim_global(var="bio", res=10, path=filedir, version="2.1")
fun <- function(x) {x/10}

bioclim_temp <- calc(bioclim[[c(1,2,5,6,7,8,9,10,11)]], fun=function(x){x/10}) #this gives an error
bioclim <- stack(bioclim_temp[[c(1:2)]], bioclim[[c(3:4)]], bioclim_temp[[c(3:9)]], bioclim[[c(12:19)]]); rm(bioclim_temp)
names(bioclim) <- paste0("bio", 1:19)

bioclim_temp <- app(bioclim[[c(1,2,5,6,7,8,9,10,11)]], fun=fun)
bioclim <- stack(bioclim_temp[[c(1:2)]], bioclim[[c(3:4)]], bioclim_temp[[c(3:9)]], bioclim[[c(12:19)]])
names(bioclim) <- paste0("bio", 1:19)

#mean body mass measure----
agedat_combined1$mean_mass <- apply(agedat_combined[c(10,15)], 1, function(row) {
  if (is.na(row[1]) & !is.na(row[2])) {
    return(row[2])
  } else if (!is.na(row[1]) & is.na(row[2])) {
    return(row[1])
  } else if (is.na(row[1]) & is.na(row[2])) {
    return(NA)
  } else {
    return(mean(c(row[1], row[2])))
  }
})






#also kind of want to add ceph data----
cephlong <- read.csv("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/cephdat.csv")
cephlong <- cephlong[c(1,14)]
cephlong <- na.omit(cephlong) #45 species, maximum age in days
cephlong$ceph_maxlong_m <- cephlong$lifespan.max/30
#likely giant squid should be changed to 24-36 months not 170
#that would make its lifespan considerably shorter than max of much smaller V. infernalis
#idk
