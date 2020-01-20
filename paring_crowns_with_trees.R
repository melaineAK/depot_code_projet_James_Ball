# 20/01/2019

# pairing algorithm described in Aubry-Kientz 2019 (remote sensing)
# to pair crowns segmented from als with trees from inventories

library(rgeos)
library(rgdal)
library(raster)

rm(list=ls())

#######################################################################################
# data:

#######################################################################################
### Trees (has to be adapted to your data)
  
# dendrot: data of Paracou for all plots. Diameter measured in 2015. Coordinates, identifiaction.
dendrotot <- readOGR(dsn="D://postdoc_Melaine/comparaison_algo/donnees/Paracou/Shapes/Dendro",layer = "dendro_paracou") 
# Keep only trees that are alive in 2015
dendrotot <- dendrotot[dendrotot$Vvn2015=="VRAI",]
# keep only trees that are in the plots
dendrotot <- dendrotot[dendrotot$Parcell == p,]

# order trees by decreasing DBH
dendrotot <- dendrotot[order(dendrotot@data$C2015_c,decreasing=T),]

# number of trees
Ntrees <- nrow(dendrotot)
  
#######################################################################################
### Segmented crowns (has to be adapted to your data)

load(file = paste0("D://postdoc_Melaine/comparaison_algo/results3/AMS3D-computree-Pollock/CrownsData_P",p,"_tree.Rdata"))
parcelles <- readOGR('D://postdoc_Melaine/comparaison_algo/donnees/Paracou/SIGParacou2.0/Parcelles',layer="parcellaire")
  
crowns <- shape

# order crowns by decreasing height
crowns <- crowns[order(crowns@data$z90,decreasing=T),]

# number of crowns
Ncrowns <- nrow(crowns)
  
# same projection
proj4string(crowns) <- dendrotot@proj4string
  
# area of each crown
area <- sapply(slot(crowns, "polygons"), slot, "area")
  
# distance between all crowns and all trees
distances <- spDists(crowns,dendrotot)
colnames(distances) <- dendrotot$idArbre
  
#############################################
## load results of distance model
  
load(file = "D://postdoc_Melaine/appariement/distance tronc_couronne/RankExponentialFit_2.Rdata")
  
#############################################
## load results of allometric model
  
load('D://postdoc_Melaine/appariement/allometrie/posterior_mean_2.Rdata')
  
##
# DBH distribution for each crown (computed using the allometric model)
  
# height
H <- crowns$z90
# diameter (crowns are supposed to be circles...)
CD <- 2 * sqrt(area / pi)
# DBH (all trees)
DBH <- dendrotot$C2015_c / pi
  
# species (corresponds to the species used to calibrate the model, see bayesian_estimation.R)
load("D://postdoc_Melaine/appariement/allometrie/taxo_ref_2.Rdata")
taxo_ref$Espece_complete <- paste(taxo_ref$Espece,taxo_ref$Genre,sep="_")
taxo_ref$name_Espece_complete <- paste(taxo_ref$name_Espece,taxo_ref$name_Genre,sep="_")
Espece_unique <- 1:length(unique(taxo_ref$Espece_complete))
names(Espece_unique) <- unique(taxo_ref$name_Espece_complete)
Gen_unique <- 1:length(unique(taxo_ref$Genre))
names(Gen_unique) <- unique(taxo_ref$name_Genre)
Fam_unique <- 1:length(unique(taxo_ref$Famille))
names(Fam_unique) <- unique(taxo_ref$name_Famille)
  
taxo_trees <- data.frame(cbind(as.character(dendrotot$Espece), as.character(dendrotot$Genre) ,as.character(dendrotot$Famille)))
colnames(taxo_trees) <- c('Espece','Genre','Famille')
taxo_trees$Espece_complete <- paste(taxo_trees$Espece,taxo_trees$Genre,sep="_")
Sp <- Espece_unique[taxo_trees$Espece_complete]

# alpha and beta: parameters of the model, at the species level
alpha_hat <- m[paste0("alpha_sp[",as.character(Sp),"]")]
beta_hat <- m[paste0("beta_sp[",as.character(Sp),"]")]
  
# some species are not in the reference data
# let's try with the genus
Gn <- Gen_unique[as.character(taxo_trees$Genre)]
alpha_hat[is.na(alpha_hat)] <- m[paste0("alpha_gen[",as.character(Gn),"]")][is.na(alpha_hat)]
beta_hat[is.na(beta_hat)] <- m[paste0("beta_gen[",as.character(Gn),"]")][is.na(beta_hat)]
  
# some genus are not in the reference data
# let's try with the family
Fm <- Fam_unique[as.character(taxo_trees$Famille)]
alpha_hat[is.na(alpha_hat)] <- m[paste0("alpha_fam[",as.character(Fm),"]")][is.na(alpha_hat)]
beta_hat[is.na(beta_hat)] <- m[paste0("beta_Fam[",as.character(Fm),"]")][is.na(beta_hat)]
  
# some families are not in the reference data, let's use the hyperparameter
alpha_hat[is.na(alpha_hat)] <- m['alpha0']
beta_hat[is.na(beta_hat)] <- m['beta0']
  
# for each tree and each crown (create matrix with alpha, beta, DBH, CD and H):
alpha_hat <- matrix(alpha_hat, ncol=Ntrees, nrow = Ncrowns, byrow = T)
beta_hat <- matrix(beta_hat, ncol=Ntrees, nrow = Ncrowns, byrow = T)
DBH <- matrix(DBH , ncol=Ntrees, nrow = Ncrowns, byrow = T)
CD <- matrix(CD, ncol=Ntrees, nrow = Ncrowns, byrow = F)
H <- matrix(H, ncol=Ntrees, nrow = Ncrowns, byrow = F)

# d_DBH (g in the paper) 
mu_hat <- alpha_hat + beta_hat * log(CD * H)
d_DBH <- dlnorm(DBH, mean = mu_hat, sd = sqrt(m["V.y"]))
colnames(d_DBH) <- dendrotot$idArbre
  
#######################################################################################
# pairing algo

# ID for crowns
rownames(distances) <- 1:nrow(distances)

# attrib: id of the attributed tree/crown
# ecart: difference between p of the attributed one and the second one
dendrotot$attrib <- dendrotot$ecart <- 0
crowns$attrib <- crowns$ecart <- 0
  
# for each tree
for (i in 1:Ntrees){
    
  # select crowns that are < 15m far from the tree trunk
  crownsnear <- which(distances[,i] < 15)
  # rank by distance
  ranks <- rank(distances[crownsnear,i])
  # compute d (in the paper)
  prob_dist <- dexp((ranks),exp.fit$estimate)
    
  # in the paper : p = d * g 
  d_tot <- d_DBH[crownsnear,i] * prob_dist
  names(d_tot) <- rownames(distances)[crownsnear]
    
  # find the best candidate 
  match <- which.max(d_tot[!names(d_tot)%in%dendrotot$attrib])
  if (length(match)!=0)
  {
    dendrotot$attrib[i] <- names(match)
    crowns$attrib[as.double(names(match))] <- as.character(dendrotot$idArbre[i])
    
    match2 <- which.max(d_tot[!names(d_tot)%in%dendrotot$attrib])
    if (length(match2)!=0)
    {
      dendrotot$ecart[i] <- d_tot[names(match)] - d_tot[names(match2)]
      crowns$ecart[as.double(names(match))] <- d_tot[names(match)] - d_tot[names(match2)]
    }
  }
}
  
  



