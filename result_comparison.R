
setwd("/media/ckerouanton/DONNEES/these/analyse_R/ST-DBSCAN_sensitivity/")

### DATA

# resultats
d = read.csv("data/results_stdbscan/r_trajectories.csv")
d = d[, -1]
colnames(d) = stringr::str_replace(colnames(d), "X", "")

# sample
samp = read.csv("data/sample/trajectory_sample.csv")

# parametres
p = read.csv("data/parameters.csv")

# resultats agreges

res = p

### COMPARAISON POINTS

res$param = paste0(res$seq_eps, "_", res$seq_eps2, "_", res$seq_minpts, "_") # id pour jointure du résultat
samp$v_exp2 = as.numeric(samp$v_exp > 0)

for (i in 1:nrow(res)){
  s = d[colnames(d) == res$param[i]][,1] # s étant la séquence de résultat pour un ensemble de paramètres
  s = as.numeric(s > 0)
  res$ntrue[i] = sum(as.numeric(samp$v_exp2 == s)) / nrow(samp) # N(eps, minpts, eps2, E(T))  je ne sais pas comment on nomme "ensemble des trajectoires" alors j'ai mis E(T)
  print(i)
}

write.csv(res, "data/results_comparison/r_comparison_points.csv")
hist(res$ntrue, 40, xlim = c(0,1))


### SUBSET 

res2 = res[res$seq_eps2 == 600,] # a priori eps2 est juste un facteur discriminant, donc je fait un subset pour regarder les effets des deux autres paramètres

lep = c(levels(as.factor(res2$seq_eps))) # x
lmp = c(levels(as.factor(res2$seq_minpts))) # y

try2 = matrix(data = res2$ntrue, nrow = length(lep), ncol = length(lmp)) # z

### REPRESENTATION RESULTAT
library(plot3D)
hist3D(x = as.numeric(lep), y = as.numeric(lmp), z = try2*100, theta = -45, phi = 40, lighting = F, scale = T,
       ticktype="detailed", shade = 0.1, resfac = 3,
       col = "white", 
       border = "black", labels = T, nticks = 5, zlab = "% SUCCESS", xlab = "EPS", ylab = "MINPTS")




############### COMPARAISON SURFACES

source("functions.R")
library(rgeos)
library(sp)
library(sf)



samp$v_exp3 = paste0(samp$idt, "_", samp$v_exp) # vexp3 = identifiant traj + identifiant cluster
samp2 = samp[samp$v_exp3 != "243_6" & !stringr::str_detect(samp$v_exp3, "_0"),] # le 243_6 doit buguer, et la détection "_0" c'est pour enlever les non-cluster
# table(samp$idt)

## Ici, je convertis les points experts en surfaces en utilisant des env convexes 
## (fonction convex_dbscan), ce qui va donner un spatial df stockant les surfaces expert

surfaces = convex_dbscan(cbind(samp2$x, samp2$y), samp2$v_exp3, percent = 100) # transfo des points en env convexe (on pourrait mettre un 95% conf)
surfaces = st_as_sf(surfaces) # conversion en format sf
surfaces$idt = stringr::str_split(surfaces$id, pattern = "_", simplify = T)[,1] # id traj aux surfaces
# plot(st_geometry(surfaces))
id_traj = c(levels(as.factor(samp$idt))) # vecteur pour la boucle


## Dans la boucle suivante, je construis des polygones pour chaque resultat stdbscan, et par trajectoire,
## je compare la surface totale expert et la surface totale stdbscan. La boucle parcourt les differentes 
## combinaisons de parametres

cl = makeCluster(20, type = "MPI") # parallelisation 
registerDoSNOW(cl)



for (i in 1:nrow(res)){
  s = d[,colnames(d) == res$param[i]] # s étant la séquence de résultat pour un ensemble de paramètres
  newd = data.frame(idt = samp$idt, x = samp$x, y = samp$y, cl = unlist(s)) # ensemble de points pour la combinaison
  newd = newd[newd$cl > 0,] # uniquement les points clusters
  
  # s'il n'y a pas de points cluster, alors surface stdbscan == 0
  if(nrow(newd)==0){ 
    res$strue[i] = 0
    next()
  }
  
  newd$cl2 = paste0(newd$idt, "_", newd$cl) # id unique pour le lot de trajectoires et clusters
  
  # enveloppes convexes pour les nuages de points
  surf2 = convex_dbscan(cbind(newd$x, newd$y), newd$cl2, percent = 100)
  surf2 = st_as_sf(surf2)
  #
  
  surf2$idt = stringr::str_split(surf2$id, pattern = "_", simplify = T)[,1] # id trajectoire pour chaque surface
  
  inter = c()
  surfst = c()
  surfex = c()
  
  # boucle pour comparer chaque trajectoire (n = 15)
  for (j in 1:15){
    surfex = c(surfex, sum(st_area(surfaces[surfaces$idt == id_traj[j],])))
    
    # s'il n'y pas de cluster pour une trajectoire, alors pas d'intersection
    if(length(surf2[surf2$idt == id_traj[j],])==0){
      inter = c(inter, 0)
      surfst = c(inter,0)
      next()
    }
    
    # intersection entre cluster expert et cluster stdbscan
    poly_int = st_intersection(subset(surf2, idt == id_traj[j]), subset(surfaces, idt == id_traj[j]))
    
    # surface intersectée, et surface stdbscan
    inter = c(inter, sum(st_area(poly_int)))
    surfst = c(surfst, sum(st_area(surf2[surf2$idt == id_traj[j],])))
  }
  
  res$strue[i] = sum(inter) # vrai == surface correctement intersectée
  res$sfalsepos[i] = sum(surfst)-sum(inter) # faux positif == surface stdbscan en plus de la surface intersectée
  res$sfalseneg[i] = sum(surfex)-sum(inter) # faux négatif == surface expert en plus de la surface intersectée
  print(i)
}
write.csv(res, "data/results_comparison/r_comparison_surface.csv")
# res2 = read.csv("data/results_comparison/r_comparison_surface.csv")

### SUBSET 

res3 = res2[res2$seq_eps2==600,] # a priori eps2 est juste un facteur discriminant, donc je fait un subset pour regarder les effets des deux autres paramètres
res3$strue_per = (res3$strue/sum(st_area(surfaces)))*100
res3$surestim = (res3$sfalsepos/sum(st_area(surfaces)))

lep = c(levels(as.factor(res3$seq_eps))) # x
lmp = c(levels(as.factor(res3$seq_minpts))) # y


### REPRESENTATION RESULTAT
library(plot3D)
library(ggplot2)

try2 = matrix(data = res3$strue_per, nrow = length(lep), ncol = length(lmp)) # z
try3 = matrix(data = res3$surestim, nrow = length(lep), ncol = length(lmp)) # z

# classes pour la couleur
try4 = round(try3, 1)
try4[try4 >= 0] = 2
try4[try3 < 1] = 1
try4[try3 > 2] = 3
try4 = as.numeric(as.factor(try4))

### figure pour comparaison de surfaces
hist3D(x = as.numeric(lep), y = as.numeric(lmp), z = try2, theta = -65, phi = 10, lighting = F, scale = T,
       ticktype="detailed", shade = 0.1,
       opaque.top = TRUE,
       # col = "white",
       col = c("grey80", "grey50", "grey20"),
       colvar = try4, # couleur correspondant à la surestimation (il est possible de ne rien mettre)
       alpha = 0.4,
       border = "black", 
       labels = T, nticks = 5, zlab = "% SUCCESS", xlab = "EPS", ylab = "MINPTS", main = "Surface intersectée")

# ajout de reperes-plan % de réussite z1 et z2 (ici un plan pour 60% de reussite et pour 80%)
z1 = 60
polygon3D(x = c(1,30,30,1,1), y = c(6,6,36,36,6), z = c(z1,z1,z1,z1,z1), alpha = 0.3, add = T, col = "navy")
lines3D(x = c(1,30,30,1,1), y = c(6,6,36,36,6), z = c(z1,z1,z1,z1,z1), lwd = 2, add = T, col = "navy")
z2 = 80
polygon3D(x = c(1,30,30,1,1), y = c(6,6,36,36,6), z = c(z2,z2,z2,z2,z2), alpha = 0.3, add = T, col = "red")
lines3D(x = c(1,30,30,1,1), y = c(6,6,36,36,6), z = c(z2,z2,z2,z2,z2), lwd = 2, add = T, col = "red")


### figure pour la surestimation
hist3D(x = as.numeric(lep), y = as.numeric(lmp), z = sqrt(try3), theta = -45, phi = 40, lighting = F, scale = T,
       ticktype="detailed", shade = 0.1,
       col = "red", 
       border = "black", labels = T, nticks = 5, zlab = "SQRT SURESTIMATION", xlab = "EPS", ylab = "MINPTS", main = "Surestimation surface (S.falsepositive/S.Vexp)")


