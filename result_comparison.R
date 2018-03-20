
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

### COMPARAISON

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

