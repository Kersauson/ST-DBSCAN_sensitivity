
setwd("/media/ckerouanton/DONNEES/these/analyse_R/ST-DBSCAN_sensitivity/")

### FONCTIONS

source("stdbscan.R")


### DATA FRAME PARAMETRES / RESULTATS
# seq_eps = seq(1,30, 1)
# seq_minpts = seq(6, 36, 1)
# seq_eps2 = seq(100, 1000, 100)
# df = expand.grid(seq_eps = seq_eps, seq_eps2 = seq_eps2, seq_minpts = seq_minpts)
# write.csv(df, "/media/ckerouanton/DONNEES/these/analyse_R/ST-DBSCAN_sensitivity/data/parameters.csv")


### TRAJECTOIRE A COMPARER

# data
d = read.csv("data/sample/trajectory_sample.csv")
p = read.csv("data/parameters.csv")

### ST-DBSCAN SUR TRAJECTOIRES
ldf = nrow(p) # nb de combinaisons
id_traj = c(levels(as.factor(d$idt)))
nbt = length(id_traj)

library(foreach)
library(doSNOW)
library(snow)


cl = makeCluster(20, type = "MPI")
registerDoSNOW(cl)

for (j in 1:nbt){
  t = d[d$idt == id_traj[j],]
  t$tcum = cumsum(t$dt)
  d2 = as.data.frame(matrix(nrow = nrow(p)*nrow(t), ncol = ncol(t)+4))
  colnames(d2) = c(colnames(t), c("cl", "eps", "eps2", "minpts"))
  
  listassign = c(seq(1, nrow(d2), (nrow(d2)/nrow(p))), nrow(d2))
  
  mat = foreach(i = 1:ldf, .combine = cbind) %dopar% {
          stdbscan(t, t$x, t$y, t$tcum, eps = p$seq_eps[i], eps2 = p$seq_eps2[i], minpts = p$seq_minpts[i])$cluster
          # d2$cl[listassign[i]:listassign[i+1]] = stdbscan(t, t$x, t$y, t$tcum, eps = p$seq_eps[i], eps2 = p$seq_eps2[i], minpts = p$seq_minpts[i])$cluster
          # d2$eps[listassign[i]:listassign[i+1]] = p$seq_eps[i]
          # d2$eps2[listassign[i]:listassign[i+1]] = p$seq_eps2[i]
          # d2$minpts[listassign[i]:listassign[i+1]] = p$seq_minpts[i]
  }
  colnames(mat) = paste0(p$seq_eps, "_", p$seq_eps2, "_", p$seq_minpts, "_")
  print(j)
  write.csv(mat, paste0("data/results_stdbscan/r_traj", j, ".csv"))
}

### fichiers en un seul
d4 = data.frame()
for(j in 1:nbt){
  dtemp = read.csv(paste0("data/results_stdbscan/r_traj", j, ".csv"))
  dtemp$X = j
  d4 = rbind(d4, dtemp)
  }

colnames(d4)[1] = c("idt")

write.csv(d4, "/media/ckerouanton/DONNEES/these/analyse_R/ST-DBSCAN_sensitivity/data/results_stdbscan/r_trajectories.csv")

