### FONCTIONS NECESSAIRES POUR LA COMPARAISON SURFACES

# SPATIAL_TRAJ93
# comme son nom l'indique, fonction pour convertir rapidos un df en spatial projeté lambert93

spatial_traj93 <- function(traji, longitude, latitude){
  traji <- subset(traji, is.na(longitude)!= TRUE | is.na(latitude)!= TRUE)
  traji2 <- data.frame(longitude, latitude) # recup coords from 1st file
  traj_sp <- SpatialPointsDataFrame(traji2, traji, proj4string = CRS("+init=epsg:2154")) # coords et df into spatial
  return(traj_sp)
}

# CONVEX_MCP
# utilisation d'une fonction du package "adehabitat LT" pour produire l'enveloppe convexe
# d'un ensemble de points, à partir d'un critère de confiance 'percent'
# PS : j'ai desactivé le critère qui permettait d'arrêter la fonction si < 5 points

# output = spatial polygon avec surface (pour la surface, à vérifier)

convex_mcp = function (xy, percent = 95, unin = c("m", "km"), unout = c("ha", "km2", "m2")) {
  if (!inherits(xy, "SpatialPoints")) 
    stop("xy should be of class SpatialPoints")
  if (ncol(coordinates(xy)) > 2) 
    stop("xy should be defined in two dimensions")
  pfs <- proj4string(xy)
  if (length(percent) > 1) 
    stop("only one value is required for percent")
  if (percent > 100) {
    warning("The MCP is estimated using all relocations (percent>100)")
    percent <- 100
  }
  unin <- match.arg(unin)
  unout <- match.arg(unout)
  if (inherits(xy, "SpatialPointsDataFrame")) {
    if (ncol(xy) != 1) {
      warning("xy should contain only one column (the id of the animals), id ignored")
      id <- factor(rep("a", nrow(as.data.frame(xy))))
    }
    else {
      id <- xy[[1]]
    }
  }
  else {
    id <- factor(rep("a", nrow(as.data.frame(xy))))
  }
  if (percent > 100) {
    warning("The MCP is estimated using all relocations (percent>100)")
    percent <- 100
  }
  # if (min(table(id)) < 5) 
  #   stop("At least 5 relocations are required to fit an home range")
  id <- factor(id)
  xy <- as.data.frame(coordinates(xy))
  r <- split(xy, id)
  est.cdg <- function(xy) apply(xy, 2, mean)
  cdg <- lapply(r, est.cdg)
  levid <- levels(id)
  res <- SpatialPolygons(lapply(1:length(r), function(i) {
    k <- levid[i]
    df.t <- r[[levid[i]]]
    cdg.t <- cdg[[levid[i]]]
    dist.cdg <- function(xyt) {
      d <- sqrt(((xyt[1] - cdg.t[1])^2) + ((xyt[2] - cdg.t[2])^2))
      return(d)
    }
    di <- apply(df.t, 1, dist.cdg)
    key <- c(1:length(di))
    acons <- key[di <= quantile(di, percent/100)]
    xy.t <- df.t[acons, ]
    coords.t <- chull(xy.t[, 1], xy.t[, 2])
    xy.bord <- xy.t[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    so <- Polygons(list(Polygon(as.matrix(xy.bord))), k)
    return(so)
  }))
  # are <- unlist(lapply(1:length(res), function(i) {
  #   .arcpspdf(res[i, ])
  # }))
  # if (unin == "m") {
  #   if (unout == "ha") 
  #     are <- are/10000
  #   if (unout == "km2") 
  #     are <- are/1e+06
  # }
  # if (unin == "km") {
  #   if (unout == "ha") 
  #     are <- are * 100
  #   if (unout == "m2") 
  #     are <- are * 1e+06
  # }
  df <- data.frame(id = unlist(lapply(1:nlevels(id), function(i) res[i]@polygons[[1]]@ID)))
  row.names(df) <- df[, 1]
  res <- SpatialPolygonsDataFrame(res, df)
  if (!is.na(pfs)) 
    proj4string(res) <- CRS(pfs)
  return(res)
}


## CONVEX_DBSCAN
# production d'une enveloppe convexe pour un cluster dbscan
# utilisation de convex_mcp, avec 95% des points pre-enregistre

# input : xy (df des x et y des points), cl (liste des id cluster des points)
# output : enveloppes convexes (sp polygons) pour chaque identifiant de cluster (chaque nuage dbscan)

convex_dbscan = function(xy, cl, percent = 95){
  xy = data.frame(x = xy[,1], y = xy[,2], cl = cl)
  xy = xy[xy$cl != 0,]
  ps = spatial_traj93(xy, xy$x, xy$y)
  ps@data = data.frame(cl = xy$cl)
  # convex = adehabitatHR::mcp(ps, percent = percent)
  convex = convex_mcp(ps, percent = percent)
  
  return(convex)
}
