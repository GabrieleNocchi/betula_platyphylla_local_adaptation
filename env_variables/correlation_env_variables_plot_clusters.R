source("lib.R", keep.source = TRUE)

# Coordinates
coo <- read.table("coordinates.txt", h = T, stringsAsFactors = F)

# Environmental matrix
env <- matrix(rep(NA, 83 * 26), 83)
colnames(env) <- paste("Bio", 1:26, sep = "")
rownames(env) <- coo$ID


# Let's fill the matrix
for (i in 1:nrow(env)) {
  bio <- raster::getData(
    "worldclim", download = T, 
    lon = coo$LON[i], lat = coo$LAT[i], 
    var = "bio", res = 0.5
    )
  env[i, 1:19] <- raster::extract(
    x = bio, y = coo[i, c("LON", "LAT")]
    )
  }
  
myalt <- raster::getData("alt", country = "CHN", mask = FALSE)
myslope <- raster::terrain(myalt, opt = "slope")
myaspect <- raster::terrain(myalt, opt = "aspect")
myTPI <- raster::terrain(myalt, opt = "TPI")
myTRI <- raster::terrain(myalt, opt = "TRI")
myroughness <- raster::terrain(myalt, opt = "roughness")
myflow <- raster::terrain(myalt, opt = "flowdir")


for (i in 1:nrow(env)) {
  env[i, 20] <- raster::extract(
    x = myalt, y = coo[i, c("LON", "LAT")]
    )
	 env[i, 21] <- raster::extract(
    x = myslope, y = coo[i, c("LON", "LAT")]
    )
	 env[i, 22] <- raster::extract(
    x = myaspect, y = coo[i, c("LON", "LAT")]
    )
	 env[i, 23] <- raster::extract(
    x = myTPI, y = coo[i, c("LON", "LAT")]
    )
	 env[i, 24] <- raster::extract(
    x = myTRI, y = coo[i, c("LON", "LAT")]
    )
	 env[i, 25] <- raster::extract(
    x = myroughness, y = coo[i, c("LON", "LAT")]
    )
	env[i, 26] <- raster::extract(
    x = myflow, y = coo[i, c("LON", "LAT")]
    )
  }
  
# Plotting as an example bio1

uga <- rgdal::readOGR("./CHN_adm0.shp")

# Raster file altitude
alt <- raster::raster("./wc2.1_30s_elev.tif")

# Crop altitude to Uganda extent
alt <- raster::crop(x = alt, raster::extent(uga))

# Mask altitude to Uganda only
alt <- raster::mask(x = alt, mask = uga)


grid <- rgdal::readOGR("./samgrid-epsg4326.shp") 
  

windows()
par(mar=c(5, 5, 1, 5))
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1)
			 axis(side = 1, las=1, col="white", at = c(70,80,90,100,110,120,130,140), 
     tick = T, col.ticks = "gray")
axis(side = 2, las=1, col="white", at = c(20,30,40,50),
     tick = T, col.ticks = "gray")


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))


prettymapr::addscalebar(plotunit = "latlon", pos = "bottomright", 
                        label.cex = 1.1, htin = 0.15,
                        padin = c(0.8, 0.3))
prettymapr::addnortharrow(pos = "topleft", padin = c(0.7, 0.7), 
                          scale = 0.7)
points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env[, "Bio1"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(min(env[, 1]/10), max(env[, 1]/10)), mycol(199), cex = 1, gradient = "v") 
box(col="gray")						  


# Managing redundant information in the environmental matrix --------------
cor <- cor(env)

windows()
corrplot::corrplot(
  cor, 
  order = "original",
  type = "upper", diag = T,
  tl.cex = 0.4,
  tl.col = c(rep("red", 11), 
             rep("blue",8), rep("brown",7)),
  tl.srt=45
  )
  
  
windows()
corrplot::corrplot(
  cor, 
  order = "original",
  type = "upper", diag = T,
  tl.cex = 0.4,
  tl.col = c(rep("red", 11), 
             rep("blue",8), rep("brown",7)),
  tl.srt=45, addCoef.col = "darkgray", addCoefasPercent = T
  )


# Cluster-based approach --------------------------------------------------

# Hierarchy
tree <- ClustOfVar::hclustvar(X.quanti = env)
windows()
plot(tree)

stab <- ClustOfVar::stability(tree, B = 40)


# Suggested number of groups
ng <- names(which(stab$meanCR == max(stab$meanCR)))
ng <- as.numeric(substr(ng, 2, nchar(ng)))
print(ng)

ng <- 11
print(ng)

# Cutting the tree
part <- ClustOfVar::cutreevar(obj = tree, k = ng, matsim = T)
summary(part)


# Central (synthetic) environmental variables for each
# clusters
env.synt <- part$scores
head(env.synt)


# Pairwise correlations among the synthetic
windows()
corrplot::corrplot(
  cor(env.synt), 
  order = "original",
  type = "upper", diag = T,
  tl.cex = 0.8, 
  tl.srt=45, addCoef.col = "black", 
  addCoefasPercent = T
)



# Plotting all clusters
windows()
par(mfrow = c(2,2),mar=c(1,1,1,5))


# Cluster 1
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster1")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster1"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 1]), digits = 2), round(max(env.synt[, 1]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  


text(105,62, "temp. seasonality + mean temp. coldest quarter + \n min temp. coldest month + mean temp. driest quarter + \n temp. annual range + annual mean temp.", cex = 0.5)


# Cluster 2
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster2")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster2"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 2]), digits = 2), round(max(env.synt[, 2]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  

text(105,62, "mean diurnal range", cex = 0.5)


# Cluster 3
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster3")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster3"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 3]), digits = 2), round(max(env.synt[, 3]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  

text(105,62, "isothermality + altitude", cex = 0.5)


# Cluster 4
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster4")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster4"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 4]), digits = 2), round(max(env.synt[, 4]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  


text(105,62, "mean temp. warmest quarter + mean temp. wettest quarter + \n max temp. warmest month", cex = 0.5)


# Second page
windows()
par(mfrow = c(2,2),mar=c(1,1,1,5))


# Cluster 5
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster5")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster5"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 5]), digits = 2), round(max(env.synt[, 5]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  

text(105,62, "prec. warmest quarter + prec. wettest quarter + \n annual prec. + prec. wettest month", cex = 0.5)


# Cluster 6
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster6")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster6"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 6]), digits = 2), round(max(env.synt[, 6]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  


text(105,62, "prec. driest quarter + prec. coldest quarter + \n prec. driest month ", cex = 0.5)


# Cluster 7
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster7")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster7"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 7]), digits = 2), round(max(env.synt[, 7]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  

text(105,62, "prec. seasonality", cex = 0.5)


# Cluster 8
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster8")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster8"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 8]), digits = 2), round(max(env.synt[, 8]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  

text(105,62, "slope + TRI + roughness", cex = 0.5)


# Third page
windows()
par(mfrow = c(2,2),mar=c(1,1,1,5))


# Cluster 9
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster9")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster9"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 9]), digits = 2), round(max(env.synt[, 9]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  

text(105,62, "aspect", cex = 0.5)


# Cluster 10
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster10")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster10"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 10]), digits = 2), round(max(env.synt[, 10]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  


text(105,62, "TPI", cex = 0.5)


# Cluster 11
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1, main = "Cluster11")
			 


raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

mycol <- colorRampPalette(c("blue", "white", "gold", "red"))
raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))



points(
  coo$LON, coo$LAT, pch=21,
  col="black", cex=1.1,
  bg=mycol(199)[
    as.numeric(
      cut(env.synt[, "cluster11"], breaks = 199)
    )
    ]
)

plotrix::color.legend(150, 25, 155, 50, c(round(min(env.synt[, 11]), digits = 2), round(max(env.synt[, 11]), digits = 2)), mycol(199), cex = 0.5, gradient = "v") 				  


text(105,62, "flowdir", cex = 0.5)
