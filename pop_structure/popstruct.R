as.data.frame(list.files())


# Source libraries and functions ------------------------------------------
source("lib.R", keep.source = TRUE)
source("fn-popstruct-win.R", echo = F, keep.source = TRUE) # for Windows

# plotrix -----------------------------------------------------------------

library(plotrix)


# Molecular dataset -------------------------------------------------------
# Import genomic file

# Name of the dataset
name <- "birch"
# .ped file
ped <- data.table::fread(paste(name, ".ped", sep = ""))
ped[1:6, 1:10]
nrow(ped)

# .map file
map <- data.table::fread(paste(name, ".map", sep = ""))
head(map)
nrow(map)

# Coordinates
coo.qc <- read.table("coordinates.txt",
                  h = T, stringsAsFactors = F)





x = "birch"


# Discriminant analysis of principal components (DAPC) --------------------

# Working dataset for DAPC
system(paste(
  "./plink.exe ", 
  " --file ", x, 
  " --recodeA --allow-extra-chr --out ", x, 
  sep = "")
  )

# Reading input file for DAPC
dapc_input <- read.PLINK(
  paste(x, ".raw", sep = ""),
  parallel = FALSE
  )

# K-means analysis
# Max K to test in K-means analysis
maxk <- 10
grp <- adegenet::find.clusters(
  dapc_input, 
  pca.select = "percVar", perc.pca = 99, 
  max.n.clust = maxk, 
  choose.n.clust = TRUE
  )
grp$grp

# Optimizing the number of PCs to keep...
dapc <- adegenet::dapc(
  dapc_input, grp$grp, 
  pca.select = "percVar", perc.pca = 99
  )

ascore <- adegenet::optim.a.score(dapc)

# DAPC with optimal parameterization
dapc <- adegenet::dapc(
  dapc_input, 
  grp$grp,
  n.pca = ascore$best,
  n.da = length(grp$size)-1
  )
  
  
print(dapc)
adegenet::scatter.dapc(
  dapc, 1, 1, 
  bg="white", 
  legend=T, cleg = 0.6, solid=.4
  )









#########################################
# Spatial distribution of DAPC groups 

# Assigning them group membership as derived by DAPC
coo.qc$assign <- dapc$assign
head(coo.qc)

# Shapefile of Ugandan border
uga <- rgdal::readOGR("./CHN_adm0.shp")

# Raster file altitude
alt <- raster::raster("./wc2.1_30s_elev.tif")

# Crop altitude to Uganda extent
alt <- raster::crop(x = alt, raster::extent(uga))

# Mask altitude to Uganda only
alt <- raster::mask(x = alt, mask = uga)

grid <- rgdal::readOGR("./samgrid-epsg4326.shp")


windows()
par(mar=c(5, 5, 1, 1))

# Empty plot to give the right extent to the map
raster::plot(grid,
             col="white", border="white",
             ylim=c(25, 55), xlim=c(70,140),
             xlab="Longitude (°E)", ylab="Latitude (°N)", 
             las=1, cex.lab=1)
			 axis(side = 1, las=1, col="white", at = c(70,80,90,100,110,120,130,140), 
     tick = T, col.ticks = "gray")
axis(side = 2, las=1, col="white", at = c(20,30,40,50),
     tick = T, col.ticks = "gray")
			 
# Add raster with altitude
raster::plot(
  alt, legend=F, add=T,
  col=rev(gray.colors(n = 1000, alpha = 1))
)

raster::plot(uga, add=T, lwd=1.5, 
             col=scales::alpha("white", 0))


			 
# Add scale bar and North arrow to the map 
prettymapr::addscalebar(plotunit = "latlon", pos = "bottomright", 
                        label.cex = 1.3, htin = 0.15,
                        padin = c(0.8, 0.3))
prettymapr::addnortharrow(pos = "topleft", padin = c(0.4, 0.5), 
                          scale = 0.8)


# Color palette
mycol <- colorRampPalette(colors = c("red", "blue"))

# Create the dataframe for plotting individuals
# with correct symbols/colors
coo.qc$pch <- NA
coo.qc$pch[which(coo.qc$assign == 1)] <- 21
coo.qc$pch[which(coo.qc$assign == 2)] <- 24
coo.qc$pch[which(coo.qc$assign == 3)] <- 22
coo.qc$LD1 <- dapc$ind.coord
coo.qc$Color <- mycol(20)[as.numeric(cut(coo.qc$LD1, breaks = 20))]

# Add individuals to the map with proper
# symbols/colors
points(coo.qc[, 2:3], pch = coo.qc$pch, 
       bg = alpha(coo.qc$Color, 0.6), col = "black", cex = 1.4)
box(col = "gray")







# Fst among DAPC clusters -------------------------------------------------

SNPRelate::snpgdsPED2GDS(
  ped.fn = paste(x, ".ped", sep = ""), 
  map.fn = paste(x, ".map", sep = ""),
  out.gdsfn = paste(x, ".gds", sep = "")
  )
fst_input <- SNPRelate::snpgdsOpen(
  filename = paste(x,".gds", sep = "") 
  )
fst <- SNPRelate::snpgdsFst(
  gdsobj = fst_input, 
  population = dapc$assign,
  method = "W&C84", 
  autosome.only = F, 
  remove.monosnp = FALSE,
  )
SNPRelate::snpgdsClose(fst_input)

fst$Fst
windows()
hist(fst$FstSNP, main="Fst distribution")

# Fst outliers: SNP name
fst.thr <- quantile(fst$FstSNP, 0.99)
fst.snps <- map[which(fst$FstSNP>fst.thr), ][, c(2, 1, 4)]
colnames(fst.snps) <- c("SNP", "CHR", "BP")

# Fst outliers: Fst value
fst.res <- fst$FstSNP[which(fst$FstSNP>fst.thr)]
fst.res <- cbind.data.frame("Fst" = fst.res, fst.snps)
fst.res <- fst.res[order(fst.res$Fst, decreasing = T), ]
fst.res









# tess3r ------------------------------------------------------------------

# Input file for tess3r
tess_input <- data.table::fread(
  paste(x, ".raw", sep = ""), na.strings = "NA"
  )
tess_input <- tess_input[, -c(1:6)]
tess_input <- as.matrix(tess_input)
tess_input[1:5, 1:5]

# Estimating ancestry coefficients and Fst outliers
tess3.obj <- tess3(
  X = tess_input, 
  coord = as.matrix(coo.qc[,2:3]), 
  K = 1:maxk, 
  method = "projected.ls", 
  ploidy = 2, 
  keep = "best"
  )

plot(tess3.obj, pch = 19, col = "blue", 
     main = paste("Cross-validation score vs.",
                  "number of ancestral populations", sep = "\n"),
     xlab = "Number of ancestral populations", 
     ylab = "Cross-validation score")



# Let's combine tess3r results with DAPC results 
# First, let's derive group membership from DAPC
dapc_assign <- as.character(dapc$assign)

# Then, let's prepare input for producing tess3r ancestry barplots 
# with individuals ordered by DAPC grouping
K2 <- as.data.frame(
  tess3r::qmatrix(tess3.obj, K = 2)[, ]
  )
head(K2)
K2 <- data.frame(
  K2, 
  Grp = dapc_assign
  )
head(K2)

# Ancestry coefficients ordered by group
K2 <- K2[order(K2$Grp), ]
head(K2)
K2 <- K2[, -3]

K3 <- as.data.frame(
  tess3r::qmatrix(tess3.obj, K = 3)[, ]
  )
K3 <- data.frame(
  K3, 
  Grp = dapc_assign
  )
K3 <- K3[order(K3$Grp), ]
K3 <- K3[, -4]

K4 <- as.data.frame(
  tess3r::qmatrix(tess3.obj, K = 4)[, ]
  )

# here we add coordinates for future need...
K4 <- data.frame(
  K4, 
  Grp = dapc_assign, 
  coo.qc$LON, 
  coo.qc$LAT)
K4 <- K4[order(K4$Grp), ]


dapc_assign <- K4$Grp
dapc_assign <- as.data.frame(dapc_assign)
colnames(dapc_assign) <- "Cl."
dapc_assign$Cl. <- as.character(dapc_assign$Cl.)
coo.qc.ord <- K4[, c(6, 7)]

K4 <- K4[, -c(5:7)]








###############
# We create a qlist object for the pophelper package
qlist <- list(K2, K3, K4)
qlist <- pophelper::as.qlist(qlist)
qlist <- pophelper::alignK(qlist)
names(qlist) <- c("K=2", "K=3", "K=4")

# Ancestry barplots
windows()
tess3r_bar <- pophelper::plotQ(
  qlist, basesize = 11,
  clustercol = c("blue", "green4", "gold", "red"),
  sortind = "Cluster1", sharedindlab = FALSE, showsubtitle = T,
  subtitlelab = "Global ancestry coefficients", showlegend = F,
  showsp = T, splab = names(qlist), grplab = dapc_assign,
  grplabsize = 4, linesize = 0.5, pointsize = 4,
  imgoutput = "join", returnplot = T, exportplot = F,
  quiet = T, showyaxis = T, panelspacer = 0.25, panelratio = c(3, 1)
  )
grid.arrange(tess3r_bar$plot[[1]])






# Interpolation map K=3
windows()
par(
  mar = c(0.5, 0.5, 0.5, 0.5),
  oma = c(0, 0, 0, 0)
  )
plot(
  as.qmatrix(qlist[[2]]),
  cbind(coo.qc.ord$coo.qc.LON,
        coo.qc.ord$coo.qc.LAT), 
  ylim = c(15, 60), method = "map.max",
  resolution = c(1000, 1000), map.polygon = uga,
  interpolation.model = FieldsKrigModel(10), cex = 0.3,
  main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
  bty = "n", 
  col.palette = CreatePalette(color.vector = c("blue", "green4", 
                                                "red"), 
                              palette.length = 100)
  )

plot(uga, col = scales::alpha("white", 0), lwd = 4, add = T)


# tess3r Fst outliers
tess3.pvalues <- pvalue(tess3.obj, K = 4)
windows()
hist(tess3.pvalues, col="darkgray", border="darkgray", 
     xlab="p-values", main="tess3r Fst ouliers") 

# Benjamini-Hochberg algorithm to corret for false discovery rate
L <- length(tess3.pvalues)
fdr.level <- 0.2
w <- which(sort(tess3.pvalues) < fdr.level * (1:L)/L)
tess3r.res <- order(tess3.pvalues)[w]
length(tess3r.res)
tess3r.res <- map[tess3r.res, c(2, 1, 4)]
colnames(tess3r.res) <- c("SNP", "CHR", "BP")
tess3r.res


# Population structure variables ------------------------------------------
popvar <- qmatrix(tess3.obj, K = 3)
popvar <- as.data.frame(popvar)

# Pairwise correlations among tess3r ancestry
# coefficients
round(cor(popvar), 4)

# Principal component analysis on tess3r ancestry
# coefficients
popvar <- prcomp(popvar)
summary(popvar)

popvar <- popvar$x
popvar <- popvar[, 1:2]
colnames(popvar) <- paste("PSV", 1:2, sep = "")
round(cor(popvar), 4)
head(popvar)

# Spatial representation population structure variables
windows()
par(mar=c(5, 5, 1, 1))

# Empty plot to give the right extent to the map
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
  uga, add = T, lwd = 1.5, 
  col = scales::alpha("white", 0)
  )

prettymapr::addscalebar(
  plotunit = "latlon", 
  pos = "bottomright", 
  label.cex = 0.8, 
  htin = 0.15, padin = c(0.9, 0.2)
  )
prettymapr::addnortharrow(
  pos = "topleft", padin = c(0.7, 0.7), 
  scale = 0.7
  )
rgb <- adegenet::colorplot(
  xy = coo.qc[, 2:3], 
  X = popvar, 
  add.plot = TRUE, 
  transp = T, 
  alpha = 0)
points(
  coo.qc[, 2:3],
  pch = coo.qc$pch,
  bg = alpha(rgb, 0.6),
  col = "black", cex = 1.4
  )
box(col = "gray")
write.table(popvar, "popvar.txt", row.names = F, col.names = T, sep = "\t")