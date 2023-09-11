# Return a data.frame reporting the expected number of removed individuals as a function of a given ICR threshold 
mind <- function(x=NULL, threshold=NULL) {
  thr <- threshold
  x <- read.table(x, header = T,
                  stringsAsFactors = F)
  y1 <- array()
  y2 <- array()
  for (i in 1:length(thr)) {
    y1[i] <- (
      length(which(x$F_MISS >= thr[i]))/nrow(x)
    )*100
    y2[i] <- length(which(x$F_MISS >= thr[i]))
  }; rm(i)
  y1 <- as.data.frame(cbind(thr, y1, y2))
  colnames(y1) <- c("Threshold", "%_rm", "#_rm")
  return(y1)
}

# Return an histogram of ICRs
mind.hist <- function(x=NULL) {
  x <- read.table(x, header = T,
                  stringsAsFactors = F)
  xr <- range(x$F_MISS)
  xm <- mean(x$F_MISS)
  # windows(width = 6, height = 6)
  hist(x$F_MISS, main="Individual call rate",
       xlab="Missingness (1-ICR)",
       sub=paste("Mean: ", round(xm, 4),
                 " [", round(xr[1], 6),
                 " - ", round(xr[2], 6), "]", sep=""),
       las=1, col="gray", cex.sub=0.8)
}

# Plot the resulting data.frame from the mind function
mind.plot <- function(x=NULL) {
  # windows(width = 6, height = 6)
  plot(x$Threshold, x$`%_rm`,
       main = "Individual call rate",
       xlab = "Threshold",
       ylab = "Removed individuals",
       ylim = c(0, 100), t="l", las=1)
  points(x$Threshold, x$`#_rm`, t="l", col="red")
  legend("topright", 
         legend = c("Percentage", "Number"), 
         lty = 1, col=c("black", "red"), 
         bg = "white")
}

# Return a data.frame reporting the expected number of removed SNPs as a function of a given GCR threshold 
geno <- function(x=NULL, threshold=NULL) {
  thr <- threshold
  x <- read.table(x, header = T,
                  stringsAsFactors = F)
  y1 <- array()
  y2 <- array()
  for (i in 1:length(thr)) {
    y1[i] <- (
      length(which(x$F_MISS >= thr[i]))/nrow(x)
    )*100
    y2[i] <- length(which(x$F_MISS >= thr[i]))
  }; rm(i)
  y1 <- as.data.frame(cbind(thr, y1, y2))
  colnames(y1) <- c("Threshold", "%_rm", "#_rm")
  return(y1)
}

# Return an histogram of the GCRs
geno.hist <- function(x=NULL) {
  x <- read.table(x, header = T,
                  stringsAsFactors = F)
  xr <- range(x$F_MISS)
  xm <- mean(x$F_MISS)
  # windows(width = 6, height = 6)
  hist(x$F_MISS, main="Genotype call rate",
       xlab="Missingness (1-GCR)",
       sub=paste("Mean: ", round(xm, 4),
                 " [", round(xr[1], 6),
                 " - ", round(xr[2], 6), "]", sep=""),
       las=1, col="gray", cex.sub=0.8)
}

# Plot the resulting data.frame from the geno function
geno.plot <- function(x=NULL) {
  # windows(width = 6, height = 6)
  plot(x$Threshold, x$`%_rm`,
       main = "Genotype call rate",
       xlab = "Threshold",
       ylab = "Removed SNPs",
       ylim = c(0, 100), t="l", las=1)
  points(x$Threshold, x$`#_rm`, t="l", col="red")
  legend("topright", 
         legend = c("Percentage", "Number"), 
         lty = 1, col=c("black", "red"), 
         bg = "white")
}

# Return a data.frame reporting the number of SNPs removed as a function of a given MAF threshold 
maf <- function(x=NULL, threshold=NULL) {
  thr <- threshold
  x <- read.table(x, header = T,
                  stringsAsFactors = F)
  y1 <- array()
  y2 <- array()
  for (i in 1:length(thr)) {
    y1[i] <- (
      length(which(x$MAF <= thr[i]))/nrow(x)
    )*100
    y2[i] <- length(which(x$MAF <= thr[i]))
  }; rm(i)
  y1 <- as.data.frame(cbind(thr, y1, y2))
  colnames(y1) <- c("MAF", "%_rm", "#_rm")
  return(y1)
}

# Return an histogram of MAFs
maf.hist <- function(x=NULL) {
  x <- read.table(x, header = T,
                  stringsAsFactors = F)
  xr <- range(x$MAF)
  xm <- mean(x$MAF)
  # windows(width = 6, height = 6)
  hist(x$MAF, main="Minor allele frequency",
       xlab="Minor allele frequency",
       sub=paste("Mean: ", round(xm, 4),
                 " [", round(xr[1], 6),
                 " - ", round(xr[2], 6), "]", sep=""),
       las=1, col="gray", cex.sub=0.8)
}

# Plot the resulting data.frame from the maf function
maf.plot <- function(x=NULL) {
  # windows(width = 6, height = 6)
  plot(x$MAF, x$`%_rm`,
       main = "Minor allele frequency",
       xlab = "Threshold",
       ylab = "Removed SNPs",
       ylim = c(0, 100), t="l", las=1)
  points(x$MAF, x$`#_rm`, t="l", col="red")
  legend("topright", 
         legend = c("Percentage", "Number"), 
         lty = 1, col=c("black", "red"), 
         bg = "white")
}

# Perform IBD analysis
ibd <- function(x=NULL, species=NULL, ibd.t=NULL, icr=NULL, gcr=NULL, maf=NULL) {
	x <- x
	species <- species
	ibd.t <- ibd.t
	icr <- icr
	gcr <- gcr
	maf <- maf
	system(
		paste(
			"plink.exe",
			" --", species,
			" --file ", x,
			" --genome rel-check --min ", ibd.t, 
			sep=""
		)
	)
	ibd <- read.table("plink.genome", 
					  h = T, 
					  stringsAsFactors = F)[, c(1:4, 10)]
	ibd1 <- ibd[, 1:2] 
	ibd1 <- unique(ibd1)
	ibd2 <- ibd[, 3:4]
	ibd2 <- unique(ibd2)
	colnames(ibd2) <- colnames(ibd1)
	ibd12 <- rbind(ibd1, ibd2)
	ibd12 <- unique(ibd12)
	write.table(ibd12, "ibd.txt", 
				col.names = F, row.names = F, 
				quote = FALSE, sep = "\t")
	system(
		paste(
			"plink.exe",
			" --", species,
			" --file ", x,
			" --keep ibd.txt",
			" --recode --out ibd",
			sep=""
		)
	)
	system(
		paste(
			"plink.exe",
			" --", species,	
			" --file ibd", 
			" --missing",
			sep=""
		)
	)
	imiss <- read.table("./plink.imiss", h = T, 
						stringsAsFactors = F)
	ibd$m1 <- rep(NA, nrow(ibd))
	ibd$m2 <- rep(NA, nrow(ibd))
	for (i in 1:nrow(ibd)) {
		ibd[i, 6] <- imiss$F_MISS[which(ibd[i, 2] == imiss$IID)]
	}; rm(i)
	for (i in 1:nrow(ibd)) {
		ibd[i, 7] <- imiss$F_MISS[which(ibd[i, 4] == imiss$IID)]
	}; rm(i)
	ibdres <- data.frame()
	pop <- unique(ibd$FID1)
	for (i in 1:length(pop)) {
		to1 <- subset(x = ibd, subset = ibd$FID1 == pop[i])
		to1res <- data.frame(FID=rep(NA, nrow(to1)), 
							 IID=rep(NA, nrow(to1)))
		to1res[, 1] <- pop[i]
		for (j in 1:nrow(to1)){
			to2 <- to1[j, 6:7]
			if (to2[1] == to2[2]) {
				to1res[j, 2] <- sample(x = c(to1[j, 2], 
								to1[j, 4]), size = 1, 
								prob = c(0.5, 0.5))
			} 
			else if (to2[1] > to2[2]) {
				to1res[j, 2] <- to1[j, 2]
				}
			else to1res[j, 2] <- to1[j, 4]
		}; rm(j)
		rm(to1, to2)
		ibdres <- rbind(ibdres, to1res)
		rm(to1res)
	}; rm(i)
	ibdres <- unique(ibdres)
	write.table(ibdres, "ibd.rm.txt", 
				col.names = F, row.names = F, 
				quote = FALSE, sep = "\t")
	system(
		paste(
			"plink.exe",
			" --", species,
			" --file ", x,
			" --remove ibd.rm.txt ",
			" --mind ", icr,
			" --geno ", gcr,
			" --maf ", maf,
			" --recode --out ",
			paste(x, "_ibd", ibd.t, sep=""),
			sep=""
		)
	)
}
