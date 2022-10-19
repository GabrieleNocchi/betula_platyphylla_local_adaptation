# LFMM2 ANALYSIS
library(raster)
library(LEA)
library(scales)
source("lib.R", keep.source = TRUE)


# Landscape genomics with LFMM --------------------------------------------
# Latent factor mixed models (LFMM)
# x must be the prefix of the raw file
# The raw format file can be generated with plink (--recodeA) from map/ped files
x = "birch"
library(data.table)
# Molecular input file
Y <- data.table::fread(paste(x, ".raw", sep = ""),
                       na.strings = "NA", header = T)
Y <- Y[, -c(1:6)]
fwrite(Y, paste(x, ".lfmm", sep = ""), col.names = F,
            row.names = F, sep = "\t", na = "9")
Y <- LEA::lfmm2geno(paste(x, ".lfmm", sep = ""))





# Environmental input file provided
X <- read.table("environmental_data.txt", h = T, stringsAsFactors = F)
X <- as.matrix(X[, c(4:14)])




# PCA
pc <- LEA::pca(paste(x, ".lfmm", sep = ""), scale = TRUE)



# Tracy-Widom tests on all eigenvalues to identify significant components
tw <- LEA::tracy.widom(pc)
print(tw$pvalues[1:10])

# Scree plot
windows()
plot(tw$eigenvalues, main = "Scree plot",
     ylab = "Eigenvalues", xlab = "PCs", t = "b")




# sNMF (sparse nonnegative matrix factorization)
obj.snmf <- LEA::snmf(Y, K = 1:10, entropy = T, ploidy = 2,
                      project = "new")



# let’inspect the values of the cross-entropy criterion for each K:
plot(obj.snmf, pch = 16, col = "blue")


K <- 3

barplot(t(Q(obj.snmf, K = K)),
        col = c("forestgreen",
                "blue","pink"), border = NA, las=2)
box(lwd=2)


# Reading imputed genotype matrix
Y <- data.table::fread(paste(x, ".lfmm",
                             sep = ""), header = F)




# Running LFMM
mod.lfmm <- lfmm::lfmm_ridge(Y = Y,
                             X = X,
                             K = K)
str(mod.lfmm)

# P-values
pv <- lfmm::lfmm_test(Y = Y, X = X, lfmm = mod.lfmm,
                      calibrate = "gif")


pvalues <- pv$calibrated.pvalue
dim(pvalues)
head(pvalues)



# P-value histograms per env. variables
# Histograms of p-values
par(mfrow = c(3, 3), mar = c(3.5, 3.5, 3, 0.5),
    mgp = c(2.5, 0.8, 0))
for (i in 1:ncol(pvalues)) {
  hist(pvalues[, i],
       breaks = 20, xlab = "p-values",
       xlim = c(0, 1),
       main = paste("Env. variable ", i, sep = ""),
       col = "darkgray", border = "darkgray")
}
# Q-Q plot
qqplot(rexp(length(pvalues),
            rate = log(10)),
       -log10(pvalues),
       xlab = "Expected quantile",
       ylab = expression(-log[10] ~ p - values),
       pch = 19, cex = 0.4)
abline(0, 1)


# False discovery rate control with q-values  -----------------------------

pvalues_env1 <- pvalues[, 1]
qobj <- qvalue::qvalue(pvalues_env1)
hist(qobj)
plot(qobj)

# Q-values need to be computed by environmental variable, so let’s define a
# function that transform p- into p-values by column
head(pvalues)
my_qvalue <- function(x) {
  q <- qvalue::qvalue(x)
  q <- q$qvalues
  return(q)
}

qvalues <- apply(pvalues, 2, my_qvalue)
head(qvalues)

# Q-value histograms per env. variables
par(mfrow = c(3, 3), mar = c(3.5, 3.5, 3, 0.5),
    mgp = c(2.5, 0.8, 0))
for (i in 1:ncol(qvalues)) {
  hist(qvalues[, i], breaks = 20, xlab = "q-values",
       xlim = c(0, 1),
       main = paste("Env. variable ", i, sep = ""),
       col = "darkgray", border = "darkgray")
}

# Let's add SNP information (name, chromosome,
# physical position)  -- use the map file --
map <- read.table(paste("birch", ".map", sep = ""))
qvalues <- as.data.frame(qvalues)
qvalues$SNP <- as.character(map$V2)
qvalues$CHR <- map$V1
qvalues$BP <- map$V4
head(qvalues)

# Let's decide a q-value cut-off of 1%
qcut <- 0.01
# LFMM results
source("fn-landgen1.R", echo = F, keep.source = TRUE)
lfmm.res <- lfmm_qvalcut(qvalues = qvalues, cutoff = qcut)
print(lfmm.res)

# Manhattan plot
pvalues <- as.data.frame(pvalues)
pvalues$SNP <- as.character(map$V2)
pvalues$CHR <- map$V1
pvalues$BP <- map$V4



# Dataframe with ordered columns for Manhattan
# plots
manh <- pvalues[, c(12:14, 1:11)]
head(manh)




#### PLOT SINGLE ####
par(mfrow = c(3, 1), oma = c(1, 1, 1, 1),
    mar = c(4, 5, 3, 1), mgp = c(2.5, 0.8, 0))
for (i in 4:ncol(manh)) {
  manh_i <- manh[, c(1:3, i)]
  SNPs <- lfmm.res$SNP[lfmm.res$ENV == colnames(manh_i)[4]]
  colnames(manh_i)[4] <- "P"
  qqman::manhattan(
    manh_i, genomewideline = F, suggestiveline = F,
    main = paste("LFMM - Env. variable ", colnames(manh)[i], sep = ""),
    highlight = SNPs,
    col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)),
    xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1
    )
  box()
}


#### Plot Splitted ####
png("results_1.png")
par(mfrow = c(3, 1), oma = c(1, 1, 1, 1),
    mar = c(4, 5, 3, 1), mgp = c(2.5, 0.8, 0))
for (i in 4:6) {
  manh_i <- manh[, c(1:3, i)]
  SNPs <- lfmm.res$SNP[lfmm.res$ENV == colnames(manh_i)[4]]
  colnames(manh_i)[4] <- "P"
  qqman::manhattan(
    manh_i, genomewideline = F, suggestiveline = F,
    main = paste("LFMM - Env. variable ", colnames(manh)[i], sep = ""),
    highlight = SNPs,
    col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)),
    xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1
    )
  box()
}
dev.off()

png("results_2.png")
par(mfrow = c(3, 1), oma = c(1, 1, 1, 1),
    mar = c(4, 5, 3, 1), mgp = c(2.5, 0.8, 0))
for (i in 7:9) {
  manh_i <- manh[, c(1:3, i)]
  SNPs <- lfmm.res$SNP[lfmm.res$ENV == colnames(manh_i)[4]]
  colnames(manh_i)[4] <- "P"
  qqman::manhattan(
    manh_i, genomewideline = F, suggestiveline = F,
    main = paste("LFMM - Env. variable ", colnames(manh)[i], sep = ""),
    highlight = SNPs,
    col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)),
    xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1
    )
  box()
}
dev.off()

png("results_3.png")
par(mfrow = c(3, 1), oma = c(1, 1, 1, 1),
    mar = c(4, 5, 3, 1), mgp = c(2.5, 0.8, 0))
for (i in 10:12) {
  manh_i <- manh[, c(1:3, i)]
  SNPs <- lfmm.res$SNP[lfmm.res$ENV == colnames(manh_i)[4]]
  colnames(manh_i)[4] <- "P"
  qqman::manhattan(
    manh_i, genomewideline = F, suggestiveline = F,
    main = paste("LFMM - Env. variable ", colnames(manh)[i], sep = ""),
    highlight = SNPs,
    col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)),
    xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1
    )
  box()
}
dev.off()

png("results_4.png")
par(mfrow = c(3, 1), oma = c(1, 1, 1, 1),
    mar = c(4, 5, 3, 1), mgp = c(2.5, 0.8, 0))
for (i in 13:14) {
  manh_i <- manh[, c(1:3, i)]
  SNPs <- lfmm.res$SNP[lfmm.res$ENV == colnames(manh_i)[4]]
  colnames(manh_i)[4] <- "P"
  qqman::manhattan(
    manh_i, genomewideline = F, suggestiveline = F,
    main = paste("LFMM - Env. variable ", colnames(manh)[i], sep = ""),
    highlight = SNPs,
    col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)),
    xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1
    )
  box()
}
dev.off()
