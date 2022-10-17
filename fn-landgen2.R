# Log-likelihood ratio test function
log.lik.test <- function(x) {
  
  res_i$snp_geno[x] <- alternative_i$Marker[x]
  t <- which(alternative_i$Marker[x] == null$Marker)
  
  if (length(t) == 1) {
    res_i$Loglik_null[x] <- null$Loglikelihood[t]
    res_i$Loglik_alternative[x] <- alternative_i$Loglikelihood[x]
    
    # Log-lik test statistic
    res_i$D[x] <- -2 * (res_i$Loglik_null[x] - res_i$Loglik_alternative[x])
    
    # P-value
    res_i$D_pval[x] <- pchisq(res_i$D[x], 1, lower = FALSE)

    # Regression coefficients
    res_i[x, 6:(6+dimmax)] <- alternative_i[x, (13+dimmax):ncol(alternative_i)]
      
  } else {
    res_i$Loglik_alternative[x] <- alternative_i$Loglikelihood[x]
    print(paste(tv, " - ", res_i$snp_geno[x], ": unconverged null model", sep = ""))
  }
  
  rm(t)
  return(res_i[x, ])
  
}

# extract qvalue slot from an object created by qvalue::qvalue
my_qvalue <- function(x) {
  q <- qvalue::qvalue(x)
  q <- q$qvalues
  return(q)
}

# Find and match MAP information
match_pp <- function(x) {
  snp <- substr(samLogLik_i$snp_geno[x], 1, 
                nchar(samLogLik_i$snp_geno[x])-3)
  w <- which(snp == as.character(map$V2))
  r <- map[w, c(1, 4)]
  colnames(r) <- c("CHR", "BP")
  return(r)
}

# function to prepare Sambada results for inspection
sam_prepres <- function(x) {
  env <- lapply(sam.res, function(x) {return(length(x$snp_geno))})
  env <- env[which(env > 0)]
  env <- rep(names(env), env)
  qvalue <- as.vector(unlist(lapply(sam.res, function(x) {return(x$qvalue)})))
  snp <- as.vector(unlist(lapply(sam.res, function(x) {return(x$snp_geno)})))
  genotype <- substr(snp, nchar(snp)-1, nchar(snp))
  snp <- substr(snp, 1, nchar(snp)-3)
  chr <- as.vector(unlist(lapply(sam.res, function(x) {return(x$CHR)})))
  bp <- as.vector(unlist(lapply(sam.res, function(x) {return(x$BP)})))
  # as.vector(unlist(lapply(res, function(x) {return(x[, 6:ncol(x)])})))
  outliers <- data.frame("ENV" = env,
                         "qval" = qvalue,
                         "SNP" = snp,
                         "GNT" = genotype,
                         "CHR" = chr,
                         "BP" = bp)
  return(outliers)
  } 

# Mid position in genome bins
middist <- function(x) {
  bin <- min(x) + ((max(x)-min(x))/2)
  return(bin)
}

# Gene annotation
anf <- function(x) {
  anfx <- getBM(attributes=c("ensembl_gene_id", "wikigene_name", 
                             "start_position", "end_position", 
                             "description"),
                filters=c("chromosome_name","start","end"),
                values=list(g_annotation$CHR[x], 
                            g_annotation$Start[x], g_annotation$End[x]),
                mart=ensembl)
  if(dim(anfx)[1] != 0){
    anfx <- cbind.data.frame(g_annotation[x, ], anfx, row.names = NULL)
    return(anfx)
  }
}
