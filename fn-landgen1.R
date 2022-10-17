lfmm_qvalcut <- function (qvalues = NULL, cutoff = NULL) {
  
  qvalues <- qvalues
  cutoff <- cutoff
  res <- list()
  e <- qvalues[, 1:(ncol(qvalues)-3)]
  s <- qvalues[, (ncol(qvalues)-2):ncol(qvalues)]
  
  for (i in 1:ncol(e)) {
    j <- which(e[, i] < cutoff)
    
    if (length(j) == 0) {
      t <- data.frame("ENV" = colnames(e)[i], 
                      "qval" = NA, 
                      "SNP" = NA, 
                      "CHR" = NA, 
                      "BP" = NA)
      res[[i]] <- t
    } else {
      t <- data.frame("ENV" = colnames(e)[i], 
                      "qval" = e[j, i], 
                      "SNP" = s[j, 1], 
                      "CHR" = s[j, 2], 
                      "BP" = s[j, 3])
      res[[i]] <- t
      }
  }
  res <- do.call(rbind.data.frame, res)
  res <- na.omit(res)
  return(res)
}
