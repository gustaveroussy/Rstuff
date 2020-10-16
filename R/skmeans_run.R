## skmeans_run
##
## DESCRIPTION : A wrapper to ease the use of the skmeans package for R. This allows
##               spherical-kmeans clustering using various methods.
##
## AUTHOR:
## Bastien JOB (bastien.job@gustaveroussy.fr)
##
## DEPENDS ON:
## . skmeans, foreach, kmndirs, cluster
## . (gmeans, CLUTO : external softwares - optional)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 3.3.1
##
## VERSION NOTES:
##
## v1.1 20160426
## . Limited the methods for comparison to six methods (others made comparison fail)
## . Changed default maxIter to 5000 (from 3000)
## . Added a help function
## . Corrected a bug where features scores and contributions to metagenes were not saved

## v1.0 20160302
## . First release


## Just gives the names of built-in skmeans methods
skmeans.methods.list <- function(echo = FALSE) {
  skmlist <- c("standard", "genetic", "pclust", "kmndirs", "gmeans", "cluto")
  if (echo) { print('<<NOTA : "kmndirs" requires the kmndirs package installed ; "gmeans" and "cluto" require the corresponding eponymic external softwares installed.>>'); print(skmlist) }
  return(skmlist)
}

## Help for main function
help.skmeans_run <- function() {
  cat("
skmeans_run <- function(data = NULL, k.test = 2:15, method = \"standard\", maxIter = 5000, gmeans = NULL, gmeans.a = \"s\", vcluster = NULL, clutomethod = \"rbr\", odir = NULL)

data            :  (matrix)             A features (rows) by samples (cols) matrix of numerical values.
k.test          :  (integer)            A vector of integers corresponding to desired clusters.
method          :  (character)          The skmeans method to use. Run skmeans.methods.list() to obtain a list of available methods.
maxIter         :  (numeric)            The maximum number of iterations before stopping.
gmeans          :  (character)          Path to the \"gmeans\" software.
gmeans.a        :  (character)          Submethod for \"gmeans\" (see its own help).
vcluster        :  (character)          Path to the \"vcluster\" software (for CLUTO method).
clutomethod     :  (character)          Submethod for \"vcluster\" (see its own help).
odir            :  (character)          The directory where results will be put (default NULL keeps current working directory).

")
}

skmeans_run <- function(data = NULL, k.test = 2:15, method = "standard", maxIter = 5000, gmeans = NULL, gmeans.a = "s", vcluster = NULL, clutomethod = "rbr", odir = getwd()) {
  
  if (!(method %in% skmeans.methods.list(echo = FALSE))) stop("Unsupported method ! Please use a method listed in skmeans.methods.list()")
  if (is.null(odir)) odir <- getwd()
  if (!dir.exists(odir)) stop(paste0("Can't find output directory : ", odir, " !"))
  if (tolower(method) == "gmeans") if (is.null(gmeans)) stop("Please provide a valid path to the \"gmeans-\" binary !")
  if (tolower(method) == "cluto") if (is.null(vcluster)) stop("Please provide a valid path to the \"vcluster\" binary !")
  
  if (is.null(method)) methodword <- "standard" else if(tolower(method == "cluto")) methodword <- paste0(method, "-", clutomethod) else methodword <- method
  
  require(foreach)
  skmeans.res <- foreach(k = k.test) %do% {
    require(skmeans)
    print(paste0("Testing k=", k, " ..."))
    if (tolower(method) == "standard") {
      skr <- skmeans(x = t(data), method = NULL, k = k, control=list(maxIter = maxIter))
    } else if (tolower(method) %in% c("genetic", "pclust")) {
      skr <- skmeans(x = t(data), method = method, k = k, control=list(maxIter = maxIter))
    } else if (tolower(method) == "kmndirs") {
      require(kmndirs)
      # odir <- paste0(method)
      # skr <- skmeans(x = t(data), method = method, k = k, control=list(maxIter = maxIter))
      skr <- skmeans(x = t(data), method = method, k = k, control=list(nstart = 100, maxIter = 10))
    } else if (tolower(method) == "cluto") {
      # odir <- paste0(method, "-", clmethod)
      skr <- skmeans(x = t(data), method = method, k = k, control=list(vcluster = vcluster, control = paste0("-clmethod=", clutomethod)))
    } else if (tolower(method) == "gmeans") {
      # odir <- paste0(method, "-", gmeans.a)
      skr <- skmeans(x = t(data), method = method, k = k, control=list(gmeans = gmeans, control = paste0(" -a ", gmeans.a)))
    }
    return(skr)
  }
  names(skmeans.res) <- k.test
  
  # oridir <- getwd()
  outdir <- paste0(odir, '/', method)
  # dir.create(methodword, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  # setwd(methodword)
  saveRDS(skmeans.res, paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_results.RDS"))
  
  ## Plotting silhouettes
  require(cluster)
  silh.mean <- silh.q25 <- silh.q75 <- vector()
  pdf(paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_silhouettes.pdf"), width = 29.7/cm(1), height = 21/cm(1))
  for (k in 1:length(skmeans.res)) {
    mysil <- silhouette(skmeans.res[[k]])
    silh.mean <- c(silh.mean, mean(mysil[,3]))
    silh.q25 <- c(silh.q25, quantile(mysil[,3], .25))
    silh.q75 <- c(silh.q75, quantile(mysil[,3], .75))
    plot(mysil, main = paste0("Silhouettes for k = ", names(skmeans.res[k])))
  }
  plot(as.numeric(names(skmeans.res)), as.numeric(silh.q75), type = "b", pch = 20, main = "All silhouettes", xlab = "k", ylab = "Silhouette (Q25, mean, Q75)", col = 2, ylim = range(c(silh.q75, silh.q25, silh.mean), finite = TRUE))
  lines(as.numeric(names(skmeans.res)), silh.mean, type = "b", pch = 20)
  lines(as.numeric(names(skmeans.res)), silh.q25, type = "b", pch = 20, col = 2)
  abline(h = 0, lty = 2)
  dev.off()
  write.table(data.frame(k = k.test, silhouette = silh.mean, stringsAsFactors = FALSE), paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_silhouettes.txt"), sep="\t", quote = FALSE, row.names = FALSE)
  
  crit.values <- vapply(1:length(k.test), function(k) { return(skmeans.res[[k]]$value) }, .1)
  best.idx <- which.max(crit.values)
  bestK <- k.test[best.idx]
  
  ## Ploting criterion values
  png(paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_best.", bestK, ".png"), 1024, 1024)
  plot(k.test, crit.values, type = "b", pch = 20, xlab = "K", ylab = "Criterion", main = paste0("skmeans criterion (to maximize)\nmethod = ", methodword))
  points(bestK, crit.values[best.idx], pch = 18, cex= 2, col = 2)
  segments(bestK, 0, bestK, crit.values[best.idx], lty = 2, col = 2)
  dev.off()
  
  ## Membership
  skmeans.membership <- data.frame(Sample = colnames(data), foreach(k = k.test, .combine = "cbind") %do% skmeans.res[[as.character(k)]]$cluster, stringsAsFactors = FALSE)
  colnames(skmeans.membership) <- c("Sample", k.test)
  saveRDS(skmeans.membership[,c(1,which(colnames(skmeans.membership) == bestK))], paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_membership", "_best.", bestK, ".RDS"))
  write.table(skmeans.membership[,c(1,which(colnames(skmeans.membership) == bestK))], paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_membership", "_best.", bestK, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(skmeans.membership, paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_membership", "_AllK.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  ## AVG plot
  bestk.class <- skmeans.membership[,which(colnames(skmeans.membership) == bestK)]
  medsampz <- foreach(k = unique(bestk.class), .combine = "cbind") %do% { return(vapply(1:nrow(data), function(x) { return(median(data[x,bestk.class == k])) }, .1)) }
  png(paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_best.", bestK, "_mediansamples.png"), width=1650, 1024)
  plot(0,0, type = "n", xlim = c(1,nrow(medsampz)), ylim = range(medsampz), xaxs = "i", xlab = "Value index", ylab = "Value", main = paste0("SKmeans (", method, ") median samples for best K (", bestK, ").\nPopulations = ", paste0(as.vector(table(bestk.class)), collapse = ", ")))
  for (k in 1:ncol(medsampz)) lines(medsampz[,k], col = k)
  dev.off()
  
  # msrange <- vapply(1:nrow(medsampz), function(x) {diff(range(medsampz[x,]))}, .1)
  # png(paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_best.", bestK, "_mediansamples_sorted.png"), width=1650, 1024)
  # plot(0,0, type = "n", xlim = c(1,nrow(medsampz)), ylim = range(medsampz), xaxs = "i", xlab = "Value index", ylab = "Value", main = paste0("SKmeans (", method, ") median samples for best K (", bestK, ").\nPopulations = ", paste0(as.vector(table(bestk.class)), collapse = ", ")))
  # plot(0,0, type = "n", xlim = c(1,100), ylim = range(medsampz), xaxs = "i", xlab = "Value index", ylab = "Value", main = paste0("SKmeans (", method, ") median samples for best K (", bestK, ").\nPopulations = ", paste0(as.vector(table(bestk.class)), collapse = ", ")))
  # for (k in 1:ncol(medsampz)) lines(medsampz[order(msrange, decreasing = TRUE),k], col = k)
  # dev.off()
}

# source("/home/job/svn/genomics/CGH/R/chrload.R")

# k.test <- 2:15
# maxIter <- 1000
# method <- NULL ### NULL corresponds to standard spherical kmeans. Other available are : genetic, pclust
# method <- "genetic"
# method <- "pclust"
# method <- "kmndirs"
# method <- "gmeans"; gmeans <- "/home/job/Tools/gmeans/gmeans-"; gmeans.a <- "s"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "rbr"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "rb"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "direct"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "agglo"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "graph"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "bagglo"

## RUN

# mydatfile <- "/mnt/data_cigogne/job/cna_joint_admixed_analysis/04_most_centered_peak_centering/multipcf_ADMIX2_sub_GC100000_l2r_125s/OS2K6_ADMIX2_sub_Gpeak18_mcpc_125s_hg19_20160621/Data/OS2K6_ADMIX2_Gpeak_mcpc_125s_hg19_20160621.reg"
# mydf <- read.table(mydatfile, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
# mymat <- as.matrix(mydf[,-c(1:9)])

