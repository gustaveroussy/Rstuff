## NMF_run
##
## DESCRIPTION : A wrapper to ease the use of the NMF package for R. This allows
##               NMF clustering using various methods, along with capturing the
##               features/variables contributing to the different clusters, and
##               assess the clusters to a null distribution by shuffling data.
##
## AUTHOR:
## Bastien JOB (bastien.job@gustaveroussy.fr)
##
## DEPENDS ON:
## . NMF
## . foreach
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 3.2.0 to 3.3.0
##
## VERSION NOTES:
##
## v1.1b 20160621
## . Added a parameter to modifiy the nmf.option("shared.memory") to inactivate shared
##   memory through big.matrix as it makes nmf crash on R3.3.0 under Ubuntu 14.04 with
##   NMF package version 0.20.6.
## . Also added a parameter to set the seed method.
## . Changed the default seed method used to "none", so that the default fixed seed
##   123456 is actually used (was not the case before this!)
##
## v1.1 20160426
## . Limited the methods for comparison to six methods (others made comparison fail)
## . Changed default maxIter to 5000 (from 3000)
## . Added a help function
## . Corrected a bug where features scores and contributions to metagenes were not saved

## v1.0 20160302
## . First release


## Just gives the names of built-in NMF methods in package NMF
## NOTA : "ls-nmf" and "pe-nmf" seems to fail with my typical datasets...
nmf.methods.list <- function(echo = FALSE) {
  require(NMF)
  nmfm <- nmfAlgorithm()
  if (echo) print(nmfm)
  return(nmfm)
}

## Help for main function
help.nmf_run <- function() {
  cat("
nmf_run <- function(data = NULL, rank = 2:3, method = \"brunet\", default.seed = \"none\", seed = 123456, nrun = 100, maxIter = 5000, shift.if.neg = TRUE, methods.compare = TRUE, methods.list = NULL, track = TRUE, classes = NULL, shuffle = TRUE, ncores = NULL, odir = NULL, shared.memory = TRUE)

data            :  (matrix)             A features (rows) by samples (cols) matrix of numerical values.
rank            :  (integer)            A vector of integers corresponding to desired clusters.
method          :  (character)          The NMF method to use. Run nmf.methods.list() to obtain a list of available methods.
default.seed    :  (character)          The default seed method.
seed            :  (numeric)            The random seed.
nrun            :  (integer)            Nomber of runs to proceed to compute a consensus.
maxIter         :  (integer)            Number of maximum iterations for convergence.
shift.if.neg    :  (logical)            Perform a shift of the matrix if negative values are found.
methods.compare :  (logical)            Perform a quick comparison of NMF methods efficiency using 10 runs for each.
methods.list    :  (character)          Character vector giving the methods to compare if methods.compare is set to true (by defaut, all 6 compatible methods).
track           :  (logical)            Store and plot tracks of residuals.
classes         :  (character/factor)   Known classes for the assessed samples, to compare with identified clusters.
shuffle         :  (logical)            Rerun NMF with shuffled data to compare robustness against a null distribution of your data.
ncores          :  (integer)            Number of threads to use (by default, all available minus one).
odir            :  (character)          The directory where results will be put.
shared.memory   :  (logical)            Activate shared memory (using big.matrix). Set it to FALSE in case of crash related to this issue.

")
}

## Main function
nmf_run <- function(data = NULL, rank = 2:3, method = "brunet", default.seed = "none", seed = 123456, nrun = 100, maxIter = 5000, shift.if.neg = TRUE, methods.compare = TRUE, methods.list = NULL, track = TRUE, classes = NULL, shuffle = TRUE, ncores = NULL, odir = getwd(), shared.memory = TRUE) {
  # data = mymat
  # rank = 2:3
  # method = "brunet"
  # default.seed = "none"
  # seed = 123456
  # nrun = 10
  # maxIter = 5000
  # shift.if.neg = TRUE
  # methods.compare = FALSE
  # methods.list = NULL
  # track = TRUE
  # classes = NULL
  # shuffle = FALSE
  # ncores = 1
  # odir = getwd()
  # shared.memory = TRUE
  
  require(NMF)
  if (is.null(methods.list)) methods.list <- nmf.methods.list()
  if (!is.matrix(data)) stop("'data' must be a numeric matrix !")
  if (!method %in% nmfAlgorithm()) stop("Given NMF algorithm is not supported ! See nmf.methods.list()")
  if (!is.null(classes)) {
    if (length(classes) != ncol(data)) stop("'classes' length and number of columns in 'data' do not match !")
    classes <- as.character(classes)
  }
  if (!is.null(ncores) & !is.numeric(ncores)) stop("'ncores' must be a positive numeric value !")
  if (!dir.exists(odir)) stop("'odir' must be an existing and accessible directory !")
  datmin <- min(data, na.rm = TRUE)
  if (shift.if.neg & datmin < 0) data <- data - datmin
  if (length(unique(rank)) < 2) stop("At least two differentranks are required !")
  if (!is.null(ncores)) nmf.options("cores" = ncores)
  # if (track) nmf.options("track" = track)
  nmf.options("maxIter" = maxIter)
  nmf.options("random.seed" = default.seed)
  nmf.options("grid.patch" = TRUE)
  nmf.options("shared.memory" = shared.memory)
  
  methodword <- sub(pattern = "/", replacement = "-", x = method)
  print("Running NMF ...")
  nmf.res <- NMF::nmf(x = data, rank = rank, method = method, seed = seed, nrun = nrun)
  
  print("Saving results ...")
  oroot <- paste0("NMF_", nrow(data), ".", ncol(data), "_", methodword, "_r", paste(range(rank), collapse = "."), "_n", nrun, "_mI", maxIter)
  odir <- paste0(odir, "/", oroot, "/")
  dir.create(odir)
  oroot <- paste0(odir, oroot)
  # oroot <- paste0(odir, "/NMF_", nrow(data), ".", ncol(data), "_", method, "_r", paste(range(rank), collapse = "."), "_n", nrun)
  
  if (length(rank) > 1) {
    require(foreach)
    mg.scores <- foreach(ra = 1:length(nmf.res$fit)) %do% {
      return(featureScore(nmf.res$fit[[ra]]))
    }
    names(mg.scores) <- names(nmf.res$fit)
    mg.class <- foreach(ra = 1:length(nmf.res$fit)) %do% {
      extf <- extractFeatures(nmf.res$fit[[ra]])
      names(extf) <- seq_along(extf)
      return(extf)
    }
    names(mg.class) <- names(nmf.res$fit)
  } else {
    mg.scores <- featureScore(nmf.res)
    mg.class <- extractFeatures(nmf.res)
  }
  nmfobj <- list(data = data, nmf.options = nmf.options(), nmf.res = nmf.res, nmf.features.score = mg.scores, nmf.features.metagene.class = mg.class)
  saveRDS(object = nmfobj, file = paste0(oroot, ".RDS"))
  
  if (length(rank) > 1) write.table(nmf.res$measures, paste0(oroot, "_measures.txt"), sep="\t", quote = FALSE, row.names = FALSE)
  
  require(foreach)
  nmf.membership <- data.frame(Samples = colnames(nmf.res$fit[[1]]@fit@H), foreach (r = 1:length(nmf.res$measures$rank), .combine = "cbind") %do% {
    as.numeric(predict(nmf.res$fit[[r]]))
  }, stringsAsFactors = FALSE)
  colnames(nmf.membership) <- c("Sample", paste0("rank.", nmf.res$measures$rank))
  write.table(nmf.membership, file = paste0(oroot, "_membership.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  print("Plotting consensus heatmap ...")
  pdf(paste0(oroot, "_consensus.hmap.pdf"), width = 29.7/cm(1), height = 21/cm(1))
  if (is.null(classes)) consensusmap(nmf.res) else consensusmap(nmf.res, annCol = classes)
  dev.off()
  
  if(length(rank) > 1) {
    
    print("Plotting ranks measures ...")
    pdf(paste0(oroot, "_ranks.survey.pdf"), width = 29.7/cm(1), height = 21/cm(1))
    rsplot <- plot(nmf.res, what = c("all", "cophenetic", "rss", "residuals", "dispersion", "evar", "sparseness", "sparseness.basis", "sparseness.coef", "silhouette", "silhouette.coef", "silhouette.basis", "silhouette.consensus"))
    plot(rsplot)
    dev.off()
    
    pdf(paste0(oroot, "_mixt.coeff.hmap.pdf"), width = 29.7/cm(1), height = 21/cm(1))
    for (my.r in 1:length(nmf.res$fit)) coefmap(nmf.res$fit[[my.r]], main = paste0("Mixture coefficients (r=", names(nmf.res$fit)[my.r], ")"))
    dev.off()
    
    if (shuffle) {
      print("Re-running with shuffled data ...")
      data.rand <- randomize(data)
      nmf.res.rand <- nmf(x = data.rand, rank = rank, method = method, seed = seed, nrun = nrun)
      
      print("Saving shuffled results ...")
      saveRDS(object = nmf.res.rand, file = paste0(oroot, "_SHUFFLED.RDS"))
      
      print("Plotting shuffled results ...")
      pdf(paste0(oroot, "_ranks.survey_SHUFFLED.pdf"), width = 29.7/cm(1), height = 21/cm(1))
      rscplot <- plot(nmf.res, nmf.res.rand)
      plot(rscplot)
      dev.off()
    }
  }
  
  if (methods.compare) {
    ## Limiting methods to those that are robust
    comp.methods.list <- methods.list[methods.list %in% c("brunet", "KL", "lee", "Frobenius", "offset", "nsNMF")]
    if (length(comp.methods.list) <2) print("Can't perform comparison with less than two methods !")
    else {
      print("Performing methods comparison ...")
      if (track) pdf(paste0(oroot, "_track.pdf"), width = 29.7/cm(1), height = 21/cm(1))
      for (my.r in rank) {
        print(paste0(" ... for rank ", my.r, "."))
        my.nrun <- ifelse(nrun < 10, nrun, 10)
        nmf.res.comp <- nmf(x = data, rank = my.r, method = as.list(comp.methods.list), seed = seed, nrun = my.nrun, .options = "t")
        nmf.res.comp.df <- compare(nmf.res.comp)
        write.table(nmf.res.comp.df, paste0(odir, "/NMF_", nrow(data), ".", ncol(data), "_methods.compare_r", my.r, "_m", length(comp.methods.list), "_n", my.nrun, ".txt"), sep="\t", row.names = FALSE, quote = FALSE)
        if (track) {
          plot(nmf.res.comp, main = paste0("NMF Residuals (r=", my.r, ")"))
        }
      }
      if (track) dev.off()
    }
  } else if (track) {
    print("Computing track ...")
    pdf(paste0(oroot, "_track.pdf"), width = 29.7/cm(1), height = 21/cm(1))
    for (my.r in rank) {
      nmf.res.t <- nmf(x = data, rank = my.r, method = method, seed = seed, nrun = 1, .options = "t")
      plot(nmf.res.t, main = paste0("NMF Residuals (r=", my.r, ")"))
    }
    dev.off()
  }
}




#######
## DEMO FROM THE NMF VIGNETTE

# data("esGolub")
# 
# 
# ## Reduce (subseting)
# # esGolub <- esGolub[1:200, ]
# 
# ## Remove the unwanted Sample variable
# esGolub$Sample <- NULL
# 
# ## Single run (default)
# system.time(res <- nmf(esGolub, 3))
# 
# ## Retrieving of the fitted model
# fit(res)
# 
# ## Retrieving of the target matrix
# V.hat <- fitted(res)
# dim(V.hat)
# 
# ## Retrieving of the basis matrix (W, giving the metagenes)
# w <- basis(res)
# dim(w)
# 
# ## Retrieving of the mixture coefficient matrix (H, giving the metagenes)
# h <- coef(res)
# dim(h)
# 
# ## Quality and performance
# summary(res)
# summary(res, target = esGolub)
# 
# ## Getting genes scores (contribution to the metagenes)
# s <- featureScore(res)
# 
# ## Getting genes classified to their metagene
# ## NOTA : Only selecting contributing genes !
# sclass <- extractFeatures(res)
# 
# 
# 
# ## SEEDING METHODS (for random seeding)
# ## NOTA : According to the vignette, if the seeding method is deterministic, no need for multiple NMF runs !!
# nmfSeed()
# 
# 
# ## Multiple ranks testing
# ncores <- 6
# estim.r <- nmf(x = esGolub, rank = 2:6, nrun = 10, seed = 123456, .opt = paste0("vp", ncores))
# plot(estim.r)
# consensusmap(estim.r)
# V.random <- randomize(esGolub)
# system.time(estim.r.rand <- nmf(x = V.random, rank = 2:6, nrun = 10, seed = 123456, .opt = paste0("vp", ncores)))
# plot(estim.r, estim.r.rand)
# 
# # estim.mono <- nmf(x = esGolub, rank = 2:6, nrun = 1, seed = 123456, .opt = paste0("vp", ncores))
#   
# ## Comparison of different algorithms (works for a single rank)
# res.multi.method <- nmf(x = esGolub, rank = 3, nrun = 10, list("brunet", "lee", "ns"), seed = 123456, .opt = paste0("vp", ncores))
