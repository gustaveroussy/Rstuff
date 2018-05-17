## Build snpp file
facets.build.snpp <- function(snp.pileup.path = NULL, normal.bam = NULL, tumor.bam = NULL, vcf.file = NULL, samplename =  NULL, min.depth = 20, min.baseQ = 20) {

  if (is.null(normal.bam)) stop("Normal BAM file is required !")
  if (!file.exists(normal.bam)) stop("Normal BAM not found !")
  if (is.null(tumor.bam)) stop("Tumor BAM file is required !")
  if (!file.exists(tumor.bam)) stop("Tumor BAM not found !")
  if (is.null(vcf.file)) stop("A VCF file is required !")
  if (!file.exists(vcf.file)) stop("VCF file not found !")
  if (is.null(samplename)) stop("A sample name is required !")
  if (!is.null(snp.pileup.path)) snp.pileup.path <- paste0(snp.pileup.path, "/")
  
  system(paste0(
    snp.pileup.path,"snp-pileup",
    " -p",
    " -g", 
    " -r ", min.depth,
    " -Q ", min.baseQ,
    " ", vcf.file, 
    " ", samplename, ".snpp",
    " ", normal.bam,
    " ", tumor.bam
  ))
}

## Run facets on snpp file
facets.run <- function(snpp.file = NULL, samplename = NULL, genome.build = NULL, ndepth = 35, hetscale = TRUE, unmatched = FALSE, cval = 100) {
  
  print(samplename)
  
  if (is.null(snpp.file)) stop("SNPP file is required !")
  if (is.null(samplename)) stop("samplename file is required !")
  if (is.null(genome.build)) stop("A genome build is required !")
  
  require(facets, quietly = TRUE)
  
  ## Load SNPmatrix data
  message("Loading SNPmatrix data ...")
  facets.snpm <- facets::readSnpMatrix(filename = paste0(samplename, ".snpp.gz"))
  
    ## TWEAK TO INVERT N and T
    # colnames(facets.snpm) <- colnames(facets.snpm)[c(1,2,5,6,3,4)]
    
  ## Pre-process
  message("Preprocessing ...")
  facets.preproc <- facets::preProcSample(rcmat = facets.snpm, ndepth = ndepth, hetscale = hetscale, gbuild = genome.build, unmatched = unmatched)
  
  ## Segment
  message("Segmenting ...")
  facets.proc <- facets::procSample(x = facets.preproc, cval = cval)
  rootname <- paste0(samplename, "_c", cval)
  saveRDS(facets.proc, file = paste0(rootname, "_FACETS_segments.RDS"), compress = "xz")
  
  ## Estimate CN and cellularity from EM
  message("Evaluating CN, ploidy and cellularity ...")
  facets.emcncf <- facets::emcncf(x = facets.proc)
  saveRDS(facets.emcncf, file = paste0(rootname, "_FACETS_CNestimates.RDS"), compress = "xz")

  ## PLOT
  message("Plotting ...")
  png(paste0(rootname, "_FACETS_plot.png"), 1500,950)
  plotSample(x = facets.proc, emfit = facets.emcncf, sname = samplename)
  dev.off()
}

## Run facets.build.snpp() from a targetfile
facets.build.snpp.tgf <- function(targetfile = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  if (is.null(targetfile)) stop("targetfile is required !")
  if (!file.exists(targetfile)) stop("targetfile not found !")
  tgf <- read.table(targetfile, header = TRUE, sep = "\t", as.is = TRUE)
  
  require(foreach)
  
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  foreach(s = seq_len(nrow(tgf)), .export = "facets.build.snpp") %dopar% {
    facets.build.snpp(normal.bam = tgf$normal.bam[s], tumor.bam = tgf$tumor.bam[s], samplename = tgf$samplename[s], ...)
  }
  parallel::stopCluster(cl)
}

## Run facets.run() in multithreading from a list of snpp.gz files
facets.run.mt <- function(snpp.list = list.files(pattern = "\\.snpp.*$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE), nthread = 1, cluster.type = "PSOCK", ...) {
  
  snames <- sub(".snpp.*", "", basename(snpp.list), ignore.case = TRUE)
  
  require(foreach)
  
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  foreach(s = seq_along(snpp.list), .export = "facets.run") %dopar% {
    facets.run(snpp.file = snpp.list[s], samplename = snames[s], ...)
  }
  parallel::stopCluster(cl)
}
###

# setwd("")
# normal.BAM <- "/home/job/WORKSPACE/sequenza_tests/TCGA-AC-A2BK-11A-13D-A21Q-09_sorted.bam"
# tumor.BAM <- "/home/job/WORKSPACE/sequenza_tests/TCGA-AC-A2BK-01A-11D-A21Q-09_sorted.bam"
# samplename <- "TCGA-AC-A2BK-01A_vs_11A"
# vcf.file <- "/home/job/WORKSPACE/facets_tests/dbsnp_144.hs37d5.vcf.gz"
# 
# snpp.list <- list.files(pattern = "\\.snpp\\.gz$", full.names = TRUE, ignore.case = TRUE, recursive = FALSE)
# snames <- sub(".snpp.*", "", basename(snpp.list))
# for (s in seq_along(snpp.list)) {
#   facets.run(snpp.file = snpp.list[s], samplename = snames[s])
# }
# 
# 
# 
