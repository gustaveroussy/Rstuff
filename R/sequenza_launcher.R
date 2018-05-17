# setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/SEQUENZA")
# setwd("/home/job/WORKSPACE/sequenza_tests")

### PLEASE BE SURE THAT modules bioinfo/shared/python/3.6.1 AND bioinfo/shared/sequenza-utils/2.1.9999b0-py3.6 HAVE BEEN PROPERLY LOADED IF RUNNING ON CIGOGNE !

# require(sequenza)

# gc.file <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/SEQUENZA/hs37d5.gc50.gz"
# fasta.file <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/SEQUENZA/hs37d5.fasta"
# normal.BAM <- "/home/job/WORKSPACE/sequenza_tests/TCGA-AC-A2BK-11A-13D-A21Q-09_sorted.bam"
# tumor.BAM <- "/home/job/WORKSPACE/sequenza_tests/TCGA-AC-A2BK-01A-11D-A21Q-09_sorted.bam"
# samplename <- "TCGA-AC-A2BK-01A_vs_11A"

## Build a GC pack
sequenza.gc.compute <- function(sequenza.utils.path = NULL, fasta.file = NULL, window.size = 50) {
  
  require(sequenza, quietly = TRUE)
  if (is.null(fasta.file)) stop("A reference genome FASTA file is required !")
  if (!file.exists(fasta.file)) stop("Reference genome FASTA file not found !")
  if (!is.null(sequenza.utils.path)) sequenza.utils.path <- paste0(sequenza.utils.path, "/")
  
  system(paste0("sequenza-utils",
                " gc_wiggle",
                " -f ", fasta.file,
                " -w ", window.size,
                " -o ", fasta.file, ".GC.w", window.size, ".wig"
                ))
}

## Build seqz file and bin it (like would say Mickael)
sequenza.build.bin.seqz <- function(sequenza.utils.path = NULL, normal.bam = NULL, tumor.bam = NULL, samplename =  NULL, gc.file = NULL, fasta.file = NULL, min.depth = 20, min.baseQ = 20, bin.size = 50) {
  
  require(sequenza, quietly = TRUE)
  if (is.null(normal.bam)) stop("Normal BAM file is required !")
  if (!file.exists(normal.bam)) stop("Normal BAM not found !")
  if (is.null(tumor.bam)) stop("Tumor BAM file is required !")
  if (!file.exists(tumor.bam)) stop("Tumor BAM not found !")
  if (is.null(gc.file)) stop("A GC file is required !")
  if (!file.exists(gc.file)) stop("GC file not found !")
  if (is.null(fasta.file)) stop("A reference genome FASTA file is required !")
  if (!file.exists(fasta.file)) stop("Reference genome FASTA file not found !")
  if (is.null(samplename)) stop("A sample name is required !")
  if (!is.null(sequenza.utils.path)) sequenza.utils.path <- paste0(sequenza.utils.path, "/")
  
  system(paste0("sequenza-utils",
                " bam2seqz",
                " -gc ", gc.file,
                " -F ", fasta.file,
                " -n ", normal.bam,
                " -t ", tumor.bam,
                " -o ", samplename, ".seqz.gz",
                " -q ", min.baseQ,
                " -N ", min.depth
                ))
  
  system(paste0("sequenza-utils",
                " seqz_binning",
                " -s ", samplename, ".seqz.gz",
                " -w ", bin.size,
                " -o ", samplename, ".bin", bin.size, ".seqz.gz"
                ))
}

## Run sequenza.build.seqz() from a targetfile
sequenza.build.bin.seqz.tgf <- function(targetfile = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  if (is.null(targetfile)) stop("targetfile is required !")
  if (!file.exists(targetfile)) stop("targetfile not found !")
  tgf <- read.table(targetfile, header = TRUE, sep = "\t", as.is = TRUE)
  
  require(foreach, quietly = TRUE)
  
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  foreach(s = seq_len(nrow(tgf)), .export = "sequenza.build.bin.seqz") %dopar% {
    sequenza.build.bin.seqz(normal.bam = tgf$normal.bam[s], tumor.bam = tgf$tumor.bam[s], samplename = tgf$samplename[s], ...)
  }
  parallel::stopCluster(cl)
}

## Run sequenza with default parameters
sequenza.run <- function(seqz.bin.file = NULL, chromosome.list = c(as.character(1:22), "X", "Y"), assembly = "hg19", samplename = NULL) {
  require(sequenza, quietly = TRUE, warn.conflicts = FALSE)
  sname <- sub("\\.seqz.*", "", basename(seqz.bin.file))
  seqz.ext <- sequenza.extract(file = seqz.bin.file, chromosome.list = chromosome.list, assembly = assembly)
  seqz.fit <- sequenza.fit(sequenza.extract = seqz.ext, chromosome.list = chromosome.list)
  sequenza.results(sequenza.extract = seqz.ext, cp.table = seqz.fit, sample.id = samplename, chromosome.list = chromosome.list)
}

## Run sequenza.run() in multithreading from a list of snpp.gz files
sequenza.run.mt <- function(seqz.bin.list = list.files(pattern = "\\.bin.*\\.seqz\\.gz$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE), nthread = 1, cluster.type = "PSOCK", ...) {
  
  snames <- sub("\\.seqz.*", "", basename(seqz.bin.list), ignore.case = TRUE)
  
  # require(foreach, quietly = TRUE)
  
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  foreach(s = seq_along(seqz.bin.list), .export = "sequenza.run") %dopar% {
    sequenza.run(seqz.bin.file = seqz.bin.list[s], samplename = snames[s], ...)
  }
  parallel::stopCluster(cl)
}
#### 

# ## Loading the seqz

