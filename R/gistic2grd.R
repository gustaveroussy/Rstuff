## gistic2grd
##
## DESCRIPTION : A wrapper to reformat GISTIC2 results and annotate it with grd
##
## AUTHOR:
## Bastien JOB (bastien.job@gustaveroussy.fr)
##
## DEPENDS ON:
## . some_useful_functions.R (svn)
## . matrixStats
## . foreach
## . doParallel
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 3.3.0
##
## VERSION NOTES:
##
## v1.0c 20161018
##      . Replaced support for the script some_useful_functions.R by its packaged form 'maelstrom"
##
## v1.0b 20160923
##      . Added support for LDB.
## v1.0  20160531
##      . This script is a R rewrite to replace the Perl5 script 'gistic2grd'.
##      . Some fix has been added to remove duplicated genomic regions in the
##        'regions' type (appeared in GISTIC versions above 2.020).
##      . Multithreading is now handled in R using foreach/doParallel.

help.gistic2grd <- function() {
  cat("
gistic2grd(gisticfile = NULL, sp = \"hs\", gb = 19, ldb = \"/mnt/data_cigogne/bioinfo/\", ncores = 3)

gisticfile    :  (character)  Path+name of the \"all_lesions.conf_XX.txt\" results file from a GISTIC2 run [NULL]
sp            :  (character)  Species (hs, mm, rn). This is required for grd [\"sp\"]
gb            :  (integer)    Genome build (17,18,19,38 for hs. 9, 10 for mm. 4, 5, 6 for rn). This is required for grd [19]
ldb           :  (character)  Path to a local data base [\"/mnt/data_cigogne/bioinfo/\"]
ncores        :  (integer)    Number of threads. Values above 3 will be capped to 3 as only 3 files will be annotated with grd [3]
")
}

gistic2grd <- function(gisticfile = NULL, sp = "hs", gb = 19, ldb = "/mnt/data_cigogne/bioinfo/", ncores = 3) {
  if (is.null(gisticfile)) stop("Please provide a GISTIC2 output file 'all_lesions.conf_NN.txt' !")
  if (!file.exists(gisticfile)) stop ("Could not find the input GISTIC2 file (should be 'all_lesions.conf_NN.txt')")
  # if (!file.exists(useful.script)) stop("some_useful_functions.R is required [more precisely, the chrConv() function] !")
  # source(useful.script)
  require(maelstrom)
  if (ncores > 3) ncores <- 3
  oridir <- getwd()
  wdir <- dirname(gisticfile)
  setwd(wdir)
  mygres <- read.table(gisticfile, sep="\t", header = TRUE, as.is = TRUE, check.names = FALSE, strip.white = TRUE)
  mygres <- mygres[1:(nrow(mygres)/2),-c(ncol(mygres))]
  nsamples <- length(10:ncol(mygres))
  require(matrixStats)
  freq <- 1- (rowCounts(x = as.matrix(mygres[,10:ncol(mygres)]), value = 0, na.rm = TRUE) / nsamples)
  mygres[["Unique Name"]][grep(pattern = "^Amplification", x = mygres[["Unique Name"]])] <- "G"
  mygres[["Unique Name"]][grep(pattern = "^Deletion", x = mygres[["Unique Name"]])] <- "L"
  region_sizes <- list(peak = "Peak Limits", wide = "Wide Peak Limits", region = "Region Limits")
  
  require(foreach)
  require(doParallel)
  
  cl <- makeCluster(spec = ncores, type = "FORK")
  registerDoParallel(cl)
  
  foreach (x = names(region_sizes)) %dopar% {
    
    locword <- gsub(pattern = "[()]", replacement = "", x = mygres[[region_sizes[[x]]]])
    loc_probes <- strsplit(x = locword, split = "(probes |-|:)")
    out.df <- as.data.frame(matrix(unlist(loc_probes), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
    colnames(out.df) <- c("Chrom", "Start", "End", "Probe.start", "Probe.end")
    for(y in c(2:5)) out.df[,y] <- as.numeric(out.df[,y])
    out.df$Loc <- paste0(out.df$Chrom, ":", out.df$Start, "-", out.df$End)
    out.df$Chr <- chrConv(chrvec = out.df$Chrom, chr.in = "chrom", chr.out = "chr")
    out.df$Probes <- out.df$Probe.end - out.df$Probe.start +1
    out.df$Status <- mygres[["Unique Name"]]
    out.df$Freq <- freq
    # out.df[["Q-value"]] <- mygres[["Residual q values after removing segments shared with higher peaks"]]
    out.df[["Q-value"]] <- mygres[["q values"]]
    out.df <- out.df[order(out.df$Chr, out.df$Start, out.df$End),]
    checkword <- paste0(out.df$Status, out.df$Loc)
    out.df <- out.df[!duplicated(checkword),]
    out.filename <- paste0(sub(pattern = ".txt$", replacement = "", x = basename(gisticfile)), "_", x, ".g2")
    write.table(out.df[,-c(1:5,7)], file = out.filename, sep = "\t", quote = FALSE, row.names = FALSE)
    cmdl <- paste0("grd --sp ", sp, " --gb ", gb, " -m g2 --ldb ", ldb, " ", out.filename)
    # print(cmdl)
    system(command = cmdl, ignore.stdout = TRUE)
    # system(command = cmdl, ignore.stdout = FALSE)
  }
  stopCluster(cl)
  setwd(oridir)
}

