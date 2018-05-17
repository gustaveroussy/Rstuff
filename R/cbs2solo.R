cbs2solo <- function(cbs.file = NULL, genome.build = "hg19") {
  if (is.null(cbs.file)) stop("Please provide a CBS file.")
  if(!file.exists(cbs.file)) stop("Please provide a valid CBS file.")
  my.cbs <- maelstrom::read.table.fast(file = cbs.file, header = TRUE)
  data(list = genome.build, package = "chromosomes", envir = environment())
  solo.df <- data.frame(loc = paste0(unlist(cs$chr2chrom[my.cbs$Chr]), ":", my.cbs$Start, "-", my.cbs$End), probes = my.cbs$Probes, status = "N", l2r = my.cbs$Log2Ratio, ratio = 2^my.cbs$Log2Ratio, stringsAsFactors = FALSE)
  solo.df$status[solo.df$l2r>0] <- "G"
  solo.df$status[solo.df$l2r<0] <- "L"
  solo.df <- solo.df[!solo.df$status == "N",]
  write.table(solo.df, file = sub(pattern = ".cbs$", replacement = ".solo", x = cbs.file, ignore.case = TRUE), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

my.lf <- list.files(pattern = "nrf.*cbs$", full.names = TRUE, recursive = TRUE)

nthread <- 7
require(foreach)
require(doParallel)
cl <- makeCluster(spec = nthread, type = "FORK")
registerDoParallel(cl)
foreach (mlf = my.lf, .inorder = FALSE) %dopar% {
  oridir <- getwd()
  setwd(dirname(mlf))
  cbs2solo(cbs.file = basename(mlf))
  system(paste0("grd -m solo ", sub(pattern = ".cbs$", replacement = ".solo", x = basename(mlf), ignore.case = TRUE)))
  setwd(oridir)
}
stopCluster(cl)