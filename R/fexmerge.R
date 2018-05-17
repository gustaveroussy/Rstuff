## fexmerge.R
##
## DESCRIPTION :  This script performs the merging of signal intensity tracks from different
## FEX files and generate an artifical FEX containing the merged signals as an output. This
## allows the use of intensity tracks of different arrays to generate "in silico" hybridizations.
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
## Clément MAZOYER
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.12.2
##
## DEPENDS ON:
## snowfall (if ncores>1)
##
## VERSION NOTES
##
## v1 20120403
##      . First release.
##

fexmerge <- function(inputfile, ncores=1) {
  inlist <- read.table(inputfile, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  gv <- ""
  if (sp == "hs") {
    gv = "hg"
  } else if (sp == "mm") {
    gv = "mm"
  }
  gv <- paste(gv, gb, sep="")
  
  ## Replacing dots
  inlist[,3] <- chartr(".", "_", inlist[,3])
  inlist[,4] <- chartr(".", "_", inlist[,4])
  
  ## Running multithreaded...
  if (ncores > 1) {
    library(snowfall)
    sfStop()
    sfInit(parallel=T,cpus=ncores)
    sfExportAll()
    outdf <- sfApply(inlist, 1, function(x) {
      #       x <- as.vector(inlist[z,])
      p2c.core(gcx1=x[1], gcx2=x[2], sample1=x[3], sample2=x[4], sp=sp, gb=gb, gv=gv, nrf=nrf, uSD=uSD, noX=noX, cor.cut=cor.cut, cor.method=cor.method, mnt.proj=mnt.proj, in.type=in.type)
    })
    sfStop()
  }
  ## Or normaly...
  else {
    outdf <- apply(inlist, 1, function(x) { p2c.core(gcx1=x[1], gcx2=x[2], sample1=x[3], sample2=x[4], sp=sp, gb=gb, gv=gv, nrf=nrf, uSD=uSD, noX=noX, cor.cut=cor.cut, cor.method=cor.method, mnt.proj=mnt.proj, in.type=in.type) })
  }
  
  outdf <- as.data.frame(matrix(unlist(outdf, recursive=T), ncol=14, byrow=T), stringsAsFactors=F)
  colnames(outdf) <- c("Sample1", "Sample2", "undo.SD", "Nrf", "noX", "Cor.cutoff", "Cor.method", "Smoothed.cor", "Noise.MAD", "S1.more.S2", "S1.less.S2", "S1.diff.S2", "Fit.Slope", "Fit.Intercept")
  write.table(outdf, paste(outroot, "_", gv, "_U", uSD, "_nrf", nrf, "_CC", cor.cut, ".txt", sep=""), sep="\t", quote=F, row.names=F)
}