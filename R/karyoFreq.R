## karyoFreq.R pour "Karyotype-like frequencies plotter in R"
##
## DESCRIPTION:
##  This script generates a PNG image representing the frequencies of aberrations across a population
##  of profiles taken from a PROB data file, in a per-chromosome, karyotypical view.
##
## AUTHOR:
## Bastien JOB (bastien.job@gustaveroussy.fr)
##
## DEPENDS ON:
## . r-base
## . chrload.R
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.12.2
##
## VERSION NOTES
## v1.1 20140313
##    . Removed internal chrload function, which must now be sourced from outer, up to date, chrload.R script.
## v1.0b 20140303
##      . Added virtual support for rattus norvegicus (untested, though!).
## v1.0 20120215
##      . First version, built upon request of Stephane VIGNOT for projet P21_SV_CGH.




cat("\nWARNING !\n*********\nScript chrload.R from svn (genomics/CGH/R/chrload.R) must be sourced before execution.\n\n")

karyo.freq.plot <- function(probname, sp="hs", gb=19, gcol="blue", acol="cyan", lcol="orangered4", dcol="orange", alim=1.5, dlim=-1.5, mnt.proj="/mnt/proj/") {
  
  ## Loading chr data
  cs <<- chrload(sp, gb, mnt.proj)
#   mnt.proj <<- mnt.proj
  
  ## Loading PROB data
  prob <- read.table(probname, header=T, sep="\t", stringsAsFactors=F)
  l2r <- prob[,-c(1:6)]
  pos <- prob[,c(1:6)]
  rm(prob)
  mat <- t(matrix(as.numeric(as.matrix(l2r)), nrow=nrow(l2r), ncol=ncol(l2r)))
  gc()
  
  ## Computing frequences
  gainv <- apply(mat, 2, function(z) { length(which(z > 0))/nrow(mat) })
  lossv <- apply(mat, 2, function(z) { length(which(z < 0))/nrow(mat) })
  ampv <- apply(mat, 2, function(z) { length(which(z >= alim))/nrow(mat) })
  delv <- apply(mat, 2, function(z) { length(which(z <= dlim))/nrow(mat) })
  
  ## Plotting
  oripar <- par()
  
  if (sp == "hs") {
    zone <- vector()
    z <- 1
    while(z <= 72) {
      zone <- c(zone, z+1, z, z+2)
      z <- z+3
    }
    zone <- matrix(zone, nrow=2, byrow=T)
    layW = rep(c(.031667/2, 0.01, .031667/2), 12)
  } else if ( (sp == "mm") || (sp == "rn") ) {
    zone <- matrix(seq(1,22*3,1), nrow=2, ncol=22*1.5, byrow=T)
    layW = rep(c(.034545/2, 0.01, .034545/2), 11)
  } else (stop( "Unknwon species ! Are supporte : hs, mm, rn."))
  
  layH = c(0.5,0.5)
  chrnames <- unique(cs$cytob$chrA)
  png(paste(probname,".kfreq.png", sep=""), width=1680, height=1050)
  par(mgp=c(0,0,0), mar=c(1,0,1,0), omi=c(0,0.2,0.1,0), xaxt="n", yaxt="n", bty="n")
  layout(zone, widths=layW, heights = layH)
  cat("Plotting...\n")
  for (k in unique(pos$Chr)) {
    plot(0, 0, xlim=c(0,1), ylim=c(-cs$lchrxx[1], 0), type="n", xlab="")
    cytochr <- which(cs$cytob$chr == k)
  	rect(cs$cytob$x1[cytochr], -cs$cytob$chromStart[cytochr], cs$cytob$x2[cytochr], -cs$cytob$chromEnd[cytochr], col=cs$cytob$gieStain[cytochr])
    text(0.5, -cs$lchrxx[k]-cs$lchrxx[1]*0.05, chrnames[k], pos=3, cex=1.5)
    kp <- which(pos$Chr == k)
    plot(0,0, type="n", ylim=c(-cs$lchrxx[1], 0), xlim=c(-1,0), xlab="")
    segments(0, -pos$Start[kp], -lossv[kp], -pos$Start[kp], col="orangered4")
    segments(0, -pos$Start[kp], -delv[kp], -pos$Start[kp], col="orange")
    segments(c(0,-0.5), -min(pos$Start[kp]), c(0,-0.5), -max(pos$Start[kp]), lty=2)
    plot(0,0, type="n", ylim=c(-cs$lchrxx[1], 0), xlim=c(0,1), xlab="")
    segments(0, -pos$Start[kp], gainv[kp], -pos$Start[kp], col="blue")
    segments(0, -pos$Start[kp], ampv[kp], -pos$Start[kp], col="cyan")
    segments(c(0,0.5), -min(pos$Start[kp]), c(0,0.5), -max(pos$Start[kp]), lty=2)
  }
  dev.off()
}

# ## IMPORT CHR DATA
# chrload <- function(sp, gb, mnt.proj=mnt.proj) {
#   
#   sptxt <- sp
#   if (sp == "hs") sptxt <- "hg"
#   gv <- paste(sptxt, gb, sep="")
#   
#   ## IMPORT BANDES CYTO ET PREPARATION POUR PLOT KARYO
#   cat("Importing chromosomes data  ...\n")
#   if (gv == "hg18") cytob <- read.table(paste(mnt.proj, "cgh/hg18/cytoBandIdeo.hg18", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
# 	else if (gv == "hg19") cytob <- read.table(paste(mnt.proj, "cgh/hg19/cytoBandIdeo.hg19", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
# 	else if (gv == "mm9") cytob <- read.table(paste(mnt.proj, "cgh/mm9/cytoBandIdeo.mm9", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
# 	
# 	cytob$chrA <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
#   cytob$chr <- cytob$chrA
#   if (sp == "hs") {
#     cytob$chr[which(cytob$chr == "X")] <- 23
#     cytob$chr[which(cytob$chr == "Y")] <- 24
#   } else if (sp == "mm") {
#     cytob$chr[which(cytob$chr == "X")] <- 20
#     cytob$chr[which(cytob$chr == "Y")] <- 21
#   }
#   
# 	cytob$chr <- as.numeric(cytob$chr)
# 	cytob <- cytob[order(cytob$chr),]
# 
# 	lchrxx <- sapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) })
# 	lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
# 	glen <- sum(as.numeric(lchrxx))
# 	lchrsum <- sapply(c(1:length(lchrtoadd)), function(x) { sum(lchrtoadd[1:x]) })
# 
# 	cytob$x1 <- 0.15
# 	cytob$x2 <- 0.85
# 	cytob$x1[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.25
# 	cytob$x2[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.75
# 	cytob$gieStain[which(cytob$gieStain == "gneg")] <- "white"
#   cytob$gieStain[which(cytob$gieStain == "gpos25")] <- "grey75"
#   cytob$gieStain[which(cytob$gieStain == "gpos33")] <- "grey66"
#   cytob$gieStain[which(cytob$gieStain == "gpos50")] <- "grey50"
#   cytob$gieStain[which(cytob$gieStain == "gpos66")] <- "grey33"
#   cytob$gieStain[which(cytob$gieStain == "gpos75")] <- "grey25"
#   cytob$gieStain[which(cytob$gieStain == "gpos100")] <- "black"
#   cytob$gieStain[which(cytob$gieStain == "acen")] <- "yellow"
#   cytob$gieStain[which(cytob$gieStain == "gvar")] <- "darkred"
#   cytob$gieStain[which(cytob$gieStain == "stalk")] <- "darkolivegreen"
#   
#   return(list(cytob=cytob, lchrxx=lchrxx, lchrsum=lchrsum, lchrtoadd=lchrtoadd, glen=glen, cur.glen=glen))
# }
