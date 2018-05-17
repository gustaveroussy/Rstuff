## DESCRIPTION
##
## This scripts aims to find genes/probes for which aberration levels and expression are
## correlated, for the same samples in the two datasets. Currently It doesn't allow to
## find correlated probes with a differential between N subpopulations.
##
## VERSION NOTES
##
##
## v2.2 20130405
##		. Change in cheese plot : change the method to calculate the radii of circles in order
##		  to compare the size of gain and loss circles
##
## v2.1 20130320
##		. Change annotation from GEX data. Instead of hgug4112a, used the annotation from earray
##
## v2.0b 20111116
##      . Introduced a "help" function to display the description of parameters for the main function, even when
##        this script is sourced and not directly read.
##
## v2.0  20110315   . Script remade :
##                    * Wraped the script into sourcable functions.
##                    * Simplified, optimized code.
##                  . More readable plots.
##                  . Support for REG/XREG formats added.
##                  . Support for ABRS3 removed.
##                  . Added the mnt.proj option for custom /proj mountpoint.
##
## v1.1  20091016   . Added ANOVA test for N-class differentiation.
##
## v1.0  20090528	. First release



###############
## FUNCTIONS ##
###############

getChr <- function(chr){
	loc=list()
	loc$chr=NULL
	loc$start=NULL
	loc$end=NULL
	
	for(i in 1:length(chr)){
		tmp = unlist(strsplit(strsplit(chr[i],"chr")[[1]][2],":"))
		tmp2 = unlist(strsplit(tmp[2],"-"))
	
		loc$chr = c(loc$chr, tmp[1])
		loc$start = c(loc$start,as.numeric(tmp2[1]))
		loc$end = c(loc$end,as.numeric(tmp2[2]))
	}
		
	return(loc)
}

## Loading chromosomal informations
chrload <- function(hg=19, mnt.proj="/mnt/proj/") {
  if (hg == 17) cytob <- read.table(paste(mnt.proj, "cgh/hg17/cytoBandIdeo.hg17", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
  if (hg == 18) cytob <- read.table(paste(mnt.proj, "cgh/hg18/cytoBandIdeo.hg18", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
  if (hg == 19) cytob <- read.table(paste(mnt.proj, "cgh/hg19/cytoBandIdeo.hg19", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
  
  cytob$chr <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
  cytob$chr[which(cytob$chr == "X")] <- 23
  cytob$chr[which(cytob$chr == "Y")] <- 24
  cytob$chr <- as.numeric(cytob$chr)
  cytob <- cytob[order(cytob$chr),]
  
  lchrxx <- sapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) })
  lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
  glen <- sum(as.numeric(lchrxx))
  lchrsum <- sapply(c(1:length(lchrtoadd)), function(x) { sum(lchrtoadd[1:x]) })
  return(list(chr.size=lchrxx, chr.toadd=lchrtoadd, chr.additive=lchrsum, glen=glen))
}
help.cghgexcor.basic <- function() {
  cat("
PARAMETERS FOR THE cghgexcor.basic() FUNCTION :
----------------------------------------
## OPTIONS :
##  CGHd            A dataframe for CGH values, directly obtained from a REG or XREG table.
##  GEd             A dataframe for GEX values. These must be normalized values, from a public Agilent design and samples must be in the same order in GEd and CGHd. First column has to be \"ProbeName\".
##  design          The \"hgug\" Agilent design name in R. By default : hgug4112a corresponding to Agilent Human 44K arrays.
##  outname         The root for output filenames.
##  isGElog         (boolean) As GEX data have to be in log for the computation of correlations, the script has to know if they were given as logN or not. If not, they will be converted to log2.
##  cor.method      Either \"pearson\" or \"spearman\".
##  cor.min.value   Minimum number of values to consider computing the correlation for each probe. If there are less values than this number, the probe will be discarded.
##  p.min           BH-adjusted minimal p-value cut-off for the correlation test.
##  mnt.proj        Pelican's /proj mountpoint.
##  hg              The human genome build to use. Either 17, 18 or (default) 19.
##  symbol          (boolean) Choose to display HGNC Hugo gene symbols on the genomic plots or not.
\n")
}

## MAIN FUNCTION
cghgexcor.basic <- function(CGHd, GEd, geAnnot, outname="toto", isGElog=T, cor.method="pearson", cor.min.value=3, p.min=0.05, hg=19, mnt.proj="/mnt/proj/", symbol=T) {
  
  # Dealing with chromosomes
  cat("Loading chr data ...\n")
  chrdata <- chrload(hg=hg, mnt.proj=mnt.proj)
  
  # Adding annotations
  cat("Getting probes' annotations ...\n")
  GEprobloc <- as.data.frame(cbind(ProbeName=as.vector(GEd$ProbeName)), stringsAsFactors=F)
  GEprobloc$Symbol <- geAnnot[as.vector(GEprobloc$ProbeName),"GeneName"]
  ChrLoc = getChr(geAnnot[as.vector(GEprobloc$ProbeName),"ChromosomalLocation"])
  GEprobloc$Chr <- ChrLoc$chr
  GEprobloc$Chr[which(GEprobloc$Chr == "X")] <- 23
  GEprobloc$Chr[which(GEprobloc$Chr == "Y")] <- 24
  GEprobloc$Chr <- as.numeric(GEprobloc$Chr)
  GEprobloc$Start <- ChrLoc$start
  GEprobloc$End <- ChrLoc$end
#  GEprobloc$CytoBand <- annotation(GEprobloc$ProbeName, coll=design, ressource="MAP", multi=F, na.action=NA)
  GEprobloc$RefSeq <-geAnnot[as.vector(GEprobloc$ProbeName),"Name"]
  
  # Removing redundant probes
  GEprobloc <- GEprobloc[!duplicated(GEprobloc$ProbeName),]
  GEd <- GEd[!duplicated(GEd$ProbeName),]
  
  # Intersecting
  intersec.probes <- intersect(GEd$ProbeName, GEprobloc$ProbeName)
  located.GEd <- GEd[is.element(GEd$ProbeName, intersec.probes),]
  exp.GEprobloc <- GEprobloc[is.element(GEprobloc$ProbeName, intersec.probes),]
    
  # Sort by ProbeName before fuse
  located.GEd <- located.GEd[order(located.GEd$ProbeName),]
  exp.GEprobloc <- exp.GEprobloc[order(exp.GEprobloc$ProbeName),]
  
  GE_loc <- exp.GEprobloc
  GE_log2.main <- located.GEd[,-c(1)]
  
  # Converting GE data to log2 if needed
  if (isGElog == F) {
    cat("Converting GE values to log2 ...\n")
    GE_log2.main <- log2(GE_log2.main)
  }
  
  # Cleaning
  rm(GEd, GEprobloc, located.GEd, exp.GEprobloc)
  
  # Preparing CGH objects : localizations and values
  CGH_loc <- CGHd[,c(2:4)]
  CGH_log2.main <- CGHd[,-c(1:9)]
  
  # Cleaning
  rm(CGHd)
  gc()
	
  # No annotation in this function : just for direct correlation of datasets
  GE_log2 <- GE_log2.main
  CGH_log2 <- CGH_log2.main
  
  # Performing the genomic mapping
  cat("Mapping CGH values to GE features ...\n")
  bigmat <- array(dim=c(nrow(GE_log2), 2*ncol(GE_log2)))    ## New version : creating a big matrix with CGH then GEX data stuck by columns
  CGHpart <- 1:ncol(CGH_log2)
  GEpart <- (ncol(GE_log2)+1):(2*ncol(GE_log2))
  for (k in 1:nrow(GE_loc)) {
    v <- GE_loc[k,]
    CGH_index <- which( (v$Chr == CGH_loc$Chr) & (v$Start >= CGH_loc$Start) & (v$End <= CGH_loc$End) )
    if (length(CGH_index) == 1) {
      bigmat[k,CGHpart] <- as.numeric(CGH_log2[CGH_index,])
      #bigmat[k,GEpart] <- as.numeric(as.vector(GE_log2[k,]))
	  bigmat[k,GEpart] <- apply(GE_log2[k,],2,as.numeric)
    }
    else if (length(CGH_index) > 1) {
      print("MULTIMATCH !")
    }
    # Returning a line counter to the stdev out
    if ((k %% 1000) == 0) {
      print(paste(" ", k, " probes mapped.", sep=""))	# Etat d'avancement sur sortie ?cran
    }
  }
  # Computing correlations
  correlation.probes <- GE_loc
  correlation.probes$cor.coef <- apply(bigmat, 1, function(x) { if (length(which(!is.na(x[GEpart]))) > cor.min.value) cor.test(x[CGHpart][which(!is.na(x[GEpart]))], x[GEpart][which(!is.na(x[GEpart]))], method=cor.method)$estimate else NA } )
  correlation.probes$cor.rawp <- apply(bigmat, 1, function(x) { if (length(which(!is.na(x[GEpart]))) > cor.min.value) cor.test(x[CGHpart][which(!is.na(x[GEpart]))], x[GEpart][which(!is.na(x[GEpart]))], method=cor.method)$p.value else NA } )
  
  non.NA.keep <- which(!is.na(correlation.probes$cor.coef))
  correlation.probes <- correlation.probes[non.NA.keep,]
  bigmat <- bigmat[non.NA.keep,]
  correlation.probes$cor.adjp <- p.adjust(correlation.probes$cor.rawp, method="BH")
  
  # Sort on genomic position (useful for barplot order)
  genomic.order <- order(correlation.probes$Chr, correlation.probes$Start, correlation.probes$End)
  correlation.probes <- correlation.probes[genomic.order,]
  bigmat <- bigmat[genomic.order,]
  
  # Dump values table
  write.table(correlation.probes, paste(outname, "_CGH-GEXcor_", cor.method, "_cmp", cor.min.value, "_hg", hg, ".txt", sep=""), sep="\t", quote=F, row.names=F)
  
  # Getting index of positive/negative correlations, and their significant subpart at a defined p-value
  pos.cor <- which(correlation.probes$cor.coef > 0)
  sig.pos.cor <- pos.cor[which(correlation.probes$cor.adjp[pos.cor] < p.min)]
  min.sig.pcoef <- correlation.probes$cor.coef[sig.pos.cor][which.max(correlation.probes$cor.adjp[sig.pos.cor])]
  neg.cor <- which(correlation.probes$cor.coef < 0)
  sig.neg.cor <- neg.cor[which(correlation.probes$cor.adjp[neg.cor] < p.min)]
  min.sig.ncoef <- correlation.probes$cor.coef[sig.neg.cor][which.max(correlation.probes$cor.adjp[sig.neg.cor])]
  
  # Pangenomic plot
  grnmd <- runmed(correlation.probes$cor.coef, nrow(correlation.probes)/100)
  pdf(file=paste(outname, "_CGH-GEXcor_", cor.method, "_p", p.min, "_cmp", cor.min.value, "_hg", hg, "_pangenomicplot", ".pdf", sep=""), height=21/cm(1),width=29.7/cm(1))
  plot(correlation.probes$Start[pos.cor]+chrdata$chr.additive[correlation.probes$Chr[pos.cor]], correlation.probes$cor.coef[pos.cor], xlim=c(0, chrdata$glen), ylim=c(-1, 1), type="h", col="grey75", xaxs="i", xlab="Genomic position", ylab="Correlation coefficient", main=paste("CGH - GEX correlation (", cor.method, ")\n", length(sig.pos.cor), " positive correlation(s) / ",length(sig.neg.cor), " negative correlation(s) at adj.p<", sprintf("%1.1E", p.min), " out of ", nrow(correlation.probes), " features", sep=""))
  if (length(sig.pos.cor) > 0) {
    abline(h=min.sig.pcoef, lty=2, col=4)
    text(0, min.sig.pcoef-0.05, sprintf("%1.2f", min.sig.pcoef), pos=4, col=4)
    segments(correlation.probes$Start[sig.pos.cor]+chrdata$chr.additive[correlation.probes$Chr[sig.pos.cor]], 0, correlation.probes$Start[sig.pos.cor]+chrdata$chr.additive[correlation.probes$Chr[sig.pos.cor]], correlation.probes$cor.coef[sig.pos.cor], col=4)
    if (symbol) text(correlation.probes$Start[sig.pos.cor]+chrdata$chr.additive[correlation.probes$Chr[sig.pos.cor]], correlation.probes$cor.coef[sig.pos.cor], labels=correlation.probes$Symbol[sig.pos.cor], col=4, pos=3)
  }
  if (length(sig.neg.cor) > 0) {
    lines(correlation.probes$Start[neg.cor]+chrdata$chr.additive[correlation.probes$Chr[neg.cor]], correlation.probes$cor.coef[neg.cor], type="h", col="grey25")
    abline(h=min.sig.ncoef, lty=2, col=2)
    text(0, min.sig.ncoef+0.05, sprintf("%1.2f", min.sig.ncoef), pos=4, col=2)
    segments(correlation.probes$Start[sig.neg.cor]+chrdata$chr.additive[correlation.probes$Chr[sig.neg.cor]], 0, correlation.probes$Start[sig.neg.cor]+chrdata$chr.additive[correlation.probes$Chr[sig.neg.cor]], correlation.probes$cor.coef[sig.neg.cor], col=2)
    if (symbol) text(correlation.probes$Start[sig.neg.cor]+chrdata$chr.additive[correlation.probes$Chr[sig.neg.cor]], correlation.probes$cor.coef[sig.neg.cor], labels=correlation.probes$Symbol[sig.neg.cor], col=2, pos=1)
  }
  lines(correlation.probes$Start+chrdata$chr.additive[correlation.probes$Chr], grnmd, col="darkgreen")
  abline(v=chrdata$chr.additive, lty=2, col=6)
  abline(h=0, lty=2)
  for (k in sort(unique(correlation.probes$Chr))) { text(mean(c(chrdata$chr.additive[k],chrdata$chr.additive[k+1])), (k%%2)*2-1, k) }
  
  for (k in sort(unique(correlation.probes$Chr))) {
    kfeatures <- which(correlation.probes$Chr == k)
    cpk <- correlation.probes[kfeatures,]
    pos.cor.k <- which(cpk$cor.coef > 0)
    sig.pos.cor.k <- pos.cor.k[which(cpk$cor.adjp[pos.cor.k] < p.min)]
    neg.cor.k <- which(cpk$cor.coef < 0)
    sig.neg.cor.k <- neg.cor.k[which(cpk$cor.adjp[neg.cor.k] < p.min)]
    
    plot(cpk$Start[pos.cor.k], cpk$cor.coef[pos.cor.k], xlim=c(0, chrdata$chr.toadd[k+1]), ylim=c(-1, 1), type="h", col="grey75", xaxs="i", xlab="Chromosomal position", ylab="Correlation coefficient", main=paste("CGH - GEX correlation (", cor.method, ") on chr", k, "\n", length(sig.pos.cor.k), " positive correlation(s) / ",length(sig.neg.cor.k), " negative correlation(s) at adj.p<", sprintf("%1.1E", p.min), " out of ", nrow(cpk), " features", sep=""))
    if (length(sig.pos.cor) > 0) {
      abline(h=min.sig.pcoef, lty=2, col=4)
      text(0, min.sig.pcoef-0.05, sprintf("%1.2f", min.sig.pcoef), pos=4, col=4)
      if (length(sig.pos.cor.k) > 0) {
        segments(cpk$Start[sig.pos.cor.k], 0, cpk$Start[sig.pos.cor.k], cpk$cor.coef[sig.pos.cor.k], col=4)
        if (symbol) text(cpk$Start[sig.pos.cor.k], cpk$cor.coef[sig.pos.cor.k], labels=cpk$Symbol[sig.pos.cor.k], col=4, pos=3)
      }
    }
    if (length(sig.neg.cor) > 0) {
      lines(cpk$Start[neg.cor.k], cpk$cor.coef[neg.cor.k], type="h", col="grey25")
      abline(h=min.sig.ncoef, lty=2, col=2)
      text(0, min.sig.ncoef+0.05, sprintf("%1.2f", min.sig.ncoef), pos=4, col=2)
      if (length(sig.neg.cor.k) > 0) {
        segments(cpk$Start[sig.neg.cor.k], 0, cpk$Start[sig.neg.cor.k], cpk$cor.coef[sig.neg.cor.k], col=2)
        if (symbol) text(cpk$Start[sig.neg.cor.k], cpk$cor.coef[sig.neg.cor.k], labels=cpk$Symbol[sig.neg.cor.k], col=2, pos=1)
      }
    }
    lines(cpk$Start, grnmd[kfeatures], col="darkgreen")
    abline(h=0, lty=2)
  }
  dev.off()
  
  # Sort on Symbol
  genomic.order <- order(correlation.probes$Symbol)
  correlation.probes <- correlation.probes[genomic.order,]
  bigmat <- bigmat[genomic.order,]
  
  # Getting index of positive/negative correlations, and their significant subpart at a defined p-value
  pos.cor <- which(correlation.probes$cor.coef > 0)
  sig.pos.cor <- pos.cor[which(correlation.probes$cor.adjp[pos.cor] < p.min)]
  min.sig.pcoef <- correlation.probes$cor.coef[sig.pos.cor][which.max(correlation.probes$cor.adjp[sig.pos.cor])]
  neg.cor <- which(correlation.probes$cor.coef < 0)
  sig.neg.cor <- neg.cor[which(correlation.probes$cor.adjp[neg.cor] < p.min)]
  min.sig.ncoef <- correlation.probes$cor.coef[sig.neg.cor][which.max(correlation.probes$cor.adjp[sig.neg.cor])]
  
  
  # Cheeseplots
  # Positive correlations
  if (length(sig.pos.cor > 0)) {
    cpp <- correlation.probes[sig.pos.cor,]
    bmp <- bigmat[sig.pos.cor,]
    pdf(file=paste(outname,"_CGH-GEXcor_", cor.method, "_p", p.min, "_cmp", cor.min.value, "_hg", hg, "_cheeseplots_pos", ".pdf", sep=""), height=21/cm(1),width=29.7/cm(1))
    par(mfrow=c(3,4), mar=c(2,2,2,2))
    for (m in 1:nrow(cpp)) {
      boxplot(bmp[m,GEpart], fg=1, main=paste(cpp$Symbol[m], ", ", cpp$ProbeName[m], ", ", cpp$CytoBand[m], sep=""), cex.main=1, col="grey75", border="grey50")
      gain <- which(bmp[m,CGHpart] > 0)
      loss <- which(bmp[m,CGHpart] < 0)
      normal <- which(bmp[m,CGHpart] == 0)
	  circlesVec = (abs(bmp[m,CGHpart]) / max(abs(bmp[m,CGHpart])) * 0.15)
      if (length(gain) > 0) symbols(rep(1, length(gain)), bmp[m,GEpart][gain], circles=circlesVec[gain], inches=F, add=T, fg=4)
      if (length(loss) > 0) symbols(rep(1, length(loss)), bmp[m,GEpart][loss], circles=circlesVec[loss], inches=F, add=T, fg=2)
      if (length(normal) > 0) points(rep(1, length(normal)), bmp[m,GEpart][normal], pch=20)
      text(1, min(bmp[m,GEpart], na.rm=T), labels=c(paste("Cor=", round(cpp$cor.coef[m], digits=2), ", rawP=", format(cpp$cor.rawp[m], digits=2, scientific=T), ", adjP=", format(cpp$cor.adjp[m], digits=2, scientific=T), sep="")), cex=0.9)
	  abline(h=0, lty=2, col="gray40")
    }
    dev.off()
  }
  
  # Negative correlations
  if (length(sig.neg.cor > 0)) {
    cpn <- correlation.probes[sig.neg.cor,]
    bmn <- bigmat[sig.neg.cor,]
    pdf(file=paste(outname, "_CGH-GEXcor_", cor.method, "_p", p.min, "_cmp", cor.min.value, "_hg", hg, "_cheeseplots_neg", ".pdf", sep=""), height=21/cm(1),width=29.7/cm(1))
    par(mfrow=c(3,4), mar=c(2,2,2,2))
    for (m in 1:nrow(cpn)) {
      boxplot(bmn[m,GEpart], fg=1, main=paste(cpn$Symbol[m], ", ", cpn$ProbeName[m], ", ", cpn$CytoBand[m], sep=""), cex.main=1, col="grey75", border="grey50")
      gain <- which(bmn[m,CGHpart] > 0)
      loss <- which(bmn[m,CGHpart] < 0)
      normal <- which(bmn[m,CGHpart] == 0)
	  circlesVec = (abs(bmp[m,CGHpart]) / max(abs(bmp[m,CGHpart])) * 0.15)
      if (length(gain) > 0) symbols(rep(1, length(gain)), bmn[m,GEpart][gain], circles=circlesVec[gain], inches=F, add=T, fg=4)
      if (length(loss) > 0) symbols(rep(1, length(loss)), bmn[m,GEpart][loss], circles=circlesVec[loss], inches=F, add=T, fg=2)
      if (length(normal) > 0) points(rep(1, length(normal)), bmn[m,GEpart][normal], pch=20)
      text(1, min(bmn[m,GEpart], na.rm=T), labels=c(paste("Cor=", round(cpn$cor.coef[m], digits=2), ", rawP=", format(cpn$cor.rawp[m], digits=2, scientific=T), ", adjP=", format(cpn$cor.adjp[m], digits=2, scientific=T), sep="")), cex=0.9)
      abline(h=0, lty=2, col="gray40")
  	}
    dev.off()
  }
  
  return(list(bigmat=bigmat,correlation.probes=correlation.probes, cpp=cpp) )
}




## CALL EXAMPLE
# setwd("/mnt/bioinfo/P21/P21_SHB_CGH/Analyses/P21_SHB_CGH_GC5_36s_hg19_20101112_noY/CGH-GEX_correlation")
# GEdata <- read.table("P21_SHB_GE_BrB_QN_int_samplenameranked.txt", header=T, sep="\t", stringsAsFactors=F)
# CGHreg <- read.table("P21_SHB_CGH_GC5_33s_hg19_noY_samplenameranked.reg", header=T, sep="\t", stringsAsFactors=F)

# GEd <- GEdata
# CGHd <- CGHreg
# isGElog <- F    ## tells if data are log-based or not (ratios or intensities)
# cor.method <- "pearson"   ## "pearson" or "spearman"
# cor.min.value <- 3    ## only needed for cor.mode = "basic"
# p.min <- 0.05
# hg <- 19

# cghgexcor.basic(CGHd=CGHreg, GEd=GEdata, isGElog=F, cor.method="pearson", cor.min.value=3, p.min=0.05, hg=19)


