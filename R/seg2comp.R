## seg2comp.R
##
## DESCRIPTION : Allows the comparison of two segmented CGH profiles, based on their segments' status.
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
##
## DEPENDS ON:
## . r-base only.
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.8.2 to 2.13.0
##
## VERSION NOTES:
##
## v2.4  20111216
##      . Split the S2C output file to 4 versions :
##          . *_all.s2c, containing all regions
##          . *._diff.s2c, containing only differing regions (GN+NG+LN+NL+GL+LG)
##          . *._id.s2c, containing only identical regions (GG+LL+NN)
##          . *._nnid.s2c, containing only non-neutral identical regions (GG+LL)
##        This is to allow the possibility of characterizing identity regions with grd, which was not
##        possible earlier as grd was filtering out identity regions.
##      . Due to this move, modified the formulas to compute aberrations rates per sample, and
##        identity rates which were wrong in the new context.
##      . Cosmetic : changed the lenged labels in the plot.
##      . Cosmetic : changed the globale table filename.
##
## v2.3d  20111130
##      . Corrected a typo which made the script hang when input file had 4 columns.
##      . Corrected the help content, which mixed up real function options and file format.
##
## v2.3c  20111116
##      . Introduced a "help" function to display the description of parameters for the main function, even when
##        this script is sourced and not directly read.
##
## v2.3b  20111005
##      . Added a filter to replace dots with underscores in sample names.
##
## v2.3   20110831
##      . Added support for mus musculus.
##      . In this purpose, dropped the -hg option for two new options : -sp for species, and -gb for the genome build.
##
## v2.2   20110719
##      . Rewrote the script quite completely. Now its intent is just limited to compare two samples
##        according to their segments' status. No fit is made.
##      . Removed the "noY" option. chrY is now removed in every case.
##
## v2.1   20110628
##      . Modified the script behavior. Now, a table is required as an input, which must contain
##        6 columns with, per line, the lst filename of the barcode1, same for barcode2, then the name
##        to find the barcode1 in the lst, same for barcode2, then the log2(ratio) cut-off for the first
##        profile to discriminate normality (can be found in the GC5.R results), same for barcode2,
##        in this specified order. Whatever the column names are, a header must exist for this table.
##        So now, this script will work in batch on the gcx files given in the first two columns,
##        which will be named according to the third and fourth columns.
##        NOTA : if the samplenames for barcode1 and/or barcode2 are not given, the barcode will be used.
##      . Since the batch mode has been established, an output of the batch's results is now
##        generated as a datatable. 
##
## v2.0   20101109
##      . Refonte complete du script.
##		  . Support du format 'lst', abandon du support 'abr0'.
##      . Support des informations de chromosome / bandes cyto directement depuis les tables UCSC.
##	      'cytoBandIdeo' pour les g?nomes hg17, hg18 et hg19.
##      . Possibilit? d'?liminer les r?gions normales (ie, de log2ratio=0 pour les deux samples).
##      . Possibilit? d'?liminer les segments des chrX et/ou chrY.
##      . Le format de sortie d?fini a pour extension *.compdf (seg2comp)

chrconv.list <<- list(
  chrom2chr = list(
    hs=list(
      "chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chr21"=21,"chr22"=22,"chrX"=23,"chrY"=24
    ),
    mm=list(
      "chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chrX"=20,"chrY"=21
    ),
    rn=list(
      "chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chrX"=21,"chrY"=22
    )
  ),
  chr2chrom=list(
    hs=list(
      "1"="chr1","2"="chr2", "3"="chr3", "4"="chr4","5"="chr5","6"="chr6","7"="chr7","8"="chr8","9"="chr9","10"="chr10","11"="chr11","12"="chr12","13"="chr13","14"="chr14","15"="chr15","16"="chr16","17"="chr17","18"="chr18","19"="chr19","20"="chr20","21"="chr21","22"="chr22","23"="chrX","24"="chrY"
    ),
    mm=list(
      "1"="chr1","2"="chr2", "3"="chr3", "4"="chr4","5"="chr5","6"="chr6","7"="chr7","8"="chr8","9"="chr9","10"="chr10","11"="chr11","12"="chr12","13"="chr13","14"="chr14","15"="chr15","16"="chr16","17"="chr17","18"="chr18","19"="chr19","20"="chrX","21"="chrY"
    ),
    rn=list(
      "1"="chr1","2"="chr2", "3"="chr3", "4"="chr4","5"="chr5","6"="chr6","7"="chr7","8"="chr8","9"="chr9","10"="chr10","11"="chr11","12"="chr12","13"="chr13","14"="chr14","15"="chr15","16"="chr16","17"="chr17","18"="chr18","19"="chr19","20"="chr20","21"="chrX","22"="chrY"
    )
  )
)

## PAREMETERS DESCRIPTOR FUNCTION
help.seg2comp <- function() {
  cat("
PARAMETERS FOR THE seg2comp() FUNCTION :
----------------------------------------
outroot   chr       Name of the rootname for output files.
inputfile chr       (path+)Name of the list of files and barcodes (see \"INPUT FILE FORMAT\" below).
sp        chr       Species.
gb		    int       Genome build. 18, 19 for homo sapiens, 9 for mus musculus
noX       bool      Discard chrX data.
mnt.proj  chr       Mount point for pelican:/proj.


INPUT FILE FORMAT :
-------------------
This input file must contain a header (column names are non-fixed), and at least 4 columns are
mandatory, in this order :
lst1  	  chr		     First LST filename (for barcode 1).
lst2		  chr		     Second LST filename (for barcode 2).
barcode1	chr		     First barcode, to get from lst1.
barcode2	chr		     Second barcode, to get from lst2.

Two optional columns can be added, giving the sample names to display instead of barcodes :
sample1   chr       Name of the sample corresponding to barcode1.
sample2   chr       Name of the sample corresponding to barcode2.
\n")
}

## MAIN FUNCTION
seg2comp <- function(outroot="seg2comp_results", inputfile, sp="hs", gb=19, noX=F, cytoBandIdeo = "/mnt/data_cigogne/bioinfo/GoldenPath/hg19/cytoBandIdeo.hg19") {
  gv <- ""
  if (sp == "hs") {
    gv = "hg"
  } else if (sp == "mm") {
    gv = "mm"
  }
  gv <- paste(gv, gb, sep="")
  suffix <- paste("_", gv, sep="")
  if (noX) suffix <- paste(suffix, "_noX", sep="")
	inlist <- read.table(inputfile, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  ## Replacing dots if samplenames exist
  if (ncol(inlist) > 4) {
    inlist[,5] <- chartr(".", "_", inlist[,5])
    inlist[,6] <- chartr(".", "_", inlist[,6])
  }
  outdf <- apply(inlist, 1, function(x) {
    if(length(x) == 4) {
      dsp1 = x[3]
      dsp2 = x[4]
    }
    else {
      dsp1=x[5]
      dsp2=x[6]
    }
    s2c.core(lst1=x[1], lst2=x[2], barcode1=x[3], barcode2=x[4], display1=dsp1, display2=dsp2, sp=sp, gb=gb, gv=gv, noX=F, cytoBandIdeo=cytoBandIdeo)
  })
  outdf <- as.data.frame(matrix(unlist(outdf, recursive=T), ncol=14, byrow=T), stringsAsFactors=F)
  colnames(outdf) <- c("Sample1", "Barcode1", "Sample2" ,"Barcode2", "ab.rate.sample1", "ab.rate.sample2", "ab.rate.union", "frac.id.gen", "frac.opp.gen", "frac.id.ab", "frac.opp.ab", "spe.gain", "spe.loss", "noX")
  write.table(outdf, paste(outroot, "_seg2comp_results", suffix, ".txt", sep=""), sep="\t", quote=F, row.names=F)
}

s2c.core <- function(lst1, lst2, barcode1, barcode2, display1, display2, sp=sp, gb=gb, gv=gv, noX=F, cytoBandIdeo = NULL) {
	
  suffix <- paste("_", gv, sep="")
	if (noX) suffix <- paste(suffix, "_noX", sep="")
	
  l2r.display1 <- paste("l2r", display1, sep=".")
  l2r.display2 <- paste("l2r", display2, sep=".")
  status.display1 <- paste("status", display1, sep=".")
  status.display2 <- paste("status", display2, sep=".")
 
	if (sp == "hs") {
    xnum <- 23
    ynum <- 24
  } else if (sp == "mm") {
    xnum <- 20
    ynum <- 21
  }
  cytob <- read.table(cytoBandIdeo, sep="\t", header=T, stringsAsFactors=F, comment.char="_", fill=T)
#   cytob <- read.table(paste(mnt.proj, "/cgh/", gv, "/cytoBandIdeo.", gv, sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
  
  ## Filtering out invalid lines (UCSC recently added incomplete lines which seem to refer to chrN_unknown sequence sets, which is quite totally illogical. BJ 20140227)
  cytob <- cytob[which(!is.na(cytob$chromStart)),]
  ## Filtering out ANY chromosome not defined in chrconvlist$chrom2chr
  cytob <- cytob[which(cytob[["X.chrom"]] %in% names(chrconv.list$chrom2chr[[sp]])),]

	cytob$chr <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
	cytob$chr[which(cytob$chr == "X")] <- xnum
	cytob$chr[which(cytob$chr == "Y")] <- ynum
	cytob$chr <- as.numeric(cytob$chr)
	cytob <- cytob[order(cytob$chr),]

	lchrxx <- sapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) })
	lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
	glen <- sum(as.numeric(lchrxx))
	lchrsum <- sapply(c(1:length(lchrtoadd)), function(x) { sum(lchrtoadd[1:x]) })
	
  vsword = paste(display1, "_vs_", display2, sep="")
	cat(vsword, "\n")
 
	## Ouverture de table1 et conservation des aberrations du barcode1
	tab1 <- read.table(lst1, sep="\t", header=T, stringsAsFactors=F)
	sel1 <- which(tab1$Sample==barcode1)
	if (length(sel1) == 0) ab.rate.barcode1 <- NA else tab1sel <- tab1[sel1,]
	
	## Si table et table2 sont le meme fichier, alors on prend juste les aberrations de barcode2 depuis table1
  tab2 <- read.table(lst2, sep="\t", header=T, stringsAsFactors=F)
	if (lst1 == lst2) sel2 <- which(tab1$Sample==barcode2) else sel2 <- which(tab2$Sample==barcode2)
  if (length(sel2) == 0) ab.rate.barcode2 <- NA else tab2sel <- tab2[sel2,]
  
  ## Case where at least one of sel1 / sel2 is empty
  if ( (length(sel1) == 0) | (length(sel2) == 0) ) return(c(display1, barcode1, display2, barcode2, ab.rate.barcode1, ab.rate.barcode2, NA, NA, NA, NA, NA, NA, NA, noX))
  
  cur.glen <- glen
  ## Conditional removing of chrX
  if (noX) {
    tab1sel <- tab1sel[which(tab1sel$Chrom != "chrX"),]
    tab2sel <- tab2sel[which(tab2sel$Chrom != "chrX"),]
    cur.glen <- cur.glen - lchrxx[xnum]
  }
  
  ## Forced removing of chrY
  tab1sel <- tab1sel[which(tab1sel$Chrom != "chrY"),]
  tab2sel <- tab2sel[which(tab2sel$Chrom != "chrY"),]
	cur.glen <- cur.glen - lchrxx[ynum]
  
  ## Computing per-sample aberration rate
  ab.index.barcode1 <- which(tab1sel[["Log2ratio"]] != 0)
  ab.rate.barcode1 <- sum(tab1sel$End[ab.index.barcode1] - tab1sel$Start[ab.index.barcode1] +1) / cur.glen
  ab.index.barcode2 <- which(tab2sel[["Log2ratio"]] != 0)
  ab.rate.barcode2 <- sum(tab2sel$End[ab.index.barcode2] - tab2sel$Start[ab.index.barcode2] +1) / cur.glen
#   ab.rate.barcode2 <- sum(tab2sel$End - tab2sel$Start +1) / cur.glen
  
	## Pour chaque chromosome present dans au moins 1 des 2 profils
	compdf <- data.frame(stringsAsFactors=F)
	
	for (k in sort(unique(c(tab1sel$Chr, tab2sel$Chr)))) {
		
		kchrom <- paste("chr", k, sep="")
		if (kchrom == paste("chr", xnum, sep="")) kchrom <- "chrX"
    if (kchrom == paste("chr", ynum, sep="")) kchrom <- "chrY"
		tab1chr <- tab1sel[which(tab1sel$Chr == k),]
		tab2chr <- tab2sel[which(tab2sel$Chr == k),]
		
		## Listing des bornes utiles
		vecpos <- sort(unique(c(1, lchrxx[k], tab1chr$Start, tab1chr$End, tab2chr$Start, tab2chr$End)))
		vecpos2 <- vecpos + sum(lchrtoadd[1:k])
		
		for (i in 1:(length(vecpos)-1)) {
			le <- vecpos[i]
			re <- vecpos[i+1]
			
			lr1 <- 0
			if (nrow(tab1chr) > 0) {
				lr1index <- which( (tab1chr$Start < re) & (tab1chr$End > le) )
				if (length(lr1index) > 1) lr1index <- lr1index[1]
				lr1 <- tab1chr$Log2ratio[lr1index]
			}
			if (length(lr1) == 0) lr1 <- 0
			lr2 <- 0
			if (nrow(tab2chr) > 0) {
				lr2index <- which( (tab2chr$Start < re) & (tab2chr$End > le) )
				if (length(lr2index) > 1) lr2index <- lr2index[1]
				lr2 <- tab2chr$Log2ratio[lr2index]
			}
			if (length(lr2) == 0) lr2 <- 0
			
      ## Filtered no-normal/cut case.
#       if ( (lr1 != 0) | (lr2 != 0) ) compdf <- rbind(compdf, data.frame(Loc= paste(kchrom, ":", le, "-", re, sep=""), Chr=k, Start=le, End=re, Width=re-le+1, l2r.display1=lr1, l2r.display2=lr2, stringsAsFactors=F))
      compdf <- rbind(compdf, data.frame(Loc= paste(kchrom, ":", le, "-", re, sep=""), Chr=k, Start=le, End=re, Width=re-le+1, l2r.display1=lr1, l2r.display2=lr2, stringsAsFactors=F))
		}
	}
  compdf <- cbind(compdf$Loc, vsword, compdf[,c(2:ncol(compdf))], stringsAsFactors=F)
  colnames(compdf) <- c("Loc", "Comparison", "Chr", "Start", "End", "Width", l2r.display1, l2r.display2)
  
  compdf$l2r.diff <- compdf[[l2r.display1]] - compdf[[l2r.display2]]
  compdf[[status.display2]] <- compdf[[status.display1]] <- NA
  compdf[[status.display1]][which(compdf[[l2r.display1]] < 0)] <- "L"
  compdf[[status.display1]][which(compdf[[l2r.display1]] == 0)] <- "N"
  compdf[[status.display1]][which(compdf[[l2r.display1]] > 0)] <- "G"
  compdf[[status.display2]][which(compdf[[l2r.display2]] < 0)] <- "L"
  compdf[[status.display2]][which(compdf[[l2r.display2]] == 0)] <- "N"
  compdf[[status.display2]][which(compdf[[l2r.display2]] > 0)] <- "G"
  compdf$status.diff <- apply(compdf, 1, function(x) { if (x[10] == x[11]) x[10] else paste(x[10], x[11], sep="")})
  compdf$status.status <- compdf$status.diff
  compdf$status.status[which(nchar(compdf$status.status) == 1)] <- "Id"
  compdf$status.status[which( (compdf$status.status == "GL") | (compdf$status.status == "LG") )] <- "Opp"
  compdf$status.status[which( (compdf$status.status == "NG") | (compdf$status.status == "GN") | (compdf$status.status == "NL") | (compdf$status.status == "LN") )] <- "Diff"
  
#   cur.ablen <- sum(compdf$Width)
  
#   index.id.gain <- which((nchar(compdf$status.diff) == 1) & (compdf[[status.display1]] == "G") & (compdf[[status.display2]] == "G"))
  index.id.gain <- which(compdf$status.diff == "G")
#   index.id.loss <- which((nchar(compdf$status.diff) == 1) & (compdf[[status.display1]] == "L") & (compdf[[status.display2]] == "L"))
  index.id.loss <- which(compdf$status.diff == "L")
#   index.id.norm <- which((nchar(compdf$status.diff) == 1) & (compdf[[status.display1]] == "N") & (compdf[[status.display2]] == "N"))
  index.id.norm <- which(compdf$status.diff == "N")
#   index.diff.gain <- which((nchar(compdf$status.diff) == 2) & (((compdf[[status.display1]] == "G") & (compdf[[status.display2]] == "N")) | ((compdf[[status.display1]] == "N") & (compdf[[status.display2]] == "G"))))
  index.diff.gain <- which((compdf$status.diff == "GN") | (compdf$status.diff == "NG"))
#   index.diff.loss <- which((nchar(compdf$status.diff) == 2) & (((compdf[[status.display1]] == "L") & (compdf[[status.display2]] == "N")) | ((compdf[[status.display1]] == "N") & (compdf[[status.display2]] == "L"))))
  index.diff.loss <- which((compdf$status.diff == "LN") | (compdf$status.diff == "NL"))
#   index.opposite <- which((nchar(compdf$status.diff) == 2) & (((compdf[[status.display1]] == "L") & (compdf[[status.display2]] == "G")) | ((compdf[[status.display1]] == "G") & (compdf[[status.display2]] == "L"))))
  index.opposite <- which((compdf$status.diff == "GL") | (compdf$status.diff == "LG"))
  
  nb.id.gain <- length(index.id.gain)
  nb.id.loss <- length(index.id.loss)
  nb.id.norm <- length(index.id.norm)
  nb.id.tot <- nb.id.gain + nb.id.loss + nb.id.norm
  nb.diff.gain <- length(index.diff.gain)
  nb.diff.loss <- length(index.diff.loss)
  nb.diff.tot <- nb.diff.gain + nb.diff.loss
  nb.opposite <- length(index.opposite)
  
  width.id.gain <- sum(compdf$Width[index.id.gain])
  width.id.loss <- sum(compdf$Width[index.id.loss])
  width.id.norm <- sum(compdf$Width[index.id.norm])
  cur.ablen <- cur.glen - width.id.norm
#   width.id.norm <- cur.glen - cur.ablen
  width.id.tot <- width.id.gain + width.id.loss + width.id.norm
  width.diff.gain <- sum(compdf$Width[index.diff.gain])
  width.diff.loss <- sum(compdf$Width[index.diff.loss])
  width.opposite <- sum(compdf$Width[index.opposite])
  width.diff.tot <- width.diff.gain + width.diff.loss + width.opposite
  
  if (nb.id.gain > 0) spe.gain <- width.id.gain / sum(width.id.gain, width.diff.gain, width.opposite) else spe.gain <- NA
  if (nb.id.loss > 0) spe.loss <- width.id.loss / sum(width.id.loss, width.diff.loss, width.opposite) else spe.loss <- NA
  
  ab.rate.union <- cur.ablen / cur.glen
  
  if (nb.id.gain > 0) frac.id.gain.gen <- width.id.gain / cur.glen else frac.id.gain.gen <- NA
  if (nb.id.loss > 0) frac.id.loss.gen <- width.id.loss / cur.glen else frac.id.loss.gen <- NA
  frac.id.norm.gen <- width.id.norm / cur.glen
  frac.id.gen <- (cur.glen - sum(width.diff.gain, width.diff.loss, width.opposite)) / cur.glen
  if (nb.diff.gain > 0) frac.diff.gain.gen <- width.diff.gain / cur.glen else frac.diff.gain.gen <- NA
  if (nb.diff.loss > 0) frac.diff.loss.gen <- width.diff.loss / cur.glen else frac.diff.loss.gen <- NA
  if (nb.opposite > 0) frac.opp.gen <- width.opposite / cur.glen else frac.opp.gen <- NA
  if ( (nb.diff.gain > 0) | (nb.diff.loss > 0) | (nb.opposite > 0) ) frac.diff.tot <- 
  
  if (nb.id.gain > 0) frac.id.gain.ab <- width.id.gain / cur.ablen else frac.id.gain.ab <- NA
  if (nb.id.loss > 0) frac.id.loss.ab <- width.id.loss / cur.ablen else frac.id.loss.ab <- NA
  frac.id.ab <- (cur.ablen - sum(width.diff.gain, width.diff.loss, width.opposite)) / cur.ablen
  if (nb.diff.gain > 0) frac.diff.gain.ab <- width.diff.gain / cur.ablen else frac.diff.gain.ab <- NA
  if (nb.diff.loss > 0) frac.diff.loss.ab <- width.diff.loss / cur.ablen else frac.diff.loss.ab <- NA
  if (nb.opposite > 0) frac.opp.ab <- width.opposite / cur.ablen else frac.opp.ab <- NA
  
  ## DUMPING THE RESULT TABLES
  ## Full
  write.table(compdf, paste(vsword, suffix, "_all.s2c", sep=""), sep="\t", quote=F, row.names=F)
  ## Diff
  diff.index <- which((compdf[["status.status"]] == "Diff") | (compdf[["status.status"]] == "Opp"))
  write.table(compdf[diff.index,], paste(vsword, suffix, "_diff.s2c", sep=""), sep="\t", quote=F, row.names=F)
  ## All Id
  id.index <- which(compdf[["status.status"]] == "Id")
  write.table(compdf[id.index,], paste(vsword, suffix, "_id.s2c", sep=""), sep="\t", quote=F, row.names=F)
  ## Non-neutral Id
  nnid.index <- which((compdf[["status.status"]] == "Id") & (compdf[["status.diff"]] != "N"))
  write.table(compdf[nnid.index,], paste(vsword, suffix, "_nnid.s2c", sep=""), sep="\t", quote=F, row.names=F)
 
  ## Plot
  png(paste(vsword, suffix, "_s2c.png", sep=""), width=1680, height=1050)
  par(xaxs="i", yaxs="i", xaxt="n", yaxt="n", cex=2)
  plot(0, 0, xlim=c(0, cur.glen), ylim=c(-1,6), type="n", xlab="Genomic position", ylab="", main=paste("Segmented genomic profiles of samples ", display1, " and ", display2, "\nand their differences in segmentation status", sep=""))
  
  col.N <- "grey75"
  col.G <- "blue"
  col.L <- "red"
  col.GN <- "cyan"
  col.LN <- "orange"
  col.GL <- "black"
  
  colvec <- c(col.N, col.G, col.L, col.GN, col.LN, col.GL)
  namvec <- c("Normal (N)", "Gain (G)", "Loss (L)", "G vs N, N vs G", "L vs N, N vs L", "G vs L, L vs G")
  
  for (k in 0:5) {
    points(k*glen/6+glen*0.02, -0.75, col=colvec[k+1], pch=20, cex=3)
    text((k*glen/6+glen*0.03), -0.75, namvec[k+1], pos=4)
  }
  
  rect(0,4,cur.glen,5, col=col.N, border=col.N)
  rect(0,2,cur.glen,3, col=col.N, border=col.N)
  rect(0,0,cur.glen,1, col=col.N, border=col.N)
  
  text(lchrsum+lchrxx/2, 5.8, 1:23)
  g1 <- which(compdf[[status.display1]] == "G")
  if(length(g1)>0) rect(compdf$Start[g1]+lchrsum[compdf$Chr[g1]], 4, compdf$End[g1]+lchrsum[compdf$Chr[g1]], 5, col=col.G, border=col.G)
  l1 <- which(compdf[[status.display1]] == "L")
  if(length(l1)>0) rect(compdf$Start[l1]+lchrsum[compdf$Chr[l1]], 4, compdf$End[l1]+lchrsum[compdf$Chr[l1]], 5, col=col.L, border=col.L)
  
  g2 <- which(compdf[[status.display2]] == "G")
  if(length(g2)>0) rect(compdf$Start[g2]+lchrsum[compdf$Chr[g2]], 2, compdf$End[g2]+lchrsum[compdf$Chr[g2]], 3, col=col.G, border=col.G)
  l2 <- which(compdf[[status.display2]] == "L")
  if(length(l2)>0) rect(compdf$Start[l2]+lchrsum[compdf$Chr[l2]], 2, compdf$End[l2]+lchrsum[compdf$Chr[l2]], 3, col=col.L, border=col.L)
  
  idg <- which(compdf$status.diff == "G")
  if(length(idg)>0) rect(compdf$Start[idg]+lchrsum[compdf$Chr[idg]], 0, compdf$End[idg]+lchrsum[compdf$Chr[idg]], 1, col=col.G, border=col.G)
  diffg <- which( (compdf$status.diff == "GN") | (compdf$status.diff == "NG") )
  if(length(diffg)>0) rect(compdf$Start[diffg]+lchrsum[compdf$Chr[diffg]], 0, compdf$End[diffg]+lchrsum[compdf$Chr[diffg]], 1, col=col.GN, border=col.GN)
  idl <- which(compdf$status.diff == "L")
  if(length(idl)>0) rect(compdf$Start[idl]+lchrsum[compdf$Chr[idl]], 0, compdf$End[idl]+lchrsum[compdf$Chr[idl]], 1, col=col.L, border=col.L)
  diffl <- which( (compdf$status.diff == "LN") | (compdf$status.diff == "NL") )
  if(length(diffl)>0) rect(compdf$Start[diffl]+lchrsum[compdf$Chr[diffl]], 0, compdf$End[diffl]+lchrsum[compdf$Chr[diffl]], 1, col=col.LN, border=col.LN)
  opp <- which( (compdf$status.diff == "LG") | (compdf$status.diff == "GL") )
  if(length(opp)>0) rect(compdf$Start[opp]+lchrsum[compdf$Chr[opp]], 0, compdf$End[opp]+lchrsum[compdf$Chr[opp]], 1, col=col.GL, border=col.GL)
  
  if(length(idg)>0) segments(compdf$Start[idg]+lchrsum[compdf$Chr[idg]], .83, compdf$End[idg]+lchrsum[compdf$Chr[idg]], .83, col=col.G, lwd=10)
  if(length(idl)>0) segments(compdf$Start[idl]+lchrsum[compdf$Chr[idl]], .17, compdf$End[idl]+lchrsum[compdf$Chr[idl]], .17, col=col.L, lwd=10)
  if(length(diffg)>0) segments(compdf$Start[diffg]+lchrsum[compdf$Chr[diffg]], .67, compdf$End[diffg]+lchrsum[compdf$Chr[diffg]], .67, col=col.GN, lwd=10)
  if(length(diffl)>0) segments(compdf$Start[diffl]+lchrsum[compdf$Chr[diffl]], .33, compdf$End[diffl]+lchrsum[compdf$Chr[diffl]], .33, col=col.LN, lwd=10)
  if(length(opp)>0) segments(compdf$Start[opp]+lchrsum[compdf$Chr[opp]], .5, compdf$End[opp]+lchrsum[compdf$Chr[opp]], .5, col=col.GL, lwd=10)
  text(cur.glen*0.01, 5.1, display1, pos=4)
  text(cur.glen*0.01, 3.1, display2, pos=4)
  text(cur.glen*0.01, 1.1, "Differences in segmentation", pos=4)
  abline(v=lchrsum, lty=2)
  dev.off()
 
  return(c(display1, barcode1, display2, barcode2, ab.rate.barcode1, ab.rate.barcode2, ab.rate.union, frac.id.gen, frac.opp.gen, frac.id.ab, frac.opp.ab, spe.gain, spe.loss, noX))

}
