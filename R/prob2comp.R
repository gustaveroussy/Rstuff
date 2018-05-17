## prob2comp.R
##
## DESCRIPTION :  Script destiné à la comparaison de deux profils de CGH-arrays. Elle commence par déterminer
##			le profil avec la plus faible dynamique (comparaison des IQR), l'ajuste sur celui de plus
##			forte dynamique (least-square fit sur une sous-population de sondes recueillie via la courbe
##			de densité de la distribution des différences absolues entre profils), puis détermine si les
##			régions génomiques observées sont relevantes ou non (via segmentation du profil des différences,
##			et utilisation du bruit intersondes pour la détermination d'un cut-off).
##			C'est une refonte from scratch (3eme version) qui cette fois prend en compte les valeurs des
##			sondes des profils normalisés (post-GC5) des fichers "*.gcx", soit directement, soit compressés
##			en bzip2 (directement lus).
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.10.0 to 2.15.1
##
## DEPENDS ON:
## DNAcopy
##
## VERSION NOTES
##
## v3.9g 20140422
##      . Just updated the help.prob2comp()
##
## v3.9f 20140303
##      . Added support for rattus norvegicus (untested, though)
##
## v3.9e 20121105
##      . Added chromosomes names on plots.
##
## v3.9d 20120328
##      . Corrected a bug in the last plot, which did not correct well the seg.med value
##        after the 2nd step of centralization.
##      . Added a CBS output.
##
## v3.9c 20120209
##      . Corrected a typo in the GCX output filename.
##
## v3.9b 20120201
##      . Added a command to stop the threads when ncores > 1 at the end of the process.
##
## v3.9 20120104
##      . Added (conditional) multithreading thanks to the 'snowfall' package.
## 
## v3.8 20111130
##      . Adding an output of the differential profile, in GCX format.
##      . Added support for GZ and XZ compression formats.
##
## v3.7b 20111116
##      . Introduced a "help" function to display the description of parameters for the main function, even when
##        this script is sourced and not directly read.
##
## v3.7 20111114
##      . Modified the CBS results to have median value of probes per segment, instead of default,
##        forced mean value, which causes problems when "smears" appear in a chr or whole profile.
##      . Changed default nrf value to 0.5 (lowered thanks to this mean to median conversion).
##      . True segments' ends are now retrieved in a more robust way, using gcdata info by invoking indexes.
##      . Coordinates are not converted to Kb anymore for segmentation, it's unneeded in current-gen
##        CBS implementation.
##      . Set blue 4 color for genomic regions with more content (kept red 1 for less).
##
## v3.6b  20111005
##      . Added a filter to replace dots with underscores in sample names.
##
## v3.6   20110831
##      . Added support for mus musculus.
##      . In this purpose, dropped the -hg option for two new options : -sp for species, and -gb for the genome build.
##
## v3.5b  20110809
##      . Just changed the output filename for the table of the differential segments from
##        *_segdiff_list.txt to *.p2c, to be more homogeneous with the *.s2c output files from
##        the seg2comp.R script.
##
## v3.5   20110728
##      . Added an option to specify the type of input data in the GCX files (log10 or log2), as
##        two versions exist (older GCXs contained log10(ratio), now they are in log2).
##      . Added a condition based on a pearson correlation score computed for the pair of smoothed
##        profiles, so that if this correlation score is lower than a given threshold value, the
##        fit won't be performed.
##      . chrY is now always discarded.
##      . Added a new option to discard the chrX (inactive by default).
##      . Modified output summary table. Now it contains the values of the used parameters, as well as
##        the value of the internal noise between arrays, used to set the differential threshold.
##      . Bugs corrected in the format of the segments' output table.
##
## v3.4c  20110721
##      . Rename the script from dynadjust2.R to prob2comp.R, to be more logical with the other
##        script dedicated to comparison of pairs of samples, seg2comp.R.
##
## v3.4b  20110627
##      . Formatted the script header according to the new requirements in the bioinfo team.
##
## v3.4   20110624
##      . Modified the script behavior. Now, a table is required as an input, which must contain
##        4 columns with, per line, the gcx filename for sample1, same for sample2, then the name
##        to display for sample1, same for sample2, in this specified order. Whatever the column
##        names are, a header must exist for this table. So now, this script will work in batch
##        on the gcx files given in the first two columns, which will be named according to the
##        third and fourth columns. NOTA : if the samplenames for sample1 and/or sample2 are not
##        given, the barcode will be used.
##      . Since the batch mode has been established, an output of the batch's results is now
##        generated as a datatable.
##
## v3.3   20110527
##      . Added an option to modify the mount point of /proj on pelican, so that the script
##        can be called from a pelican session or a local one, whichever the user's /proj mount
##        point is.
##      . Removed (commented) the lines for conversion of GCX log10(ratio) to log2(ratio), as
##        of now GCX files contain log2(ratio) and no more log10(ratio).
##
## v3.2   20110321
##      . Modification de la routine de détection de pic pour l'ajustement de dynamique.
##        Désormais, le pic retenu ne sera plus le second après le principal, mais le second
##        le plus eleve après le principal (toujours si au moins 2 pics sont identifiables).
##
## v3.1   20110209
##      . Ajout d'une centralisation des segments (sur la médiane des segments normaux)
##
## v3.0b  20110208
##      . Ajout des coefficients du fit dans un des fichiers de sortie.
##
## v3.0a	20101109
##      . Refonte complète du script.
##      . Support des données sous la forme de fichiers gcx, compressés bz2 ou non.

## PAREMETERS DESCRIPTOR FUNCTION
help.prob2comp <- function() {
  cat("
      PARAMETERS FOR THE prob2comp() FUNCTION :
      -----------------------------------------
      outroot	  	  chr		  Root of the output filename for the summary table.
      inputfile		  chr		  (path+)Name of the 4-col table giving the names of the GCX files and the samplenames.
      sp    		    chr		  Species. Currently hs and mm are supported
      gb            int     Genome build (19 or 18 for hs, 9 for mm, 4 for rn).
      nrf		        real+   Multiplicative factor for calling the differential segments.
      uSD		        real+   CBS parameter for merging similar segments.
      noX           bool    Discard chrX (default to F).
      cor.cut       float   Correlation threshold to perform the fit or not. Recommanded values are 0.6 for Spearman, 0.5 for Pearson.
      cor.method    float   Correlation method, \"spearman\" or \"pearson\".
      mnt.proj      chr    Mount point for the pelican:/proj position.
      in.type       chr     Specify the type of data : \"log10\" or \"log2\".
      ncores        int     Number of threads/cores to use.
      
      NOTA :
      ------
      1. The script can handle direct GCX files, or compressed versions (.gcx.bz2, .gcx.gz, .gcx.xz).
      2. The script always discards data for chrY.
      \n")
}

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


## MAIN FUNCTION
prob2comp <- function(outroot="prob2comp_results", inputfile, sp="hs", gb=19, nrf=0.5, uSD=0, noX=F, cor.cut=0.5, cor.method="spearman", ldb = "/mnt/data_cigogne/bioinfo/", in.type="log2", ncores=2) {
  inlist <- read.table(inputfile, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  
  if (!(sp %in% names(chrconv.list[["chrom2chr"]]))) stop ("Unknown species ! Are supported : hs, mm or rn.")
  
  sptxt <- sp
  if (sp == "hs") sptxt <- "hg"
  gv <- paste(sptxt, gb, sep="")
  
  # mnt.proj <<- mnt.proj
  
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
      p2c.core(gcx1=x[1], gcx2=x[2], sample1=x[3], sample2=x[4], sp=sp, gb=gb, nrf=nrf, uSD=uSD, noX=noX, cor.cut=cor.cut, cor.method=cor.method, ldb=ldb, in.type=in.type)
    })
    sfStop()
  }
  ## Or normaly...
  else {
    outdf <- apply(inlist, 1, function(x) { p2c.core(gcx1=x[1], gcx2=x[2], sample1=x[3], sample2=x[4], sp=sp, gb=gb, nrf=nrf, uSD=uSD, noX=noX, cor.cut=cor.cut, cor.method=cor.method, ldb=ldb, in.type=in.type) })
  }
  
  outdf <- as.data.frame(matrix(unlist(outdf, recursive=T), ncol=14, byrow=T), stringsAsFactors=F)
  colnames(outdf) <- c("Sample1", "Sample2", "undo.SD", "Nrf", "noX", "Cor.cutoff", "Cor.method", "Smoothed.cor", "Noise.MAD", "S1.more.S2", "S1.less.S2", "S1.diff.S2", "Fit.Slope", "Fit.Intercept")
  write.table(outdf, paste(outroot, "_", gv, "_U", uSD, "_nrf", nrf, "_CC", cor.cut, ".txt", sep=""), sep="\t", quote=F, row.names=F)
}



## IMPORT CHR DATA
chrload <- function(sp, gb, ldb = "/mnt/data_cigogne/bioinfo/") {
  
  sptxt <- sp
  if (sp == "hs") sptxt <- "hg"
  gv <- paste(sptxt, gb, sep="")
  
  ## IMPORT BANDES CYTO ET PREPARATION POUR PLOT KARYO
  cat("Importing chromosomes data  ...\n")
  cytob <- read.table(paste(ldb, "/GoldenPath/", gv, "/cytoBandIdeo.", gv, sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_", fill=T)
  
  ## Filtering out invalid lines (UCSC recently added incomplete lines which seem to refer to chrN_unknown sequence sets, which is quite totally illogical. BJ 20140227)
  cytob <- cytob[which(!is.na(cytob$chromStart)),]
  ## Filtering out ANY chromosome not defined in chrconvlist$chrom2chr
  cytob <- cytob[which(cytob[["X.chrom"]] %in% names(chrconv.list$chrom2chr[[sp]])),]
  
  cytob$chrA <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
  cytob$chr <- cytob$chrA
  
  ##
  cytob$chr[which(cytob$chr == "X")] <- chrconv.list[["chrom2chr"]][[sp]]["chrX"][[1]]
  cytob$chr[which(cytob$chr == "Y")] <- chrconv.list[["chrom2chr"]][[sp]]["chrY"][[1]]
  
  cytob$chr <- as.numeric(cytob$chr)
  cytob <- cytob[order(cytob$chr),]
  
  lchrxx <- sapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) })
  lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
  glen <- sum(as.numeric(lchrxx))
  lchrsum <- sapply(c(1:length(lchrtoadd)), function(x) { sum(lchrtoadd[1:x]) })
  
  cytob$x1 <- 0.15
  cytob$x2 <- 0.85
  cytob$x1[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.25
  cytob$x2[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.75
  cytob$gieStain[which(cytob$gieStain == "gneg")] <- "white"
  cytob$gieStain[which(cytob$gieStain == "gpos25")] <- "grey75"
  cytob$gieStain[which(cytob$gieStain == "gpos33")] <- "grey66"
  cytob$gieStain[which(cytob$gieStain == "gpos50")] <- "grey50"
  cytob$gieStain[which(cytob$gieStain == "gpos66")] <- "grey33"
  cytob$gieStain[which(cytob$gieStain == "gpos75")] <- "grey25"
  cytob$gieStain[which(cytob$gieStain == "gpos100")] <- "black"
  cytob$gieStain[which(cytob$gieStain == "acen")] <- "yellow"
  cytob$gieStain[which(cytob$gieStain == "gvar")] <- "darkred"
  cytob$gieStain[which(cytob$gieStain == "stalk")] <- "darkolivegreen"
  
  return(list(cytob=cytob, lchrxx=lchrxx, lchrsum=lchrsum, lchrtoadd=lchrtoadd, glen=glen, cur.glen=glen, gv=gv))
}

chrconvlist <- list(hs=list("X"))

## CODE
p2c.core <- function(gcx1, gcx2, sample1=NULL, sample2=NULL, sp="hs", gb=19, nrf, uSD, noX, cor.cut, cor.method, ldb = "/mnt/data_cigogne/bioinfo/", in.type="log2") {
  
  library(DNAcopy)
  
  ## Chromosomes and genomic coordinates
  cs <<- chrload(sp, gb)
  
  if (sp == "hs") {
    xnum <- 23
    ynum <- 24
  } else if (sp == "mm") {
    xnum <- 20
    ynum <- 21
  } else if (sp == "rn") {
    xnum <- 21
    ynum <- 22
  } else stop( "Unkwnon species ! Are supported : hs, mm or rn.")
  
  gv <- cs$gv
  
  ## Recreating samples names from input gcx.bz2 filenames if needed
  if ((is.null(sample1)) | (is.na(sample1))) sample1 <- unlist(strsplit(gcx1, "_GC[0-9]{1}", perl=T))[1]
  if ((is.null(sample2)) | (is.na(sample2))) sample2 <- unlist(strsplit(gcx2, "_GC[0-9]{1}", perl=T))[1]
  
  samplename <- paste(sample1, "_vs_", sample2, sep="")
  cat(samplename,"\n")
  
  ## Importing data from GCX datatables
  extens <- rev(unlist(strsplit(gcx1, "_(hg|mm)[0-9]+", perl=T)))[1]
  if (extens == ".gcx.bz2") {
    gcx1.df <- read.table(bzfile(gcx1), header=T, sep="\t", check.names=F, stringsAsFactors=F)
    gcx2.df <- read.table(bzfile(gcx2), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  } else if (extens == ".gcx.gz") {
    gcx1.df <- read.table(gzfile(gcx1), header=T, sep="\t", check.names=F, stringsAsFactors=F)
    gcx2.df <- read.table(gzfile(gcx2), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  } else if (extens == ".gcx.xz") {
    gcx1.df <- read.table(xzfile(gcx1), header=T, sep="\t", check.names=F, stringsAsFactors=F)
    gcx2.df <- read.table(xzfile(gcx2), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  } else if (extens == ".gcx") {
    gcx1.df <- read.table(gcx1, header=T, sep="\t", check.names=F, stringsAsFactors=F)
    gcx2.df <- read.table(gcx2, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  } else {
    stop("Input files are not .gcx, nor .gcx.bz2, nor .gcx.gz, nor gcx.xz.")
  }
  
  ## Crossing tables based on common ProbeName
  gcx1.df <- gcx1.df[gcx1.df[["ProbeName"]] %in% gcx2.df[["ProbeName"]],]
  gcx2.df <- gcx2.df[gcx2.df[["ProbeName"]] %in% gcx1.df[["ProbeName"]],]
  
  ## Converting log10(ratio) to log2(ratio)
  if(in.type == "log10") {
    cat("Converting log10 to log2...\n")
    gcx1.df[["LogRatio"]] <- log2(10^gcx1.df[["LogRatio"]])
    gcx2.df[["LogRatio"]] <- log2(10^gcx2.df[["LogRatio"]])
  }
  
  ## Combining both dataframes
  gcx.df <- cbind(gcx1.df, sample2=gcx2.df[["LogRatio"]])
  colnames(gcx.df) <- c("ProbeName", "Chr", "Start", "End", sample1, sample2)
  
  ## Removing chrY
  gcx.df <- gcx.df[which(gcx.df[["Chr"]] != ynum),]
  cs$cur.glen <- cs$cur.glen - cs$lchrxx[ynum]
  
  ## Optional removing of chrX
  if (noX) {
    gcx.df <- gcx.df[which(gcx.df[["Chr"]] != xnum),]
    cs$cur.glen <- cs$cur.glen - cs$lchrxx[xnum]
  }
  
  ## Computing differences
  gcx.df[["Diff"]] <- gcx.df[[sample1]] - gcx.df[[sample2]]
  gcx.df[["AbsDiff"]] <- abs(gcx.df[["Diff"]])
  
  ## Converting chromosomal coordinates' starts to genomic ones
  genstart <- gcx.df[["Start"]] + cs$lchrsum[gcx.df[["Chr"]]]
  
  ## Calculating the smoothing window (1/500th of the total probes)
  smo <- as.integer(nrow(gcx.df)/500)
  
  rm.df <- gcx.df
  rm.df[[sample1]] <- as.numeric(runmed(gcx.df[[sample1]], smo))
  rm.df[[sample2]] <- as.numeric(runmed(gcx.df[[sample2]], smo))
  rm.df[["Diff"]] <- as.numeric(runmed(gcx.df[["Diff"]], smo))
  rm.df[["AbsDiff"]] <- as.numeric(runmed(abs(gcx.df[["Diff"]]), smo))
  
  
  ## Computing the correlation score
  gcx.cor <- cor(rm.df[[sample1]], rm.df[[sample2]], method="spearman")
  
  if (gcx.cor > cor.cut) {
    
    ## Determining which profile has the higher or lower dynamics, based on their runmed IQR
    compqdif <- abs(IQR(rm.df[[sample1]])) - abs(IQR(rm.df[[sample2]]))
    
    h.sample <- ''
    l.sample <- ''
    
    if (compqdif > 0) {
      h.sample <- sample1
      l.sample <- sample2
      cat(sample1,"is high,",sample2,"is low.\n")
    } else if (compqdif < 0) {
      h.sample <- sample2
      l.sample <- sample1
      cat(sample2,"is high,",sample1,"is low.\n")
    }
    
    ## Computing the AbsDiff density
    den.ad <- density(abs(rm.df[[sample1]] - rm.df[[sample2]]))	## Difference on runmeds, not on real profiles. Less accurate, more efficient !
    
    ## Getting the second local minimum AFTER THE HIGHEST PEAK as a threshold
    peak.x <- which.max(den.ad$y)
    sloperide <- rle(sign(diff(den.ad$y[which(den.ad$x >= den.ad$x[peak.x])])))
    
    sloperide$deni <- rep(NA, length(sloperide$lengths))
    for (k in 1:length(sloperide$deni)) {
      sloperide$deni[k] <- sum(sloperide$lengths[1:k]) + peak.x
    }
    sloperide$denx <- den.ad$x[sloperide$deni]
    sloperide$deny <- den.ad$y[sloperide$deni]
    
    if (length(sloperide$lengths) >= 4) {
      sr.peak2 <- which.max(sloperide$deny)
      xcut <- sloperide$denx[sr.peak2+1]
    } else {
      stop("couldn't find the 2nd local minimum")
    }
    
    ## Performing the fit
    sel.prob <- which(abs(rm.df[[sample1]] - rm.df[[sample2]]) < xcut)
    rmfit.coefs <- lsfit(rm.df[[l.sample]][sel.prob], rm.df[[h.sample]][sel.prob])$coefficients
    
    ## Creating the fitted profile, the Diff and AbsDiff
    if (rmfit.coefs[2] < 1) rmfit.coefs[2] <- 1/rmfit.coefs[2]
    adj.gcx.df <- gcx.df
    adj.gcx.df[[l.sample]] <- gcx.df[[l.sample]] * rmfit.coefs[2] + rmfit.coefs[1]
    adj.gcx.df[["Diff"]] <- adj.gcx.df[[sample1]] - adj.gcx.df[[sample2]]
    adj.rm.df <- adj.gcx.df
    adj.rm.df[[sample1]] <- as.numeric(runmed(adj.gcx.df[[sample1]], smo))
    adj.rm.df[[sample2]] <- as.numeric(runmed(adj.gcx.df[[sample2]], smo))
    adj.rm.df[["Diff"]] <- as.numeric(runmed(adj.gcx.df[["Diff"]], smo))
    adj.rm.df[["AbsDiff"]] <- as.numeric(runmed(abs(adj.gcx.df[["Diff"]]), smo))
  }
  else {
    adj.gcx.df <- gcx.df
    adj.rm.df <- rm.df
    rmfit.coefs <- c(NA, NA)
    cat("Skipping fit :", round(gcx.cor, digits=3), "\n")
  }
  
  ## Segmenting the differences profile
  SEG.df <- cbind(adj.gcx.df[,c(1:3)], adj.gcx.df[["Diff"]])
  colnames(SEG.df) <- c("Clone", "Chromosome", "Position", "Log2Rat")
  SEG.CNA <- CNA(as.matrix(SEG.df$Log2Rat), SEG.df$Chromosome, SEG.df$Position, data.type = "logratio", sampleid = "Diff")
  if (uSD == 0) { cbs.SEG.CNA <- segment(SEG.CNA, min.width=3, nperm=20000) } else { cbs.SEG.CNA <- segment(SEG.CNA, undo.splits = "sdundo", undo.SD = uSD, min.width=3, nperm=20000) }
  
  ## Getting real segments ends
  ## Probes index for start and end
  cbs.SEG.CNA$output$mark.start.i <- c(1, sapply(c(1:(nrow(cbs.SEG.CNA$output)-1)), function(x) { (sum(cbs.SEG.CNA$output$num.mark[1:x])+1) }))
  cbs.SEG.CNA$output$mark.end.i <- c(cbs.SEG.CNA$output$mark.start.i[2:length(cbs.SEG.CNA$output$mark.start.i)]-1, sum(cbs.SEG.CNA$output$num.mark))
  
  ## Getting new ends from gcdata, thanks to this index
  cbs.SEG.CNA$output$loc.end <- adj.gcx.df$End[cbs.SEG.CNA$output$mark.end.i]
  cbs.SEG.CNA$output$seg.med <- apply(cbs.SEG.CNA$output, 1, function(x) { median(adj.gcx.df$Diff[x[7]:x[8]])})
  
  ## Creating the output segment table with bp coordinates back
  cbs.SEG.out <- cbs.SEG.CNA$output
  cbs.SEG.out$loc.start <- cbs.SEG.out$loc.start
  cbs.SEG.out$loc.end <- cbs.SEG.out$loc.end
  cbs.SEG.out$genstart <- cbs.SEG.out$loc.start + cs$lchrsum[cbs.SEG.out$chrom]
  cbs.SEG.out$genend <- cbs.SEG.out$loc.end + cs$lchrsum[cbs.SEG.out$chrom]
  
  ## Defining the differences significances by the MedAbsDiff criterion used for classical CGH profiles in GC5
  noiseMAD <- median(abs(diff(adj.gcx.df[["Diff"]])))
  adj.xcut <- noiseMAD * nrf
  
  ## Temporary assessing normal, "gained" or "lost" segments
  temp.norm.index <- which( (cbs.SEG.out$seg.med > -adj.xcut) & (cbs.SEG.out$seg.med < adj.xcut) )
  
  ## Centralisation finale sur la médiane des valeurs de segments considérés normaux
  centz <- median(cbs.SEG.out$seg.med[temp.norm.index])
  adj.gcx.df[["Centered.Diff"]] <- adj.gcx.df[["Diff"]] - centz
  cbs.SEG.out$seg.med <- cbs.SEG.out$seg.med - centz
  cbs.gain.index <- which(cbs.SEG.out$seg.med >= adj.xcut)
  cbs.loss.index <- which(cbs.SEG.out$seg.med <= -adj.xcut)
  cbs.norm.index <- which( (cbs.SEG.out$seg.med > -adj.xcut) & (cbs.SEG.out$seg.med < adj.xcut) )
  
  ## Selecting the probes inside the different genomic regions
  if(length(cbs.gain.index) > 0) gain.probes <- sort(unlist(apply(cbs.SEG.out[cbs.gain.index,], 1, function(x) { which( (adj.gcx.df[["Chr"]] == as.numeric(x[["chrom"]])) & (adj.gcx.df[["Start"]] <= as.numeric(x[["loc.end"]])) & (adj.gcx.df[["End"]] >= as.numeric(x[["loc.start"]])) ) })))
  if(length(cbs.loss.index) > 0) loss.probes <- sort(unlist(apply(cbs.SEG.out[cbs.loss.index,], 1, function(x) { which( (adj.gcx.df[["Chr"]] == as.numeric(x[["chrom"]])) & (adj.gcx.df[["Start"]] <= as.numeric(x[["loc.end"]])) & (adj.gcx.df[["End"]] >= as.numeric(x[["loc.start"]])) ) })))
  
  ## FINAL PLOTS
  ## Initial, unscaled
  ymax <- 1.5
  png(paste(samplename, "_", gv, "_U", uSD, "_nrf", nrf, "_CC", cor.cut, ".png", sep=""), width=1680, height=1050)
  par(mgp=c(1,0,0), mar=c(2,2,2,2), xaxs="i")
  par(mfrow=c(3,1))
  plot(genstart, rm.df[[sample1]], type="l", xlim=c(1,cs$cur.glen), ylim=c(-ymax,ymax), col=4, main=paste(sample1, " (blue) vs ", sample2, " (green) raw profiles. Cor = ", round(gcx.cor, digits=3), sep=""), xlab="Genomic position", ylab="log2(ratio)", cex.main=2)
  lines(genstart, rm.df[[sample2]], col=3)
  abline(h=0, col="grey50", lty=2)
  abline(v=cs$lchrsum, lty=2, col=1)
  text(cs$lchrsum+cs$lchrxx/2, ymax*(-2*(unique(cs$cytob$chr)%%2)+1), unique(cs$cytob$chrA), cex=2)
  
  ## Final, scaled (or not!) with colored differences
  if (gcx.cor > cor.cut) scaleword <- "scaled" else scaleword <- "unscaled"
  plot(0, 0, type="n", xlim=c(1,cs$cur.glen), ylim=c(-1.5,1.5), main=paste(sample1, " (blue) vs ", sample2, " (green) ", scaleword, " profiles", sep=""), xlab="Genomic position", ylab="log2(ratio)", cex.main=2)
  if (length(cbs.gain.index) > 0) segments(genstart[gain.probes], adj.rm.df[[sample1]][gain.probes], genstart[gain.probes], adj.rm.df[[sample2]][gain.probes], col=2)
  if (length(cbs.loss.index) > 0) segments(genstart[loss.probes], adj.rm.df[[sample1]][loss.probes], genstart[loss.probes], adj.rm.df[[sample2]][loss.probes], col=2)
  lines(genstart, adj.rm.df[[sample1]], col=4)
  lines(genstart, adj.rm.df[[sample2]], col=3)
  abline(h=0, col="grey50", lty=2)
  if (length(cbs.gain.index) > 0) segments(cbs.SEG.out$genstart[cbs.gain.index], 1.5, cbs.SEG.out$genend[cbs.gain.index], 1.5, col=4, lwd=4)
  if (length(cbs.loss.index) > 0) segments(cbs.SEG.out$genstart[cbs.loss.index], -1.5, cbs.SEG.out$genend[cbs.loss.index], -1.5, col=2, lwd=4)
  abline(v=cs$lchrsum, lty=2, col=1)
  text(cs$lchrsum+cs$lchrxx/2, ymax*(-2*(unique(cs$cytob$chr)%%2)+1), unique(cs$cytob$chrA), cex=2)
  
  ## Differences profile
  plot(0, 0, type="n", xlim=c(1,cs$cur.glen), ylim=c(-1.5,1.5), main="Significant differences between profiles", xlab="Genomic position", ylab="log2(ratio)", cex.main=2)
  points(genstart, adj.gcx.df[["Diff"]]-centz, pch=20, cex=0.5, col="grey75")
  lines(genstart, adj.rm.df[["Diff"]]-centz, col=4)
  if (length(cbs.norm.index) > 0) segments(cbs.SEG.out$genstart[cbs.norm.index], cbs.SEG.out$seg.med[cbs.norm.index], cbs.SEG.out$genend[cbs.norm.index], cbs.SEG.out$seg.med[cbs.norm.index], col=1, lwd=4)
  if (length(cbs.gain.index) > 0) segments(cbs.SEG.out$genstart[cbs.gain.index], cbs.SEG.out$seg.med[cbs.gain.index], cbs.SEG.out$genend[cbs.gain.index], cbs.SEG.out$seg.med[cbs.gain.index], col=4, lwd=4)
  if (length(cbs.loss.index) > 0) segments(cbs.SEG.out$genstart[cbs.loss.index], cbs.SEG.out$seg.med[cbs.loss.index], cbs.SEG.out$genend[cbs.loss.index], cbs.SEG.out$seg.med[cbs.loss.index], col=2, lwd=4)
  abline(h=0, col="grey25", lty=2)
  abline(h=c(-adj.xcut, adj.xcut), col="grey25", lty=2)
  abline(v=cs$lchrsum, lty=2, col=1)
  text(cs$lchrsum+cs$lchrxx/2, ymax*(-2*(unique(cs$cytob$chr)%%2)+1), unique(cs$cytob$chrA), cex=2)
  dev.off()
  
  ## Dumping the output GCX for the differential profile
  out.gcx <- adj.gcx.df[,c(1:4,7)]
  colnames(out.gcx) <- c("ProbeName", "Chr", "Start", "End", "LogRatio")
  write.table(out.gcx, paste(samplename, "_", gv, "_U", uSD, "_nrf", nrf, "_CC", cor.cut, ".gcx", sep=""), sep="\t", quote=F, row.names=F)
  
  ## Dumping out CBS (table of segmented differences, with a significance column)
  cbs.CBS.out <- cbind(samplename, cbs.SEG.out[,c(2:5,9)], stringsAsFactors=F)
  colnames(cbs.CBS.out) <- c(samplename, "Chr", "Start", "End", "Probes", "Log2Ratio")
  cbs.CBS.out[["Log2Ratio"]][cbs.norm.index] <- 0
  write.table(cbs.CBS.out, paste(samplename, "_", gv, "_nrf", nrf, "_uSD", uSD, ".cbs", sep=""), sep="\t", quote=F, row.names=F)
  
  ## Dumping out P2C (table of segmented differences, with a significance column)
  cbs.SEG.out <- cbind(cbs.SEG.out[,2], samplename, cbs.SEG.out[,c(2:5,9)], stringsAsFactors=F)
  colnames(cbs.SEG.out) <- c("Loc", "Comparison", "Chr", "Start", "End", "Probes", "Log2RatioDiff")
  cbs.SEG.out$Loc[which(cbs.SEG.out$Loc == 23)] <- "X"
  cbs.SEG.out$Loc[which(cbs.SEG.out$Loc == 24)] <- "Y"
  cbs.SEG.out$Loc <- paste("chr", cbs.SEG.out$Loc, ":", cbs.SEG.out$Start, "-", cbs.SEG.out$End, sep="")
  cbs.SEG.out[["SigDif"]] <- 0
  cbs.SEG.out[["SigDif"]][cbs.gain.index] <- 1
  cbs.SEG.out[["SigDif"]][cbs.loss.index] <- -1
  write.table(cbs.SEG.out, paste(samplename, "_", gv, "_nrf", nrf, "_uSD", uSD, ".p2c", sep=""), sep="\t", quote=F, row.names=F)
  
  ## Computing and writing the sum of different genomic lengthes
  less.genlen <- more.genlen <- 0
  if(length(cbs.gain.index) > 0) more.genlen <- sum(unlist(apply(cbs.SEG.out[cbs.gain.index,], 1, function(x) { as.numeric(x[["End"]]) - as.numeric(x[["Start"]]) }))) / cs$cur.glen
  if(length(cbs.loss.index) > 0) less.genlen <- sum(unlist(apply(cbs.SEG.out[cbs.loss.index,], 1, function(x) { as.numeric(x[["End"]]) - as.numeric(x[["Start"]]) }))) / cs$cur.glen
  
  return(c(sample1, sample2, uSD, nrf, noX, cor.cut, cor.method, gcx.cor, noiseMAD, more.genlen, less.genlen, more.genlen+less.genlen, rmfit.coefs[2], rmfit.coefs[1]))
}
