## clicoM, the Clinical Comparator
##
## DESCRIPTION:
## The aim of the script is to perform automaticaly the search for differential genomic regions between defined
## subpopulations of aCGH profiles, based on their aberrations. Currently, the statistical tests performed are a
## T-test for a 2-classes comparison, and an ANOVA for N-classes. Later, other tests should be added, like Wilcoxon,
## and limma tests, among others.
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
##
## DEPENDS ON:
## . r-base only.
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.8.2 to 3.3.0
##
## VERSION NOTES:
##
## v2.8   20161019
##      . Now an IGV XML session will be saved when significant regions will be found (in pcut mode)
##        to ease data loading into IGV browser.
##
## v2.7b  20160930
##      . Ameliorated output plot (PNG) : now last track reflects -log10(p) much more accurately.
##      . Also small blue dots on the last track in pcut mode.
##      . Now the pangenomic plots show chromosome names (1:n,X,Y) instead of number.
##
## v2.7   20160930
##      . Added violin plots for all significant regions in pcut mode.
##      . Now the -log10(pvalue) plots show the requested pcut value instead of default 5E-02 (still default)
##      . Corrected the dashed lines for p-value thresholds which were wrong for all .5 value at all orders.
##
## v2.6   20160728
##      . Added pcut mode that generates a bed and bedgraph, then runs cc2grd, for a given p-value threshold.
##      . Replaced the obsolete 'mnt.proj' option to a more generic 'ldb' specifying a path to a
##        local mirror of public databases (by default '/mnt/data_cigogne/bioinfo/').
##
## v2.5   20160609
##      . Added support for continuous data (like, patients age), using linear regression or correlation tests.
##      . Added multithreading !
##
## v2.4   20140306
##      . Corrected an error in the help function
##      . Added support for newer cytoBandsIdeo table from UCSC
##
## v2.3c  20111116
##      . Introduced a "help" function to display the description of parameters for the main function, even when
##        this script is sourced and not directly read.
##
## v2.3b  20111005
##      . Added a filter to replace dots in colnames with underscores. Needed for genrdesc to work well in -cc mode.
##
## v2.3   20110831
##      . Added support for mus musculus.
##      . In this purpose, dropped the -hg option for two new options : -sp for species, and -gb for the genome build.
##
## v2.2  20110623
##      . Corrected some bugs (handling of input files).
##      . Changed the colors to the new recommanded red/blue set.
##      . Added 4 new options to control the colours of gains, losses, amps and dels on the plots
##
## v2.1b 20110430
##      . Modified structure for generated files : A folder is automaticaly created for each analyzed annotation
##        according to the colname, the test type and a timestamp, then files are generated
##        within this folder.
##
## v2.1
##      . Added an option to specify the mount point of pelican /proj.
##      . Wraped the script as a function so that it can be directly sourced from anywhere.
##
## v2.0
##      . Script modified to add support for REG/XREG formats. Support of ABRS3 is abandonned.
##
## v1.0
##      . First release.
##

## PAREMETERS DESCRIPTOR FUNCTION
help.clicom <- function() {
  cat("
      PARAMETERS FOR THE clicom() FUNCTION :
      ----------------------------------------
      reg.table   : path+filename to the REG/XREG table
      annot.table : path+filename to the annotation table
      colvec      : a vector of integers corresponding to the index of columns from the annotation table that must be analyzed.
      test.type.2 : type of statistical test for 2-classes variable (W : Wilcoxon, or T : Student's T-test).
      test.type.N : type of statistical test for N-classes (ANOVA : Analysis of variance,  or KW : Kruskal-Wallis).
      test.type.continuous : type of statistical test for continuous variable (COR.P : Pearson correlation, or COR.S : Spearman correlation).
      numeric.as.continuous : Treat numeric columns (integer, double) as continuous data. If FALSE, these columns are converted to factors.
      gaincol     : color to use to plot gains.
      losscol     : color to use to plot losses.
      ampcol      : color to use to plot amplifications.
      delcol      : color to use to plot deletions.
      sp          : Species (hs, mm).
      gb          : genome build. (18, 19 for homo sapiens, 9 for mus musculus).
      amplim      : the log2ratio cut-off to plot amplifications' and deletions' frequencies in a darker color.
      ldb         : Path to a local mirror of public databases (here for the GoldenPath tables).
      pcut        : If non-NULL, subtract regions with raw p-value lower than pcut.
      grd         : If Call 'grd' on the pcut-selected regions (if any). Requires the 'grd' in your PATH, and grd set to pcut set to non-NULL.
      nthread      : number of threads to use.
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
clicom <-function(reg.table, annot.table, colvec, test.type.2="W", test.type.N="KW", test.type.continuous="COR.S", numeric.as.continuous = TRUE, gaincol="blue", losscol="red", ampcol="cyan4", delcol="orangered4", sp="hs", gb=19, amplim=1, ldb = "/mnt/data_cigogne/bioinfo/", width = 1440, height = 920, pcut = 1.0E-03, grd = TRUE, out.dir = getwd(), nthread = 1) {
  
  if (!is.null(pcut)) {
    message("PCUT mode activated!")
    if (!is.numeric(pcut)) stop("pcut must be a numeric !")
    if ( (pcut > 1) | (pcut <= 0) )stop("pcut must be in the ]0,1[ interval !")
  }
  
  message("Reading data table ...")
  datatable <- read.table(reg.table, header=T, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  message("Reading annotation table ...")
  annot <- read.table(annot.table, header=T, sep="\t", stringsAsFactors=TRUE, check.names=FALSE, comment.char = "", quote = "")
  
  ## removing dots
  message("Cleaning ...")
  # colnames(annot) <- chartr(".", "_", colnames(annot))
  
  gv <- ""
  if (sp == "hs") {
    gv = "hg"
    xnum = 23
    ynum = 24
  } else if (sp == "mm") {
    gv = "mm"
    xnum = 20
    ynum = 21
  }
  gv <- paste(gv, gb, sep="")
  
  ## IMPORT BANDES CYTO ET PREPARATION POUR PLOT KARYO
  cat("Importing chromosomes data  ...\n")
  # cytob <- read.table(paste(ldb, "/GoldenPath/", gv, "/cytoBandIdeo.", gv, sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_", fill=T)
  cytob <- read.table(paste(ldb, "/GoldenPath/", gv, "/cytoBandIdeo.", gv, sep=""), sep="\t", header=F, stringsAsFactors=F, comment.char="_", fill=T)
  colnames(cytob) <- c('X.chrom', 'chromStart', 'chromEnd', 'cytoband', 'gieStain')
  
  ## Filtering out invalid lines (UCSC recently added incomplete lines which seem to refer to chrN_unknown sequence sets, which is quite totally illogical. BJ 20140227)
  cytob <- cytob[which(!is.na(cytob$chromStart)),]
  ## Filtering out ANY chromosome not defined in chrconvlist$chrom2chr
  cytob <- cytob[which(cytob[["X.chrom"]] %in% names(chrconv.list$chrom2chr[[sp]])),]
  
  cytob$chrA <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
  cytob$chr <- cytob$chrA
  cytob$chr[which(cytob$chr == "X")] <- xnum
  cytob$chr[which(cytob$chr == "Y")] <- ynum
  cytob$chr <- as.numeric(cytob$chr)
  cytob <- cytob[order(cytob$chr),]
  
  lchrxx <- sapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) })
  lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
  glen <- sum(as.numeric(lchrxx))
  lchrsum <- sapply(c(1:length(lchrtoadd)), function(x) { sum(lchrtoadd[1:x]) })
  
  cghdata <- as.matrix(datatable[, -c(1:9)])
  
  require(foreach)
  require(doParallel)
  cl <- makeCluster(spec = nthread, type = "PSOCK")
  registerDoParallel(cl)
  foreach (colsel = colvec) %dopar% {
    
    print(paste("Working on column : ",colnames(annot)[colsel]," ...", sep=""))
    
    outable <- datatable[,1:9]
    outable$GStart <- outable$Start + lchrsum[outable$Chr]
    outable$GEnd <- outable$End + lchrsum[outable$Chr]
    colz <- colnames(outable)
    
    test.name <- "NULL"
    ## Handling contiuous variable case
    if (numeric.as.continuous & (is.numeric(annot[,colsel]))) {
      
      print(" Continuous data!")
      test.name <- test.type.continuous
      testype <- "conti"
      
      require(matrixStats)
      outable$MedianL2R <- rowMedians(cghdata, na.rm = TRUE)
      
      ## PEARSON correlation test
      if (test.type.continuous == "COR.P") {
        froot <- paste0(colnames(annot)[colsel], "_COR.P")
        outable$Cor.P <- NA
        outable$RawP <- NA
        for (k in 1:nrow(cghdata)) {
          ct.res <- cor.test(cghdata[k,], annot[,colsel], method = "pearson")
          outable$Cor.P[k] <- ct.res$estimate
          outable$RawP[k] <- ct.res$p.value
        }
      }
      
      ## SPEARMAN correlation test
      if (test.type.continuous == "COR.S") {
        froot <- paste0(colnames(annot)[colsel], "_COR.S")
        outable$Cor.S <- NA
        outable$RawP <- NA
        for (k in 1:nrow(cghdata)) {
          ct.res <- cor.test(cghdata[k,], annot[,colsel], method = "spearman")
          outable$Cor.S[k] <- ct.res$estimate
          outable$RawP[k] <- ct.res$p.value
        }
      }
      
      outable$AdjP.BH <- p.adjust(outable$RawP, method="BH")
      
      ### Plot
      q95cn <- rowQuantiles(cghdata, probs = .95, na.rm = TRUE)
      q05cn <- rowQuantiles(cghdata, probs = .05, na.rm = TRUE)
      q75cn <- rowQuantiles(cghdata, probs = .75, na.rm = TRUE)
      q25cn <- rowQuantiles(cghdata, probs = .25, na.rm = TRUE)
      globmedian <- rowMedians(cghdata, na.rm = TRUE)
      
      g9top <- g9low <- q95cn
      g9top[g9top < 0] <- 0
      g9low[g9low > 0] <- 0
      l9top <- l9low <- q05cn
      l9top[l9top < 0] <- 0
      l9low[l9low > 0] <- 0
      
      g7top <- g7low <- q75cn
      g7top[g7top < 0] <- 0
      g7low[g7low > 0] <- 0
      l7top <- l7low <- q25cn
      l7top[l7top < 0] <- 0
      l7low[l7low > 0] <- 0
      
      dir.create(paste0(out.dir, '/', froot), mode="0775")
      
      xmlcon <- file(paste0(out.dir, '/', froot, "/", froot, "_", pcut, "_igv_session.xml"), open = "w")
      writeLines(
        text = c(
          '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
          paste0('<Session genome="', gv, '" hasGeneTrack="true" hasSequenceTrack="true" locus="All" path="igv_session.xml" version="8">'),
          '<Resources>'
        ), con = xmlcon
      )
        
      png(paste0(out.dir, '/', froot, "/", froot, ".png"), width=width, height=height)
      par(mfrow=c(3,1), mgp=c(1,0,0), mar=c(2,2,3,2), xaxs="i")
      
      lidx <- which(globmedian < 0)
      gidx <- which(globmedian > 0)

      my.ylim <- c(-1,1)
      if (length(lidx)>0) my.ylim[1] <- median(globmedian[lidx], na.rm = TRUE) * 2
      if (length(gidx)>0) my.ylim[2] <- median(globmedian[gidx], na.rm = TRUE) * 2

      plot(outable$GStart, globmedian, type = "n", ylim=my.ylim, main=paste0("Whole population median L2R profile (", ncol(cghdata), ")"), xlab="Genomic position", ylab="Median profile", xaxt="n", cex.main=2)
      if (length(gidx)>0) {
        rect(outable$GStart[gidx], 0, outable$GEnd[gidx], globmedian[gidx], col=adjustcolor(gaincol, alpha.f = .25), lwd=1, border = adjustcolor(gaincol, alpha.f = .25))
        rect(outable$GStart[gidx], q75cn[gidx], outable$GEnd[gidx], q25cn[gidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
        rect(outable$GStart[gidx], q95cn[gidx], outable$GEnd[gidx], q05cn[gidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
      }
      if (length(lidx)>0) {
        rect(outable$GStart[lidx], 0, outable$GEnd[lidx], globmedian[lidx], col=adjustcolor(losscol, alpha.f = .25), lwd=1, border = adjustcolor(losscol, alpha.f = .25))
        rect(outable$GStart[lidx], q75cn[lidx], outable$GEnd[lidx], q25cn[lidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
        rect(outable$GStart[lidx], q95cn[lidx], outable$GEnd[lidx], q05cn[lidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
      }
      abline(h = 0)

      text(lchrsum+lchrxx/2, my.ylim[1], unique(cytob$chrA), cex=2)
      abline(v=lchrsum)
      
      mbgdf <- data.frame(chrom = outable$Chrom, start = outable$Start, end = outable$End, value = outable$MedianL2R)
      BGhead <- paste0("track type=bedGraph name=\"", froot," Median L2R\" description=\"", froot, " Median log2(ratio)\" color=0,0,255 altColor=255,0,0 graphType=bar viewLimits=-0.75:0.75 autoScale=off windowingFunction=none smoothingWindow=off visibility=full")
      bgfile <- paste0(out.dir, '/', froot, "/", froot, "_MedianL2R.bedgraph")
      write.table(BGhead, bgfile, col.names=F, row.names=F, sep="\t", quote=F)
      write.table(mbgdf, bgfile, col.names=F, row.names=F, sep="\t", quote=F, append = TRUE)
      
      reslist <- list(l2r=paste0('<Resource path="', basename(bgfile), '"/>'))
      tracklist <- list(
        l2r = c(
          paste0('<Track altColor="255,0,0" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" color="0,0,255" colorScale="ContinuousColorScale;0.0;-0.75;0.0;0.75;255,0,0;255,255,255;0,0,255" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="100" id="', basename(bgfile), '" name="Population median L2R" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="none">'),
          paste0('<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="0.75" minimum="-0.75" type="LINEAR"/>'),
          paste0('</Track>')
        )
      )

      if (test.type.continuous == "COR.P") {

        ctop <- clow <- outable$Cor.P
        ctop[ctop < 0] <- 0
        clow[clow > 0] <- 0

        my.ylim <- c(-1,1)
        plot(outable$GStart, outable$Cor.P, type = "n", ylim = c(-1,1), main = paste0("Pearson correlation with ", colnames(annot)[colsel]), xlab = "Genomic position", ylab = "Pearson correlation coefficient", xaxt = "n", cex.main=1.8)
        rect(outable$GStart, 0, outable$GEnd, ctop, col = "blue", lwd=1, border = "blue")
        rect(outable$GStart, 0, outable$GEnd, clow, col = "red", lwd=1, border = "blue")
        abline(h = 0)
        text(lchrsum+lchrxx/2, my.ylim[1], unique(cytob$chrA), cex=2)
        abline(v=lchrsum)

        cbgdf <- data.frame(chrom = outable$Chrom, start = outable$Start, end = outable$End, value = outable$Cor.P)
        BGhead <- paste0("track type=bedGraph name=\"", froot," Pearson Correlation\" description=\"", froot, " Pearson Correlation\" color=0,0,255 altColor=255,0,0 graphType=bar viewLimits=-1:1 autoScale=off windowingFunction=none smoothingWindow=off visibility=full")
        bgfile <- paste0(out.dir, '/', froot, "/", froot, "_correlation.bedgraph")
        write.table(BGhead, bgfile, col.names=F, row.names=F, sep="\t", quote=F)
        write.table(cbgdf, bgfile, col.names=F, row.names=F, sep="\t", quote=F, append = TRUE)
        
        reslist <- c(reslist, cor=paste0('<Resource path="', basename(bgfile), '"/>'))
        tracklist <- c(
          tracklist,
          cor = list(
            c(
              paste0('<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" color="100,0,0" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="100" id="', basename(bgfile), '" name="', colnames(annot)[colsel], ' Corr Coeff" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="none">'),
              paste0('<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1.0" minimum="-1.0" type="LINEAR"/>'),
              paste0('</Track>')
            )
          )
        )
        
        # print(tracklist)
      }
      if (test.type.continuous == "COR.S") {
        ctop <- clow <- outable$Cor.S
        ctop[ctop < 0] <- 0
        clow[clow > 0] <- 0

        my.ylim <- c(-1,1)
        plot(outable$GStart, outable$Cor.S, type = "n", ylim = c(-1,1), main = paste0("Pearson correlation with ", colnames(annot)[colsel]), xlab = "Genomic position", ylab = "Pearson correlation coefficient", xaxt = "n", cex.main=1.8)
        rect(outable$GStart, 0, outable$GEnd, ctop, col = "blue", lwd=1, border = "blue")
        rect(outable$GStart, 0, outable$GEnd, clow, col = "red", lwd=1, border = "red")
        abline(h = 0)
        text(lchrsum+lchrxx/2, my.ylim[1], unique(cytob$chrA), cex=2)
        abline(v=lchrsum)

        cbgdf <- data.frame(chrom = outable$Chrom, start = outable$Start, end = outable$End, value = outable$Cor.S)
        BGhead <- paste0("track type=bedGraph name=\"", froot," Spearman Correlation\" description=\"", froot, " Spearman Correlation\" color=0,0,255 altColor=255,0,0 graphType=bar viewLimits=-1:1 autoScale=off windowingFunction=none smoothingWindow=off visibility=full")
        bgfile <- paste0(out.dir, '/', froot, "/", froot, "_correlation.bedgraph")
        write.table(BGhead, bgfile, col.names=F, row.names=F, sep="\t", quote=F)
        write.table(cbgdf, bgfile, col.names=F, row.names=F, sep="\t", quote=F, append = TRUE)
        
        reslist <- c(reslist, cor=paste0('<Resource path="', basename(bgfile), '"/>'))
        tracklist <- c(
          tracklist,
          cor = list( 
            c(
              paste0('<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" color="100,0,0" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="100" id="', basename(bgfile), '" name="', colnames(annot)[colsel], ' Corr Coeff" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="none">'),
              paste0('<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1.0" minimum="-1.0" type="LINEAR"/>'),
              paste0('</Track>')
            )
          )
        )
        
        # print(tracklist)
      }
      abline(h=0)
      abline(v=lchrsum)
      
      outable$MedianL2R <- NULL
      
      
    }
    
    ## Handling categorial variable case
    else {
      
      testype <- "categ"
      
      # Transforming to a factor, if needed
      if (is.numeric(annot[,colsel])) annot[,colsel] <- as.factor(annot[,colsel])
      
      ncateg <- nlevels(annot[,colsel])
      
      # 1-class : no test !
      # if (ncateg == 1) next("Single category !")
      if (ncateg == 1) return()
      
      # 2-classes : T-test or Wilcoxon
      if (ncateg == 2) testname <- paste0("_", toupper(test.type.2))
      
      # N>2 -classes : ANOVA
      if (ncateg > 2) testname <- paste0("_", toupper(test.type.N))
      
      # Creating output filenames' root
      froot <- paste0(colnames(annot)[colsel], testname)
      
      # Creating outdir
      dir.create(paste0(out.dir, '/', froot), mode="0775")
      
      xmlcon <- file(paste0(out.dir, '/', froot, "/", froot, "_", pcut, "_igv_session.xml"), open = "w")
      writeLines(
        text = c(
          '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
          paste0('<Session genome="', gv, '" hasGeneTrack="true" hasSequenceTrack="true" locus="All" path="igv_session.xml" version="8">'),
          '<Resources>'
        ), con = xmlcon
      )
      reslist <- tracklist <- list()
      
      png(paste0(out.dir, '/', froot, "/", froot, ".png"), width=width, height=height)
      par(mfrow=c(ncateg+1,1), mgp=c(1,0,0), mar=c(2,2,3,2), xaxs="i")
      
      # Computing median profiles for each class
      for (n in 1:ncateg) {
        
        catname <- paste(colnames(annot)[colsel], levels(annot[,colsel])[n], sep=".")
        print(paste(catname, " ...", sep=""))
        
        nlist <- which(as.numeric(annot[,colsel]) == n)
        medcn <- paste0(catname, ".MedianL2R")
        require(matrixStats)
        outable[[medcn]] <- rowMedians(cghdata[,nlist], na.rm = TRUE)
        q95cn <- rowQuantiles(cghdata[,nlist], probs = .95, na.rm = TRUE)
        q05cn <- rowQuantiles(cghdata[,nlist], probs = .05, na.rm = TRUE)
        q75cn <- rowQuantiles(cghdata[,nlist], probs = .75, na.rm = TRUE)
        q25cn <- rowQuantiles(cghdata[,nlist], probs = .25, na.rm = TRUE)
        
        g9top <- g9low <- q95cn
        g9top[g9top < 0] <- 0
        g9low[g9low > 0] <- 0
        l9top <- l9low <- q05cn
        l9top[l9top < 0] <- 0
        l9low[l9low > 0] <- 0
        
        g7top <- g7low <- q75cn
        g7top[g7top < 0] <- 0
        g7low[g7low > 0] <- 0
        l7top <- l7low <- q25cn
        l7top[l7top < 0] <- 0
        l7low[l7low > 0] <- 0
        
        lidx <- which(outable[[medcn]] < 0)
        gidx <- which(outable[[medcn]] > 0)
        
        my.ylim <- c(-1,1)
        if (length(lidx)>0) my.ylim[1] <- median(outable[[medcn]][lidx], na.rm = TRUE) * 2
        if (length(gidx)>0) my.ylim[2] <- median(outable[[medcn]][gidx], na.rm = TRUE) * 2
        
        plot(outable$GStart, outable[[medcn]], type = "n", ylim=my.ylim, main=paste(catname, " median L2R profile (",length(nlist), ")", sep=""), xlab="Genomic position", ylab="Median profile", xaxt="n", cex.main=2)
        if (length(gidx)>0) {
          rect(outable$GStart[gidx], 0, outable$GEnd[gidx], outable[[medcn]][gidx], col=adjustcolor(gaincol, alpha.f = .25), lwd=1, border = adjustcolor(gaincol, alpha.f = .25))
          rect(outable$GStart[gidx], q75cn[gidx], outable$GEnd[gidx], q25cn[gidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
          rect(outable$GStart[gidx], q95cn[gidx], outable$GEnd[gidx], q05cn[gidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
        }
        if (length(lidx)>0) {
          rect(outable$GStart[lidx], 0, outable$GEnd[lidx], outable[[medcn]][lidx], col=adjustcolor(losscol, alpha.f = .25), lwd=1, border = adjustcolor(losscol, alpha.f = .25))
          rect(outable$GStart[lidx], q75cn[lidx], outable$GEnd[lidx], q25cn[lidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
          rect(outable$GStart[lidx], q95cn[lidx], outable$GEnd[lidx], q05cn[lidx], col=adjustcolor("black", alpha.f = .25), lwd=0)
        }
        abline(h = 0)
        
        text(lchrsum+lchrxx/2, my.ylim[1], unique(cytob$chrA), cex=2)
        abline(v=lchrsum)
        
        cbgdf <- data.frame(chrom = outable$Chrom, start = outable$Start, end = outable$End, value = outable[[medcn]])
        BGhead <- paste0("track type=bedGraph name=\"", catname," median L2R\" description=\"", catname, " median log2ratio\" color=0,0,255, altColor=255,0,0 graphType=bar viewLimits=-0.75:0.75 autoScale=off windowingFunction=none smoothingWindow=off visibility=full")
        bgfile <- paste0(out.dir, '/', froot, "/", catname, "_medianL2R.bedgraph")
        write.table(BGhead, bgfile, col.names=F, row.names=F, sep="\t", quote=F)
        write.table(cbgdf, bgfile, col.names=F, row.names=F, sep="\t", quote=F, append = TRUE)
        
        reslist <- c(reslist, l2r = paste0('<Resource path="', basename(bgfile), '"/>'))
        tracklist <- c(
          tracklist,
          l2r = list(
            c(
              paste0('<Track altColor="255,0,0" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" color="0,0,255" colorScale="ContinuousColorScale;0.0;-0.75;0.0;0.75;255,0,0;255,255,255;0,0,255" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="100" id="', basename(bgfile), '" name="', catname, ' median L2R" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="none">'),
              paste0('<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="0.75" minimum="-0.75" type="LINEAR"/>'),
              paste0('</Track>')
            )
          )
        )
        colz <- colnames(outable)
      }
      
      ## T-test
      if ((ncateg == 2) & (test.type.2 == "T")) {
        test.name <- test.type.2
        outable$Tscore <- NA
        outable$RawP <- NA
        for (k in 1:dim(cghdata)[1]) {
          test.T <- t.test(as.numeric(cghdata[ k, which(annot[,colsel] == levels(annot[,colsel])[1]) ]), as.numeric(cghdata[ k, which(annot[,colsel] == levels(annot[,colsel])[2]) ]))
          outable$Tscore[k] <- test.T$statistic
          outable$RawP[k] <- test.T$p.value
        }
      }
      
      ## Wilcoxon
      if ((ncateg == 2) & (test.type.2 == "W")) {
        test.name <- test.type.2
        outable$Wscore <- NA
        outable$RawP <- NA
        for (k in 1:dim(cghdata)[1]) {
          test.T <- wilcox.test(as.numeric(cghdata[ k, which(annot[,colsel] == levels(annot[,colsel])[1]) ]), as.numeric(cghdata[ k, which(annot[,colsel] == levels(annot[,colsel])[2]) ]))
          outable$Wscore[k] <- test.T$statistic
          outable$RawP[k] <- test.T$p.value
        }
      }
      
      ## ANOVA
      if ((ncateg > 2) & (test.type.N == "ANOVA")) {
        test.name <- test.type.N
        outable$Fscore <- NA
        outable$RawP <- NA
        for (k in 1:nrow(cghdata)) {
          test.lm <- lm(as.numeric(cghdata[k,]) ~ annot[,colsel])
          test.anova <- anova(test.lm)
          outable$Fscore[k] <- test.anova$"F value"[1]
          outable$RawP[k] <- test.anova$"Pr(>F)"[1]
        }
      }
      
      ## Kruskal-Wallis
      if ((ncateg > 2) & (test.type.N == "KW")) {
        test.name <- test.type.N
        outable$KWscore <- NA
        outable$RawP <- NA
        for (k in 1:dim(cghdata)[1]) {
          test.ks <- kruskal.test(as.numeric(cghdata[k,]) ~ annot[,colsel])
          outable$KWscore[k] <- test.ks$statistic
          outable$RawP[k] <- test.ks$p.value
        }
      }
      outable$AdjP.BH <- p.adjust(outable$RawP, method="BH")
    }
    
    ## Identifying significant regions
    if (!is.null(pcut)) p.ok <- which(outable$RawP < pcut)
    
    plot(outable$GStart, -log10(outable$RawP), type = "n", main="P-value", xlab="Genomic position", ylab="-log10(Pval)", xaxt="n", cex.main=1.8)
    rect(outable$GStart, 0, outable$GEnd, -log10(outable$RawP), border = "red4", col = "red4")
    rect(outable$GStart, 0, outable$GEnd, -log10(outable$AdjP.BH), border = "orangered", col = "orangered")
    if (length(p.ok) > 0) points(outable$GStart[p.ok], rep(0, length(p.ok)), col = 4, pch = 20, cex = 3)
    abline(v=lchrsum)

    abline(h=-log10(c(1 %o% 10^(-1:min(log10(outable$RawP), na.rm = TRUE)))), lty=2, col="grey25")
    abline(h=-log10(c(5 %o% 10^(-1:min(log10(outable$RawP), na.rm = TRUE)))), lty=2, col="grey75")
    abline(h=ifelse(is.null(pcut), -log10(5E-02), -log10(pcut)), lty=3, col=2)
    text(lchrsum+lchrxx/2, max(-log10(outable$RawP), na.rm = TRUE), unique(cytob$chrA), cex=2)
    dev.off()
    
    ## Dumping results table
    write.table(outable, file=paste0(out.dir, '/', froot, "/", froot, ".txt"), row.names=F, sep="\t", quote=F)
    
    if (!is.null(pcut)) {
      if (length(p.ok) > 0) {
        
        if (testype == "categ") {
          ## Violinplots
          require(ggplot2)
          pdf(paste0(out.dir, '/', froot, "/", froot, "_", pcut, "_violin.pdf"), width = 29.7/1.5/cm(1), height = 21/1.5/cm(1))
          for (myp in p.ok) {
            ymin <- min(cghdata[myp,])
            print(p <- qplot(factor(annot[,colsel]), as.numeric(cghdata[myp,]), geom = "violin", trim = FALSE, scale = "count", fill = factor(annot[,colsel]), xlab = froot, ylab = "Log2(ratio)", main = paste0(datatable$Chrom[myp], ":", datatable$Start[myp], "-", datatable$End[myp], " (", sprintf("%.2E", outable$RawP[myp]), ")")) + geom_boxplot(width=.2) + geom_jitter(width = 0) + guides(fill=FALSE) + geom_hline(yintercept = 0, linetype = "dashed"))
          }
          dev.off()
        }
        
        ## Bedgraph
        bgdf <- data.frame(chrom = outable$Chrom, start = outable$Start, end = outable$End, value = -log10(outable$RawP))
        bnam <- paste0(colnames(annot)[colsel], " (", test.name, ")")
        bnam <- gsub(pattern = "\\_", replacement = " ", x = bnam)
        BGhead <- paste0("track type=bedGraph name=\"", bnam,"\" description=\"", bnam, "\" color=100,0,0 graphType=bar autoScale=off windowingFunction=none smoothingWindow=off visibility=full")
        bgfile <- paste0(out.dir, '/', froot, "/", froot, ".bedgraph")
        write.table(BGhead, bgfile, col.names=F, row.names=F, sep="\t", quote=F)
        write.table(bgdf, bgfile, col.names=F, row.names=F, sep="\t", quote=F, append = TRUE)
        
        reslist <- c(reslist, pval = paste0('<Resource path="', basename(bgfile), '"/>'))
        tracklist <- c(
          tracklist,
          pval = list(
            c(
              paste0('<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" color="100,0,0" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="100" id="', basename(bgfile), '" name="', bnam, ' -log10(p)" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="none">'),
              paste0('<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="', max(bgdf$value, na.rm = TRUE), '" minimum="0.0" type="LINEAR"/>'),
              paste0('</Track>')
            )
          )
        )

        ## Bed
        subdf <- outable[p.ok,]
        bdf <- data.frame(chrom = subdf$Chrom, start = subdf$Start, end = subdf$End, name = paste0("Sig", 1:nrow(subdf)), value = subdf$RawP)
        Bhead <- paste0("track name=\"", bnam, " p ", pcut, "\" description=\"", bnam, " p ", pcut, "\" visibility = 1 color=255,0,0")
        bgfile <- paste0(out.dir, '/', froot, "/", froot, "_", pcut, ".bed")
        write.table(Bhead, bgfile, col.names=F, row.names=F, sep="\t", quote=F)
        write.table(bdf, bgfile, col.names=F, row.names=F, sep="\t", quote=F, append = TRUE)
        
        reslist <- c(reslist, pcut = paste0('<Resource path="', basename(bgfile), '"/>'))
        tracklist <- c(
          tracklist,
          pcut = list(
            c(
              paste0('<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="255,0,0" colorScale="ContinuousColorScale;0.0;0.0;255,255,255;255,0,0" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="100" id="', basename(bgfile),'" name="', bnam, ' p0.001" renderer="BASIC_FEATURE" sortable="true" visible="true" windowFunction="count">'),
              paste0('<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="0.0" minimum="0.0" type="LINEAR"/>'),
              paste0('</Track>')
            )
          )
        )

        pblock <- c('<Regions>', paste0('<Region chromosome="', bdf$chrom, '" description="', bdf$name, '" end="', bdf$end, '" start="', bdf$start, '"/>'), '</Regions>')

        ## GRD
        if (grd) {
          oridir <- getwd()
          setwd(paste0(out.dir, '/', froot))
          tmpfile <- paste0(out.dir, '/', froot, "_", pcut, ".txt")
          grd.df <- data.frame(loc = paste0(outable$Chrom, ":", outable$Start, "-", outable$End), width = outable$End - outable$Start + 1, outable[,c(8,9,12:ncol(outable))])
          grd.df <- grd.df[p.ok,]
          colnames(grd.df)[1:4] <- c("#Loc", "Width", "Band1", "Band2")
          write.table(grd.df, tmpfile, sep = "\t", quote = FALSE, row.names = FALSE)
          ccword <- paste0("cc_", testype)
          mycmd <- paste0("grd --sp ", sp, " --gb ", gb, " -m ", ccword, " --ldb ", ldb, " ", tmpfile)
          try(system(mycmd))
          setwd(oridir)
        }
        
        file.rename(from = paste0(out.dir, '/', froot), to = paste0(out.dir, '/', "_", froot))
      }
    }

    ## XML RESOURCES part
    for (r in c("l2r", "cor", "pval", "pcut")) {
      r.res.exists <- which(names(reslist) == r)
      if (length(r.res.exists) > 0) for (l in r.res.exists) writeLines(text = reslist[[l]], con = xmlcon)
    }
    writeLines(text = '</Resources>', con = xmlcon)
    
    ## XML DATAPANEL part
    writeLines(text = '<Panel name="DataPanel">', con = xmlcon)
    for (r in c("l2r", "cor", "pval", "pcut")) {
      r.track.exists <- which(names(tracklist) == r)
      if (length(r.track.exists) > 0) for (l in r.track.exists) writeLines(text = tracklist[[l]], con = xmlcon)
    }
    writeLines(text = '</Panel>', con = xmlcon)
    
    ## XML FIXEDPANEL part
    writeLines(text = c(
      '<Panel name="FeaturePanel">',
      '<Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/>',
      '<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;308.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="hg19_genes" name="RefSeq Genes" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">',
      '<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="308.0" minimum="0.0" type="LINEAR"/>',
      '</Track>',
      '</Panel>',
      '<PanelLayout dividerFractions="0.8710059171597633"/>'
    ), con = xmlcon)
    
    ## XML REGIONS part
    if (!is.null(pcut)) {
      if (length(p.ok) > 0) {
        writeLines(text = pblock, con = xmlcon)
      }
    }
    
    ## XML end
    writeLines(text = c(
      '<HiddenAttributes>',
      '<Attribute name="DATA FILE"/>',
      '<Attribute name="DATA TYPE"/>',
      '<Attribute name="NAME"/>',
      '</HiddenAttributes>',
      '</Session>'
    ), con = xmlcon)
    
    close(xmlcon)
  }
  stopCluster(cl)
}






