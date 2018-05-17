## NHM2 pour "New HeatMap, version 2"
##
## DESCRIPTION:
##  This script generates a clusterized heatmap + dendrogram from a PROB data file. It was created to
##  replace the included heatmap() function which used too much RAM for very big matrices. The script
##  has been enhanced with support for hs or mm, and several CGH-speicific plots were added.
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
##
## DEPENDS ON:
## . limma (for QN)
## . cluster (base)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.8.2 to 2.15.1
##
## VERSION NOTES
##
## v2.4 20170503
##      - Added support for missing values when using correlation distances.
##
## v2.3d 20141117
##      . Adaptated the construction method to newer name "ward.D" (method is still the very same).
##      . Added support for direct TIFF output.
##      . Switched to TIFF output by default
##
## v2.3c 20140623
##      . Added by default "useRaster=T" to the 'image()' command, as it generates simpler and
##        smaller output for PDF format, without altering quality so much.
##
## v2.3b 20140331
##      . Modified default graphical device for the heatmap in png to 'Xlib'. This was required
##        as the image() function was not rendering the heatmap (empty plot) for high-dimensional
##        data with the default device 'cairo' under linux environments (bug obtained with
##        R v3.0.2 on a Fedora 19 x86_64 kernel 3.13.7.100). If this also fails one day, it still
##        seems to work with 'quartz'.
##
## v2.3 20140303
##      . Added support for rattus norvegicus (untested, though).
##
## v2.2 20120824
##      . Added a new output : a text file containing the matrix of 'cutree's for clusters from 2 to n(samples)-1.
##        This is to be used in annotation tables for CGH projects. Samples order is the same as input.
##      . The 'FRAC' output table (giving global genomic instability), is not sorted by sample names anymore.
##        The output is now in the same order as samples in input.
##
## v2.1e 20120626
##      . Modified the handling of graphics output. There are now 2 separated options : one for the file type
##        (png, bmp or pdf), and a second one for the dimensions (width for bitmap, L or P for vectorial).
##        This allows to create bitmap graphs with custom size (but a fixed ratio of 1.6), and also allows
##        the script to handle cex according to the requested graphics type.
##      . Corrected a bug where PNG were produced when requesting BMP (even though the files had a .bmp extension).
##
## v2.1d 20120106
##      . Fixed a bug that made the frequency and breakpoint plots for samples asynchronous with the heatmap
##        and dendrogram, when samples in the PROB file were not ordered (typically when columns names were
##        replaced from barcode to known sample names).
##
## v2.1c 20111116
##      . Introduced a "help" function to display the description of parameters for the main function, even when
##        this script is sourced and not directly read.
##
## v2.1b 20110916
##      . Corrected a minor typo preventing the script to be directly sourced.
##
## v2.1 20110831
##      . Added support for mus musculus genome.
##      . Translated the description of the main function.
##
## v2.3 20110623
##      . Added SGL plot (new graph) and values (in the table formerly containing frequencies only).
##      . Modified output filenames suffix to include the values of amp and del.
##      . Corrected the support of different image output format for the frequency plot (which was nroken since version 2.2).
##      . Added new output formats : bmp1680, bmp3200, bmp4800.
##      . Changed default loss color to blue (to fit with the new official recommandations for color-blind people)
##      . Added customization for the amplification and deletion dots colors.
##      . Corrected a typo (unneeded backspace in the bottom-left corner text).
##
## v2.2  20110610
##      . Changed header to the official common internal structure.
##      . Switched new comments to English.
##      . Added customization for lower ad higher values' colors for the heatmap
##      . Computation of FRAC is now made before the QN (when one is asked) so that the FRAC plots won't be
##        messed up by the QN.
##      . Renamed the default imgout "pdf" to "pdfL" (corresponding to PDF in a landscape disposition).
##      . Added a new imgout "pdfP" for a PDF in a portrait disposition.
##
## v2.1b 20110528
##      . Set 'path' option to "." by default.
##      . Added a description of options for the main function.
##
## v2.1  20100917
##      . Ajout d'une sortie texte des taux de genome aberrant et nb breakpoints par échantillon.
##    	. Ajout d'une sortie texte des fréquences de gain et perte (et ampli, del) par probe.
##	  	. Ajout d'une sortie graphique du profil moyen (graphe de fréquences) dans un autre PDF, car la version
##	  	  générée avec la heatmap sans annotation est trop petite pour être autre chose qu'une description rapide.
##
## v2.0b
##      . Gestion "simplifiee" de l'annotation clinique (fichier + colonne). Un header est obligatoire pour le fichier
##		    d'annotation.
##	  	. Correction d'un bug dans l'affichage et l'ordonancement des chromosomes sur la heatmap, dans le cas de
##	  	  regions issues de chr non consecutifs.
##
## v2.0
##      . Ajout de la possibilite d'afficher des annotations cliniques (avec légende), grace a Guillaume.
##		  . Ajout d'autres informations sous forme de representations graphiques :
##			  _ Taux de genome aberrant en gain, perte (et par deduction : total).
##		  	_ Taux de breakpoints (relatif a l'echantillon en contenant le plus).
##	  		_ Si aucune annotation clinique n'est specifiee, la frequence d'apparition
##		  	  des aberrations, en gain comme en perte.
##			  _ Un cartouche de texte rappelant :
##				  * Le fichier de donnees sources
##				  * Les distances et methodes utilisees pour le clustering
##		  		* Si une QN a ete effectuee
##		  		* le coefficient agglomeratif du clustering obtenu
##		  		* Si annotation clinique : le nom de colonne correspondant.
##
## v1.3
##      . Ajout de la prise en charge de distances Spearman comme modification mineure du support existant
##		    des distances Pearson.
##	  	. Ajout d'une option pour gerer le cut-off d'amplitude pour les couleurs maxi (utile pour P20_IALT_CGH
##	  	  par exemple, pour lequel les amplitudes sont tres faibles).
##	  	. Idem pour parametrer le nombre de couleurs.
##
## v1.2
##      . Ajout des options 'noX' et 'noY' (pour eliminer les chrX et chrY des clusterings, en particulier
##		    pour les projets externes pour lesquels le sexe n'est pas matché. (NOTA : ne fonctionne que pour
##		    l'humain X=23 et Y=24 !!)
##		  . Ajout d'un "autotweak" : pour les distances Pearson, l'ajout d'une valeur non-nulle pour tout
##		    profil dont l'ecart-type est nul. (NOTA : ecart-type et pas somme, car avec la normalisation
##		    quantile, on peut avoir un profil excentre mais d'ecart-type nul).
##		  . Ajout d'une option de normalisation par quantiles, pour faire joli.
##
## v1.1
##      . Ajout des distances Pearson grace a Cathy.
##		  . Inversion des couleurs des ampli/del pour un meilleur contraste.
##
## v1.0
##      . 1ere version.



library(limma)
library(cluster)

# cat("\nWARNING !\n*********\nGraphical device 'cairo' fails to properly display the heatmap\non high-dimensial data (typically 1M arrays for few samples),\nso this script requires you to have 'Xlib' available as graphical device.\n\n")


## GENERATION DU VECTEUR DE COULEURS
colz <- function(low = lowcol, high = highcol, mid=midcol) {
  n <- 100
  if (is.character(low)) low <- col2rgb(low)/255
  if (is.character(high)) high <- col2rgb(high)/255
  if (is.character(mid)) mid <- col2rgb(mid)/255
  col <- rgb(c(seq(low[1], mid[1], len=n/2), seq(mid[1], high[1], len = n/2)), c(seq(low[2], mid[2], len=n/2), seq(mid[2], high[2], len = n/2)), c(seq(low[3], mid[3], len=n/2), seq(mid[3], high[3], len = n/2)))
  return(col)
}

## CREATION DU LAYOUT
displayLayout = function(x){
  if (is.null(x)) {
    zone <- matrix(c(1,2,3,4,6,5,0,0), nrow=2, ncol=4, byrow=T)
    layW = c(.25,.675,.05,.025)
    layH = c(.85,.15)
  } else {
    zone <- matrix(c(1,2,3,4,5,7,6,0,0,0), nrow=2, ncol=5, byrow=T)
    layW = c(.25,.625,.05,.025,.05)
    layH = c(.85,.15)
  }
  layout(zone, widths=layW, heights = layH)
}

## PARAMETERS DESCRIPTOR FUNCTION
help.nhm.prob <- function() {
  cat("
PARAMETERS FOR THE nhm.prob() FUNCTION :
----------------------------------------
probname   = PROB filename.
distance   = Distance measurement type for the clustering. See dist() included with R. [euclidean]
method     = Construction method for the clustering. See hclust() included with R. [ward.D]
sp         = Species (hs, mm or rn). Needed to know which number have chrX and chrY. [hs]
amp        = Minimal log2(ratio) value to consider an amplification. [1.5]
del        = Maximal negative log2(ratio) value to consider deletion. [-1.5]
lowcol     = Color for losses (log2ratio < 0). [red]
midcol     = Color for normality (log2ratio == 0). [grey75]
highcol    = Color for gains (log2ratio > 0). [blue]
ampcol     = Color for amplifications (see 'amp'). [cyan4]
delcol     = Color for deletions (see 'del'). [orangered4]
colcut     = Maximal absolute log2(ratio) value for the heatmap saturation. [0.7]
ncol       = Number of different colors to use to generate the heatmap. [100]
labwidth   = Width (in characters) for the dendrogram labels. [7]
colscale   = Perform a scaling for heatmap coloring. If set to TRUE, it is recommended to increase colcut to a higher value. [FALSE]
noX        = Discard chrX. [FALSE]
noY        = Discard chrY. [TRUE]
qn         = Perform a quantiles normalization on the dataset, for the heatmap image only (does not impac clustering). See normalizeQuantile() from the limma package. [FALSE]
annotfile  = (optional) Clinical annotations filename.
annotcol   = Column number to use from the annotation file. Only needed if an annotation file name is also given.
gtype      = Graphical device ('cairo', 'cairo-png', 'Xlib', 'quartz'). Depends on the host system. [Xlib]
gformat    = Graphics output format : 'png', 'bmp', 'tif' (or 'tiff') and 'pdf'. [png]
gdim       = Dimension to use for the output device. An integer for PNG, BMP, and TIF, which corresponds to its width. The char 'L' or 'P' for PDF, corresponding to Landscape or Portrait. [1440]
useCairo   = Use cairo as plotting device. [FALSE]
\n")
}

## MAIN FUNCTION
nhm.prob <- function(probname, distance="euclidean", method="ward.D", sp="hs", amp=1.5, del=-1.5, lowcol="red", midcol="grey75", highcol="blue", ampcol="cyan4", delcol="orangered4", colcut=0.7, ncol=100, labwidth=7, colscale=FALSE, noX=FALSE, noY=TRUE, qn=FALSE, annotfile=NULL, annotcol=NULL, gtype = "Xlib", gformat="png", gdim = 1440, useCairo=FALSE) {
  
  require(data.table)
  # setwd(path)
  
  annot <- NULL
  annotname <- NULL
  prefix <- suffix <- xnum <- ynum <- ""
  chrnames <- vector()
  fnam <- basename(probname)
  fnam <- sub(pattern = ".prob", replacement = "", x = fnam)
  
  ## Specifying the chrX and chrY numbers
  if (sp == "hs") {
    xnum = 23
    ynum = 24
    chrnames <- c(1:22, "X", "Y")
  } else if (sp == "mm") {
    xnum = 20
    ynum = 21
    chrnames <- c(1:19, "X", "Y")
  } else if (sp == "rn") {
    xnum = 21
    ynum = 22
    chrnames <- c(1:20, "X", "Y")
  }else stop("Unknown species ! Are supported : hs, mm or rn.")
  
  ## Si l'annotation existe, on en tire ses nom et valeur
  if (!is.null(annotfile)) {
    annotall <- fread(input = annotfile, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
    annot <- as.factor(annotall[,annotcol])
    suffix <- paste0("_", annotname)
  }
  
  ## Loading PROB file
  prob <- fread(input = probname, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
  ## Discarding chrX data if asked
  if (noX == T) {
    prob <- prob[prob$Chr != xnum,]
    suffix <- paste0(suffix, "_noX")
  }
  ## Discarding chrY data if asked
  if (noY == T) {
    prob <- prob[prob$Chr != ynum,]
    suffix <- paste0(suffix, "_noY")
  }
  l2r <- prob[,-c(1:6)]
  samples <- colnames(l2r)
  
  suffix <- paste0(suffix, "_A", amp, "_D", del)
  suffix <- paste0(suffix, "_", gdim)
  
  pos <- prob[,c(1:6)]
  rm(prob)
  gc()
  
  mat <- matrix(as.numeric(as.matrix(l2r)), nrow=nrow(l2r), ncol=ncol(l2r))
  gc()
  
  
  ## Computing frequencies of aberrations (and number ofbreakpoints) for every sample
  ## NOTE : This has to be done BEFORE the QN (if one is asked), because it can mess up
  ##        these frequencies
  gainv <- lossv <- totv <- nbrk <- vector()
  gainv <- apply(mat, 2, function(z) { length(which(z > 0))/nrow(mat) })
  lossv <- apply(mat, 2, function(z) { length(which(z < 0))/nrow(mat) })
  totv <- gainv + abs(lossv)
  nbrk_int <- apply(mat, 2, function(z) { length(which(diff(z) != 0)) })
  nbrk <- nbrk_int / max(nbrk_int)
  fga.out <- data.frame(Sample=samples, Gain.Frac=gainv, Loss.Frac=lossv, Aberrant.Frac=totv, Breakpoints=nbrk_int, stringsAsFactors=F)
  write.table(fga.out, paste0(prefix, fnam, suffix, "_FRAC.txt"), quote=F, sep="\t", row.names=F)
  
  ## Performing limma's quantile normalization if asked
  if (qn == T) {
    mat <- normalizeQuantiles(mat, ties=T)
    prefix <- "QN_"
  }
  
  ## Transposing matrix for clustering
  mat <- t(mat)
  
  ## Computing distances
  if ( (distance == "pearson") | (distance == "spearman") ) {
    ## Autotweak for totaly flat profiles.
    for (k in 1:nrow(mat)) {
      if (sd(mat[k,], na.rm = TRUE) == 0) {
        mat[k, ncol(mat)] <- 0.01
      }
    }
    if (any(is.na(mat))) dd <- as.dist((1-cor(t(mat), method=distance, use = "na.or.complete"))/2) else dd <- as.dist((1-cor(t(mat), method=distance, use = "everything"))/2)
  } else {
    dd <- dist(mat, method=distance)
  }
  
  ## Clustering samples
  hc <- hclust(dd, method=method)
  hc$labels <- samples
  
  ## Dummping cuts if asked
  cut.df <- as.data.frame(matrix(data=NA, nrow=nrow(mat), ncol=nrow(mat)-1), stringsAsFactors=F)
  cut.df[,1] <- hc$labels
  colnames(cut.df)[1] <- "Sample"
  for (c in 2:ncol(cut.df)) {
    cut.df[,c] <- cutree(hc, c)
    colnames(cut.df)[c] <- paste0("Cut_", c)
  }
  write.table(cut.df, paste0(prefix, fnam, "_", distance, suffix, "_cuts.txt"), sep="\t", quote=F, row.names=F)
  
  ## Computing the cophenetic coefficient
  hcoef <- coef.hclust(hc)
  
  ## Retransposing for the heatmap
  mat <- t(mat)
  mat <- mat[,hc$order]
  
  # Color cutoff for saturation on heatmap
  pmat <- mat
  pmat[which(pmat > colcut)] <- colcut
  pmat[which(pmat < -colcut)] <- -colcut
  
  cpvec <- c(0, sapply(unique(pos$Chrom), function(x) { max(which(pos$Chrom == x)) }))
  cpdif <- diff(cpvec)
  cpvec2 <- cpvec[1:length(cpdif)] + cpdif/2
  
  if (tolower(gformat) == "png") {
    if (useCairo) {
      require(Cairo)
      Cairo(paste0(prefix, fnam, "_", distance, suffix, ".png"), width=gdim, height=gdim/1.6, bg = "white", type = "png", onefile = FALSE)
    }
    else {
      png(paste0(prefix, fnam, "_", distance, suffix, ".png"), width=gdim, height=gdim/1.6, type=gtype, res=96)
    }
    par(cex = gdim/840)
  }
  if (tolower(gformat) %in% c("tif", "tiff")) {
    tiff(paste0(prefix, fnam, "_", distance, suffix, ".tif"), width=gdim, height=gdim/1.6, type=gtype, compression="lzw", res=96)
    par(cex = gdim/840)
  }
  if (tolower(gformat) %in% c("jpg", "jpeg")) {
    jpeg(paste0(prefix, fnam, "_", distance, suffix, ".jpg"), width=gdim, height=gdim/1.6, type=gtype, quality = 85, res=96)
    par(cex = gdim/840)
  }
  if ((tolower(gformat) == "pdf") & (gdim == "L")) pdf(paste0(prefix, fnam, "_", distance, suffix, ".pdf"), height=21/cm(1), width=29.7/cm(1))
  if ((tolower(gformat) == "pdf") & (gdim == "P")) pdf(paste0(prefix, fnam, "_", distance, suffix, ".pdf"), height=29.7/cm(1), width=21/cm(1))
  
  displayLayout(annot)
  
  par(mar=c(1,1,2,labwidth))
  plot(as.dendrogram(hc, cex=.2), axes=F, horiz=T, yaxs="i", xaxs="i")

  par(mar=c(1,1,2,1))
  if (colscale) pmat <- scale(pmat, center = FALSE, scale = TRUE)
  image(x=1:nrow(l2r), y=1:ncol(l2r), z=pmat, col=colz(low=lowcol, high=highcol, mid=midcol), zlim=c(-colcut,colcut), axes=F, xaxs="i", yaxs="i", useRaster = TRUE)
  abline(v=cpvec, col="white")
  mtext(chrnames[sort(unique(pos$Chr))], side=1, at = cpvec2)
  for (n in 1:dim(pmat)[2]) {
    yind <- which(mat[,n] >= amp)
    bind <- which(mat[,n] <= del)
    points(yind, rep(n, length(yind)), pch = 20, col=ampcol, cex=1.5)
    points(bind, rep(n, length(bind)), pch = 20, col=delcol, cex=1.5)
  }
  
  # Vecteur de couleurs
  myColors = c("aquamarine2", "chocolate2", "darkgoldenrod1", "grey75", "cyan2", "darkorange", "darkorchid", "grey25")
  
  ## Plot frequences internes
  par(mar=c(1,0,2,1))
  fga.out <- fga.out[hc$order,]
  midbar <- (fga.out[["Gain.Frac"]] - fga.out[["Loss.Frac"]]) / 2
  plot(0,0, xlim=c(-1,1), ylim=c(0,ncol(mat)), axes=F, xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n")
  for (j in seq(-0.75,0.75,0.25)) abline(v=j, lty=2, col = "grey25")
  rect(0, 0:(ncol(mat)-1), fga.out[["Gain.Frac"]], 1:ncol(mat), col=highcol, border="black")
  rect(0, 0:(ncol(mat)-1), -fga.out[["Loss.Frac"]], 1:ncol(mat), col=lowcol, border="black")
  segments(midbar, 0:(ncol(mat)-1), midbar, 1:ncol(mat), col="yellow")
  rect(-1, 0:(ncol(mat)-1), 1, 1:ncol(mat), border="black")
  for (u in c(-1,0,1)) mtext(u, side=3, at = u)
  
  par(mar=c(1,0,2,1))
  plot(0,0, xlim=c(0,1), ylim=c(0,ncol(mat)), axes=F, xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n")
  rect(0, 0:(ncol(mat)-1), fga.out[["Breakpoints"]]/max(fga.out[["Breakpoints"]]), 1:ncol(mat), col="grey75", border="black")
  rect(0, 0:(ncol(mat)-1), 1, 1:ncol(mat), border="black")
  for (u in c(0,1)) mtext(u, side=3, at = u)
  
  if (!is.null(annotfile)) {
    ## Plot ANNOTATIONS
    par(mar=c(1,0,2,1))
    plot(0,0, xlim=c(0,1), ylim=c(0,length(annot)), axes=F, xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n")
    rect(0,0,1, length(annot), col="grey50", border="black") ## sous-couche grise pour les NA.
    rect(0, 0:(length(annot)-1), 1, 1:length(annot), col=myColors[as.numeric(annot[hc$order])], border="black")
    
    ## Plot LEGENDS
    xall <- seq(0,100,30)
    xmax <- max(xall)*1.4
    boxw <- xmax/length(xall)/3
    levlen <- length(levels(annot))
    par(mar=c(1,0,2,1))
    plot(0, 0, xlim=c(0,xmax), ylim=c(0,levlen+1), axes=F, xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n")
    
    Xcoord <- matrix(c(rep(xall,2), rep(levlen,length(xall)), rep(levlen-1,length(xall))), byrow=T, ncol=length(xall)*2)
    for (g in 1:levlen) {
      rect(Xcoord[1,g],Xcoord[2,g],Xcoord[1,g]+boxw,Xcoord[2,g]+0.80, col=myColors[g], border="black")
      text(Xcoord[1,g]+boxw,Xcoord[2,g]+0.4, levels(annot)[g], cex=2, pos=4)
    }
    if(length(which(is.na(annot))) > 0) {
      rect(Xcoord[1,levlen+1],Xcoord[2,levlen+1],Xcoord[1,levlen+1]+boxw,Xcoord[2,levlen+1]+0.8, col="grey50", border="black")
      text(Xcoord[1,levlen+1]+boxw,Xcoord[2,levlen+1]+0.4, "NA", cex=2, pos=4)
    }
    
  } else {
    par(mar=c(1,1,1,1))
    plot(0,0, xlim=c(0,nrow(mat)), ylim=c(-1,1), axes=F, xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n")
    for (j in seq(-0.75,0.75,0.25)) abline(h=j, lty=2)
    Gfreq <- apply(mat, 1, function(q) { length(which(q > 0)) / ncol(mat) })
    Lfreq <- apply(mat, 1, function(q) { length(which(q < 0)) / ncol(mat) })
    Afreq <- apply(mat, 1, function(q) { length(which(q >= amp)) / ncol(mat) })
    Dfreq <- apply(mat, 1, function(q) { length(which(q <= del)) / ncol(mat) })
    lines(Gfreq, type="h", col=highcol)
    lines(-Lfreq, type="h", col=lowcol)
    lines(Afreq, type="h", col=ampcol)
    lines(-Dfreq, type="h", col=delcol)
    for (k in unique(pos$Chrom)) {
      chrpos <- max(which(pos$Chrom==k))
      abline(v=chrpos, col="black")
    }
    rect(0,-1, nrow(mat), 1, border="black")
    for (u in c(-1,0,1)) mtext(u, side=2, at = u)
  }
  
  plot(0,0, xlim=c(0,100), ylim=c(-6,0), axes=F, xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n")
  const <- c(paste0(fnam, "\nDistances: ", toupper(distance), "\nMethod: ", toupper(method), "\nQuantileNorm: ", qn, "\nAgglomerative coef: ", format(hcoef, digits=4)))
  text(0, -2, const, pos=4, cex=1, font = 2)
  text(0, -5, paste0("Amp (", amp, ")"), pos=4, col=ampcol, cex=1, font = 2)
  text(40, -5, paste0("Del (", del, ")"), pos=4, col=delcol, cex=1, font = 2)
  dev.off()
  
  gc()
  
  ## Plot fréquences externe
  if (tolower(gformat) == "png")  {
    png(paste0(prefix, fnam, suffix, "_Freqs.png"), width=gdim, height=gdim/1.6, type = gtype)
    par(cex = gdim/840)
  }
  if (tolower(gformat) %in% c("tif","tiff")) {
    tiff(paste0(prefix, fnam, suffix, "_Freqs.tif"), width=gdim, height=gdim/1.6, compression="lzw", type = gtype)
    par(cex = gdim/840)
  }
  if (tolower(gformat) %in% c("jpg","jpeg")) {
    jpeg(paste0(prefix, fnam, suffix, "_Freqs.jpg"), width=gdim, height=gdim/1.6, quality = 85, type = gtype)
    par(cex = gdim/840)
  }
  if ((tolower(gformat) == "pdf") & (gdim == "L")) pdf(paste0(prefix, fnam, suffix, "_Freqs.pdf"),height=21/cm(1),width=29.7/cm(1))
  if ((tolower(gformat) == "pdf") & (gdim == "P")) pdf(paste0(prefix, fnam, suffix, "_Freqs.pdf"),height=29.7/cm(1),width=21/cm(1))
  
  plot(0,0, xlim=c(0,nrow(mat)), ylim=c(-1,1), xaxs="i", yaxs="i", xaxt="n", type="n", main=paste0("Average profile of the population (", ncol(mat), " samples)"), xlab="Genomic position", ylab="Frequency of losses | gains")
  
  Gfreq <- apply(mat, 1, function(q) { length(which(q > 0)) / ncol(mat) })
  Lfreq <- apply(mat, 1, function(q) { length(which(q < 0)) / ncol(mat) })
  Afreq <- apply(mat, 1, function(q) { length(which(q >= amp)) / ncol(mat) })
  Dfreq <- apply(mat, 1, function(q) { length(which(q <= del)) / ncol(mat) })
  lines(Gfreq, type="h", col=highcol)
  lines(-Lfreq, type="h", col=lowcol)
  lines(Afreq, type="h", col=ampcol)
  lines(-Dfreq, type="h", col=delcol)
  for (j in seq(-0.75,0.75,0.25)) abline(h=j, lty=3)
  abline(h=0)
  for (k in unique(pos$Chrom)) {
    chrpos <- max(which(pos$Chrom==k))
    abline(v=chrpos, col="black", lty=2)
  }
  if (tolower(gtype) == "pdf") mtext(chrnames[sort(unique(pos$Chr))], side=1, at = cpvec2)
  if (tolower(gtype) != "pdf") mtext(chrnames[sort(unique(pos$Chr))], side=1, at = cpvec2, cex=gdim/840)
  
  dev.off()
  
  ## Plot SGL
  if (tolower(gformat) == "png") {
    png(paste0(prefix, fnam, suffix, "_Sums.png"), width=gdim, height=gdim/1.6, type = gtype)
    par(cex = gdim/840)
  }
  if (tolower(gformat) %in% c("tif", "tiff")) {
    tiff(paste0(prefix, fnam, suffix, "_Sums.tif"), width=gdim, height=gdim/1.6, compression = "lzw", type = gtype)
    par(cex = gdim/840)
  }
  if (tolower(gformat) %in% c("jpg", "jpeg")) {
    jpeg(paste0(prefix, fnam, suffix, "_Sums.jpg"), width=gdim, height=gdim/1.6, quality = 85, type = gtype)
    par(cex = gdim/840)
  }
  if ((tolower(gformat) == "pdf") & (gdim == "L")) pdf(paste0(prefix, fnam, suffix, "_Sums.pdf"),height=21/cm(1),width=29.7/cm(1))
  if ((tolower(gformat) == "pdf") & (gdim == "P")) pdf(paste0(prefix, fnam, suffix, "_Sums.pdf"),height=29.7/cm(1),width=21/cm(1))
  
  Gsum <- apply(mat, 1, function(q) { sum(q[which(q > 0)]) })
  Lsum <- apply(mat, 1, function(q) { sum(q[which(q < 0)]) })
  Asum <- apply(mat, 1, function(q) { sum(q[which(q >= amp)]) })
  Dsum <- apply(mat, 1, function(q) { sum(q[which(q <= del)]) })
  
  plot(0,0, xlim=c(0,nrow(mat)), ylim=range(c(Gsum, Lsum)), xaxs="i", yaxs="i", xaxt="n", type="n", main=paste0("Sums of gains and losses for the population (", ncol(mat), " samples)"), xlab="Genomic position", ylab="Sum of log2(ratio) in losses | gains")
  
  lines(Gsum, type="h", col=highcol)
  lines(Lsum, type="h", col=lowcol)
  lines(Asum, type="h", col=ampcol)
  lines(Dsum, type="h", col=delcol)
  for (k in unique(pos$Chrom)) {
    chrpos <- max(which(pos$Chrom==k))
    abline(v=chrpos, col="black", lty=2)
  }
  if (tolower(gtype) == "pdf") mtext(chrnames[sort(unique(pos$Chr))], side=1, at = cpvec2)
  if (tolower(gtype) != "pdf") mtext(chrnames[sort(unique(pos$Chr))], side=1, at = cpvec2, cex= gdim/840)
  dev.off()
  
  ## SORTIE FREQUENCES PAR PROBE
  sfreq <- data.frame(pos, GainFreq=Gfreq, LossFreq=Lfreq, AmpFreq=Afreq, DelFreq=Dfreq, GainSum=Gsum, LossSum=Lsum, AmpSum=Asum, DelSum=Dsum, stringsAsFactors=F)
  write.table(sfreq, paste0(prefix, fnam, suffix, "_FREQSUM.txt"), quote=F, sep="\t", row.names=F)
  gc()
  
}
