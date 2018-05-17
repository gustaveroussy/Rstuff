## This function reads a cytoBandIdeo file from the UCSC genome browser
## and generates an object (a list) with various informations about chromosomes.
## This object should be usefull for all and any genomics analysis project here.
##
## AUTHOR : Bastien Job (bastien.job@inserm.fr)
## VERSION : 1.0 (20161007)

read_ucsc_cytobandideo <- function(file = "cytoBandIdeo.hg19", genomebuild = "hg19", species = "Homo sapiens") {
  
  if(any(is.null(c(file, genomebuild, species)))) stop("Please provide information for all of the 3 parameters (file, genomebuild, species) !")
  
  ## IMPORT BANDES CYTO ET PREPARATION POUR PLOT KARYO
  cat("Importing chromosomes data  ...\n")
  cytob <- read.table(file = file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE, check.names = FALSE)
  colnames(cytob) <- c("chrom", "start", "end", "cytoband", "giestain")
  
  ## Filtering out invalid lines (UCSC recently added incomplete lines which seem to refer to chrN_unknown sequence sets, which is quite totally illogical. BJ 20140227)
  cytob <- cytob[!is.na(cytob$start),]
  
  ## Filtering out ANY chromosome not defined in chrconvlist$chrom2chr
  cytob <- cytob[grep(pattern = "^chr([0-9]+|X|Y|M)$", x = cytob$chrom),]
  
  ## Adding alternative versions of chr names
  cytob$chrN <- cytob$chrA <- sub(pattern = "^chr", replacement = "", x = cytob$chrom)
  non.auto <- c("X", "Y", "M")
  num.chr.max <- max(as.numeric(cytob$chrN[!cytob$chrN %in% non.auto]))
  for (cx in 1:length(non.auto)) cytob$chrN[cytob$chrN == non.auto[cx]] <- num.chr.max + cx
  cytob$chrN <- as.integer(cytob$chrN)
  
  ## Ordering by chromosomal number
  cytob <- cytob[order(cytob$chrN, cytob$start, cytob$end),]
  
  ## Additional data for chromosomes
  chrom <- data.frame(chrom = unique(cytob$chrom), chrA = unique(cytob$chrA), chrN = unique(cytob$chrN), stringsAsFactors = FALSE)
  chrom$chr.length <- vapply(unique(cytob$chrN), function(x) { max((cytob[which(cytob$chrN == x),])$end) }, 1)
  lchrtoadd <- c(0, chrom$chr.length[1:length(chrom$chr.length)-1])
  chrom$chr.length.sum <- cumsum(chrom$chr.length)
  chrom$chr.length.toadd <- c(0, chrom$chr.length.sum[-c(length(chrom$chr.length.sum))])
  chrom$mid.chr <- round(chrom$chr.length / 2)
  chrom$mid.chr.geno <- chrom$mid.chr + chrom$chr.length.toadd[chrom$chrN]
  chrom$centromere <- cytob$end[vapply(unique(cytob$chrN), function(k) {
    centro.pos <- which(cytob$chrN == k & cytob$giestain == "acen" & strtrim(cytob$cytoband, 1) == "p")
    return(ifelse(length(centro.pos) == 0, NA, centro.pos))
  }, 1L)]
  names(chrom$centromere) <- unique(cytob$chrN)
  glen <- sum(chrom$chr.length)
  
  ## Adding genomic coordinates for cytob
  cytob$start.geno <- cytob$start + chrom$chr.length.toadd[cytob$chrN]
  cytob$end.geno <- cytob$end + chrom$chr.length.toadd[cytob$chrN]
  
  ## Adding genomic coordinates for centromeres too
  chrom$centromere.geno <- chrom$centromere + chrom$chr.length.toadd[unique(cytob$chrN)]
  
  ## Adding converters
  chrom2chr <- as.list(unique(cytob$chrN))
  names(chrom2chr) = unique(cytob$chrom)
  chr2chrom <- as.list(unique(cytob$chrom))
  names(chr2chrom) = unique(cytob$chrN)
  
  ## Adding bands widths for plots
  cytob$x1 <- 0.15
  cytob$x2 <- 0.85
  atypical.bands <- cytob$giestain %in% c("acen", "gvar", "stalk")
  cytob$x1[atypical.bands] <- .25
  cytob$x2[atypical.bands] <- .75
  
  ## Adding bands color for plots
  bandtype2col <- list(gneg="white", gpos25="grey75", gpos33="grey66", gpos50="grey50", gpos66="grey33", gpos75="grey25", gpos100="black", gpos="black", acen="yellow", gvar="darkred", stalk="darkolivegreen")
  cytob$giestaincol <- unlist(bandtype2col[cytob$giestain])
  
  cs <- list(species = species,
             genomebuild = genomebuild,
             cytobands = cytob,
             chromosomes = chrom,
             chrom2chr = chrom2chr,
             chr2chrom = chr2chrom,
             # chrom = unique(cytob$chrom),
             # chrA = unique(cytob$chrA),
             # chr = unique(cytob$chrN),
             # centromere = centromere,
             # centromere.geno = centromere.geno,
             # midchr = midchr,
             # midchr.geno = midchrsum,
             # lchrxx=lchrxx,
             # lchrsum=lchrsum,
             # lchrtoadd=lchrtoadd,
             genome.length=glen)
  return(cs)
}

