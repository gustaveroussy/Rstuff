## chrload.R
##
## DESCRIPTION : Load chromosomal informations from cytoBandIdeo for a given species and genome build
##
## AUTHOR:
## Bastien JOB (bastien.job@gustaveroussy.fr)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 3.2.3
##
## DEPENDS ON:
## R-base
##
## VERSION NOTES
##
## v1.2b 20160830
##      . Removed support for repository on pelican
##      . Replaced it with the "ldb" parameter (local data base), currently pointing to cigogne/data/bioinfo
##
## v1.2 20160407
##      . Added vector of chrom, chrN, chrA and midpoint directly in cs
##
## v1.1 20160125
##      . Cleaned some code
##      . Divided code into two functions to allow the use of custom sources (code would earlierly work only using pelican's data structure path)
##
## v1.0 20140307
##      . Added chrom2chr and chr2chrom converter lists, species long name, and genome version to the output object.

## Build a cs object from a cytOBandIdeo table from the UCSC genome browser
read_ucsc_cytobandideo <- function(file = NULL, sp = NULL, gb = NULL) {
  
  gv <- paste0(sub("hs", "hg", sp), gb)
  
  chrom2chr = list(
    hs = list("chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chr21"=21,"chr22"=22,"chrX"=23,"chrY"=24),
    mm = list("chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chrX"=20,"chrY"=21),
    rn = list("chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chrX"=21,"chrY"=22)
  )
  
  chr2chrom=list(
    hs = list("1"="chr1", "2"="chr2",  "3"="chr3",  "4"="chr4", "5"="chr5", "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10", "11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", "16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19", "20"="chr20", "21"="chr21", "22"="chr22", "23"="chrX", "24"="chrY"), 
    mm = list("1"="chr1", "2"="chr2",  "3"="chr3",  "4"="chr4", "5"="chr5", "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10", "11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", "16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19", "20"="chrX", "21"="chrY"), 
    rn = list("1"="chr1", "2"="chr2",  "3"="chr3",  "4"="chr4", "5"="chr5", "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10", "11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", "16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19", "20"="chr20", "21"="chrX", "22"="chrY")
  )
  
  ## Converting sp to text
  sptxtlong <- list(hs="homo_sapiens", mm="mus_musculus", rn="rattus_norvegicus")
  
  ## IMPORT BANDES CYTO ET PREPARATION POUR PLOT KARYO
  cat("Importing chromosomes data  ...\n")
  cytob <- read.table(file = file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "_", fill = TRUE)
  
  ## Filtering out invalid lines (UCSC recently added incomplete lines which seem to refer to chrN_unknown sequence sets, which is quite totally illogical. BJ 20140227)
  cytob <- cytob[!is.na(cytob$chromStart),]
  
  ## Filtering out ANY chromosome not defined in chrconvlist$chrom2chr
  cytob <- cytob[cytob[["X.chrom"]] %in% names(chrom2chr[[sp]]),]
  
  ## Modified versions of chr names
  cytob$chrA <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
  cytob$chr <- as.numeric(chrom2chr[[sp]][cytob$X.chrom])
  
  ## Ordering by chromosomal number
  cytob <- cytob[order(cytob$chr),]
  
  ## Additional data for chromosomes
  lchrxx <- vapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) }, 1)
  lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
  lchrsum <- cumsum(lchrtoadd)
  centromere <- cytob$chromEnd[cytob$gieStain == "acen" & strtrim(cytob$name, 1) == "p"]
  glen <- sum(lchrxx)
  
  ## Adding genomic coordinates for cytob
  cytob$genoStart <- cytob$chromStart + lchrsum[cytob$chr]
  cytob$genoEnd <- cytob$chromEnd + lchrsum[cytob$chr]
  
  ## Adding genomic coordinates for centromeres too
  centromere.geno <- centromere + lchrsum[unique(cytob$chr)]
  
  ## Adding bands widths for plots
  cytob$x1 <- 0.15
  cytob$x2 <- 0.85
  atypical.bands <- cytob$gieStain %in% c("acen", "gvar", "stalk")
  cytob$x1[atypical.bands] <- .25
  cytob$x2[atypical.bands] <- .75
  
  ## Adding bands color for plots
  bandtype2col <- list(gneg="white", gpos25="grey75", gpos33="grey66", gpos50="grey50", gpos66="grey33", gpos75="grey25", gpos100="black", gpos="black", acen="yellow", gvar="darkred", stalk="darkolivegreen")
  cytob$gieStainCol <- unlist(bandtype2col[cytob$gieStain])
  
  
  cs <- list(cytob = cytob,
             species = sptxtlong[[sp]],
             chrom2chr = chrom2chr[[sp]],
             chr2chrom = chr2chrom[[sp]],
             chrom = unique(cytob$X.chrom),
             chrA = unique(cytob$chrA),
             chr = unique(cytob$chr),
             centromere = centromere,
             centromere.geno = centromere.geno,
             lchrxx=lchrxx,
             lchrsum=lchrsum,
             lchrtoadd=lchrtoadd,
             glen=glen,
             gv=gv)
  return(cs)
}

## IMPORT CHR DATA FOR A GIVEN SPECIES AND GENOME BUILD, USING PELICAN DATA STRUCTURE
chrload <- function(sp = NULL, gb = NULL, ldb = "/mnt/data_cigogne/bioinfo") {
  gv <- paste0(sub("hs", "hg", sp), gb)
  read_ucsc_cytobandideo(file = paste(ldb, "/GoldenPath/", gv, "/cytoBandIdeo.", gv, sep=""), sp = sp, gb = gb)
}
