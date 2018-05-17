## This function converts a GEML design file provided by Agilent to three files needed
## for our Agilent CGH analysis pipeline :
## - DNABack_BCLeft info table
## - FASTA file containing sequence for all available probes
## - BED file (needed to generate GC/CnGs tracks)
##
## AUTHOR : Bastien Job (bastien.job@inserm.fr)
## VERSION : 1.0 (20161007)

Agilent_CGH_GEML_converter <- function(xml.file = NULL, ncores = 6) {
  
  if (is.null(xml.file)) stop("Please provide a valid Agilent CGH GEML design file !")
  
  ## Reading header
  print("Reading header ...")
  xh <- readLines(xml.file, n = 12)
  xh <- gsub(pattern = "(<!--|-->|  <|<|>)", replacement = "", xh)
  strucline <- unlist(strsplit(x = xh[10], split = "(, |=)"))
  NumRows <- as.numeric(strucline[which(strucline == "NumFeaturesRows")+1])
  NumCols <- as.numeric(strucline[which(strucline == "NumFeaturesColumns")+1])/2
  hline <- unlist(strsplit(x = xh[12], split = "(\" |=)"))
  species <- gsub(pattern = "\"", replacement = "", hline[which(hline == "species_database")+1])
  dline <- unlist(strsplit(x = xh[11], split = "(\" |=)"))
  date <- gsub(pattern = "(\"|/)", replacement = "", dline[which(dline == "date")+1])
  barcode <- gsub(pattern = "(\"|Agilent-)", replacement = "", dline[which(dline == "name")+1])
  
  ## Getting XML data
  print("Scanning XML ...")
  myXML <- XML::xmlParse(xml.file)
  xmltop <- XML::xmlRoot(myXML)[[1]]
  
  ## Creating the DNABack_BCLeft file
  print("Retrieving and formating needed data ...")
  
  block.list <- c("reporter", "reporter", "feature", "reporter", "gene", "gene")
  attrib.list <- c("systematic_name", "name", "number", "control_type", "primary_name", "description")
  type.list <- c("a", "a", 1L, "a", "a", "a")
  
  `%do%` <- foreach::"%do%"
  xmldf <- foreach::foreach(xk = 1:length(block.list), .combine = "cbind", .export = "xmltop") %do% {
    print(attrib.list[xk])
    vapply(XML::getNodeSet(xmltop, paste0("//", block.list[xk])), function(gns) {
      xres <- XML::xmlGetAttr(gns, attrib.list[xk])
      if(length(xres) > 0) return(xres) else return("")
    }, type.list[xk])
  }

  colnames(xmldf) <- c("Name", "ID", "RefNumber", "ControlType", "GeneName", "Description")
  xmldf <- as.data.frame(xmldf, stringsAsFactors = FALSE)
  xmldf$RefNumber <- as.numeric(xmldf$RefNumber)
  xmldf$Description <- gsub(pattern = "ref\\|", replacement = "", xmldf$Description)
  xmldf$ControlType[xmldf$ControlType == ""] <- "false"
  xmldf <- xmldf[!is.na(xmldf$RefNumber),]
  xmldf <- xmldf[order(xmldf$RefNumber),]
  ### Generating column and row position from RefNumber
  xmldf$Column <- xmldf$RefNumber %% NumCols
  xmldf$Column[xmldf$Column == 0] <- NumCols
  xmldf$Row <- (xmldf$RefNumber %/% NumCols) + 1
  xmldf$Row[xmldf$Column == NumCols] <- xmldf$Row[xmldf$Column == NumCols] - 1
  xmldf$ChromosomalLocation <- xmldf$Name
  xmldf$ChromosomalLocation[grep(pattern = "chr([0-9]+|X|Y|M):", x = xmldf$ChromosomalLocation, invert = TRUE)] <- "Unknown"
  ### Generating dummy columns
  xmldf$TopHit <- xmldf$Go <- xmldf$EntrezGeneID <- xmldf$PerformanceScore <- ""
  ### Reordering to the fixed column order of DNABack_BCLeft file
  xmldf <- xmldf[,c(7,8,1:5,13,6,12,9,11,10)]
  print("Writing DNABack_BCLeft file ...")
  ### Writing
  write.table(xmldf, paste0(barcode, "_D_DNABack_BCLeft_", date, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  ## Creating the FASTA file
  print("Writing FASTA file ...")
  myseq <- vapply(XML::getNodeSet(xmltop, "//reporter"), function(gns) { xres <- XML::xmlGetAttr(gns, "active_sequence") ; if(length(xres) > 0) return(xres) else return("") }, "a")
  okseq <- which(myseq != "")
  facon <- file(paste0(barcode, "_D_Fasta_", date, ".txt"), open = "w")
  for(xk in okseq) writeLines(c(paste0(">", xmldf$ID[xk]), myseq[xk]), con = facon)
  close(facon)
  
  ## Creating the BED file
  print("Writing BED file ...")
  ### Removing non-CGH probes
  xmldf <- xmldf[xmldf$ControlType == "false",]
  ### Removing replicated probes
  xmldf <- xmldf[!duplicated(xmldf$ID),]
  locprobes <- which(xmldf$ChromosomalLocation != "Unknown")
  chrstartend <- unlist(strsplit(x = xmldf$ChromosomalLocation[locprobes], split = "(:|-)"))
  mybed <- as.data.frame(matrix(chrstartend, ncol = 3, byrow = TRUE), stringsAsFactors = FALSE)
  colnames(mybed) <- c("Chrom", "Start", "End")
  mybed$Start <- as.numeric(mybed$Start)
  mybed$End <- as.numeric(mybed$End)
  mybed$ProbeName <- xmldf$ID[locprobes]
  ### Generating the numerical chrN vector to allow sorting on chr > start > end
  bedchrN <- mybed[,1]
  bedchrN <- sub(pattern = "chr", replacement = "", x = bedchrN)
  bedchrN[bedchrN == "X"] <- 23
  bedchrN[bedchrN == "Y"] <- 24
  bedchrN <- as.numeric(bedchrN)
  ### Sorting
  mybed <- mybed[order(bedchrN, mybed$Start, mybed$End),]
  ### Writing
  bedcon <- file(paste0(barcode, "_D_BED_", date, ".bed"), open = "w")
  writeLines(paste0("browser position", xmldf$ChromosomalLocation[locprobes[1]]), con = bedcon)
  writeLines(paste0("track name=\"Agilent-", barcode, "\" description=\"Agilent-", barcode, "\" visibility=2 color=0,128,0 useScore=1"), con = bedcon)
  close(bedcon)
  write.table(mybed, file = paste0(barcode, "_D_BED_", date, ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
  print("Done.")
}


setwd("/mnt/toucan_seq3/B16_AUVI/design/GEML_parser")
xml.file <- "../046984_D_F_20130131.xml"
