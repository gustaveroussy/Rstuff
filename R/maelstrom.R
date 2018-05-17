## Set of useful functions


#####################
### GRAPH DEVICES ###
#####################

## Opens a "connection" to output plots to a multipage TIFF with rasters
## It invisibly writes those plots to individual PNGs
## Further ... argments will be passed to png()
rasterpdf.open <- function(file = "test.pdf", pattern = "abcxyz", ...) {
  pngname <- paste0(pattern, "%04d.png")
  ## Checking if temporary files corresponding to this pattern exist
  toclean <- list.files(pattern = paste0(pattern, "[0-9]{4}.png"), full.names = TRUE)
  ## If so, delete them
  if (length(toclean) > 0) file.remove(toclean)
  png(pngname, ...)
  return(list(pdfname = file, pattern = pattern))
}

## Closes the "connexion" opened with rasterpdf.open()
## It needs output values from rasterpdf (a list containing the pdf output
## name and the pattern to find png files).
## Further ... argments will be passed to pdf()
rasterpdf.close <- function(pdf.optionslist = NULL, clean = TRUE, width = 29.7/cm(1), height = 21/cm(1), ...) {
  try(dev.off())
  multipng2pdf(pdfname = pdf.optionslist$pdfname, pattern = pdf.optionslist$pattern, clean = clean, width = width, height = height, ...)
}

## Function to concatenate PNG files to a multipage (raster) pdf
multipng2pdf <- function(pdfname = NULL, pattern = "abcxyz", clean = TRUE, width = 29.7/cm(1), height = 21/cm(1), ... ) {
  pnglist <- list.files(path = dirname(pdfname), pattern = paste0(pattern, "[0-9]+.png"), full.names = TRUE)
  pnglist <- sort(vapply(pnglist, tools::file_path_as_absolute, "a", USE.NAMES = FALSE))
  if (length(pnglist) == 0) stop("MULTIPNG2PDF : No PNG file found to aggregate (missing, or bad pattern?) !")
  pdf(file = pdfname, width = width, height = height, ...)
  par(mai=c(0,0,0,0))
  for (pl in pnglist) {
    mypng <- png::readPNG(pl)
    graphics::plot(c(0,1),c(0,1),type="n")
    rasterImage(mypng,0,0,1,1)
  }
  dev.off()
  if (clean) zzz <- file.remove(pnglist)
}



####################
### FILE PARSING ###
####################

## Fast file writer using iotools::write.csv.raw but corrected for header handling (~ 5x faster)
write.table.fast <- function(x, file = NULL, header = TRUE, sep = "\t", fileEncoding="", row.names = FALSE, ...) {
  if (header) write.table(x = x[NULL,], file = file, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = fileEncoding)
  if(!row.names) rownames(x) <- NULL
  trychk <- try(iotools::write.csv.raw(x = x, file = file, sep = sep, col.names=FALSE, fileEncoding=fileEncoding, append = header, ...))
  if (!is.null(trychk)) {
    print("Fast write failed, using canonical write.table ...")
    write.table(x = x, file = file, sep = sep, row.names = row.names, quote = FALSE)
  }
  gc()
}

## Fast file reader using data.table::fread
read.table.fast <- function(file = NULL, header = TRUE, sep= "\t", row.names = FALSE, ...) {
  if (row.names) {
    if (header) h.df <- read.table(file = file, sep = sep, header = header, nrows = 1, check.names = FALSE)
    data.df <- data.table::fread(input = file, sep = sep, header = FALSE, skip = 1, data.table = FALSE, ...)
    rownames(data.df) <- data.df[,1]
    data.df[,1] <- NULL
    if (header) colnames(data.df) <- colnames(h.df)
  } else {
    data.df <- data.table::fread(input = file, sep = sep, header = header, data.table = FALSE, ...)
  }
  return(data.df)
}

## Function to load data from a SQLite ".db" file, using DBI and RSQLite. Returns a list of tables.
db.load <- function(db.file = NULL) {
  if (is.null(db.file)) stop("Please provide a .db file !")
  if(!file.exists(db.file)) stop("Provided .db file does not exist !")
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db.file)

  tablist <- DBI::dbListTables(con)
  tablist.ln <- length(tablist)
  if (tablist.ln < 1) stop("No table could be found in the provided db file !")
  print(paste0("Found ", tablist.ln, " tables :"))
  print(tablist)

  myDBdata <- sapply(tablist, function(x) {
    return(DBI::dbFetch(DBI::dbSendQuery(con, paste0("SELECT * FROM ", x)), n=-1))
  }, simplify = FALSE)

  DBI::dbDisconnect(conn = con)
  return(myDBdata)
}

## Function to load data from a HDF5 file, using rhdf5. Returns a list of tables.
# hdf5.load <- function(file = NULL, recursive = TRUE, all = TRUE, ...) {
#   if (is.null(file)) stop("Please provide a HDF5 file !")
#   if(!file.exists(file)) stop("Provided HDF5 file does not exist !")
#   hdf5.ls <- rhdf5::h5ls(file = file, recursive = recursive, all = all, ...)
#   getidx <- which(hdf5.ls$group == "/")
#   hdf5.data <- sapply(getidx, function(x) { return(rhdf5::h5read(file = file, name = hdf5.ls$name[x])) })
#   names(hdf5.data) <- hdf5.ls$name[getidx]
#   return(hdf5.data)
# }
hdf5.load <- function (file = NULL) {
  if (is.null(file)) stop("Please provide a HDF5 file !")
  if (!file.exists(file)) stop("Provided HDF5 file does not exist !")
  return(rhdf5::h5read(file = file, name = "/"))
}

## Function to read a csv/tsv file with data ordered HORIZONTALLY (to a df)
read.horiz.csv = function(file = NULL, header=TRUE, sep="\t", ...) {
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
}



##################
### CONVERSION ###
##################

## This function converts any factor column to a character column in a dataframe
factors2char.df <- function(df = NULL) {
  return(data.frame(rapply(as.data.frame(df), as.character, classes = "factor", how = "replace"), stringsAsFactors = FALSE))
}

## This function converts any factor column to a numeric column in a dataframe
factors2num.df <- function(df = NULL) {
  return(data.frame(rapply(as.data.frame(df), as.numeric, classes = "factor", how = "replace"), stringsAsFactors = FALSE))
}

## Converts chrom <-> chr (ie, "chr1" <-> 1) with support to alphanumerical output
## If alpha = TRUE, "chrX" -> "X", else "chrX" -> 23 (for homo sapiens)
chrConv <- function(chrvec = NULL, chr.in = "chrom", chr.out = "chr", alphanumeric = FALSE, sp = "hs") {
  if (chr.in == chr.out) stop ("Input and output types are identical !")
  if (!(chr.in %in% c("chrom", "chr"))) stop("Unsupported type for chr.in !")
  if (!(chr.out %in% c("chrom", "chr"))) stop("Unsupported type for chr.out !")
  if (length(chrvec) == 0) stop("chrvec is empty !")

  chrom2chr = list(
    hs = list("chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chr21"=21,"chr22"=22,"chrX"=23,"chrY"=24, "chrM"=25),
    mm = list("chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chrX"=20,"chrY"=21, "chrM"=22),
    rn = list("chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chrX"=21,"chrY"=22, "chrM"=23)
  )

  chr2chrom = list(
    hs = list("1"="chr1", "2"="chr2",  "3"="chr3",  "4"="chr4", "5"="chr5", "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10", "11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", "16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19", "20"="chr20", "21"="chr21", "22"="chr22", "23"="chrX", "24"="chrY", "25"="chrM"),
    mm = list("1"="chr1", "2"="chr2",  "3"="chr3",  "4"="chr4", "5"="chr5", "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10", "11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", "16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19", "20"="chrX", "21"="chrY", "22"="chrM"),
    rn = list("1"="chr1", "2"="chr2",  "3"="chr3",  "4"="chr4", "5"="chr5", "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10", "11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", "16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19", "20"="chr20", "21"="chrX", "22"="chrY", "23"="chrM")
  )

  alphaconv = list(
    hs = list("23"="X", "24"="Y", "25"="M"),
    mm = list("20"="X", "21"="Y", "22"="M"),
    rn = list("21"="X", "22"="Y", "23"="M")
  )

  if (!(sp %in% names(chrom2chr))) stop("Unsupported species !")

  if (chr.in == "chrom") {
    for (k in unique(chrvec)) {
      if (!(k %in% names(chrom2chr[[sp]]))) chrvec[chrvec == k] <- NA
      else {
        chrvec[chrvec == k] <- unlist(chrom2chr[[sp]][[k]])
      }
    }
    if (alphanumeric) {
      for (k in names(alphaconv[[sp]])) {
        chrvec[chrvec == k] <- unlist(alphaconv[[sp]][[k]])
      }
    } else chrvec <- as.numeric(chrvec)
  } else if (chr.in == "chr") {
    if (alphanumeric) {
      chrvec <- paste0("chr", chrvec)
    } else {
      for (k in unique(chrvec)) {
        if (!(k %in% names(chr2chrom[[sp]]))) chrvec[chrvec == k] <- NA
        else {
          chrvec[chrvec == k] <- unlist(chr2chrom[[sp]][[k]])
        }
      }
    }
  }
  return(chrvec)
}


#########################
### MATRIX OPERATIONS ###
#########################

## Rotates a matrix (90°, clockwise)
rotate.matrix.clockwise <- function(mat = NULL) {
  if ( !is.matrix(mat)) stop("Supplied object is not a matrix !")
  return(t(mat)[,nrow(mat):1])
}

## Scales the columns of a matrix using a given vector of values.
## This is used to scale additional samples to another quantiles-normalized matrix, using its ranks.
## Mainly used for supplemental data projection in PCA.
matrix.ranks.scaler <- function(mymatrix = NULL, ranks=NULL) {
  ranks <- sort(ranks)
  for(x in 1:ncol(mymatrix)) mymatrix[,x] <- ranks[rank(mymatrix[,x])]
  return(mymatrix)
}

## This function merges values from multiple lines of a numerical matrix according to their same rowname
## Supported methods are : median, mean, min, max
matrix.rows.merger <- function(mymatrix, method="median") {
  mymatrix <- mymatrix[order(rownames(mymatrix)),]
  mtbl <- as.data.frame(table(rownames(mymatrix), dnn = "rowname"), stringsAsFactors = FALSE)
  mtbl <- mtbl[order(mtbl$Freq),]
  outmat <- matrix(NA, nrow=nrow(mtbl), ncol=ncol(mymatrix))
  dimnames(outmat) <- list(mtbl$rowname, colnames(mymatrix))
  pnone <- mtbl$rowname[which(mtbl$Freq == 1)]
  outmat[1:length(pnone),] <- mymatrix[rownames(mymatrix) %in% pnone,]
  for (rn in (length(pnone)+1):nrow(outmat)) {
    if (method == "median") outmat[rn,] <- as.numeric(matrixStats::colMedians(mymatrix[rownames(mymatrix) == rownames(outmat)[rn], ], na.rm = TRUE))
    if (method == "mean") outmat[rn,] <- as.numeric(colMeans(mymatrix[rownames(mymatrix) == rownames(outmat)[rn], ], na.rm = TRUE))
    if (method == "min") outmat[rn,] <- as.numeric(matrixStats::colMins(mymatrix[rownames(mymatrix) == rownames(outmat)[rn], ], na.rm = TRUE))
    if (method == "max") outmat[rn,] <- as.numeric(matrixStats::colMaxs(mymatrix[rownames(mymatrix) == rownames(outmat)[rn], ], na.rm = TRUE))
  }
  outmat <- outmat[order(rownames(outmat)),]
  return(outmat)
}

## This function merges values from multiple lines of a numerical matrix according to their same rowname (aggregate version)
## Supported methods are any R function that can coerce a vector to a single value (such as) median, mean, min, max, ...)
matrix.rows.aggregator <- function(mymatrix = NULL, fun = "median") {
  if(length(duplicated(rownames(mymatrix))) == 0) {
    cat("Input matrix had no repeated rowname !\n")
    return(mymatrix)
  }
  myagg <- aggregate(x = mymatrix, by = list(rownames(mymatrix)), fun)
  myagg.mat <- as.matrix(myagg[,-1, drop = F])
  rownames(myagg.mat) <- myagg[,1]
  return(myagg.mat)
}

## Replace missing values in a numerical matrix using the median / mean / min or max of corresponding line
## Please only use this if knn imputation failed (see package "impute")
matrix.na.replace.byrow <- function(mymatrix = NULL, method = "median") {
  outmat <- t(as.matrix(vapply(1:nrow(mymatrix), function(f) {
    fout <- mymatrix[f,]
    repval <- 0
    if (method == "median") repval <- median(fout, na.rm=TRUE)
    if (method == "min") repval <- min(fout, na.rm=TRUE)
    if (method == "max") repval <- max(fout, na.rm=TRUE)
    if (method == "mean") repval <- mean(fout, na.rm=TRUE)
    fout[is.na(fout)] <- repval
    return(fout)
  }, rep(.1, ncol(mymatrix)))))
  dimnames(outmat) <- dimnames(mymatrix)
  return(outmat)
}

## COUNTS ON COLUMNS / ROWS (in complement to matrixStats)
### On columns
colWhichMax <- function(mymatrix = NULL) apply(mymatrix, 2, function(x) {which.max(x)})
colWhichMaxMult <- function(mymatrix = NULL) apply(mymatrix, 2, function(x) {which(x== max(x))})
colWhichMin <- function(mymatrix = NULL) apply(mymatrix, 2, function(x) {which.min(x)})
colWhichMinMult <- function(mymatrix = NULL) apply(mymatrix, 2, function(x) {which(x == min(x))})
colNbNA <- function(mymatrix = NULL) apply(mymatrix, 2, function(x) {length(which(is.na(x)))})
colIn <- function(mymatrix = NULL, term = NULL) apply(mymatrix, 2, function(x) {
  isin <- which(x %in% term)
  if (length(isin) == 0) isin <- 0
  return(isin)
})

### On rows
rowWhichMax <- function(mymatrix = NULL) apply(mymatrix, 1, function(x) {which.max(x)})
rowWhichMaxMult <- function(mymatrix = NULL) apply(mymatrix, 1, function(x) {which(x == max(x))})
rowWhichMin <- function(mymatrix = NULL) apply(mymatrix, 1, function(x) {which.min(x)})
rowWhichMinMult <- function(mymatrix = NULL) apply(mymatrix, 1, function(x) {which(x == min(x))})
rowNbNA <- function(mymatrix = NULL) apply(mymatrix, 1, function(x) {length(which(is.na(x)))})
rowIn <- function(mymatrix = NULL, term = NULL) apply(mymatrix, 1, function(x) {
  isin <- which(x %in% term)
  if (length(isin) == 0) isin <- 0
  return(isin)
})


#########################
### VECTOR OPERATIONS ###
#########################

## Interleave two NUMERIC vectors into one
interleave.numeric <- function(a, b) {
  mNrow <- nrow(cbind(a, b))
  m <- as.data.frame(cbind(a[1:mNrow], b[1:mNrow]))
  as.numeric(na.exclude(unlist(lapply(split(m, 1:mNrow), as.numeric))))
}

################
### PLOTTING ###
################

## Biplots for PCA results
## pca.res is an output of prcomp()
## comp is a numerical vector giving the components to biplot
pca2d <- function(pca.res=NULL, compo=1:3, class=NULL, toproj=NULL) {
  plotmat <- pca.res$x
  oripar <- par(no.readonly = TRUE)
  colvec <- 1
  rb1 <- rainbow(n=length(unique(class)))
  if (!is.null(class)) {
    class.centers <- t(sapply(sort(unique(class)), function(p) { return( colMeans(plotmat[class==p,]) ) }, simplify = T))
    colvec <- rb1[as.numeric(as.factor(class))]
  }
  par(mfrow=c(length(compo), length(compo)), mar=c(1,1,1,1), xaxt="n", yaxt="n")
  for(x in compo) {
    for(y in compo) {
      if (!is.null(class)) {
        plot(plotmat[,x], plotmat[,y], pch=20, col=colvec, xlab=paste("PC ", x, sep=""), ylab=paste("PC ", y, sep=""))
        points(class.centers[,x], class.centers[,y], pch = 23, col = 1, bg = rb1, cex = 2)
        text(class.centers[,x], class.centers[,y], labels = sort(unique(class)), col = 1, pos = 3, offset = 1, cex = 2)
      }
      else plot(plotmat[,x], plotmat[,y], pch=20, xlab=paste("PC ", x, sep=""), ylab=paste("PC ", y, sep=""))
      if (!is.null(toproj)) {
        projmat <- toproj$data
        if (is.numeric(pca.res$center)) projmat <- projmat - pca.res$center
        if (is.numeric(pca.res$scale)) projmat <- projmat / pca.res$scale
        projall <- t(projmat) %*% pca.res$rotation
        colvec2 <- vapply(toproj$class, function(x) {
          isin <- which(sort(unique(class)) %in% x)
          if (length(isin) > 0) return(isin)
          else return(1)
        }, 1)
        nprojclass <- length(unique(toproj$class))
        if (nprojclass == 1) colvec2 <- 1 else colvec2 <- rainbow(n=nprojclass)[as.numeric(as.factor(toproj$class))]
        points(projall[,x], projall[,y], pch=20, col=colvec2)
        text(projall[,x]+abs(diff(range(projall[,x])))*.05, projall[,y], labels=rownames(projall), col=colvec2, adj=0)
      }
      mtext(side = 1, text = paste("PC ", x, sep=""), line = 0.5)
      mtext(side = 2, text = paste("PC ", y, sep=""), line = 0)
      abline(v=0, lty=2)
      abline(h=0, lty=2)
    }
  }
  par(oripar)
}


## 3D plot of PCA with optional classes and projections
## pca.res is an output of prcomp()
## compo is a vector of 3 different digits, corresponding to which components to display
## toproj is a list of observations to project, with $data for data and $class for class...
pca3d <- function(pca.res=NULL, compo=1:3, class=NULL, toproj=NULL, adj=1) {
  colvec <- 1
  rb1 <- rainbow(n=length(unique(class)))
  if (!is.null(class)) colvec <- rb1[as.numeric(as.factor(class))]
  plotmat <- pca.res$x
  plot3d(plotmat[,compo[1]], plotmat[,compo[2]], plotmat[,compo[3]], type="s", size=.5, col=colvec)
  if (!is.null(class)) {
    class.centers <- t(sapply(sort(unique(class)), function(p) { return( colMeans(plotmat[class==p,]) ) }, simplify = T))
    rgl::plot3d(class.centers[,compo[1]], class.centers[,compo[2]], class.centers[,compo[3]], type="s", size=1, col=rb1, add=T)
    rgl::text3d(class.centers[,compo[1]], class.centers[,compo[2]], class.centers[,compo[3]], texts = sort(unique(class)), col=rb1, adj = adj)
  }
  if (!is.null(toproj)) {
    projmat <- toproj$data
    if (is.numeric(pca.res$center)) projmat <- projmat - pca.res$center
    if (is.numeric(pca.res$scale)) projmat <- projmat / pca.res$scale
    projall <- t(projmat) %*% pca.res$rotation
    colvec2 <- vapply(toproj$class, function(x) {
      isin <- which(sort(unique(class)) %in% x)
      if (length(isin) > 0) return(isin)
      else return(1)
    }, 1)
    colvec2 <- rainbow(n=length(unique(toproj$class)))[as.numeric(as.factor(toproj$class))]
    rgl::plot3d(projall[,compo[1]], projall[,compo[2]], projall[,compo[3]], type="p", size=10, col=colvec2, add=T)
  }
}



############
## SYSTEM ##
############

## A more robust way to get machine OS type
get.os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(tolower(os))
}



##############
## GENOMICS ##
##############

bed2gc <- function(bed = NULL, genome.fasta = NULL, species = "human", genome.build = "hg19", nt.add = 0, keep = "^chr([0-9]+|X|Y)$", ncores = 1) {

  ## Checks
  if(is.null(bed)) stop("Please provide a valid path to a BED file !")
  if(is.null(genome.fasta)) stop("Please provide a valid path to a genome FASTA file !")
  if(!is.numeric(nt.add)) stop("'nt.add' must be a numeric vector (actually an integer vector) !")
  if(any(nt.add < 0)) stop("'nt.add' must be a POSITIVE numeric vector !")
  if(any(is.infinite(nt.add))) stop("'nt.add' must be a FINITE positive numeric vector !")
  if(!(tolower(species) %in% c("human", "mouse"))) stop("'species' must be 'human' or 'mouse' !")
  if(!file.exists(bed)) stop(paste0("Could not open BED file ", bed, " !"))
  if(!file.exists(genome.fasta)) stop(paste0("Could not open genome FASTA file ", genome.fasta, " !"))
  logicores <- parallel::detectCores(logical = TRUE)
  if(ncores > logicores) print(paste0("WARNING ! ", ncores, "threads requested but only ", logicores, " logical cores available !"))

  ## Adjusting nt.add
  nt.add <- unique(nt.add)

  ## Reading BED
  print("Loading BED file ...")
  my.bed <- maelstrom::read.table.fast(file = bed, header = FALSE, sep = "\t", data.table = FALSE)
  ncol.bed <- ncol(my.bed)
  if (ncol.bed > 4) my.bed <- my.bed[,1:4] else if (ncol.bed == 4) colnames(my.bed) <- c("chr", "start", "end", "name") else colnames(my.bed) <- c("chr", "start", "end")

  ## Keeping desired chromosomes
  keep.idx <- grep(pattern = keep, x = my.bed[,1])
  if (length(keep.idx) > 0) my.bed <-my.bed[keep.idx,] else stop("No region to keep !")

  ## Preparing the BED data
  print("Sorting ...")
  my.bed.sorted <- bedr::bedr(engine = "bedtools", input = list(i = my.bed), method = "sort", params = "", check.zero.based = FALSE, check.sort = FALSE, check.merge = FALSE)

  if (ncores > 1) my.cidx <- cut(x = seq_len(nrow(my.bed.sorted)), breaks = ncores, labels = FALSE) else my.cidx <- rep.int(1, nrow(my.bed.sorted))

  ## Preparing MT
  print("Loading cluster ...")
  cl <- parallel::makeCluster(spec = ncores, type = "FORK")
  doParallel::registerDoParallel(cl)

  print("Computing GC tracks ...")
  nta <- 0
  myXgc <- foreach::foreach(nta = nt.add, .combine = "cbind") %do% {
    if (nta > 0)  my.bed.grown <- bedr::grow.region(x = my.bed.sorted, n.add = nta, species = species, build = genome.build, check.zero.based = FALSE, check.sort = FALSE, check.merge = FALSE, verbose = FALSE) else my.bed.grown <- my.bed.sorted
    pos.valid.chk <- bedr::is.valid.region(my.bed.grown, check.zero.based = FALSE)
    if (!all(pos.valid.chk)) {
      if (!any(pos.valid.chk)) stop("No valid region found !")
        my.bed.grown <- my.bed.grown[pos.valid.chk,]
    }
    cidx <- 0
    mygc <- foreach::foreach(cidx = unique(my.cidx), .combine = "c") %dopar% {
      my.bed.grown.c <- my.bed.grown[my.cidx == cidx,]
      my.gc <- bedr::bedr(engine = "bedtools", input = list(bed = my.bed.grown.c), method = "nuc", params = paste0("-fi ", genome.fasta), check.merge = FALSE, check.sort = FALSE, check.valid = FALSE, verbose = FALSE)
      my.gc <- my.gc[-1,]
      gc.out <- as.numeric(my.gc[,5])
      na.idx <- which(as.numeric(my.gc[,11]) == my.gc[,3] - my.gc[,2])
      if (length(na.idx) > 0) gc.out[na.idx] <- NA
      return(gc.out)
    }
    return(mygc)
  }

  print("Stopping cluster ...")
  stopCluster(cl)

  colnames(myXgc) <- paste0(nt.add, "bp")

  mypos <-cbind(my.bed.sorted, myXgc)

  ## Output
  print("Writing GC data ...")
  maelstrom::write.table.fast(x = mypos, file = sub(pattern = "\\.bed$", replacement = "_GC.txt", x = bed, ignore.case = TRUE), header = TRUE, sep = "\t")
  print("Done.")
}
