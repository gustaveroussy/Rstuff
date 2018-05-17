chrLoad <- function(genome = "hg19") {
  if(is.null(genome)) stop("Please provide a genome build (by example : 'hg19') !")
  try(data(genome, envir = environment()))
}

chrom2chr <- function(cs = NULL, chromvalues=NULL) {
  if (is.null(cs)) stop("Please provide a 'cs' object !")
  if (is.null(chromvalues)) stop("Please provide a vector of Chrom values !")
  if(any(!is.character(chromvalues))) stop("the vector of Chrom values should only contain character data !")
  return(unlist(cs$chrom2chr[chromvalues]))
}

chr2chrom <- function(cs = NULL, chrvalues=NULL) {
  if (is.null(cs)) stop("Please provide a 'cs' object !")
  if (is.null(chrvalues)) stop("Please provide a vector of Chr (numeric) values !")
  if(any(!is.numeric(chrvalues))) stop("the vector of Chr values should only contain numeric data !")
  return(unlist(cs$chr2chrom[chrvalues]))
}