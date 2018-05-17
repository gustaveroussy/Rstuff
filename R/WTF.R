## WTF.R
##
## DESCRIPTION :  This script performs a Wilcoxon rank test (W), a students' T-test (T),
##                and/or a Fisher's exact test (F) on a continuous variable and a class.
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.12.0 to 2.13.0
##
## DEPENDS ON:
## -
##
## VERSION NOTES
##
## v1.3 20120824
##      . Added more controls : Check if query or target columns do exist, if there are no duplicate.
##      . Now the script can handle
##      . Added a help.chifisher function.
##
## v1.2b 20120201
##      . Corrected a bug (annot variable badly named) that prevented the script to run at all.
##
## 1.2  . Wraped the script in a function that can be easily sourced.
##      . Translated descriptions and comments to US English.
##
## 1.1  . Wrapped the code into an commitable script
##
## 1.0 20110303
##      . First release for P19_AS_GE.


wtf <- function(annot.table, continues, classes) {

  ## Loading annotations
  annot <- read.table(annot.table, header=T, sep="\t", as.is=T, check.names=F)
  ## Removing dupes in query and target vectors
  continues <- unique(continues)
  classes <- unique(classes)
  final.W <- final.T <- final.F <- matrix(NA, ncol=length(continues), nrow=length(classes))
  
  for (k1 in 1:length(continues)) {
    colclass <- continues[k1]
    ccnum <- which(colnames(annot) == colclass)
    if (length(ccnum) == 0) {
      cat("Resquested continue column [", colclass, "] could not be found !\n", sep="")
      next
    }
    if (length(ccnum) > 1) {
      cat("Requested continue column ", colclass ," exists ", length(ccnum), " times !\nCan't choose ...\n", sep="")
      next
    } else {
      
      tpv.W <- tpv.T <- tpv.F <- vector()
      for (k2 in 1:length(classes)) {
        ## Dismiss if query and target are the same...
        if (continues[k1] == classes[k2]) {
          tpv.W <- c(tpv.W, NA)
          tpv.T <- c(tpv.T, NA)
          tpv.F <- c(tpv.F, NA)
          next
        }
        ## Else :
        cat("Crossing [", continues[k1], "] with [", classes[k2], "] ...\n", sep="")
        
        ct <- classes[k2]
        ctnum <- which(colnames(annot) == ct)
        if (length(ctnum) == 0) {
          cat("Resquested class column [", ct, "] could not be found !\n", sep="")
          tpv.W <- c(tpv.W, NA)
          tpv.T <- c(tpv.T, NA)
          tpv.F <- c(tpv.F, NA)
          next
        }
        if (length(ctnum) > 1) {
          cat("Requested class column ", ct, "exists ", length(ctnum), " times !\nCan't choose ...\n", sep="")
          tpv.W <- c(tpv.W, NA)
          tpv.T <- c(tpv.T, NA)
          tpv.F <- c(tpv.F, NA)
          next
        } else {
          if (length(unique(annot[,ctnum])) < 2) {
            cat("Requested class only has ", length(unique(annot[ctnum])), "exists ", length(ctnum), " times !\nCan't choose ...\n", sep="")
            tpv.W <- c(tpv.W, NA)
            tpv.T <- c(tpv.T, NA)
            tpv.F <- c(tpv.F, NA)
          }
          if (length(unique(annot[,ctnum])) == 2) {
            wt <- wilcox.test(annot[,ccnum] ~ annot[,ctnum])
            tpv.W <- c(tpv.W, wt$p.value)
            tt <- t.test(annot[,ccnum] ~ annot[,ctnum])
            tpv.T <- c(tpv.T, tt$p.value)
            slmt <- summary(lm(annot[,ccnum] ~ annot[,ctnum]))
            tpv.F <- c(tpv.F, slmt$coefficients[2,4])
          }
          else {
            tpv.W <- c(tpv.W, NA)
            tpv.T <- c(tpv.T, NA)
            tpv.F <- c(tpv.F, NA)
          }
        }
      }
      cat("FINALW:", dim(final.W), ", TPVW:", length(tpv.W),"\n")
      final.W[,k1] <- tpv.W
      final.T[,k1] <- tpv.T
      final.F[,k1] <- tpv.F
    }
  }
  write.table(rbind(c("W-tests", continues), cbind(classes, final.W)), paste("W-TEST_globaltable.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
  write.table(rbind(c("T-tests", continues), cbind(classes, final.T)), paste("T-TEST_globaltable.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
  write.table(rbind(c("F-tests", continues), cbind(classes, final.F)), paste("F-TEST_globaltable.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
}
