## CHIFISHER.R
##
## DESCRIPTION :  This script performs a series of Fisher / X2  or Wilcoxon / Kruskal-Wallis tests
##                on selected query / target columns from a tab-separated annotation file.
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.12.0 to 3.3.2
##
## DEPENDS ON:
## -
##
## VERSION NOTES
##
## v1.6 20171129
##      . Added classes summaries for T/W/ANOVA/KW tests (to have numerical data represented in the violin plots).
##
## v1.5b 20161104
##      . Removed the "simulate" option introduced in v1.4 : now the use of a simulated p-value
##        computation for a Fisher test is automatically triggered in case of error in the
##        normal mode, thanks to a try() call.
##
## v1.5 20160615
##      . Added support for 2-class and multi-class continous variables (T/W and ANOVA/K-W)
##
## v1.4 20130613
##      . Added a "simulate" option for cases giving a "FEXACT error 7". Default B is set to 1E+06.
##        Set to FALSE by default.
##      . Added a "star code" to the STDOUT to quickly show if some positive results were found.
##        5 "*" for p<1E-05, 3 "*" for p<1E-03, 1 "*" for p<5E-02. Also 1 "~" for p<1E-01.
##
## v1.3b 20121005
##      . Just added a "check.names=F" to the read.table(annotfile), to handle column names with
##        classical special characters (mainly "-").
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
## 1.1	. Fixed a bug in NA handling.
##
## 1.0	. First release.


## CHIFISHER
##
## Options :
##    annot.table : path+filename to the annotation table
##    query.vec : a vector of column names to use as queries (typicaly the annotations from the unsupervised analyses)
##    target.vec : a vector of column names to use as targets (typicaly the clinical annotations)
##    test.type : the type of the statistical test to perform. Either "F" or "X2"
##     very disequilibrated contingency tables, when obtained a "FEXACT error 7".


## PARAMETERS DESCRIPTOR FUNCTION
help.chifisher <- function() {
  cat("
PARAMETERS FOR THE chifisher() FUNCTION :
-----------------------------------------
annot.table   = Annotations table (tab-separated plain text).
query.vec     = Vector of query column names (ex: Hclust classes at different cuts).
target.vec    = Vector of test column names (ex: clinical annotation classes).
test.type     = Which test to use ('F' for Fisher, 'X2' for X-squared).
test.type.2 = type of statistical test for 2-classes variable (W : Wilcoxon, or T : Student's T-test).
test.type.N = type of statistical test for N-classes (ANOVA : Analysis of variance,  or KW : Kruskal-Wallis).
out.dir     = where to write results (tables, plots)
make.plot   = try to perform plots ? (in some R installations, plots fail so ...)
sim.p = Perform permutation-based p-value approximation
n.permut = if sim.p == TRUE, number of permutations to perform
\n")
}

## MAIN FUNCTION
chifisher <- function(annot.table = NULL, query.vec = NULL, target.vec = NULL, test.type = "F", test.type.2 = "W", test.type.N = "KW", numeric.as.continuous = FALSE, out.dir = getwd(), make.plot = TRUE, sim.p = TRUE, n.permut = 1E+6) {
  ## Loading annotations
  sfix <- ""
  # if (simulate) sfix <- "_SIMU"

  oridir <- getwd()

  annot <- read.table(annot.table, header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE, comment.char = "", quote = "")
  ## Removing dupes in query and target vectors
  query.vec <- unique(query.vec)
  target.vec <- unique(target.vec)
  final <- matrix(NA, nrow=length(target.vec), ncol=length(query.vec))
  for (k1 in 1:length(query.vec)) {
    colclass <- query.vec[k1]
    ccnum <- which(colnames(annot) == colclass)
    if (length(ccnum) == 0) {
      cat("Resquested query column [", colclass, "] could not be found !\n", sep="")
      next
    }
    if (length(ccnum) > 1) {
      cat("Requested query column ", colclass ," exists ", length(ccnum), " times !\nCan't choose ...\n", sep="")
      next
    } else {
      classes <- as.factor(annot[,ccnum])
      ncateg <- nlevels(classes)
      tpv <- vector()
      # outname <- paste(test.type, "-TEST_", colclass, sfix,".txt", sep="")
      for (k2 in 1:length(target.vec)) {

        ## Dismiss if query and target are the same...
        if (query.vec[k1] == target.vec[k2]) {
          tpv <- c(tpv, NA)
          next
        }
        ## Else :
        message(paste0("Crossing [", query.vec[k1], "] with [", target.vec[k2], "] ..."))
        ct <- target.vec[k2]
        ctnum <- which(colnames(annot) == ct)
        # message(ctnum)
        if (length(ctnum) == 0) {
          message(paste0("Resquested target column [", ct, "] could not be found !"))
          tpv <- c(tpv, NA)
          next
        }
        if (length(ctnum) > 1) {
          message(paste0("Requested target column ", ct, "exists ", length(ctnum), " times !\nCan't choose ..."))
          tpv <- c(tpv, NA)
          next
        } else {

          wdir <- paste0(out.dir, '/', query.vec[k1], "/", target.vec[k2])
          dir.create(path = wdir, recursive = TRUE, showWarnings = FALSE)

          if (numeric.as.continuous & (is.numeric(annot[,ctnum]))) {

            test.type.X <- ""

            if (ncateg == 2) {
              test.type.X <- test.type.2

              ## T-test
              if (test.type.2 == "T") {
                Tres <- t.test(annot[,ctnum] ~ classes)
              }
              ## Wilcoxon
              if ((ncateg == 2) & (test.type.2 == "W")) {
                Tres <- wilcox.test(annot[,ctnum] ~ classes)
              }
            }
            if (ncateg > 2) {
              # message("here")
              test.type.X <- test.type.N

              if (test.type.N == "ANOVA") {
                ## ANOVA
                test.lm <- lm(annot[,ctnum] ~ classes)
                Tres <- anova(test.lm)
                Tres$p.value <- Tres$"Pr(>F)"[1]
              }
              ## Kruskal-Wallis
              if ((ncateg > 2) & (test.type.N == "KW")) {
                Tres <- kruskal.test(annot[,ctnum] ~ classes)
              }
            }

            # message("here")
            outname <- paste(wdir, "/", test.type.X, "-TEST_", colclass, "_", target.vec[k2], sfix, sep="")

            # class.summary <- as.data.frame(t(sapply(unique(classes[!is.na(classes)]), function(x) { c(summary(annot[classes == x, ctnum]), length(which(classes == x)), sd(annot[which(classes == x), ctnum]), Tres$statistic, Tres$p.value) }, simplify = TRUE)))
            class.summary <- as.data.frame(t(sapply(unique(classes[!is.na(classes)]), function(x) {
              smry.res <- summary(annot[classes == x, ctnum])
              smry.res <- smry.res[names(smry.res) != "NA's"]
              return(c(smry.res, length(which(classes == x)), sd(annot[which(classes == x), ctnum]), Tres$statistic, Tres$p.value) )
              }, simplify = TRUE)))
            # return(class.summary)
            # colnames(class.summary) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N", "Sd", paste0(test.type.X, c(".score", ".p.value")))
            colnames(class.summary)[1:6] <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")
            csncol <- ncol(class.summary)
            colnames(class.summary)[(csncol-3):csncol] <- c("N", "Sd", paste0(test.type.X, c(".score", ".p.value")))
            rownames(class.summary) <- paste0(colclass, ".", unique(classes[!is.na(classes)]))
            # rownames(class.summary) <- unique(classes[!is.na(classes)])
            # rownames(class.summary) <- paste0(colnames(annot)[ctnum], ".", rownames(class.summary))


            rcon <- file(paste0(outname, ".txt"), "w")
            writeLines(paste(c(colnames(annot)[ctnum], colnames(class.summary)), collapse = "\t"), rcon)
            close(rcon)
            write.table(class.summary, file = paste0(outname, ".txt"), quote = FALSE, row.names = TRUE, col.names = FALSE, sep="\t", append = TRUE)

            if (make.plot) {
              mycounts <- vapply(sort(unique(annot[[query.vec[k1]]])), function(x) { length(which(annot[[query.vec[k1]]] == x))}, 1)
              library(ggplot2)
              ymin <- min(annot[[target.vec[k2]]], na.rm = TRUE) - (diff(range(annot[[target.vec[k2]]], na.rm = TRUE))*.05)
              p <- ggplot2::qplot(factor(annot[[query.vec[k1]]]), annot[[target.vec[k2]]], geom = "violin", trim = FALSE, scale = "count", fill = factor(annot[[query.vec[k1]]]), xlab = query.vec[k1], ylab = target.vec[k2]) + geom_boxplot(width=.2) + ggplot2::geom_jitter(width = 0) + ggplot2::guides(fill=FALSE) + annotate("text", x = sort(unique(annot[[query.vec[k1]]])), y = ymin, label = mycounts, vjust = 2)
              ggplot2::ggsave(filename = paste0(outname, "_gg.png"), width = 15, height = 15, units = "cm", plot = p)
            }

            saveRDS(Tres, file = paste0(outname, ".RDS"))

          } else {

            outname <- paste(wdir, "/", test.type, "-TEST_", colclass, sfix, sep="")
            groups <- as.factor(annot[,ctnum])

            vtm <- clab <- tlab <- vector()
            # clab <- paste(colclass, levels(classes), sep=".")
            # tlab <- paste(ct, levels(groups), sep=".")
            for (cclev in levels(classes)) {
              ind <- which(classes == cclev)
              for (ctlev in levels(groups)) {
                vtm <- c(vtm, length(which(groups[ind] == ctlev)))
              }
            }
            mymat <- matrix(vtm, nr=length(levels(groups)), dimnames=list(c(levels(groups)), c(levels(classes))))
      			names(dimnames(mymat)) <- c(target.vec[k2], colclass)

      			if (test.type == "F") {
      			  Tres <- try(fisher.test(mymat), silent = TRUE)
      			  if (is(Tres) == "try-error") {
      			    message(" Direct p-value computation failed. Trying in simulated mode (MCMC).")
      			    Tres <- fisher.test(mymat, simulate.p.value=sim.p, B=n.permut)
      			  }
      			} else if (test.type == "X2") Tres <- chisq.test(mymat, simulate.p.value = sim.p, B = n.permut) else stop("Unkown test!")

      			if(make.plot) {
      			library(vcd)
        			pdf(paste0(outname, "_", target.vec[k2], ".pdf"), width = 21/cm(1), height = 21/cm(1))
        			try(vcd::mosaic(mymat, shade = TRUE))
        			try(vcd::assoc(mymat, shade = TRUE))
        			dev.off()
      			}
      			clab <- paste(colclass, levels(classes), sep=".")
      			tlab <- paste(ct, levels(groups), sep=".")
      			dimnames(mymat) <- list(c(tlab), c(clab))

      			write.table( rbind(c(colclass, levels(classes)), cbind(rownames(mymat), mymat)), paste0(outname, ".txt"), quote=F, row.names=F, col.names=F, sep="\t", append=T)
      			if (test.type == "F") write.table(paste("\n",Tres$method, ":\np-value = ", Tres$p.value, "\n\n\n", sep=""), paste0(outname, ".txt"), append=T, sep="\t", quote=F, col.names=F, row.names=F)
      			if (test.type == "X2") write.table(paste("\n",Tres$method, ":\nX2 statistic = ", Tres$statistic, ", p-value = ", Tres$p.value, "\n\n\n", sep=""), paste0(outname, ".txt"), append=T, sep="\t", quote=F, col.names=F, row.names=F)
      		}

    			if (Tres$p.value < 1E-05) message("*****") else if (Tres$p.value < 1E-03) message("***") else if (Tres$p.value < 5E-02) message("*") else if (Tres$p.value < 1E-01) message("~")

    			tpv <- c(tpv, Tres$p.value)
    		}
    	}
  	  final[,k1] <- tpv
  	}
  }
  write.table(rbind(c(paste(test.type, "tests", sep="-"), query.vec), cbind(target.vec, final)), file = paste0(out.dir, '/', test.type, "-", test.type.2, "-", test.type.N, "-TEST_globaltable", sfix, ".txt"), quote=F, row.names=F, col.names=F, sep="\t")
}

