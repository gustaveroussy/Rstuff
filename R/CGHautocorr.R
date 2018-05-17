setwd("/mnt/toucan_seq3/B16_AUVI/GC5/cbs_cs_a1E-04/B16_AUVI_cs_a1-04_50s_hg19_20161017_noX_noY/")
require(maelstrom)
myreg <- read.table.fast(file = "/mnt/toucan_seq3/B16_AUVI/GC5/cbs_cs_a1E-04/B16_AUVI_cs_a1-04_50s_hg19_20161017_noX_noY/Data/B16_AUVI_cs_a1-04_50s_hg19_20161017_noX_noY.xreg")


## PARAMETERS
cor.method <- "pearson"
adjust.method <- "BH"
p.min <- 1E-03


## Preparing data
mymat <- as.matrix(myreg[,-c(1:9)])
rownames(mymat) <- paste0(myreg$Chrom, ":", myreg$Start, "-", myreg$End)

## Computing correlation and p-values
require(Hmisc)
mycorp <- rcorr(t(mymat), type = cor.method)

## Adjusting
mycorp$Pa <- matrix(data = p.adjust(mycorp$P, method = "BH"), nrow = nrow(mycorp$P))
dimnames(mycorp$Pa) <- dimnames(mycorp$P)


require(RColorBrewer)
mypal <- colorRampPalette(colors = c("red", "white", "blue"))
png("test2.png", width = 5000, height = 5000)
image(x = 1:nrow(mymat), y = 1:nrow(mymat), z = mycorp$r, zlim = c(-1,1), col = mypal(100), xlab = '', ylab = "")
non.sig <- which(mycorp$Pa > p.min)
if (length(non.sig) > 0) points(non.sig %/% nrow(mymat), non.sig %% nrow(mymat), pch = 20, col = "white", cex = .1)
chrbrk <- cumsum(rle(myreg$Chrom)$lengths)
abline(v = chrbrk, col = 1, lwd = 5)
abline(h = chrbrk, col = 1, lwd = 5)
dev.off()

for (x in 1:ncol(mycorp$Pa)) {
  overp <- which(mycorp$Pa[,x] > p.min)
  if (length(overp) > 0) points(overp, x, pch = 20, col = 2)
}
chrbrk <- cumsum(rle(myreg$Chrom)$lengths)
abline(v = chrbrk, col = 1, lwd = 5)
abline(h = chrbrk, col = 1, lwd = 5)
dev.off()


# require(corrplot)
# png("cortest.png", width = 8192, height = 8192)
# corrplot(corr = mycorp$r, type = "lower", method = "color", tl.pos = "n", p.mat = mycorp$Pa, sig.level = 1E-03, pch.cex = .1, pch = 20)
# dev.off()

