setwd("/mnt/data_colibri/pelican/volume03/ws/ugf/Projets/Osteosarcome/Projets_Microarrays/P26_OSTE/04_Bioinformatique/CGH/MIXED_JOINT_GC5/OS2K6_ADMIX4/OS2K6_ADMIX4_164s_hg19_20171201_noY/MCR")

source("/home/job/svn/genomics/CGH/R/chrload.R")
cs <- chrload(sp = "hs", gb = 19)

require(data.table)
ingistic <- fread(input = "scores.gistic", data.table = FALSE)
# ingistic$ChrN <- unlist(cs$chrom2chr[ingistic$Chromosome])
ingistic$genostart <- cs$lchrsum[ingistic$Chromosome] + ingistic$Start
ingistic$genoend <- cs$lchrsum[ingistic$Chromosome] + ingistic$End

midchr <- cs$lchrsum[-length(cs$lchrsum)] + (diff(cs$lchrsum)/2)

chr.max <- length(unique(ingistic$Chromosome))
amp.idx <- ingistic$Type == "Amp"
midchr <- cs$lchrsum + diff(c(cs$lchrsum, cs$glen)) / 2

png("GISTIC2.png", 1500, 950)
par(cex = 1.5, mai = c(1.5,1.5,1,0.25))
my.ylim <- 30
plot(0,0, type = "n", ylim = c(-1,1)*my.ylim, xlim = c(0, c(cs$lchrsum, cs$genome.length)[chr.max+1]), xaxs = "i", yaxs = "i", xaxt = "n", xlab = "Genomic position", ylab = "(-/+) log10 of GISTIC2 Q-value for losses / gains", main = "GISTIC2 RESULTS")
rect(cs$lchrsum[seq.int(2, length(cs$lchrsum), 2)], -my.ylim, c(cs$lchrsum[seq.int(3, length(cs$lchrsum), 2)], cs$glen), my.ylim, col = "grey75", border = "transparent")
abline(v = cs$centromere.geno, lty = 2, col = "grey50")
lines(ingistic$genostart[amp.idx], ingistic[["-log10(q-value)"]][amp.idx], type = "l", col = 4, lwd = 2)
lines(ingistic$genostart[!amp.idx], -ingistic[["-log10(q-value)"]][!amp.idx], type = "l", col = 2, lwd = 2)
abline(h = 0, lwd = 2, col = "grey50")
abline(h = c(1,-1)*log10(25), col = 3, lwd = 2, lty = 3)
# text(x = cs$m, y = my.ylim * .9 * (seq_len(length(unique(ingistic$Chromosome)))%%2 * 2 - 1), labels = sub(pattern = "chr", replacement = "", x = unlist(cs$chr2chrom[unique(ingistic$Chromosome)])))
text(x = midchr, y = my.ylim * .9 * (seq_len(length(unique(ingistic$Chromosome)))%%2 * 2 - 1), labels = cs$chrom[1:chr.max], cex = .5)
dev.off()
