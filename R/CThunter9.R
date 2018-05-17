setwd("/home/job/svn/genomics/CGH/R/CThunter")
require(foreach)

my.version <- 9
### V8 *Vroom!* :   Gap/Edge min size cutter introduced.
### V7 *Up!* :      Elongated chromosomes, per-chr segmentation (fixed from v6 !), inversion rate introduced


## Params
print("Defining parameters ...")
datafile <- "/mnt/data_cigogne/job/cna_joint_admixed_analysis/04_most_centered_peak_centering/multipcf_ADMIX3_sub_GC100000_l2r_124s/copynumber_segmentation_centered_mpc_mcpc.RDS"
datafile.withflat <- "/mnt/data_cigogne/job/cna_joint_admixed_analysis/04_most_centered_peak_centering/multipcf_ADMIXED_l2r_140s/copynumber_segmentation_centered_mpc_mcpc.RDS"
flat.samplenames <- c("53730_1_2", "53730_1_3", "53731_1_4", "53732_1_3", "53957_1_4", "54011_1_3", "54032_1_1", "54032_1_4", "53727_1_1", "53903_1_3", "54006_1_1")
cutfactor <- 2
my.width <- 1E+08
my.adjust <- .5
n <- 2^10
my.pad <- .1
inv.rate <- .5
nseg.min <- 10
width.min <- 0
gapedge.max <- 1E+07

## Data loading
print("Loading data ...")
data("hg19", package = "chromosomes", envir = environment())
cgh.data <- readRDS(datafile)
cgh.data.withflat <- readRDS(datafile.withflat)

chrend.idx <- vapply(unique(cgh.data$cn.data.ini$Chr), function(k) { max(which(cgh.data$cn.data.ini$Chr == k)) }, 1)

l2r.mat <- t(t(as.matrix(cgh.data$cn.data.wins[,-c(1:2)])) - cgh.data$most.centered.peak)
flat.idx <- which(colnames(cgh.data.withflat$cn.data.ini) %in% flat.samplenames)
l2r.mat.flats <- t(t(as.matrix(cgh.data.withflat$cn.data.wins[,flat.idx])) - cgh.data.withflat$most.centered.peak[flat.idx - 2])
l2r.mat <- cbind(l2r.mat, l2r.mat.flats)
l2r.mat <- l2r.mat[,order(colnames(l2r.mat))]

coords.df <- data.frame(Chr = cgh.data$cn.data.ini$Chr, Start = cgh.data$cn.data.ini$Start, CumStart = cgh.data$cn.data.ini$Start + cs$chromosomes$chr.length.toadd[cgh.data$cn.data.ini$Chr], stringsAsFactors = FALSE)
medprobspace <- round(median(abs(diff(coords.df$CumStart))))

kbrk <- vapply(unique(coords.df$Chr), function(x) { coords.df$CumStart[max(which(coords.df$Chr == x))] }, 1)
klen <- vapply(unique(coords.df$Chr), function(x) { coords.df$Start[max(which(coords.df$Chr == x))] }, 1)
ktoadd <- c(0, cumsum(klen)[-length(klen)])

## SEGMENTING with per-chr with pseudo-circularization
seg.list2 <- sapply(1:ncol(l2r.mat), function(s) {
  print(s)
  k <- 1
  kseg.list <- sapply(unique(coords.df$Chr), function(k) {
    kseg.idx <- which(coords.df$Chr == k)
    kidx.range <- range(kseg.idx)
    pad.add <- round(length(kseg.idx) * my.pad)
    lefty <- (max(kseg.idx)-pad.add):max(kseg.idx)
    righty <- min(kseg.idx):(min(kseg.idx)+pad.add)
    new.idx <- c(lefty, kseg.idx, righty)
    kloc.middle <- coords.df$CumStart[kseg.idx]
    kloc.left <- coords.df$CumStart[lefty] - max(coords.df$CumStart[lefty]) + min(kloc.middle) - min(coords.df$Start[coords.df$Chr == k])
    kloc.right <- coords.df$Start[righty] + max(kloc.middle)
    kloc.new <- c(kloc.left, kloc.middle, kloc.right)
    
    pad.l2r <- l2r.mat[new.idx,s]
    
    peltres <- changepoint::cpt.meanvar(data = pad.l2r , penalty = "MBIC", method = "PELT", param.estimates = FALSE, minseglen = 3)@cpts
    
    segdenDEF <- density(kloc.new[peltres], adjust = my.adjust, width = my.width, n = n, from = min(coords.df$CumStart[coords.df$Chr == k]), to = max(coords.df$CumStart[coords.df$Chr == k]), cut = 3)
    segdenDEF.df <- data.frame(x = segdenDEF$x, chr = k, y = segdenDEF$y, stringsAsFactors = FALSE)

    trueseg <- peltres - pad.add
    trueseg <- trueseg[trueseg > 0 & trueseg < length(kseg.idx)]
    trueseg <- kseg.idx[trueseg]
    trueseg <- c(trueseg, max(kseg.idx))
    
    return(list(seg = trueseg, den = segdenDEF.df, nseg.ori = length(peltres)))

  }, simplify = FALSE)
  names(kseg.list) <- unique(coords.df$Chr)
  return(kseg.list)
}, simplify = FALSE)
names(seg.list2) <- colnames(l2r.mat)

segz <- sapply(1:length(seg.list2), function(s) {
  foreach(k = 1:length(seg.list2[[s]]), .combine = "c") %do% {
    return(seg.list2[[s]][[k]]$seg)
  }
}, simplify = FALSE)
names(segz) <- colnames(l2r.mat)

seg.dfz <- sapply(1:length(segz), function(s) {
  seg.df <- data.frame(end.idx = segz[[s]], stringsAsFactors = FALSE)
  seg.df$start.idx <- c(1, seg.df$end.idx[-nrow(seg.df)] +1)
  seg.df$sd <- vapply(1:nrow(seg.df), function(x){ sd(l2r.mat[seg.df$start.idx[x]:seg.df$end.idx[x],s], na.rm = TRUE)}, .1)
  seg.df$median <- vapply(1:nrow(seg.df), function(x){ median(l2r.mat[seg.df$start.idx[x]:seg.df$end.idx[x],s], na.rm = TRUE)}, .1)
  seg.df$chr <- coords.df$Chr[seg.df$start.idx]
  seg.df$start <- coords.df$Start[seg.df$start.idx]
  seg.df$end <- coords.df$Start[seg.df$end.idx]
  seg.df$cumstart <- coords.df$CumStart[seg.df$start.idx]
  seg.df$cumend <- coords.df$CumStart[seg.df$end.idx]
  seg.df$width <- seg.df$end - seg.df$start +1
  return(seg.df)
}, simplify = FALSE)
names(seg.dfz) <- colnames(l2r.mat)

nseg.oriz <- vapply(1:length(seg.list2), function(s) {
  nall <- foreach(k = 1:length(seg.list2[[s]]), .combine = "c") %do% {
    return(seg.list2[[s]][[k]]$nseg.ori)
  }
  return(sum(nall))
}, 1)
names(nseg.oriz) <- colnames(l2r.mat)

denz <- sapply(1:length(seg.list2), function(s) {
  foreach(k = 1:length(seg.list2[[s]]), .combine = "rbind") %do% {
    return(seg.list2[[s]][[k]]$den)
  }
}, simplify = FALSE)
names(denz) <- colnames(l2r.mat)


## Per-chr adjustment of densities heights
denz.adj <- sapply(1:length(seg.list2), function(s) {
  foreach(k = 1:length(seg.list2[[s]]), .combine = "rbind") %do% {
    temp.df <- seg.list2[[s]][[k]]$den
    temp.df$y.adj <- temp.df$y * length(seg.list2[[s]][[k]]$seg)
    temp.df$y.adj2 <- temp.df$y * seg.list2[[s]][[k]]$nseg.ori
    return(return(temp.df))
  }
}, simplify = FALSE)
names(denz.adj) <- colnames(l2r.mat)

## Adding chromosomal affiliation to densities' breaks
den.kend <- vapply(1:length(kbrk), function(k) { max(which(denz.adj[[1]]$x <= kbrk[k])) }, 1L)
den.kstart <- c(1,c(den.kend[-length(den.kend)] +1))

expstep <- max(coords.df$CumStart) / length(unlist(segz))
expvec <- seq(1, max(coords.df$CumStart), expstep)
expavg <- length(expvec) / ncol(l2r.mat)
expden <- density(expvec, adjust = my.adjust, width = my.width, n = n, from = 0, to = max(coords.df$CumStart), cut = 3)
expden.y.adj <- expden$y*expavg*cutfactor
cutlim <- max(expden.y.adj)

plop <- maelstrom::rasterpdf.open(file = paste0("CThunter", my.version, "v1x_cf", cutfactor, "_invRate", inv.rate, "_nsmin", nseg.min, "_wmin", width.min, "_gem", gapedge.max, "_", ncol(l2r.mat), "s.pdf"), width = 2000, height = 1000)
par(mfrow = c(2,1), oma = c(2,2,2,2), mai = c(.1,.1,.2,0))

ct.res <- foreach(samp.idx = 1:ncol(l2r.mat), .combine = "rbind") %do% {
  
  samp.name <- names(seg.dfz)[samp.idx]
  print(paste0(samp.name, " [", samp.idx, "]"))
  
  
  seg.df <- seg.dfz[[samp.idx]]
  
  pos.idx <- t(sapply(1:length(kbrk), function(k) { as.numeric(denz.adj[[samp.idx]]$y.adj[den.kstart[k]:den.kend[k]] >= cutlim) }))
  
  ## Handling smmall segments (for gaps/edges)
  for(k in unique(seg.df$chr)) {
    seg.df.k <- seg.df[seg.df$chr == k,]
    smallies.idx <- which(seg.df.k$width > gapedge.max)
    if(length(smallies.idx) > 0) {
      denz.adj.k <- denz.adj[[samp.idx]][denz.adj[[samp.idx]]$chr == k,]
      den.idx <- foreach(d = 1:length(smallies.idx), .combine = "c") %do% { which(denz.adj.k$x >= seg.df.k$cumstart[smallies.idx][d] & denz.adj.k$x <= seg.df.k$cumend[smallies.idx][d]) }
      pos.idx[k,den.idx] <- 0
    }
  }
  
  posreg.flag <- FALSE
  posreg.nb <- 0
  if (sum(pos.idx) > 0) {
    # print("    *")
    posreg.flag <- TRUE
    
    pos.df <- foreach(k = 1:nrow(pos.idx), .combine = "rbind") %do% {
      pos.df.k <- data.frame(value = numeric(), idx.length = numeric(), idx.end = numeric(), idx.start = numeric(), idx.cumstart = numeric(), idx.cumend = numeric(), cumstart = numeric(), cumend = numeric(), width = numeric(), cumstart = numeric(), cumend = numeric(), width = numeric(), stringsAsFactors = NA)
      if(all(pos.idx[k,] == 0)) return(pos.df.k)
      # print(k)
      pos.rle <- rle(pos.idx[k,])
      pos.df.k <- data.frame(samplename = samp.name, chr = k, value = pos.rle$values, idx.length = pos.rle$lengths, idx.end = cumsum(pos.rle$lengths), stringsAsFactors = FALSE)
      if (nrow(pos.df.k) == 1) pos.df.k$idx.start <- 1 else pos.df.k$idx.start <- c(1, pos.df.k$idx.end[-nrow(pos.df.k)]+1)
      pos.df.k <- pos.df.k[pos.df.k$value == 1,, drop = FALSE]
      pos.df.k$value <- NULL
      pos.df.k$idx.cumstart <- pos.df.k$idx.start + ((k-1) * n)
      pos.df.k$idx.cumend <- pos.df.k$idx.end + ((k-1) * n)
      pos.df.k$cumstart <- denz.adj[[samp.idx]]$x[pos.df.k$idx.cumstart]
      pos.df.k$cumend <- denz.adj[[samp.idx]]$x[pos.df.k$idx.cumend]
      pos.df.k$start <- pos.df.k$cumstart - ktoadd[pos.df.k$chr]
      pos.df.k$end <- pos.df.k$cumend - ktoadd[pos.df.k$chr]
      pos.df.k$width <- pos.df.k$cumend - pos.df.k$cumstart + 1
      return(pos.df.k)
    }
    posreg.nb <- nrow(pos.df)
  }
  
  seg.df$segpos <- 0
  # segpos.idx <- foreach(p = 1:nrow(pos.df), .combine = "c") %do% { which(seg.df$cumstart <= pos.df$cumend[p] & seg.df$cumend >= pos.df$cumstart[p]) }
  
  midpoints <- (seg.df$cumstart + seg.df$cumend) / 2
  segdiff <- data.frame(diff = diff(seg.df$median), stringsAsFactors = FALSE)
  segdiff$sign <- sign(segdiff$diff)
  segdiff$cumstart <- midpoints[-length(midpoints)]
  segdiff$cumend <- midpoints[-1]
  segdiff$y1 <- seg.df$median[-nrow(seg.df)]
  segdiff$y2 <- seg.df$median[-1]
  colvec <- list("0" = "black", "1" = "blue", "-1" = "red")
  segdiff$col <- unlist(colvec[as.character(segdiff$sign)])
  segdiff$chr <- seg.df$chr[-nrow(seg.df)]
  
  if (posreg.flag) {
    pos.df$nseg <- vapply(1:nrow(pos.df), function(pid) { length(which(seg.df$cumstart <= pos.df$cumend[pid] & seg.df$cumend >= pos.df$cumstart[pid])) }, 1)
    pos.df$inv.ratio <- vapply(1:nrow(pos.df), function(pid) {
      posreg.idx <- which(seg.df$cumstart <= pos.df$cumend[pid] & seg.df$cumend >= pos.df$cumstart[pid])
      if (length(posreg.idx) <= 1) return(0)
      pos.absdiff <- abs(diff(segdiff$sign[posreg.idx]))
      return(length(which(pos.absdiff == 2)) / length(pos.absdiff))
    }, .1)
    
    ### Filters
    # 0 = valid
    # 2 = inversion rate too low
    # 3 = number of segments too low
    # 4 = minimum size too low
    # 5 = 2 + 3
    # 6 = 2 + 4
    # 7 = 3 + 4
    # 9 = 2 + 3 + 4
    
    pos.df$valid <- 0
    pos.df$valid[pos.df$inv.ratio < inv.rate] <- pos.df$valid[pos.df$inv.ratio < inv.rate] +2
    pos.df$valid[pos.df$nseg < nseg.min] <- pos.df$valid[pos.df$nseg < nseg.min] +3
    pos.df$valid[pos.df$width < width.min] <- pos.df$valid[pos.df$nseg < width.min] +4
  }
  
  blockcols <- list("0" = "cyan", "2" = "purple", "3" = "brown", "4" = "grey50", "5" = "tomato3", "6" = "wheat4", "7" = "mistyrose3", "9" = "red")
  
  # png(paste0(samp.name, ".png"), 1600, 950)
  # par(mfrow = c(2,1), oma = c(2,2,2,2), mai = c(.1,.1,.2,0))
  my.ylim.l2r <- c(-2,2)
  plot(coords.df$CumStart, l2r.mat[,samp.idx], pch = ".", cex = 3, type = "b", xaxs = "i", yaxs = "i", xaxt = "n", ylim = my.ylim.l2r, xlim = c(0, max(coords.df$CumStart)), main = paste0(samp.name, " PELT-segmented profile (nseg = ", nrow(seg.df), ", nCT = ", length(which(pos.df$valid == 0)), ")"), xlab = "genomic position", ylab = "L2R", col = "grey50")
  lines(coords.df$CumStart, runmed(l2r.mat[,samp.idx], 99), col = "black")
  if (posreg.flag) rect(pos.df$cumstart, my.ylim.l2r[1], pos.df$cumend, my.ylim.l2r[2], col = adjustcolor(unlist(blockcols[as.character(pos.df$valid)]), alpha.f = .35), border = unlist(blockcols[as.character(pos.df$valid)]))
  abline(h = 0, col = 3 )
  abline(v = kbrk, col = 3, lwd = 2)
  segments(seg.dfz[[samp.idx]]$cumstart, seg.dfz[[samp.idx]]$median, seg.dfz[[samp.idx]]$cumend, seg.dfz[[samp.idx]]$median, col = c("black", segdiff$col), lwd = 5)
  my.ylim.den <- range(c(expden.y.adj, denz.adj[[samp.idx]]$y.adj)) * 1.01
  plot(expden$x, expden.y.adj, lty = 2, lwd = 2, xaxs = "i", yaxs = "i", xlim = c(0, cs$genome.length), ylim = my.ylim.den, main = paste0("Local breakpoint density and HD regions (", posreg.nb, ")"), xlab = "genomic position", ylab = "Density", type = "n")
  abline(h = cutlim, col = 4, lty = 2, lwd = 2)
  abline(h = my.ylim.den[1] + (inv.rate * diff(my.ylim.den)), col = "gold", lty = 2, lwd = 2)
  if (posreg.flag) rect(pos.df$cumstart, my.ylim.den[1], pos.df$cumend, my.ylim.den[2] * pos.df$inv.ratio, col = adjustcolor("yellow", alpha.f = .35), border = NA)
  lines(denz.adj[[samp.idx]]$x, denz.adj[[samp.idx]]$y.adj, col = 1)
  abline(v = kbrk, col = 3, lwd = 2)
  # dev.off()
  
  ct.pos <- data.frame(SampleName = character(), chr = integer(), start = integer(), end = integer(), width = integer(), nseg = integer(), inv.ratio = numeric(), stringsAsFactors = FALSE)
  
  if (posreg.flag) {
    npos <- which(pos.df$valid == 0)
    
    if (length(npos) > 0) {
      ct.pos <- data.frame(SampleName = samp.name, chr = pos.df$chr, start = pos.df$start, end = pos.df$end, width = pos.df$width, nseg = pos.df$nseg, inv.ratio = pos.df$inv.ratio, stringsAsFactors = FALSE)
      ct.pos <- ct.pos[npos,]
      
      cbs.status.df <- ct.pos[,c(1:5)]
      colnames(cbs.status.df) <- c(samp.name, "Chr", "Start", "End", "Probes")
      cbs.status.df$Log2Ratio <- 1
      maelstrom::write.table.fast(x = cbs.status.df, file = paste0(samp.name, "_CTH", my.version, "_STATUS.cbs"))
      
      ## Generating the CBS file
      samp.pos.idx <- unlist(sapply(1:nrow(ct.pos), function(n) { in.idx <- which(seg.df$chr == ct.pos$chr[n] & seg.df$start >= ct.pos$start[n] & seg.df$end <= ct.pos$end[n]) }))
      seg.df$median[!((1:nrow(seg.df)) %in% samp.pos.idx)] <- 0
      
      cbs.df <- data.frame(samplename = samp.name, Chr = seg.df$chr, Start = seg.df$start, End = seg.df$end, Probes = seg.df$width, Log2Ratio = seg.df$median, stringsAsFactors = FALSE)
      colnames(cbs.df)[1] <- samp.name
      
      # cbs.df <- cbs.df[samp.pos.idx,]
      
      maelstrom::write.table.fast(x = cbs.df, file = paste0(samp.name, "_CTH", my.version, ".cbs"))
      seg.df$median[!((1:nrow(seg.df)) %in% samp.pos.idx)] <- 0
    
    }
  }
  return(ct.pos)
}
maelstrom::rasterpdf.close(pdf.optionslist = plop, width = 40/cm(1), height = 20/cm(1))

str(ct.res)
ct.res$start <- round(ct.res$start)
ct.res$end <- round(ct.res$end)
ct.res$width <- round(ct.res$width)
maelstrom::write.table.fast(ct.res, "CT_results.txt")




###################################""

plop2 <- maelstrom::rasterpdf.open(file = paste0("CThunter5v2_cf", cutfactor, ".pdf"), width = 2000, height = 1000)
par(mfrow = c(2,1), oma = c(2,2,2,2), mai = c(.1,.1,.2,0))
for (samp.idx in 1:ncol(l2r.mat)) {
  print(samp.idx)
  pos.idx <- as.numeric(denz.adj[[samp.idx]]$y.adj2 >= cutlim2)
  posreg.flag <- FALSE
  posreg.nb <- 0
  if (any(pos.idx == 1)) {
    posreg.flag <- TRUE
    pos.rle <- rle(pos.idx)
    pos.df <- data.frame(idx.start = c(0, cumsum(pos.rle$lengths[-length(pos.rle$lengths)]))+1, idx.end = cumsum(pos.rle$lengths), value = pos.rle$values, stringsAsFactors = FALSE)
    pos.df <- pos.df[pos.df$value == 1,]
    pos.df$cumstart <- denz.adj[[samp.idx]]$x[pos.df$idx.start]
    pos.df$cumend <- denz.adj[[samp.idx]]$x[pos.df$idx.end]
    posreg.nb <- nrow(pos.df)
  }
  
  seg.df <- seg.dfz[[samp.idx]]
  
  my.ylim.l2r <- c(-2,2)
  plot(coords.df$CumStart, l2r.mat[,samp.idx], pch = ".", cex = 3, type = "b", xaxs = "i", yaxs = "i", xaxt = "n", ylim = my.ylim.l2r, xlim = c(0, max(coords.df$CumStart)), main = paste0("PELT-segmented profile (", nrow(seg.df), ")"), xlab = "genomic position", ylab = "L2R", col = "grey50")
  lines(coords.df$CumStart, runmed(l2r.mat[,samp.idx], 99), col = "black")
  rect(pos.df$cumstart, my.ylim.l2r[1], pos.df$cumend, my.ylim.l2r[2], col = adjustcolor("yellow", alpha.f = .35), border = NA)
  abline(h = 0, col ="3")
  abline(v = kbrk, col = "3", lwd = 1)
  segments(seg.list[[samp.idx]]$cumstart, seg.list[[samp.idx]]$median, seg.list[[samp.idx]]$cumend, seg.list[[samp.idx]]$median, col = c("red", "blue"), lwd = 5)
  my.ylim.den <- range(c(expden.y.adj, denz.adj[[samp.idx]]$y.adj2)) * 1.01
  plot(expden$x, expden.y.adj, lty = 2, lwd = 2, xaxs = "i", yaxs = "i", xlim = c(0, cs$genome.length), ylim = my.ylim.den, main = paste0("Local breakpoint density and HD regions (", posreg.nb, ")"), xlab = "genomic position", ylab = "Density", type = "n")
  abline(h = cutlim2, col = 4, lty = 2, lwd = 2)
  rect(pos.df$cumstart, my.ylim.den[1], pos.df$cumend, my.ylim.den[2], col = adjustcolor("yellow", alpha.f = .35), border = NA)
  lines(denz.adj[[samp.idx]]$x, denz.adj[[samp.idx]]$y.adj2, col = 1)
  abline(v = kbrk, col = "3", lwd = 1)
}
maelstrom::rasterpdf.close(pdf.optionslist = plop2, width = 40/cm(1), height = 20/cm(1))









####




  seg.df <- data.frame(end.idx = sort(unique(c(changepoint::cpt.meanvar(data = l2r.mat[,s] , penalty = "MBIC", method = "PELT", param.estimates = FALSE, minseglen = 3)@cpts, chrend.idx))), stringsAsFactors = FALSE)
  seg.df$start.idx <- c(1, seg.df$end.idx[-nrow(seg.df)] +1)
  seg.df$sd <- vapply(1:nrow(seg.df), function(x){ sd(l2r.mat[seg.df$start.idx[x]:seg.df$end.idx[x],s], na.rm = TRUE)}, .1)
  seg.df$median <- vapply(1:nrow(seg.df), function(x){ median(l2r.mat[seg.df$start.idx[x]:seg.df$end.idx[x],s], na.rm = TRUE)}, .1)
  # seg.df$uSD <- c(vapply(2:nrow(seg.df), function(x) {
  #   if (any(is.na(seg.df$sd[(x-1):x]))) return(0)
  #   if (abs(diff(seg.df$median[(x-1):x])) < min(seg.df$sd[(x-1):x])) return(1) else return(0)
  # }, 1), 0)
  # seg.df <- seg.df[seg.df$uSD == 0,]
  seg.df$start.idx <- c(1, seg.df$end.idx[-nrow(seg.df)] +1)
  seg.df$sd <- vapply(1:nrow(seg.df), function(x){ sd(l2r.mat[seg.df$start.idx[x]:seg.df$end.idx[x],s], na.rm = TRUE)}, .1)
  seg.df$median <- vapply(1:nrow(seg.df), function(x){ median(l2r.mat[seg.df$start.idx[x]:seg.df$end.idx[x],s], na.rm = TRUE)}, .1)
  seg.df$chr <- coords.df$Chr[seg.df$start.idx]
  seg.df$start <- coords.df$Start[seg.df$start.idx]
  seg.df$end <- coords.df$Start[seg.df$end.idx]
  seg.df$cumstart <- coords.df$CumStart[seg.df$start.idx]
  seg.df$cumend <- coords.df$CumStart[seg.df$end.idx]

  return(seg.df)
}, simplify = FALSE)

breakidx.list <- sapply(1:length(seg.list), function(s) {
  return(seg.list[[s]]$end.idx)
})

breakpop.list <- sapply(1:length(seg.list), function(s) {
  return(coords.df$CumStart[seg.list[[s]]$end.idx])
})


expstep <- cs$genome.length / length(unlist(breakidx.list))
expvec <- seq(1, cs$genome.length, expstep)
expavg <- length(expvec) / ncol(l2r.mat)
expden <- density(expvec, adjust = my.adjust, width = my.width, n = n, from = 0, to = cs$genome.length, cut = 3)
expLOCden <- density(unlist(breakpop.list), adjust = my.adjust, width = my.width, n = n, from = 0, to = cs$genome.length, cut = 3)

expden.y.adj <- expden$y*expavg*cutfactor
cutlim <- max(expden.y.adj)
expLOCden.y.adj <- expLOCden$y*expavg

# require(maelstrom)
plop <- maelstrom::rasterpdf.open(file = "CThunterZ.pdf", width = 2000, height = 1000)
# png(file = "test.png", width = 1500, height = 1000)
par(mfrow = c(2,1), oma = c(2,2,2,2), mai = c(.1,.1,.2,0))
for (samp.idx in 1:ncol(l2r.mat)) {
  print(samp.idx)
  segdenDEF <- density(coords.df$CumStart[seg.list[[samp.idx]]$end.idx], adjust = my.adjust, width = my.width, n = n, from = 0, to = cs$genome.length, cut = 3)
  segdenDEF.y.adj <- segdenDEF$y*length(seg.list[[samp.idx]]$end.idx)
  
  pos.idx <- as.numeric(segdenDEF.y.adj >= cutlim)
  if (any(pos.idx) == 1) {
    pos.rle <- rle(pos.idx)
    pos.df <- data.frame(idx.start = c(0, cumsum(pos.rle$lengths[-length(pos.rle$lengths)]))+1, idx.end = cumsum(pos.rle$lengths), value = pos.rle$values, stringsAsFactors = FALSE)
    pos.df <- pos.df[pos.df$value == 1,]
    pos.df$cumstart <- segdenDEF$x[pos.df$idx.start]
    pos.df$cumend <- segdenDEF$x[pos.df$idx.end]
  }
  
  seg.df <- seg.list[[samp.idx]]
  
  my.ylim.l2r <- c(-2,2)
  plot(coords.df$CumStart, l2r.mat[,samp.idx], pch = ".", cex = 3, type = "b", xaxs = "i", yaxs = "i", xaxt = "n", ylim = my.ylim.l2r, xlim = c(0, cs$genome.length), main = "PELT-segmented profile (low-variance-filtered)", xlab = "genomic position", ylab = "L2R", col = "grey50")
  lines(coords.df$CumStart, runmed(l2r.mat[,samp.idx], 99), col = "black")
  rect(pos.df$cumstart, my.ylim.l2r[1], pos.df$cumend, my.ylim.l2r[2], col = adjustcolor("yellow", alpha.f = .5), border = NA)
  abline(h = 0, col ="3")
  abline(v = kbrk, col = "3", lwd = 2)
  segments(seg.list[[samp.idx]]$cumstart, seg.list[[samp.idx]]$median, seg.list[[samp.idx]]$cumend, seg.list[[samp.idx]]$median, col = 2, lwd = 5)
  # abline(v = breakpop.list[[samp.idx]], col = 2)
  # plot(segdenDEF, xaxs = "i", xlim = c(0, cs$genome.length), main = "Local breakpoint density", xlab = "genomic position", ylab = "Density")
  # lines(expden$x, expden$y*3, col = 2, lty = 2, lwd = 2)
  # lines(expLOCden, col = 2)
  my.ylim.den <- range(c(expden.y.adj, segdenDEF.y.adj, expLOCden.y.adj)) * 1.01
  plot(expden$x, expden.y.adj, lty = 2, lwd = 2, xaxs = "i", yaxs = "i", xlim = c(0, cs$genome.length), ylim = my.ylim.den, main = "Local breakpoint density", xlab = "genomic position", ylab = "Density", type = "n")
  abline(h = cutlim, col = 4, lty = 2, lwd = 2)
  rect(pos.df$cumstart, my.ylim.den[1], pos.df$cumend, my.ylim.den[2], col = adjustcolor("yellow", alpha.f = .5), border = NA)
  lines(segdenDEF$x, segdenDEF.y.adj, col = 1)
  lines(expLOCden$x, expLOCden.y.adj, col = "grey70")
  abline(v = kbrk, col = "3", lwd = 2)
}
maelstrom::rasterpdf.close(pdf.optionslist = plop, width = 40/cm(1), height = 20/cm(1))





## Refining breakpoints based on variance
seg.df$start.idx <- c(1, seg.df$end.idx[-nrow(seg.df)] +1)
seg.df$sd <- vapply(1:nrow(seg.df), function(x){ sd(samp.df$L2R[seg.df$start.idx[x]:seg.df$end.idx[x]], na.rm = TRUE)}, .1)
seg.df$median <- vapply(1:nrow(seg.df), function(x){ median(samp.df$L2R[seg.df$start.idx[x]:seg.df$end.idx[x]], na.rm = TRUE)}, .1)
seg.df$uSD <- c(vapply(2:nrow(seg.df), function(x) {
  if (any(is.na(seg.df$sd[(x-1):x]))) return(0)
  if (abs(diff(seg.df$median[(x-1):x])) < min(seg.df$sd[(x-1):x])) return(1) else return(0)
}, 1), 0)

seg.df <- seg.df[seg.df$uSD == 0,]

segdenDEF <- density(samp.df$CumStart[seg.df$end.idx], adjust = .05)
expstep <- cs$genome.length / length(seg.df$end.idx)
expvec <- seq(1, cs$genome.length, expstep)
expden <- density(expvec, adjust = .1)

par(mfrow = c(2,1), oma = c(2,2,2,2))
plot(samp.df$CumStart, samp.df$L2R, pch = ".", cex = 3, type = "b", xaxs = "i", ylim = c(-1.5,1.5), xlim = c(0, cs$genome.length))
lines(samp.df$CumStart, runmed(samp.df$L2R, 99), col = "cyan")
abline(v = samp.df$CumStart[seg.df$end.idx], col = 2)
plot(segdenDEF, xaxs = "i", xlim = c(0, cs$genome.length))
lines(expden$x, expden$y*3, col = 2, lty = 2, lwd = 2)




#############################################################################


setwd("/home/job/svn/genomics/CGH/R/CThunter")

data("hg19", package = "chromosomes", envir = environment())
## Reseg using PELT
cgh.data <- readRDS("/mnt/data_cigogne/job/cna_joint_admixed_analysis/04_most_centered_peak_centering/multipcf_ADMIX3_sub_GC100000_l2r_124s/copynumber_segmentation_centered_mpc_mcpc.RDS")
chrend <- vapply(unique(cgh.data$cn.data.ini$Chr), function(k) { which.max(cgh.data$cn.data.ini$Chr == k) }, 1)
chrlencum <- c(0,cumsum(vapply(unique(cgh.data$cn.data.ini$Chr), function(kr) { max(cgh.data$cn.data.ini$Start[cgh.data$cn.data.ini$Chr == kr])}, 1)))

samp.idx <- 1


samp.df <- data.frame(Chr = cgh.data$cn.data.ini$Chr, Start = cgh.data$cn.data.ini$Start, CumStart = cgh.data$cn.data.ini$Start + chrlencum[cgh.data$cn.data.ini$Chr], L2R = cgh.data$cn.data.wins[,samp.idx+2] - cgh.data$most.centered.peak[samp.idx], stringsAsFactors = FALSE)

# seg.L2R <- sort(unique(c(changepoint::cpt.mean(data = samp.df$L2R, penalty = "AIC0", method = "PELT", param.estimates = FALSE, minseglen = 3)@cpts, chrend)))
seg.L2R <- sort(unique(c(changepoint::cpt.meanvar(data = samp.df$L2R, penalty = "MBIC", method = "PELT", param.estimates = FALSE, minseglen = 3)@cpts, chrend)))
# seg.L2R <- sort(unique(c(changepoint::cpt.var(data = samp.df$L2R, penalty = "BIC0", method = "PELT", param.estimates = FALSE, minseglen = 3)@cpts, chrend)))

length(seg.L2R)

seg.df <- data.frame(Idx.start = c(1, seg.L2R[-length(seg.L2R)]+1), Idx.end = seg.L2R, stringsAsFactors = FALSE)
seg.df$Start <- samp.df$Start[seg.df$Idx.start]
seg.df$End <-  samp.df$Start[seg.df$Idx.end]
seg.df$CumStart <- samp.df$CumStart[seg.df$Idx.start]
seg.df$CumEnd <- samp.df$CumStart[seg.df$Idx.end]
seg.df$Width <- seg.df$End - seg.df$Start + 1
seg.df$Nprobes <- diff(c(0, seg.L2R))
seg.df$L2R <- vapply(1:nrow(seg.df), function(s) { return(median(samp.df$L2R[seg.df$Idx.start[s]:seg.df$Idx.end[s]]))}, .1)

## re-binning for variance
# probeden <- round(max(samp.df$CumStart) / nrow(samp.df))
# probeden <- median(seg.df$Width, na.rm = TRUE)
# probeden <- 5E+06
probeden <- 1E+06
require(foreach)
segbin.df <- foreach(s = 1:nrow(seg.df), .combine = "rbind") %do% {
  if(seg.df$Width[s] <= probeden) return(seg.df[s,]) else {
    # return(rep(seg.df[s,], ceiling(seg.df$Width[s]/probeden)))
    temp.df <- seg.df[rep(s, ceiling(seg.df$Width[s]/probeden)),]
    temp.df$Start <- c(temp.df$Start[1], vapply(2:nrow(temp.df), function(x) { return(temp.df$Start[x-1]+(probeden*(x-1))) }, 1))
    temp.df$End <- c(vapply(1:(nrow(temp.df)-1), function(x) { return(temp.df$Start[x]+probeden-1) }, 1), temp.df$End[nrow(temp.df)])
    temp.df$CumStart <- c(temp.df$CumStart[1], vapply(2:nrow(temp.df), function(x) { return(temp.df$CumStart[x-1]+(probeden*(x-1))) }, 1))
    temp.df$CumEnd <- c(vapply(1:(nrow(temp.df)-1), function(x) { return(temp.df$CumStart[x]+probeden-1) }, 1), temp.df$CumEnd[nrow(temp.df)])
    temp.df$Width <- temp.df$End - temp.df$Start +1
    # return(seg.df[rep(s, ceiling(seg.df$Width[s]/probeden)),])
    return(temp.df)
  }
}

## Computing the DIFF values
# segbin.df$Diff <- c(0, abs(diff(segbin.df$L2R)))
segbin.df$Diff <- c( abs(diff(segbin.df$L2R)), 0)
# segbin.df$Diff <- c(0, smooth(abs(diff(segbin.df$L2R)), twiceit = TRUE))

## Segmenting DIFF
# p2.res <- changepoint::cpt.var(data = segbin.df$Diff, penalty = "BIC0", method = "PELT", param.estimates = FALSE)@cpts
p2.res <- changepoint::cpt.var(data = segbin.df$Diff, penalty = "MBIC", method = "PELT", param.estimates = FALSE)@cpts
str(p2.res)

par(mfrow = c(2,1), oma = c(2,2,2,2))
# plot(segbin.df$Diff, pch = ".", cex = 5, type = "b", xaxs = "i", ylim = c(0,.5))
plot(samp.df$CumStart, samp.df$L2R, pch = ".", cex = 3, type = "b", xaxs = "i", ylim = c(-1.5,1.5), xlim = c(0, cs$genome.length))
lines(samp.df$CumStart, runmed(samp.df$L2R, 99), col = "cyan")
abline(v = samp.df$CumStart[seg.L2R], col = 2)

smorletest <- as.numeric(smooth(segbin.df$Diff))
rletest <- rle(smorletest)
rle.df <- data.frame(End = cumsum(rletest$lengths), Length = rletest$lengths, Value = rletest$values, stringsAsFactors = FALSE)
rle.df$Start <- c(1, rle.df$End[-nrow(rle.df)]+1)
for(s in 2:(nrow(rle.df)-1)) {
  # if(rle.df$Value[s-1] == 0 & rle.df$Value[s] != 0 & rle.df$Value[s+1] == 0 & rle.df$Length[s] < 3) smorletest[rle.df$Start[s]:rle.df$End[s]] <- 0
  if(rle.df$Value[s-1] == 0 & rle.df$Value[s] != 0 & rle.df$Value[s+1] == 0) smorletest[rle.df$Start[s]:rle.df$End[s]] <- 0
}
# plot(smorletest, pch = ".", cex = 5, type = "b", xaxs = "i", main = probeden)
# # abline(v = changepoint::cpt.var(data = smorletest, penalty = "MBIC", method = "PELT", param.estimates = FALSE)@cpts, col = 2)
# abline(v = changepoint::cpt.var(data = smorletest, penalty = "BIC", method = "PELT", param.estimates = FALSE)@cpts, col = 2)
# par(mfrow = c(1,1))

plot(0,0, type = "n", xaxs = "i", main = probeden, xlim = c(0, cs$genome.length), ylim = c(0,max(smorletest)), xlab = "Genomic position")
segments(segbin.df$CumStart, smorletest, segbin.df$CumEnd, smorletest, col = 1, lwd = 4, main = probeden)
# abline(v = changepoint::cpt.var(data = smorletest, penalty = "MBIC", method = "PELT", param.estimates = FALSE)@cpts, col = 2)
abline(v = segbin.df$CumEnd[changepoint::cpt.var(data = smorletest, penalty = "BIC", method = "PELT", param.estimates = FALSE)@cpts], col = 2)
par(mfrow = c(1,1))








##############################################"

plot(segbin.df$Diff, pch = ".", cex = 5, type = "b", xaxs = "i")
# points(segbin.df$Diff, pch = ".", cex = 5, type = "b", col = 3)
abline(v = p2.res, col = 2)



plot(abs(segbin.df$Diff), pch = ".", cex = 5, type = "b", xaxs = "i")
points(smooth(abs(segbin.df$Diff)), pch = ".", cex = 5, type = "b", col = 3)
abline(v = p2.res, col = 2)



plot(samp.df$L2R, pch = ".", xaxs = "i", col = "grey50")
abline(v = seg.L2R, col = 2)

lines(runmed(L2R, 199), col = 4, lwd = 2)
abline(h = 0, col = 3)




nseg.mini <- 5

seg.func <- "var"
seg.func <- "meanvar"

penalty.method <- "BIC"
penalty.method <- "MBIC"
penalty.method <- "AIC"
penalty.method <- "Hannan-Quinn"

data.type <- "DL2R"
data.type <- "DRM"
data.type <- "DSMOT"

metric.type <- "VAR"
metric.type <- "MAPD"
metric.type <- "SAPD"
metric.type <- "IQR"

require(foreach)
foreach(seg.func = c("meanvar", "var")) %do% {
  foreach(data.type = c("DL2R")) %do% {
    # foreach(metric.type = c("IQR", "VAR", "MAPD")) %do% {
    foreach(metric.type = c("VAR")) %do% {
      foreach(penalty.method = c("AIC0", "BIC0", "Hannan-Quinn0", "MBIC")) %do% {
        # foreach(penalty.method = c("Hannan-Quinn0")) %do% {
        pdf(paste0("CTscore_N", nseg.mini, "_", seg.func, "_", penalty.method, "_", data.type, "_", metric.type, ".pdf"), width = 35/cm(1), height = 21/cm(1))
        for (k in cbsfiles) {
          kname <- basename(k)
          # print(kname)
          cbs.df <- maelstrom::read.table.fast(k)
          cbs.df <- cbs.df[cbs.df$Chr != 24,]
          chrend <- vapply(unique(cbs.df$Chr), function(kr) { max(which(cbs.df$Chr == kr))}, 1)
          chrlencum <- c(0,cumsum(vapply(unique(cbs.df$Chr), function(kr) { max(cbs.df$End[cbs.df$Chr == kr])}, 1)))
          cbs.df$CumStart <- cbs.df$Start + chrlencum[cbs.df$Chr]
          cbs.df$CumEnd <- cbs.df$End + chrlencum[cbs.df$Chr]
          
          myrm <- runmed(cbs.df$Log2Ratio, 3)
          smooth.twice <- smooth(cbs.df$Log2Ratio, twiceit = TRUE)
          
          if (data.type == "DL2R") cbs.df$diff <- c(0, diff(cbs.df$Log2Ratio))
          if (data.type == "DRM") cbs.df$diff <- c(0, diff(myrm))
          if (data.type == "DSMOT") cbs.df$diff <- c(0, diff(smooth.twice))
          
          if (seg.func == "var") dl2r.seg <- changepoint::cpt.var(data = cbs.df$diff, penalty = penalty.method, method = "PELT", param.estimates = FALSE, minseglen = nseg.mini)@cpts
          if (seg.func == "meanvar") dl2r.seg <- changepoint::cpt.meanvar(data = cbs.df$diff, penalty = penalty.method, method = "PELT", param.estimates = FALSE, minseglen = nseg.mini)@cpts
          
          dl2r.seg <- sort(unique(c(dl2r.seg, chrend)))
          seg.df <- data.frame(Idx.start = c(1, dl2r.seg[-length(dl2r.seg)]+1), Idx.end = dl2r.seg, stringsAsFactors = FALSE)
          seg.df$Start <- cbs.df$Start[seg.df$Idx.start]
          seg.df$End <-  cbs.df$End[seg.df$Idx.end]
          seg.df$CumStart <- cbs.df$CumStart[seg.df$Idx.start]
          seg.df$CumEnd <- cbs.df$CumEnd[seg.df$Idx.end]
          seg.df$Width <- seg.df$End - seg.df$Start + 1
          seg.df$Nseg <- c(dl2r.seg[1], diff(dl2r.seg))
          
          seg.df$Var <- vapply(1:nrow(seg.df), function(s) {
            in.idx <- (seg.df$Idx.start[s]+1):seg.df$Idx.end[s]
            if (length(in.idx) < nseg.mini) return(NA)
            # if (metric.type == "VAR") return(var(cbs.df$diff[in.idx]))
            if (metric.type == "VAR") return(var(abs(cbs.df$diff[in.idx])))
            # if (metric.type == "VAR") return(var(cbs.df$Log2Ratio[in.idx]))
            if (metric.type == "MAPD") return(median(abs(diff(cbs.df$diff[in.idx]))))
            # if (metric.type == "MAPD") return(median(abs(diff(cbs.df$Log2Ratio[in.idx]))))
            # if (metric.type == "SAPD") return(log10(sum(abs(diff(cbs.df$diff[in.idx])))))
            if (metric.type == "IQR") return(IQR(cbs.df$diff[in.idx]))
            # if (metric.type == "IQR") return(IQR(cbs.df$Log2Ratio[in.idx]))
            # if (metric.type == "KURT") return(modes::kurtosis(abs(cbs.df$diff[in.idx]), finite = TRUE))
          }, .1)
          
          ## Plot
          # plot(cbs.df$Log2Ratio, ylim = range(cbs.df$Log2Ratio), xlim = c(0,max(chrlencum)), type = "n", xaxs = "i", main = kname, xlab = "Genomic position", ylab = "L2R // CTscore")
          # abline(v = chrlencum, col = "grey50", lty = 3)
          # rect(seg.df$CumStart, 0, seg.df$CumEnd, seg.df$Var, col = 3, border = NA)
          # segments(cbs.df$CumEnd[-nrow(cbs.df)], cbs.df$Log2Ratio[-nrow(cbs.df)], cbs.df$CumStart[-1], cbs.df$Log2Ratio[-1], col = 1)
          # segments(cbs.df$CumStart, cbs.df$Log2Ratio, cbs.df$CumEnd, cbs.df$Log2Ratio, lwd = 3, col = 2)
          
          ## Plot (2)
          # oripar <- par()
          par(mfrow = c(1,1))
          plot(cbs.df$Log2Ratio, xlim = c(0,max(chrlencum)), ylim = c(-1.5,1.5), yaxs = "i", type = "n", xaxs = "i", main = kname, xlab = "Genomic position", ylab = "L2R // CTscore")
          rect(seg.df$CumStart, -1.5, seg.df$CumEnd, 1.5, col = c("pink", "lightblue")[((1:nrow(seg.df))%%2)+1], border = NA)
          abline(v = chrlencum, col = 1, lty = 2, lwd = 2)
          abline(h = 0, lty = 3, lwd = 2)
          rect(seg.df$CumStart, -1.5, seg.df$CumEnd, seg.df$Var-1.5, col = 4, border = 4, lwd = 3)
          segments(cbs.df$CumEnd[-nrow(cbs.df)], cbs.df$Log2Ratio[-nrow(cbs.df)], cbs.df$CumStart[-1], cbs.df$Log2Ratio[-1], col = 1)
          segments(cbs.df$CumStart, cbs.df$Log2Ratio, cbs.df$CumEnd, cbs.df$Log2Ratio, lwd = 3, col = 2)
          
          # par(mfrow = c(8,5), oma = c(1,1,1,1))
          # zzz <- vapply(1:nrow(seg.df), function(s) {
          #   in.idx <- (seg.df$Idx.start[s]+1):seg.df$Idx.end[s]
          #   if (length(in.idx) > nseg.mini) plot(density(cbs.df$Log2Ratio[in.idx]), main = s)
          #   return(1)
          # }, 1)
          # par(oripar)
        }
        dev.off()
      }
    }
  }
}
