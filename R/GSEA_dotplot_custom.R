workdir <- '/home/job/WORKSPACE/B20038_PASA_01/ANALYSIS/PANEL_ONCO/dotplot_custom_subtypes'

setwd(workdir)

## Collecting filenames
infiles <- list.files(path = workdir, pattern = 'CANCER.HALLMARKS_GSEA.res.RDS', recursive = TRUE, full.names = TRUE)

## Getting types
my.types <- unname(vapply(infiles, function(x) unlist(strsplit(x = x, split = '/'))[11], 'a'))
## Loading GSEA results
fgsea.res.list <- lapply(infiles, readRDS)
names(fgsea.res.list) <- my.types
## Merging into a single object
fgsea.res.merged <- clusterProfiler::merge_result(enrichResultList = fgsea.res.list)

## Corrected generatio
dot_df <- fgsea.res.merged@compareClusterResult
dot_df$gene_count <- stringr::str_count(dot_df$core_enrichment, "/")
dot_df$GeneRatio <- dot_df$gene_count / dot_df$setSize

## Adding sign column (to split plot)
dot_df$sign = "activated"
dot_df$sign[dot_df$NES < 0] = "suppressed"

## Reordering the 

## Reshaping term names (to split longest ones)
dot_df$Description <- gsub(pattern = '_', ' ', dot_df$ID, fixed = TRUE)
## Reshaping clusers
levels(dot_df$Cluster) <- sub(pattern = '_vs_Other$', replacement = '', x = levels(dot_df$Cluster))
## Building plot
desc.tbl <- table(dot_df$Description)
desc.occ <- unname(desc.tbl[dot_df$Description])
# p <- ggplot2::ggplot(dot_df, ggplot2::aes(x = GeneRatio, y = Description)) + 
p <- ggplot2::ggplot(dot_df, ggplot2::aes(x = GeneRatio, y = forcats::fct_reorder(Description, desc.occ, .desc = FALSE))) + 
  ggplot2::geom_point(ggplot2::aes(size = GeneRatio, color = p.adjust)) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ggplot2::ylab(NULL) +
  ggplot2::ggtitle("Subtype GSEA (MSigDB CANCER HALLMARKS) terms enrichment")

library(forcats)
p2 <- p + ggplot2::facet_grid(. ~ Cluster + sign) + ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=40)) + ggplot2::scale_x_continuous(breaks = seq(0, 1, .5), limits = c(0,1))

png(filename = paste0(workdir, '/dot_test.png'), width = 2000, height = 1000)
print(p2)
dev.off()

