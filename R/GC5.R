## GC5.R
##
## DESCRIPTION : script performing GC% debiaising / segmentation oligo Agilent CGH/CHIP profiles (v5),
##
## AUTHOR:
## Bastien JOB (bastien.job@igr.fr)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 2.12.0 to 3.0.2
##
## DEPENDS ON:
## package :  DNAcopy (for the CBS segmentation)
## package :  limma (for the loessfits)
## package :  snowfall (for multithreading) *ONLY IF ncores > 1 !*
## script :   gc5psconv (to convert Postscript output to multipage TIF)
## program :  bzip2 (to convert GCX to GCX.BZ2)
## file :     (pelican) /proj/cgh/Agilent/design_list.txt (TSV with informations about design versions)
## file :     *.gc (table with precomputed GC% tracks per species+techno+genomebuild+AMADID+designdate)
## file :     *.cngs (table with precomputed squared quantity of CnG triplets tracks per species+techno+genomebuild+AMADID+designdate)
##
## VERSION NOTES
##
## v5.15b 20160830
##      . Replaced the mnt.proj and mnt.db options by ldb (local data base), to get rid of pelican.
##	. Introduced a new parameter : alpha (see CBS help).
##
## v5.15  20140505
##      . [RCMT] build.target return name of target file.
##
## v5.14  20140311
##      . Exported chrload() to its own script, which would be sourced by scripts requiring it, such as this one.
##      . Added a message at sourcing of GC5.R explaining that chrload() must be sourced prior gc5fit() run.
##      . Consequently adaptated GC5.R to changes made to chrload() : gv, chrom2chr, chr2chrom and species full
##        name are now included into the chrload() output object.
##      . Added automatic copying of the targets file given as input to the results directory.
##      . Changed the first centralization routine for CGH : now center of the distribution is computed thanks to
##        the ratio of left/right areas under curve for each peak instead of its position to abscisse edges.
##        Consequently, the definition of the nrf parameter changed to a fraction of the population (area) of
##        the most populated peak, instead of the fraction of the height of the highest peak. Default value
##        for nrf is currently unchanged but needs further testing to assess its stability.
##      . Added by default bzip2 compression of wiggles (halves their size).
##
## v5.13b 20140303
##      . Added "gv" to the chrload() output.
##      . Modified handling of species for the "species long name" required to access the GC and CNGS data files.
##
## v5.13 20140227
##      . Slight modification to the chrload() to adapt to the recent UCSC stupid habit to include incomplete
##        lines corresponding to chrN_unknown sequence sets (in cytoBandIdeo).
##      . Also had to remove case of chrM which has appeared in the same file. Now the script checks
##        and only keeps chromosomes that are present in the list chrconvlist$chrom2chr for the requested
##        species.
##      . Added proliminary (untested) support for rattus norvegicus (limited to rn4).
##
## v5.12 20130318
##      . Corrected bug : now, proper support of gcfit=F and cngfit=F is working (formerly, gcfit=F
##        deactivated both, and cngfit did not do anything). Early tests seem to show that using cngfit
##        only and deactivating gcfit performs a little better normalization (on ChIP/CH3 experiments).
##
## v5.11b 20121121
##      . Corrected a bug making the generation of the global PQC impossible.
##
## v5.11 20121001
##      . Rewrote the chrX/chrY tweak required after a dyefit, which appeared to have been missing since
##        many versions !
##
## v5.10b 20120821
##      . Minor correction : now the script returns the true value of arrays to process when some are set to
##        FALSE for the "Perform" column in the targets file. Precedently, it was giving the total amount of
##        arrays in the targets file. Strikingly, this bug seemed to prevent snowfall to effectively work
##        with all requested cores. Seems to be fixed now.
##      . Accordingly, better set the global declaration of 'tlist'.
##
## v5.10 20120416
##      . Added the possibility to use dyes from different slides (different FEX files), for an "in silico
##        hybridization" scheme.
##      . Accordingly modified the targetslist structure to handle the existance of dye track and source FEX
##        file for each channel independantly. Thus, by default, the build.targets() function will create a
##        targetlist in which all found FEX files will be present once, and for which the script will use both
##        channels from this FEX file. To perform such "in silico hybridization", the user should manually
##        modify the targetslist, writing the new source FEX file(s).
##      . Consequently, removed the "inv" option, not needed anymore.
##      . Added the amadid and design date in the PQC output files, for tracking. This was mainly lacking
##        because the default (and quite always used) value of "latest" did not correspond to any real date...
##      . Added a "mode" argument to the build.targets() function, which will modify the type of columns that
##        will be displayed in the generated targetlist. Consequently, the "mode" parameter won't be available
##        anymore in the calling of gc5fit(), as read in the targetlist. This allows an easier handling of
##        parameters which are CGH- or CHIP- specific. This also allows to process arrays for different modes in
##        a same batch (even if this will probably never happen!).
##      . To ease the mode separation, added a new gc5.run intermediary function.
##      . Fixed a bug that prevented the right behaviour of the "process" option.
##
## v5.9.2c 20120209
##      . Corrected  bug which prevented the generating of the global PQC table when only one slide processing
##        was requested (again!).
##      . Added a control to check the adequation of ncores value and the number of processes requested. If
##        this number is lower than the nomber of cores requested, ncores is lowered to the number of
##        processes.
##
## v5.9.2b 20120206
##      . Changed the build.target() function structure (lapply to sapply, changed data types, ...) which
##        prevented the function to work when asked for very few fex files (2 or less). Function input and
##        output did not change.
##
## v5.9.2 20120120
##      . Added support for a new normalization scheme, based on CnG counts in the sequence (C-whatever-G,
##        and the inverse G-whatever-C), taken from the *.cngs generated by the 'CnGbed' script, for each
##        Agilent BED file, using the 's' value for the 'smod' option (meaning that repetitions of the triplet
##        have a squared weight in the count). These precomputed tracks add to those already available for the
##        direct GC%, and do not replace these. It increases the computation time but adds visibly more
##        quality on the final profile.
##      . Added a stop for cores use, when ncores>1, to free some memory when the user works with other scripts
##        in the same session after using this GC5.R script.
##
## v5.9.1b 20111213
##      . Added a conditional loading of the snowfall package, only when ncores >1. This is to restore
##        support for computers with a single core and/or which do not have the snowfall package installed.
##
## v5.9.1 20111208
##      . Modified the structure of the "cut" version of CBS files : now normal segments are present, with a
##        log2(ratio) value of 0 (they formerly were removed).
##      . Removed the generation of the global PQC table when only one sample has to be processed in the batch.
##        This is mainly because it was bugged in this case (no values for this single sample), and it's easier
##        to correct it this way (=.
##
## v5.9.0 20111202
##      . Changed the structure of precomputed GC% data. Now, these GC% tracks are available in one file per
##        species+techno+genomebuild+AMADID+designdate, whereas this source was a single, big file in which were
##        concatenated all probes from different AMADID for a same species+techno+genomebuild (ie : hs, CGH, hg19).
##        This, was required for an easier and more convenient update of arrays annotations : anytime Agilent releases
##        a new annotation file (bed + TDT downloaded from eArray) for any design, the 'GCbed' perl script will be
##        used on the new BED file(s) without touching tracks for non-updated designs.
##      . As GC5 won't deal with a single GC% table, it has to know where to find which GC table. In this purpose
##        the script now needs a small TSV table in which it can find the required informations : Technology,
##        Format, AMADID, Species, Genome.build, and Design.date. This table named 'designs.list' is stored
##        in /proj/cgh/Agilent. This file is also required for the 'tabgen' perl script. This sctructure will also
##        allow to choose for a same species+techno+genomebuild+AMADID between different design dates, when
##        available. To get the right GC table for the right format, the input FEX files are required to have the
##        AMADID at the beginning of their filename (new name format).
##      . Consequently, the gc5load() function has been pulled out from the code, as deprecated.
##      . To simplify the re-processing of FEX files, the main function will now require a 'targets' file, with one
##        line per array to process, and several columns for all of the parameters (plus some additional info
##        like techno, species, build, amadid, valid, ...). A new function 'build.targets()' has been created to
##        generate a default, pre-filled targets file.
##      . GC5 is now multithreaded, thanks to the 'snowfall' package ! Fur this purpose, an option 'ncore' has been
##        added to the gc5fit() function (default=2). So, it can basically process several arrays at the same time.
##      . Added an output of the normalized log2(ratio) profile as a WIG, to export the profile to IGB/IGV/UCSC(...).
##      . Added an output of the segmentation results as BEDGRAPH files (Cut & NoCut), to export to IGB/IGV/UCSC(...).
##      . Added an output of the final segmented profile in pangenomic and karyotypic views, in two separated PNGs.
##      . Added automatic compression of GCX files to BZ2, if the "bzip2" software is available on the host system.
##      . Added chromosomal views in the graphical report.
##      . Added colored marks on the horizontal final views (pangenomic, and chromosomal) to show where segments are
##        present but not displayed due to their value out of the fixed Y limits.
##      . "ymax" parameter can now be configured (not an option, it has to be modified within the script). Default
##        value is still 1.5.
##
## v5.8c 20111118
##      . Corrected a bug in the karyotype view which still showed the segment mean values instead of the now used
##        median values, making this view incoherent with the linear genomic profile (especially for chr19).
##      . Deleted some deprecated commented code lines corresponding to segment means computations and nucleotides to
##        kilobases conversions.
##
## v5.8b 20111116
##      . Introduced a "help.gc5fit" function to display the description of parameters for the main function, even when
##        this script is sourced and not directly read.
##      . Renamed the main "GC5fit()" function to "gc5fit()".
##
## v5.8 20111104
##      . Modified the CBS results to have median value of probes per segment, instead of default,
##        forced mean value, which causes problems when "smears" appear in a chr or whole profile.
##      . Segments' ends are now retrieved in a more robust way, using gcdata info by invoking indexes.
##      . Coordinates are not converted to Kb anymore for segmentation, it's unneeded in current-gen
##        CBS implementation.
##      . Added the "tif" (or "tiff", both equally accepted) value to the "of" option for the graphical
##        output format. With this value, the script will write a .ps, then invoke the gc5psconv bash
##        script, then delete the .ps file.
##
## v5.7 20110819
##      . Rewrote the script in a more modular form.
##      . Support for chip / ch3 arrays is present, in its current form (optional quantiles
##        normalization, different centering method, no segmentation).
##      . Preliminary support for other species : modification of the structure. mm (mus musculs)
##        is currently in test. Update : mm fully supported now.
##      . Reverted mod classes : now only "cgh" & "chip" subsist : "snp" is merged with "cgh", as their
##        GC% sources were also merged. "ch3" and "chip" were merged as well, as they share the same
##        array design for the current generation, and while the specific computations of 
##        CpG composition were not as efficient as GC%, neither added some value to the results.
##      . New GC% tracks were pre-computed for lower sizes. So, now these GC tracks are available :
##        [probe, 100b, 250b, 500b, 750b, 1k, 2.5k, 5k, 10k, 50k, 100k] for "cgh", and [probe, 100b,
##        200b, 300b, 400b, 500b, 600b, 1k, 5k, 10k, 50k, 100k] for "chip". These tracks are available
##        for hg19 and mm9 ONLY, support for hg18 will be computed later.
##      . Modified output filenames, way more clear now.
##      . A qq-norm plot has been added to the final segmented profile for CGH, as it is an
##        interesting graphical representation for both the centralization efficiency, and the
##        amount of instability.
##      . Default output format for the graphical report is now postscript : file size is smaller than
##        pdf, and conversion to 24-bit lzw-compressed multipage tif is ~10x faster (using the new
##        "gc5psconv" bash script, based on ghostscript only).
##
## v5.6 20110721
##      . Added a "mod" option to support other technologies than pure-CGH. First (functional)
##        intent is to support CGH+SNP arrays. Later, support for other mods (notably "ch3" and
##        "chip") is present but very preliminary.
##      . Added the "mod" value to the pqc table.
##
## v5.5  20110630
##      . Rewrote the header to fit the team's new requirements.
##      . Changed the names of the main, most used options. Now:
##          - 'scut' is now 'nrf', standing for Normality Range Factor.
##          - 'cfac' is now 'mpf', standing for Main Peak Factor.
##          - 'colfit' is new 'dyefit', more rigorous.
##          - ('uSD' remains unchanged.)
##      . Added chromosome names to plots.
##      . Added a new output file containing all the parameters used and the quality control
##        values that are generated. One file is generated for each analyzed array, which
##        could be concatenated by a perl or R script afterwards (not written yet!). Its
##        extension is ".pqc" for "parameters and quality controls". This file will be
##        usefull for global tracking and the opportunity to fill the QC tables way more
##        easily.
##      . Added a new option to skip the GCfit loop.
##      . Added automatic support for compressed FEX files (.fex.bz2, .fex.xz, .fex.gz).
##      . Changed the dimensions of the PDF output file, to fit more with widescreen, also increased
##        its resolution.
##      . Changed color for gains, according to the new official recommanditions.
##      . Thus, heavily modified the plots.
##      . Purged some commented lines, not needed anymore.
##
## 5.4b3  20110506
##      . Changed data format in GCX output : now log2(ratio) values instead of log10(ratio).
##      . Added an option to modify the mount point of /proj on pelican, so that the script
##        can be called from a pelican session or a local one, whichever the user's /proj mount
##        point is.
##      . Shifted to English comments... will have to translate earlier ones (*some day*)
##
## 5.4b2
##      . Cosmétique : ajout du numéro d'array pour toutes les figures, pour captures
##        d'écran en vue de badigeonner le PPT...
##      . Correction de typos...
##      . Elimination du terme '_profiles' de certains noms de fichiers.
##
## 5.4b
##      . Ajout d'un tweak chrY pour les mêmes raisons que chrX.
##
## 5.4
##      . Prise en charge des véritables positions de fin de sonde de fin de segment.
##      . Support réel (dans les fichiers de sortie) du hg. Désormais, seuls les
##  	    ProbeName, rMedianSignal et gMedianSignal sont utilisés à partir des
##        fichiers cbs. Les positions génomiques utilisées seront systématiquement
##		    issues du fichier contenant les calculs de GC. On peut donc faire un shift
##		    de génome directement depuis GC5, désormais !
##		  . Changement du paramètre nperm pour DNAcopy, mis à 20000 pour plus de robustesse.
##
## 5.3b2
##      . Typos...
##
## 5.3b
##      . Réelle prise en charge du build hg18 en plus du build hg19.
##
## 5.3
##      . Ajout d'un tweak pour le chrX pour réajuster son niveau souvent modifié par
##		    la loess color (surtout pour les profils masculins).
##		  . Modification des appels à la fonction density, remise aux valeurs par défaut
##		    car l'utilisation du mod 'bcv' prend énormément de temps, surtout sur les
##		    lames 1M.
##
## 5.2
##      . Modifications pour support du GC% en hg19 calculé pour toutes les sondes de
##		    l'humain. Comme Agilent ne fournit pas de coordonnées chromosomiques de leurs
##		    nouveaux formats d'array (G3) pour les anciens builds du génome, le support
##		    des hg17 et hg18 est arrêté.
##
## 5.1c
##      . Début d'une prise en charge du build génomique selon argument "hg")
##
## 5.1b
##      . Ajout d'une option de controle "mpf" pour la centralisation, pour contrôler
##        la population (fractionnelle) minimale pour détecter le pic de normalité.
##
## 5.1
##      . Ajout d'une vue de profil en karyo.
##		  . Ajout d'un premier fit couleur par défaut (amélioration globale notable
##		    sur plusieurs plans : SumRMspread en baisse, segments normaux plus
##		    cohérents).
##		  . Ajout d'un plot des intensités pour contrôler que le test est bien en cy5
##		    et la ref en cy3.
##
## 5.0
##      . Nouvelle méthode par boucle avec sélection automatique du meilleur fit,
##		    avec arrêt dès qu'aucun fit n'améliore le profil.
##      . Remplacement du fit linéaire de la piste GC puis soustraction au profil,
##        par un loessfit, bien plus performant.
##
## 4.0
##      . Nouvelle version qui teste à chaque étape lequel des deux types de fits
##		    possibles (entre fit pour chaque dye, et fit sur le log2ratio) est le
##		    meilleur (ça peut aussi être aucun des deux!).
##		  . Nouvelle mesure de débiaisement : la somme  des différences absolues sur
##		    un lissage (via runmed) du profil.
##
## 1.1
##      . Modification de la routine de centralisation pour ignorer les gonosomes.
##		  . Ajout d'une option pour inverser le profil final (en cas d'hybridation inversee).
##
## 1.0
##      . 1ere version.

## CREATING THE CHROM2CHR AND CHR2CHROM LIST
# chrconv.list <<- list(
#   chrom2chr = list(
#     hs=list(
#       "chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chr21"=21,"chr22"=22,"chrX"=23,"chrY"=24
#     ),
#     mm=list(
#       "chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chrX"=20,"chrY"=21
#     ),
#     rn=list(
#       "chr1"=1,"chr2"=2, "chr3"=3, "chr4"=4,"chr5"=5,"chr6"=6,"chr7"=7,"chr8"=8,"chr9"=9,"chr10"=10,"chr11"=11,"chr12"=12,"chr13"=13,"chr14"=14,"chr15"=15,"chr16"=16,"chr17"=17,"chr18"=18,"chr19"=19,"chr20"=20,"chrX"=21,"chrY"=22
#     )
#   ),
#   chr2chrom=list(
#     hs=list(
#       "1"="chr1","2"="chr2", "3"="chr3", "4"="chr4","5"="chr5","6"="chr6","7"="chr7","8"="chr8","9"="chr9","10"="chr10","11"="chr11","12"="chr12","13"="chr13","14"="chr14","15"="chr15","16"="chr16","17"="chr17","18"="chr18","19"="chr19","20"="chr20","21"="chr21","22"="chr22","23"="chrX","24"="chrY"
#     ),
#     mm=list(
#       "1"="chr1","2"="chr2", "3"="chr3", "4"="chr4","5"="chr5","6"="chr6","7"="chr7","8"="chr8","9"="chr9","10"="chr10","11"="chr11","12"="chr12","13"="chr13","14"="chr14","15"="chr15","16"="chr16","17"="chr17","18"="chr18","19"="chr19","20"="chrX","21"="chrY"
#     ),
#     rn=list(
#       "1"="chr1","2"="chr2", "3"="chr3", "4"="chr4","5"="chr5","6"="chr6","7"="chr7","8"="chr8","9"="chr9","10"="chr10","11"="chr11","12"="chr12","13"="chr13","14"="chr14","15"="chr15","16"="chr16","17"="chr17","18"="chr18","19"="chr19","20"="chr20","21"="chrX","22"="chrY"
#     )
#   )
# )


cat("\nWARNING !\n*********\n1) Script chrload.R from svn (genomics/CGH/R/chrload.R) must be sourced before execution.\n\n2) If you want to generate a TIF output file, please check that the\ngc5psconv script is located in an executable path in your system !\n\n3) If bzip2 is installed on your system, GCX files will automatically\nget compressed to bz2.\n\n")
mysystem <- .Platform$OS.type

## IMPORT CHR DATA
# chrload <- function(sp, gb) {
#   
#   sptxt <- sp
#   if (sp == "hs") sptxt <- "hg"
#   gv <- paste(sptxt, gb, sep="")
#   
#   ## IMPORT BANDES CYTO ET PREPARATION POUR PLOT KARYO
#   cat("Importing chromosomes data  ...\n")
#   cytob <- read.table(paste(mnt.proj, "cgh/", gv, "/cytoBandIdeo.", gv, sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_", fill=T)
#   
#   ## Filtering out invalid lines (UCSC recently added incomplete lines which seem to refer to chrN_unknown sequence sets, which is quite totally illogical. BJ 20140227)
#   cytob <- cytob[which(!is.na(cytob$chromStart)),]
#   ## Filtering out ANY chromosome not defined in chrconvlist$chrom2chr
#   cytob <- cytob[which(cytob[["X.chrom"]] %in% names(chrconv.list$chrom2chr[[sp]])),]
#   
#   cytob$chrA <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
#   cytob$chr <- cytob$chrA
#   if (sp == "hs") {
#     cytob$chr[which(cytob$chr == "X")] <- 23
#     cytob$chr[which(cytob$chr == "Y")] <- 24
#   } else if (sp == "mm") {
#     cytob$chr[which(cytob$chr == "X")] <- 20
#     cytob$chr[which(cytob$chr == "Y")] <- 21
#   } else if (sp == "rn") {
#     cytob$chr[which(cytob$chr == "X")] <- 21
#     cytob$chr[which(cytob$chr == "Y")] <- 22
#   }
#   
#   cytob$chr <- as.numeric(cytob$chr)
#   cytob <- cytob[order(cytob$chr),]
#   
#   lchrxx <- sapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) })
#   lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
#   glen <- sum(as.numeric(lchrxx))
#   lchrsum <- sapply(c(1:length(lchrtoadd)), function(x) { sum(lchrtoadd[1:x]) })
#   
#   cytob$x1 <- 0.15
#   cytob$x2 <- 0.85
#   cytob$x1[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.25
#   cytob$x2[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.75
#   cytob$gieStain[which(cytob$gieStain == "gneg")] <- "white"
#   cytob$gieStain[which(cytob$gieStain == "gpos25")] <- "grey75"
#   cytob$gieStain[which(cytob$gieStain == "gpos33")] <- "grey66"
#   cytob$gieStain[which(cytob$gieStain == "gpos50")] <- "grey50"
#   cytob$gieStain[which(cytob$gieStain == "gpos66")] <- "grey33"
#   cytob$gieStain[which(cytob$gieStain == "gpos75")] <- "grey25"
#   cytob$gieStain[which(cytob$gieStain == "gpos100")] <- "black"
#   cytob$gieStain[which(cytob$gieStain == "acen")] <- "yellow"
#   cytob$gieStain[which(cytob$gieStain == "gvar")] <- "darkred"
#   cytob$gieStain[which(cytob$gieStain == "stalk")] <- "darkolivegreen"
#   
#   return(list(cytob=cytob, lchrxx=lchrxx, lchrsum=lchrsum, lchrtoadd=lchrtoadd, glen=glen, cur.glen=glen, gv=gv))
# }


## FUNCTION TO CONVERT CHROM TO CHR (dir="chrom2chr") OR CHR TO CHROM (dir="chr2chrom")
# chrconv <- function(dir, sp, c) {
#   return(chrconv.list[[dir]][[sp]][[c]])
# }


# ## FUNCTION TO PERFORM THE GC-FIT LOOPS ON AN INTENSITY DISTRIBUTION
# int.gcloop <- function(current, dye) {
#   
#   ### FITLOOP
# 	cat("Initiating fitloop for Cy", dye, " (", current$rm.mad,") ...\n", sep="")
# 	
# 	minigc <- gcdata[,c(8:ncol(gcdata))]
# 	gcheads <- colnames(gcdata)[c(8:ncol(gcdata))]
# 	b <- 1000
# 	
# 	while ( (b != 1) & (ncol(minigc) != 0) ) {
# 		
# 		biggy <- list()
# 		biggy <- append(biggy, list(current))
# 		rmtest <- current$rm.mad
# 		for (z in 1:length(minigc)) {
# 			tempfit <- int.fit(biggy[[1]]$int, minigc[,z])
# 			biggy <- append(biggy, list(tempfit))
# 			rmtest <- c(rmtest, tempfit$rm.mad)
# 		}
# 		b <- which.min(rmtest)
# 		if (b > 1) {
# 			
# 			current <- biggy[[b]]
# 			cat(b, "\t", gcheads[b-1], current$rm.mad, "\n")
# 			
#       int.profile.plot(current, tit=paste(gcheads[b-1], "-fitted, median-centered", sep=""), dye=dye)
#       
# 			minigc <- as.data.frame(minigc[,-c(b-1)])
# 			gcheads <- gcheads[-c(b-1)]
# 		}
#   }
#   return(current)
# }

## FUNCTION TO PERFORM THE GC-FIT LOOPS ON A LOG2(RATIO) DISTRIBUTION
l2r.gcloop <- function(current, gcd) {
  
  ### FITLOOP
  minigc <- gcd[,c(8:ncol(gcd))]
  gcheads <- colnames(gcd)[c(8:ncol(gcd))]
  b <- 1000
  
  while ( (b != 1) & (ncol(minigc) != 0) ) {
    
    biggy <- list()
    biggy <- append(biggy, list(current))
    rmtest <- current$rm.mad
    for (z in 1:length(minigc)) {
      tempfit <- l2r.fit(biggy[[1]]$l2r, minigc[,z])
      biggy <- append(biggy, list(tempfit))
      rmtest <- c(rmtest, tempfit$rm.mad)
    }
    b <- which.min(rmtest)
    if (b > 1) {
      
      current <- biggy[[b]]
      cat(b, "\t", gcheads[b-1], current$rm.mad, "\n")
      l2r.profile.plot(current, gcd, tit=paste(samplename, " ", gcheads[b-1], "-fitted, median-centered", sep=""))
      
      ## Adding values to PQC
      pqc.header.full <- c(pqc.header.full, paste("L2R.", gcheads[b-1], "-fit.MAD", sep=""), paste("L2R.", gcheads[b-1], "-fit.SSAD", sep=""))
      pqc.value.full <- c(pqc.value.full, current$l2r.mad, current$rm.mad)
      
      minigc <- as.data.frame(minigc[,-c(b-1)])
      gcheads <- gcheads[-c(b-1)]
    }
  }
  return(current)
}

## FONCTION LOESSFIT GC vs INT
int.fit <- function(int, gc) {
  intfN <- loessFit(int, gc)
  intN <- int-intfN$fitted
  sdN <- sd(intN)
  Nspread <- median(abs(diff(intN)))
  rmed <- as.numeric(runmed(intN, smo))
  Nrmspread <- sum(abs(diff(rmed)))
  return(list(int=intN, int.sd=sdN, int.mad=Nspread, rm=rmed, rm.mad=Nrmspread))
}

## FONCTION LOESSFIT GC vs LOG2(RATIO)
l2r.fit <- function(l2r, gc) {
  l2fN <- loessFit(l2r, gc)
  l2N <- l2r-l2fN$fitted
  sdN <- sd(l2N)
  Nspread <- median(abs(diff(l2N)))
  rmed <- as.numeric(runmed(l2N, smo))
  Nrmspread <- sum(abs(diff(rmed)))
  return(list(l2r=l2N, l2r.sd=sdN, l2r.mad=Nspread, rm=rmed, rm.mad=Nrmspread))
}

## FUNCTION TO DRAW GENOMIC INTENSITY PROFILE
int.profile.plot <- function(intpack, gcd, trackname, tit, track.type) {
  if (track.type == "Test") lcol <- 2 else if (track.type == "Ref") lcol <- 3
  yl <- c(0,0)
  if (mod == "cgh") yl = c(-1,1)
  else yl = c(-1.5,1.5)
  plot(gcd[["Genomic.start"]], intpack[["int"]], pch=20, cex=0.25, ylim=yl, col="grey80", main=paste(track.type, " ", trackname, " ", tit, " ", " profile | MAD = ", round(intpack$int.mad, digits=4), " | SSAD = ", round(intpack$rm.mad, digits=4), sep=""), xlab="genomic position", ylab=paste("log10(int)", sep=""))
  lines(gcd[["Genomic.start"]], intpack[["rm"]], col=lcol)
  abline(h=0, lty=2)
  for (k in cs$lchrsum) abline(v=k, lty=3, col=1)
  text(cs$lchrsum+cs$lchrxx/2, yl[2], unique(cs$cytob$chrA))
}

## FUNCTION TO DRAW GENOMIC LOG2RATIO PROFILE
l2r.profile.plot <- function(l2rpack, gcd, tit) {
  plot(gcd[["Genomic.start"]], l2rpack[["l2r"]], pch=20, cex=0.25, col="grey80", main=paste(tit, " profile | MAD = ", format(l2rpack$l2r.mad, digits=4), " | SSAD = ", format(l2rpack$rm.mad, digits=4), sep=""), ylim=c(-ymax,ymax), xlab="genomic position", ylab="log2(ratio)")
  lines(gcd[["Genomic.start"]], l2rpack[["rm"]], col=4)
  abline(h=0, lty=2)
  for (k in cs$lchrsum) abline(v=k, lty=3, col=1)
  text(cs$lchrsum+cs$lchrxx/2, ymax*(-2*(unique(cs$cytob$chr)%%2)+1), unique(cs$cytob$chrA))
}

## FUNCTION TO ADD SEGMENTS TO THE GENOMIC LOG2RATIO PROFILE
l2r.seg.add <- function (seg.obj, gcd, cut) {
  amp.index <- which(seg.obj$seg.med > ymax)
  gain.index <- which(seg.obj$seg.med >= cut)
  loss.index <- which(seg.obj$seg.med <= -cut)
  del.index <- which(seg.obj$seg.med < -ymax)
  normal.index <- which( (seg.obj$seg.med > -cut) & (seg.obj$seg.med < cut) )
  segments(seg.obj$genstart[normal.index], seg.obj$seg.med[normal.index], seg.obj$genend[normal.index], seg.obj$seg.med[normal.index], col=1, lwd=4)
  segments(seg.obj$genstart[gain.index], seg.obj$seg.med[gain.index], seg.obj$genend[gain.index], seg.obj$seg.med[gain.index], col=4, lwd=4)
  segments(seg.obj$genstart[loss.index], seg.obj$seg.med[loss.index], seg.obj$genend[loss.index], seg.obj$seg.med[loss.index], col=2, lwd=4)
  for (ai in amp.index) {
    over.probes <- which((gcd[["Genomic.end"]] >= seg.obj$genstart[ai]) & (gcd[["Genomic.start"]] <= seg.obj$genend[ai]))
    points(gcd[["Genomic.start"]][over.probes], rep(ymax, length(over.probes)), pch=2, col="cyan4")
  }
  for (di in del.index) {
    under.probes <- which((gcd[["Genomic.end"]] >= seg.obj$genstart[di]) & (gcd[["Genomic.start"]] <= seg.obj$genend[di]))
    points(gcd[["Genomic.start"]][under.probes], rep(-ymax, length(under.probes)), pch=6, col="orangered4")
  }
  abline(h=scut, col="cyan4", lty=2)
  abline(h=-scut, col="orangered4", lty=2)
}

## FUNCTION TO DRAW THE GENOMIC KARYOTYPIC LOG2RATIO PROFILE
l2r.karyo.plot <- function(l2rpack, gcd, seg.obj, cut, sp) {
  oripar <- par(no.readonly = TRUE)
  par(mgp=c(0,0,0), mar=c(1,0,1,0), omi=c(0,0.2,0.1,0), xaxt="n", yaxt="n", bty="n")
  if (sp == "hs") {
    zone <- matrix(seq(1,48,1), nrow=2, ncol=24, byrow=T)
    layW = rep(c(0.01, .031667), 12)
  } else if (sp == "mm") {
    zone <- matrix(seq(1,44,1), nrow=2, ncol=22, byrow=T)
    layW = rep(c(0.01, .034545), 11)
  }
  
  layH = c(0.5,0.5)
  chrnames <- unique(cs$cytob$chrA)
  layout(zone, widths=layW, heights = layH)
  for (k in unique(gcd$Chr)) {
    
    plot(0, 0, xlim=c(0,1), ylim=c(-cs$lchrxx[1], 0), type="n", xlab="")
    cytochr <- which(cs$cytob$chr == k)
    rect(cs$cytob$x1[cytochr], -cs$cytob$chromStart[cytochr], cs$cytob$x2[cytochr], -cs$cytob$chromEnd[cytochr], col=cs$cytob$gieStainCol[cytochr])
    segchr.loss <- which((seg.obj$chrom == k) & (seg.obj$seg.med < -cut))
    segchr.gain <- which((seg.obj$chrom == k) & (seg.obj$seg.med > cut))
    if (length(segchr.loss) > 0) segments(0, -seg.obj$loc.start[segchr.loss], 0, -seg.obj$loc.end[segchr.loss], col=2, lwd=3)
    if (length(segchr.gain) > 0) segments(1, -seg.obj$loc.start[segchr.gain], 1, -seg.obj$loc.end[segchr.gain], col=4, lwd=3)
    
    chrprob <- which(gcd$Chr == k)
    plot(0, 0, xlim=c(-1.5,1.5), ylim=c(-cs$lchrxx[1], 0), type="n", xlab="")
    points(l2rpack$l2r[chrprob], -gcd$Start[chrprob], cex=0.2, col="grey80")
    if (length(segchr.loss) > 0) rect(0, -seg.obj$loc.start[segchr.loss], seg.obj$seg.med[segchr.loss], -seg.obj$loc.end[segchr.loss], col=2, border=2)
    if (length(segchr.gain) > 0) rect(0, -seg.obj$loc.start[segchr.gain], seg.obj$seg.med[segchr.gain], -seg.obj$loc.end[segchr.gain], col=4, border=4)
    
    abline(v=0, col="grey50", lty=2)
    text(-1, -5e+06, chrnames[k], pos=3, cex=1.5)
    lines(l2rpack$rm[chrprob], -gcd$Start[chrprob], col="grey30")
  }
  layout(matrix(c(1), 1, 1, byrow=T))
  par(oripar)
}

## FUNCTION TO DRAW CHROMOSOMAL LOG2RATIO PROFILES
l2r.profile.chr.plot <- function(l2rpack, gcd, kc) {
  kc.probes <- which(gcd[["Chr"]] == kc)
  plot(gcd[["Start"]][kc.probes], l2rpack[["l2r"]][kc.probes], pch=20, cex=0.5, col="grey80", main=paste("Chromosome ", kc, sep=""), xlim=c(0, cs$lchrxx[1]), ylim=c(-ymax,ymax), xlab="chromosomal position", ylab="log2(ratio)")
  lines(gcd[["Start"]][kc.probes], l2rpack[["rm"]][kc.probes], col=4)
  abline(h=0, lty=2)
}

## FUNCTION TO ADD SEGMENTS TO THE CHROMOSOMAL LOG2RATIO PROFILES
l2r.seg.chr.add <- function(seg.obj, gcd, cut, kc) {
  kc.seg <- which(seg.obj[["chrom"]] == kc)
  kc.gcd <- which(gcd[["Chr"]] == kc)
  seg.obj <- seg.obj[kc.seg,]
  amp.index <- which(seg.obj$seg.med > ymax)
  gain.index <- which(seg.obj$seg.med >= cut)
  loss.index <- which(seg.obj$seg.med <= -cut)
  del.index <- which(seg.obj$seg.med < -ymax)
  normal.index <- which( (seg.obj$seg.med > -cut) & (seg.obj$seg.med < cut) )
  segments(seg.obj$loc.start[normal.index], seg.obj$seg.med[normal.index], seg.obj$loc.end[normal.index], seg.obj$seg.med[normal.index], col=1, lwd=4)
  segments(seg.obj$loc.start[gain.index], seg.obj$seg.med[gain.index], seg.obj$loc.end[gain.index], seg.obj$seg.med[gain.index], col=4, lwd=4)
  segments(seg.obj$loc.start[loss.index], seg.obj$seg.med[loss.index], seg.obj$loc.end[loss.index], seg.obj$seg.med[loss.index], col=2, lwd=4)
  for (ai in amp.index) {
    over.probes <- which((gcd[["End"]][kc.gcd] >= seg.obj$loc.start[ai]) & (gcd[["Start"]][kc.gcd] <= seg.obj$loc.end[ai]))
    points(gcd[["Start"]][kc.gcd][over.probes], rep(ymax, length(over.probes)), pch=2, col="cyan4")
  }
  for (di in del.index) {
    under.probes <- which((gcd[["End"]] >= seg.obj$loc.start[di]) & (gcd[["Start"]][kc.gcd] <= seg.obj$loc.end[di]))
    points(gcd[["Genomic.start"]][kc.gcd][under.probes], rep(-ymax, length(under.probes)), pch=6, col="orangered4")
  }
  abline(h=scut, col="cyan4", lty=2)
  abline(h=-scut, col="orangered4", lty=2)
}

## FUNCTION TO DRAW INT-INT PLOTS
int.plot <- function (rfX, gfX, zval=3, tit, plot.it=T) {
  A <- (log2(rfX) + log2(gfX))/2
  M <- log2(rfX/gfX)
  r.100 <- which(rfX < 100)
  g.100 <- which(gfX < 100)
  rg.100.and <- intersect(r.100, g.100)
  rg.100.or <- union(r.100, g.100)
  Z <- (M-median(M)) / sd(M)
  Zind <- which(abs(Z) > zval)
  
  if (plot.it) {
    plot(log2(gfX), log2(rfX), pch=20, cex=0.5, col="grey60", main=paste("Intensity plot (", tit, ")", sep=""), xlab="log2(Cy3-int)", ylab="log2(Cy5-int)")
    points(log2(gfX)[rg.100.or], log2(rfX)[rg.100.or], pch=20, cex=0.5, col="grey30")
    points(log2(gfX)[rg.100.and], log2(rfX)[rg.100.and], pch=20, cex=0.5, col=1)
    points(log2(gfX)[Zind], log2(rfX)[Zind], pch=20, cex=0.5, col=4)
    abline(a=0, b=1, lty=2)
  }
  
  return(data.frame(M=M, A=A))
}

## FUNCTION TO DRAW MA-PLOTS
MA.plot <- function (M, A, zval, tit) {
  order.A <- order(A)
  A <- A[order.A]
  M <- M[order.A]
  Z <- (M-median(M)) / sd(M)
  Zind <- which(abs(Z) > zval)
  plot(A, M, pch=20, cex=0.5, col="grey60", main=paste("M-A plot for ", tit, " intensities", sep=""), xlab="A", ylab="M")
  points(A[Zind], M[Zind], pch=20, cex=0.5, col=4)
  lines(A, runmed(M, smo), col=2)
  abline(h=0, lty=2)
}

## PAREMETERS DESCRIPTOR FOR BUILD.TARGETS FUNCTION
help.build.targets <- function() {
  cat("
      PARAMETERS FOR THE build.targets() FUNCTION :
      ---------------------------------------------
      design.date = Date of the design for the corresponding AMADID, as downloadable on Agilent eArray.
      mode     = Type of technology. \"cgh\" stands for CGH or CGH+SNP, \"chip\" stands for ChIP or CH3.
      rif      = \"Red intensity filter\", intensity threshold for low-intensity spot values in cy5. ONLY FOR CHIP MODE.
      gif      = \"Green intensity filter\", intensity threshold for low-intensity spot values in cy3. ONLY FOR CHIP MODE.
      nrf      = Multiplier to determine the width of the normality band. It's a multiplier for the internal noise value. ONLY FOR CGH MODE.
      mpf      = Cut-off to consider higher peaks in the log2(ratio) distribution to detect the normality. It is a % of the highest peak, so value can go from 0 to 1. ONLY FOR CGH MODE.
      uSD      = The \"undo.sd\" parameter for the CBS implementation in DNAcopy. ONLY FOR CGH MODE.
      alpha    = The \"alpha\" parameter for the CBS implementation in DNAcopy. ONLY FOR CGH MODE.
      hyb      = Hybridation used : \"d\" for direct (Cy5/Cy3), \"r\" for reverse (Cy3/Cy5).
      inv      = Has the log2(ratio) profile to be be inversed ? Useful when the slide has been inadvertidly been hybridized with dyes swaped.
      qn       = Perform a normalization by quantiles. ONLY FOR CHIP MODE.
      dyefit   = Perform a loessfit on dye channels. Recommanded for CGH ! NOT recommanded for CHIP.
      gcfit    = Perform a loessfit for GC composition. Recommanded !
      cngfit   = Perform a loessfit for the local amount of CnG trinucleotides. Recommanded !
      center   = Perform a cgh-speficic centralization. If set to F, the profile will be simply centered on its median value.
      gb       = Genome build (18 or 19 for hs, 9 for mm).
      of       = Output file format for the report : tif (default), tiff (equally accepted), ps (postscript) or pdf (Adobe portable document format).
      \n")
}

## FUNCTION TO BUILD A DEFAULT, PRE-FILLED TARGETS FILE
build.targets <- function(design.date="latest", mode="cgh", rif=0, gif=0, nrf=0.25, mpf=0.5, uSD=0, alpha = 0.01, gb=19, hyb="d", dyefit=T, qn=F, gcfit=T, cngfit=T, center=T, of="tif") {
  
  cat("\nIMPORTANT : Please check you requested the required mode for your technology!\n\n")
  cat("Requested mode :", mode, "\n")
  fexlist <-  sort(c(dir(pattern= ".fex$"), dir(pattern= ".fex.bz2$"), dir(pattern= ".fex.gz$"), dir(pattern= ".fex.xz$")))
  
  if (length(fexlist) == 0) {
    return("No FEX file found.\n")
  } else {
    cat(length(fexlist), "FEX file(s) found.\n")
    fxlist <- data.frame(t(sapply(fexlist, function(x) { c(unlist(strsplit(unlist(strsplit(x, ".", fixed=T))[1], "_", fixed=T)), x) })), stringsAsFactors=F, check.names=F)
    colnames(fxlist) <- c("AMADID", "Slide", "Row", "Col", "Filename")
    targetlist <- data.frame(AMADID=fxlist$AMADID,
                             Test.fex=fexlist,
                             Test.barcode=paste(fxlist$Slide, fxlist$Row, fxlist$Col, sep="_"),
                             Test.track="R", stringsAsFactors=F)
    targetlist$Ref.fex <- targetlist$Test.fex
    targetlist$Ref.barcode <- targetlist$Test.barcode
    targetlist$Ref.track <- "G"
    targetlist$AMADID <- as.numeric(targetlist$AMADID)
    targetlist$Genome.build <- gb
    targetlist$Design.date <- "latest"
    if (mode == "chip") {
      targetlist$Test.int.filter <- rif
      targetlist$Ref.int.filter <- gif
      targetlist$Quantiles.norm <- qn
    }
    if (mode == "cgh") {
      targetlist$Normal.range.factor <- nrf
      targetlist$Main.peak.factor <- mpf
      targetlist$undo.SD <- uSD
      targetlist$alpha <- alpha
    }
    targetlist$Centralization <- center
    targetlist$Hybridization <- hyb
    targetlist$Dye.norm <- dyefit
    targetlist$GC.norm <- gcfit
    targetlist$CnG.norm <- cngfit
    targetlist$Output.format <- of
    targetlist$Perform <- T
    targetlist$Add.title <- ""
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    toutname <- paste("targets_", mode, "_", timestamp, ".txt", sep="")
    write.table(targetlist, toutname, sep="\t", quote=F, row.names=F)
    cat(toutname, "created.\n")
    return( list.files(path=getwd(), pattern=toutname, full.names=T) )
  }
}

## DISPLAY THE DESIGNS TABLE
# designs.list <- function(mnt.db="/mnt/db/") {
#   d <- read.table(paste(mnt.db, "/Agilent/array_design/designs.list", sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
#   d <- d[order(d$AMADID, d$Genome.build, d$Design.date),]
#   d
# }
designs.list <- function(ldb="/mnt/data_cigogne/bioinfo/") {
  d <- read.table(paste(ldb, "/Agilent/array_design/designs.list", sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  d <- d[order(d$AMADID, d$Genome.build, d$Design.date),]
  return(d)
}

## PAREMETERS DESCRIPTOR FOR GC5FIT FUNCTION
help.gc5fit <- function() {
  cat("
      PARAMETERS FOR THE gc5fit() FUNCTION :
      --------------------------------------
      targets.list     = (path+)Name of the targets file, either manually made, or made using the build.targets() function.
      ldb              = Path to a local data base.
      ncores           = Number of cores (threads) to use for a multithreaded run (thanks to snowfall).
      \n")
}

## GC5 RUN FUNCTION : SNOWFALL MULTITHREADED
gc5fit <- function(targets.list=NULL, ldb="/mnt/data_cigogne/bioinfo/", ncores=2) {
  
  if (is.null(targets.list)) return("Targets file required !")
  
  ## Globalizing
  # mnt.proj <<- mnt.proj
  # mnt.db <<- mnt.db
  ldb <<- ldb
  ymax <<- 1.5
  
  ## Load targetslist
  tlist <<- read.table(targets.list, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  
  ## Checking mode by the presence of the "undo.SD" column
  if ("undo.SD" %in% colnames(tlist)) xmod <- "cgh" else xmod <- "chip"
  
  ## Requesting timestamp
  tstamp <<- format(Sys.time(), "%Y%m%d%H%M%S")
  
  ## Creating outdir
  dir.create(tstamp)
  ## Copying targets file
  file.copy(targets.list, tstamp)
  
  ## Removing arrays that won't be processed
  perf.check <- which(tlist[["Perform"]] == TRUE)
  
  cat("Found ", nrow(tlist[perf.check,]), " arrays to process : ", tlist[perf.check,][["Test.barcode"]],".\n")
  
  ## Checking the compatibility of cores and requested processes
  if (nrow(tlist[perf.check,]) < ncores) {
    ncores <- nrow(tlist[perf.check,])
    cat("\nNumber of arrays to process has been found lower than the requested number of cores.\nLowering ncores to ", ncores, "\n\n", sep="")
  }
  if (ncores > 1) {
    ## Preparing multithread
    library(snowfall)
    sfStop()
    sfInit(parallel=T,cpus=ncores)
    sfExportAll()
    
    ## Running on the targets list
    pqc.gc5 <- sfSapply(perf.check, gc5.run)
    sfStop()
  } else {
    pqc.gc5 <- sapply(c(perf.check), gc5.run)
  }
  
  if (xmod == "cgh") {
    rownames(pqc.gc5) <- c("Barcode", "AMADID", "Design.date", "Mode", "Species", "Genome.build", "Polarity", "Dye.norm", "GC.norm", "Normal.range.factor", "Main.peak.factor", "Centralization", "undo.SD", "alpha", "Nb.probes", "L10Test.raw.SD", "L10Test.raw.SSAD", "L10Ref.raw.SD", "L10Ref.raw.SSAD", "L2R.raw.MAD", "L2R.raw.SSAD", "L2R.final.MAD", "L2R.final.SSAD", "Aberration.L2R.cutoff")
  } else if (xmod == "chip") {
    rownames(pqc.gc5) <- c("Barcode", "AMADID", "Design.date", "Mode", "Species", "Genome.build", "Polarity", "Dye.norm", "GC.norm", "Quantiles.norm", "Centralization", "Nb.probes", "L10Test.raw.SD", "L10Test.raw.SSAD", "L10Ref.raw.SD", "L10Ref.raw.SSAD", "L2R.raw.MAD", "L2R.raw.SSAD", "L2R.final.MAD", "L2R.final.SSAD")
  }
  write.table(pqc.gc5, paste(tstamp, "/", tstamp, "_", ncol(pqc.gc5), "s.pqc", sep=""), quote=F, col.names=F, sep="\t")
}

## INTERMEDIATE CALLING FUNCTION TO SEPARATE CGH- or CHIP- MODES
gc5.run <- function(z) {
  x <- as.vector(tlist[z,])
  if ("undo.SD" %in% colnames(x)) {
    gc5.core(amadid=x[["AMADID"]], test.fex=x[["Test.fex"]], test.barcode=x[["Test.barcode"]], ref.fex=x[["Ref.fex"]], ref.barcode=x[["Ref.barcode"]], test.track=x[["Test.track"]], ref.track=x[["Ref.track"]], gb=x[["Genome.build"]], design.date=x[["Design.date"]], nrf=x[["Normal.range.factor"]], mpf=x[["Main.peak.factor"]], uSD=x[["undo.SD"]], alpha=x[["alpha"]], hyb=x[["Hybridization"]], dyefit=x[["Dye.norm"]], gcfit=x[["GC.norm"]], cngfit=x[["CnG.norm"]], center=x[["Centralization"]], of=x[["Output.format"]], add.title=x[["Add.title"]])
  } else {
    gc5.core(amadid=x[["AMADID"]], test.fex=x[["Test.fex"]], test.barcode=x[["Test.barcode"]], ref.fex=x[["Ref.fex"]], ref.barcode=x[["Ref.barcode"]], test.track=x[["Test.track"]], ref.track=x[["Ref.track"]], gb=x[["Genome.build"]], design.date=x[["Design.date"]], test.if=x[["Test.int.filter"]], ref.if=x[["Ref.int.filter"]], hyb=x[["Hybridization"]], dyefit=x[["Dye.norm"]], qn=x[["Quantiles.norm"]], gcfit=x[["GC.norm"]], cngfit=x[["CnG.norm"]], center=x[["Centralization"]], of=x[["Output.format"]], add.title=x[["Add.title"]])
  }
}

## CORE GC5 FUNCTION
gc5.core <<- function(amadid=amadid, test.fex=test.fex, test.barcode=test.barcode, test.track=test.track, ref.fex=ref.fex, ref.barcode=ref.barcode, ref.track=ref.track, gb=gb, design.date=design.date, test.if=test.if, ref.if=ref.if, nrf=nrf, mpf=mpf, uSD=uSD, alpha = alpha, hyb=hyb, dyefit=dyefit, qn=qn, gcfit=gcfit, cngfit=cngfit, center=center, of=of, add.title=add.title) {
  
  library(DNAcopy)
  library(limma)
  
  cat(test.barcode, test.track, test.fex, ref.barcode, ref.track, ref.fex, "\n")
  ## Creating and globalizing samplename
  if (test.barcode == ref.barcode) {
    if (test.track == ref.track) {
      stop("ERROR : same track selected as Test and Ref from the same array !")
    } else if ((test.track == "R") & (ref.track == "G")) {
      samplename <<- test.barcode
    } else {
      samplename <<- paste(test.track, test.barcode, '.', ref.track, ref.barcode, sep="")
    }
  } else {
    samplename <<- paste(test.track, test.barcode, '.', ref.track, ref.barcode, sep="")
  }
  
  if(is.na(add.title)) add.title <- NULL
  ## Loading designs list
  # dl <- read.table(paste(mnt.db, "/Agilent/array_design/designs.list", sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  dl <- read.table(paste(ldb, "/Agilent/array_design/designs.list", sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  
  ## Getting design info
  select.d.index <- 0
  if (design.date == "latest") {
    d.index <- which(dl[["AMADID"]] == amadid & dl[["Genome.build"]] == gb)
    select.d.index <- d.index[which.max(dl[["Design.date"]][d.index])]
    
  } else {
    select.d.index <- d.index <- which(dl[["AMADID"]] == amadid & dl[["Genome.build"]] == gb & dl[["Design.date"]] == design.date)
    if (length(select.d.index) == 0) stop(paste("No design found for AMADID ", amadid, ", genome build ", gb, " and design date ", design.date, " !", sep=""))
    if (length(select.d.index) > 1) stop(paste("More than 1 design found for AMADID ", amadid, ", genome build ", gb, " and design date ", design.date, " !", sep=""))
  }
  sp <- dl[["Species"]][select.d.index]
  design.date <- dl[["Design.date"]][select.d.index]
  techno <- dl[["Technology"]][select.d.index]
  slide.format <- dl[["Format"]][select.d.index]
  design.name <- dl[["Name"]][select.d.index]
  
  ## Setting mod global variables
  mod <<- ""
  if ((techno == "CGH") | (techno == "SNP")) mod <<- "cgh"
  if ((techno == "CHIP") | (techno == "CH3")) mod <<- "chip"
  
  ## Import chromosomes data
  cs <<- chrload(sp, gb, ldb)
  gv <- cs$gv
  
  cat("Working on ", gv, "\n", sep="")
  mod <- tolower(mod)
  cat("Called mode : ", mod, "\n", sep="")

  ## Import GC data for the selected design
  gc.file <- paste("0", amadid, "_D_BED_", design.date, "_", gv, ".gc", sep="")
  # gcload <- read.table(paste(mnt.proj, "/cgh/Agilent/", cs$species, "/", gv, "/", gc.file, sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  gcload <- read.table(paste(ldb, "/Agilent/GCdata/", cs$species, "/", gv, "/", gc.file, sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
  gcdata <- gcload[,1:5]
  
  if (gcfit) {
    gcdata <- cbind(gcdata, gcload[,6:ncol(gcload)])
  }
  rm(gcload)
  
  if (cngfit) {
    ## Loading cngs
    cngs.file <- paste("0", amadid, "_D_BED_", design.date, "_", gv, ".cngs", sep="")
    # cngsdata <- read.table(paste(mnt.proj, "/cgh/Agilent/", cs$species, "/", gv, "/", cngs.file, sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
    cngsdata <- read.table(paste(ldb, "/Agilent/GCdata/", cs$species, "/", gv, "/", cngs.file, sep=""), header=T, sep="\t", check.names=F, stringsAsFactors=F)
    
    ## Merging the datasets
    gcdata <- cbind(gcdata, cngsdata[,-c(1:5)])
    rm(cngsdata)
  }
  gc()
  
  ## Elimination des sondes tripliquees de gcdata
  gcdata <- gcdata[!duplicated(gcdata$Identifier),]
  cat("Duped probes removed.\n")
  
  ## Initiating PQC vectors
  pqc.header.full <<- vector()
  pqc.value.full <<- vector()
  
  ## PQC grouping table
  msg <- cat("Working on array ", samplename," ...\n", sep="")
  
  #   pqc.header.full <- c("Sample.name", "Mode", "Species", "Genome.build", "Polarity", "Dye.norm", "GC.norm")
  pqc.header.full <- c("Sample.name", "AMADID", "Design.date", "Mode", "Species", "Genome.build", "Polarity", "Dye.norm", "GC.norm")
  #   pqc.value.full <- c(samplename, mod, sp, gv, hyb, dyefit, gcfit)
  pqc.value.full <- c(samplename, amadid, design.date, mod, sp, gv, hyb, dyefit, gcfit)
  
  if (mod == "cgh") {
    pqc.header.full <- c(pqc.header.full, "Normal.range.factor", "Main.peak.factor", "Centralization", "undo.SD", "alpha")
    pqc.value.full <- c(pqc.value.full, nrf, mpf, center, uSD, alpha)
  } else if (mod == "chip") {
    pqc.header.full <- c(pqc.header.full, "Quantiles.norm", "Centralization")
    pqc.value.full <- c(pqc.value.full, qn, center)
  }
  
  ## Loading TEST FEX file, depending on its extension
  test.extens <- rev(unlist(strsplit(test.fex, ".", fixed=T)))
  if (test.extens[1] == "fex") {
    test.fexdata <- read.table(test.fex, header=T, sep="\t", stringsAsFactors=F)
  } else if ((test.extens[1] == "bz2") & (test.extens[2] == "fex")) {
    test.fexdata <- read.table(bzfile(test.fex), header=T, sep="\t", stringsAsFactors=F)
  } else if ((test.extens[1] == "gz") & (test.extens[2] == "fex")) {
    test.fexdata <- read.table(gzfile(test.fex), header=T, sep="\t", stringsAsFactors=F)
  } else if ((test.extens[1] == "xz") & (test.extens[2] == "fex")) {
    test.fexdata <- read.table(xzfile(test.fex), header=T, sep="\t", stringsAsFactors=F)
  } else {
    stop(cat("Test input file ", test.fex, " is not .fex, or compressed in an unspported format\n"))
  }
  cat("Test track read.\n")
  
  ## Loading REF FEX file, depending on its extension
  ref.extens <- rev(unlist(strsplit(ref.fex, ".", fixed=T)))
  if (ref.extens[1] == "fex") {
    ref.fexdata <- read.table(ref.fex, header=T, sep="\t", stringsAsFactors=F)
  } else if ((ref.extens[1] == "bz2") & (ref.extens[2] == "fex")) {
    ref.fexdata <- read.table(bzfile(ref.fex), header=T, sep="\t", stringsAsFactors=F)
  } else if ((ref.extens[1] == "gz") & (ref.extens[2] == "fex")) {
    ref.fexdata <- read.table(gzfile(ref.fex), header=T, sep="\t", stringsAsFactors=F)
  } else if ((ref.extens[1] == "xz") & (ref.extens[2] == "fex")) {
    ref.fexdata <- read.table(xzfile(ref.fex), header=T, sep="\t", stringsAsFactors=F)
  } else {
    stop(cat("Ref input file ", ref.fex, " is not .fex, or compressed in an unspported format\n"))
  }
  cat("Ref track read.\n")
  
  ## Crossing TEST and REF fexdataz
  test.fexdata <- test.fexdata[test.fexdata[["ProbeName"]] %in% ref.fexdata[["ProbeName"]],]
  ref.fexdata <- ref.fexdata[ref.fexdata[["ProbeName"]] %in% test.fexdata[["ProbeName"]],]
  
  #   cat(nrow(test.fexdata), nrow(ref.fexdata), "\n")
  ## Check if NO array probes are found in the GC data
  if (nrow(test.fexdata) == 0) {
    cat("\nNO COMMON PROBE FOUND BETWEEN TEST AND REF TRACKS !\n\n")
    next()
  }
  
  fexdata <- test.fexdata[,c(1:4)]
  if (test.track == "R") fexdata[["TestSignal"]] <- test.fexdata[["rMedianSignal"]] 
  if (test.track == "G") fexdata[["TestSignal"]] <- test.fexdata[["gMedianSignal"]]
  if (ref.track == "R") fexdata[["RefSignal"]] <- ref.fexdata[["rMedianSignal"]] 
  if (ref.track == "G") fexdata[["RefSignal"]] <- ref.fexdata[["gMedianSignal"]]
  
  cat("Test & ref tracks merged.\n")
  
  ## Elimination des sondes inexistantes sur l'array
  gcdata <- gcdata[which(is.element(gcdata$Identifier, fexdata$ProbeName)),]
  fexdata <- fexdata[is.element(fexdata$ProbeName, gcdata$Identifier),]
  
  ## Check if NO array probes are found in the GC data
  if (nrow(gcdata) == 0) {
    cat("\nNO COMMON PROBE FOUND BETWEEN GC DATA AND ARRAY DATA !\n\n")
    next()
  }
  
  ## Ajout des coordonnées génomiques
  genstart <- gcdata[["Start"]] + cs$lchrsum[gcdata[["Chr"]]]
  genend <- gcdata[["End"]] + cs$lchrsum[gcdata[["Chr"]]]
  
  if (ncol(gcdata) > 5) {
    gcdata <- data.frame(gcdata[,c(1:5)], Genomic.start=genstart, Genomic.end=genend, gcdata[,c(6:ncol(gcdata))], check.names=F)
  }
  else {
    gcdata <- data.frame(gcdata, Genomic.start=genstart, Genomic.end=genend, check.names=F)
  }
  gc()
  
  ## Sorting data according to the positions given by gcdata and not those from fexdata !
  gcdata <- gcdata[order(gcdata[["Chr"]], gcdata[["Start"]]),]
  fexdata <- fexdata[order(gcdata[["Chr"]], gcdata[["Start"]]),]
  
  ## Determining the smoothing factor
  smo <<- as.integer(nrow(gcdata)/400)
  smo <<- smo +1 -smo%%2
  
  #################################################
  ## PERFORMING FILTERING ON PROBES' INTENSITIES ##
  #################################################
  if (mod == "chip") {
    if (sum(c(test.if, ref.if)) > 0) {
      cat("Filtering low intensities\n")
      ikeep <- which( (fexdata[["TestSignal"]] > test.if) | (fexdata[["RefSignal"]] > ref.if) )
      fexdata <- fexdata[ikeep,]
      gcdata <- gcdata[ikeep,]
      cat(length(ikeep), " probes kept.\n", sep="")
    }
  }
  
  ## Adding probes quantity in PQC
  pqc.header.full <- c(pqc.header.full, "Nb.probes")
  pqc.value.full <- c(pqc.value.full, nrow(fexdata))
  
  cat("All data ready.\n")
  
  ## Building output filename suffix
  suffix <- ""
  if (!dyefit) suffix <- paste(suffix, "_NOdyefit", sep="")
  if (mod == "chip") {
    if (qn) suffix <- paste(suffix, "_QN", sep="")
    if (test.if > 0) suffix <- paste(suffix, "_tif", test.if, sep="")
    if (ref.if > 0) suffix <- paste(suffix, "_rif", ref.if, sep="")
  }
  if (!gcfit) suffix <- paste(suffix, "_NOGCfit", sep="")
  if (!cngfit) suffix <- paste(suffix, "_NOCnGfit", sep="")
  if (!center) suffix <- paste(suffix, "_NOcenter", sep="") else if (mod == "cgh") { suffix <- paste(suffix, "_mpf", mpf, sep="") }
  
  
  ## Opening output file
  outname <- ""
  if (mod == "cgh") outname = paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_nrf", nrf, "_", gv, sep="") else if (mod == "chip") outname = paste(tstamp, "/", samplename, "_GC5", suffix, "_", gv, sep="")
  if ((of == "ps") | (of == "tif") | (of == "tiff")) postscript(paste(outname, ".ps", sep=""), height=21/cm(1), width=37/cm(1), paper="special", horizontal=F) else if (of == "pdf") pdf(paste(outname, ".pdf", sep=""), height=21/cm(1), width=37/cm(1))
  
  par(mgp=c(1,0,0), mar=c(2,2,2,2), xaxs="i")
  
  #############
  ##   RAW   ##
  #############
  
  ## BUILDING THE CY5 OBJECT
  rfX10 <- log10(fexdata[["TestSignal"]])
  rfX10 <- rfX10 - median(rfX10)
  rfX10.mad <- median(abs(diff(rfX10)))
  rfX10med <- as.numeric(runmed(rfX10, smo))
  rfX10med.noise <- sum(abs(diff(rfX10med)))
  currentT <- list(int=rfX10, int.sd=sd(rfX10), int.mad=rfX10.mad, rm=rfX10med, rm.mad=rfX10med.noise)
  
  pqc.header.full <- c(pqc.header.full, "L10Test.raw.SD", "L10Test.raw.SSAD")
  pqc.value.full <- c(pqc.value.full, currentT$int.sd, currentT$rm.mad)
  
  ## BUILDING THE CY3 OBJECT
  gfX10 <- log10(fexdata[["RefSignal"]])
  gfX10 <- gfX10 - median(gfX10)
  gfX10.mad <- median(abs(diff(gfX10)))
  gfX10med <- as.numeric(runmed(gfX10, smo))
  gfX10med.noise <- sum(abs(diff(gfX10med)))
  currentR <- list(int=gfX10, int.sd=sd(gfX10), int.mad=gfX10.mad, rm=gfX10med, rm.mad=gfX10med.noise)
  
  pqc.header.full <- c(pqc.header.full, "L10Ref.raw.SD", "L10Ref.raw.SSAD")
  pqc.value.full <- c(pqc.value.full, currentR$int.sd, currentR$rm.mad)
  
  ## BUILDING THE L2R OBJECT
  l2rX <- log2(10^currentT$int/10^currentR$int)
  l2rX.mad <- median(abs(diff(l2rX)))
  l2rXmed <- as.numeric(runmed(l2rX, smo))
  l2rXmed.noise <- sum(abs(diff(l2rXmed)))
  current <- list(l2r=l2rX, l2r.sd=sd(l2rX), l2r.mad=l2rX.mad, rm=l2rXmed, rm.mad=l2rXmed.noise)
  
  pqc.header.full <- c(pqc.header.full, "L2R.raw.MAD", "L2R.raw.SSAD")
  pqc.value.full <- c(pqc.value.full, current$l2r.mad, current$rm.mad)
  pqc.header.out <- pqc.header.full
  pqc.value.out <- pqc.value.full
  
  ## RAW PLOTS
  cat("Plotting the raw profiles (", current$rm.mad, ") ...\n", sep="")
  par(mfrow=c(2,1))
  int.profile.plot(currentT, gcdata, trackname=paste(test.track, test.barcode, sep="."), tit="raw, median-centered", track.type="Test")
  int.profile.plot(currentR, gcdata, trackname=paste(ref.track, ref.barcode, sep="."), tit="raw, median-centered", track.type="Ref")
  par(mfrow=c(1,1))
  if (mod == "chip") {
    layout(matrix(c(1,1,2,3), 2, 2, byrow=T))
    l2r.profile.plot(current, gcdata, tit=paste(samplename, " raw, median-centered", sep=""))
    MA.df <- int.plot(10^currentT$int, 10^currentR$int, zval=3, tit="raw")
    MA.plot(MA.df$M, MA.df$A, zval=3, tit="raw")
    layout(matrix(c(1), 1, 1, byrow=T))
  } else if (mod == "cgh") l2r.profile.plot(current, gcdata, tit=paste(samplename, " raw, median-centered", sep=""))
  
  
  #####################
  ## Optional dyefit ##
  #####################
  
  if (dyefit) {
    cat("Performing the dyes loessfit as requested...\n")
    
    if (hyb == "d") {
      glo <- loessFit(currentT$int, currentR$int)
      gfX10 <- glo$fitted
      
      ## REBUILDING THE CY3 OBJECT
      gfX10 <- gfX10 - median(gfX10)
      gfX10.mad <- median(abs(diff(gfX10)))
      gfX10med <- as.numeric(runmed(gfX10, smo))
      gfX10med.noise <- sum(abs(diff(gfX10med)))
      currentR <- list(int=gfX10, int.sd=sd(gfX10), int.mad=gfX10.mad, rm=gfX10med, rm.mad=gfX10med.noise)
      
      pqc.header.full <- c(pqc.header.full, "L10Ref.dyefit.SD", "L10Ref.dyefit.SSAD")
      pqc.value.full <- c(pqc.value.full, currentR$int.sd, currentR$rm.mad)
      
    } else if (hyb == "r") {
      glo <- loessFit(currentR$int, currentT$int)
      rfX10 <- glo$fitted
      
      ## REBUILDING THE CY5 OBJECT
      rfX10 <- rfX10 - median(gfX10)
      rfX10.mad <- median(abs(diff(rfX10)))
      rfX10med <- as.numeric(runmed(rfX10, smo))
      rfX10med.noise <- sum(abs(diff(rfX10med)))
      currentT <- list(int=rfX10, int.sd=sd(rfX10), int.mad=rfX10.mad, rm=rfX10med, rm.mad=rfX10med.noise)
      
      pqc.header.full <- c(pqc.header.full, "L10Test.dyefit.SD", "L10Test.dyefit.SSAD")
      pqc.value.full <- c(pqc.value.full, currentT$int.sd, currentT$rm.mad)
    }
    
    ## X TWEAK : GETTING INITIAL CHRX -> AUTOSOMES DISTANCE (ON RUNMED)
    print("Tweaking XY ...")
    
    Xnum <- cs$chrom2chr[["chrX"]]
    Ynum <- cs$chrom2chr[["chrY"]]
    
    X.index <- which(gcdata$Chr == Xnum)
    Y.index <- which(gcdata$Chr == Ynum)
    A.index <- which(gcdata$Chr < Xnum)
    XA.diff <- median(current$rm[A.index]) - median(current$rm[X.index])
    YA.diff <- median(current$rm[A.index]) - median(current$rm[Y.index])
    
    ## X TWEAK : MEASURING NEW DISTANCE AFTER DYEFIT
    l2rX <- log2(10^currentT$int/10^currentR$int)
    nrmed <- as.numeric(runmed(l2rX, smo))
    XA.diff2 <- median(nrmed[A.index]) - median(nrmed[X.index])
    YA.diff2 <- median(nrmed[A.index]) - median(nrmed[Y.index])
    
    ## X TWEAK : correcting the chrX
    l2rX[X.index] <- l2rX[X.index] - XA.diff + XA.diff2
    l2rX[Y.index] <- l2rX[Y.index] - YA.diff + YA.diff2
    
    ## REBUILDING THE L2R OBJECT (1)
    l2rXmed <- as.numeric(runmed(l2rX, smo))
    l2rX.mad <- median(abs(diff(l2rX)))
    l2rXmed.noise <- sum(abs(diff(l2rXmed)))
    current <- list(l2r=l2rX, l2r.sd=sd(l2rX), l2r.mad=l2rX.mad, rm=l2rXmed, rm.mad=l2rXmed.noise)
    
    pqc.header.full <- c(pqc.header.full, "L2R.dyefit.MAD", "L2R.dyefit.SSAD")
    pqc.value.full <- c(pqc.value.full, current$l2r.mad, current$rm.mad)
    
    ## DYEFITTED PLOTS
    cat("Plotting the dye-fitted profiles (", current$rm.mad, ") ...\n", sep="")
    par(mfrow=c(2,1))
    int.profile.plot(currentT, gcdata, trackname=paste(test.track, test.barcode, sep="."), tit="dye-fitted", track.type="Test")
    int.profile.plot(currentR, gcdata, trackname=paste(ref.track, ref.barcode, sep="."), tit="dye-fitted", track.type="Ref")
    par(mfrow=c(1,1))
    if (mod == "chip") {
      layout(matrix(c(1,1,2,3), 2, 2, byrow=T))
      l2r.profile.plot(current, gcdata, tit=paste(samplename, " dye-fitted", sep=""))
      MA.df <- int.plot(10^currentT$int, 10^currentR$int, zval=3, tit="dye-fitted")
      MA.plot(MA.df$M, MA.df$A, zval=3, tit="dye-fitted")
      layout(matrix(c(1), 1, 1, byrow=T))
    } else if (mod == "cgh") l2r.profile.plot(current, gcdata, tit=paste( samplename, " dye-fitted", sep=""))
  }
  
  #################
  ## Optional QN ##
  #################
  if (mod != "cgh") {
    if (qn) {
      cat("Performing the QN as requested...\n")
      qn.out <- normalizeQuantiles(matrix(c(rfX10, gfX10), ncol=2))
      
      ## REBUILDING THE CY5 OBJECT
      rfX10 <- qn.out[,1]
      rfX10 < rfX10 - median(gfX10)
      rfX10.mad <- median(abs(diff(rfX10)))
      rfX10med <- as.numeric(runmed(rfX10, smo))
      rfX10med.noise <- sum(abs(diff(rfX10med)))
      currentT <- list(int=rfX10, int.sd=sd(rfX10), int.mad=rfX10.mad, rm=rfX10med, rm.mad=rfX10med.noise)
      
      pqc.header.full <- c(pqc.header.full, "L10Test.qn.SD", "L10Test.qn.SSAD")
      pqc.value.full <- c(pqc.value.full, currentT$int.sd, currentT$rm.mad)
      
      ## REBUILDING THE CY3 OBJECT
      gfX10 <- qn.out[,2]
      gfX10 < gfX10 - median(gfX10)
      gfX10.mad <- median(abs(diff(gfX10)))
      gfX10med <- as.numeric(runmed(gfX10, smo))
      gfX10med.noise <- sum(abs(diff(gfX10med)))
      currentR <- list(int=gfX10, int.sd=sd(gfX10), int.mad=gfX10.mad, rm=gfX10med, rm.mad=gfX10med.noise)
      
      pqc.header.full <- c(pqc.header.full, "L10Ref.qn.SD", "L10Ref.qn.SSAD")
      pqc.value.full <- c(pqc.value.full, currentR$int.sd, currentR$rm.mad)
      
      ## REBUILDING THE L2R OBJECT
      l2rX <- log2(10^currentT$int/10^currentR$int)
      l2rX.mad <- median(abs(diff(l2rX)))
      l2rXmed <- as.numeric(runmed(l2rX, smo))
      l2rXmed.noise <- sum(abs(diff(l2rXmed)))
      current <- list(l2r=l2rX, l2r.sd=sd(l2rX), l2r.mad=l2rX.mad, rm=l2rXmed, rm.mad=l2rXmed.noise)
      
      pqc.header.full <- c(pqc.header.full, "L2R.qn.MAD", "L2R.qn.SSAD")
      pqc.value.full <- c(pqc.value.full, current$l2r.mad, current$rm.mad)
      
      ## QN PLOTS
      cat("Plotting the quantiles-normalized profiles (", current$rm.mad, ") ...\n", sep="")
      par(mfrow=c(2,1))
      int.profile.plot(currentT, gcdata, trackname=paste(test.track, test.barcode, sep="."), tit="quantiles-normalized", track.type="Test")
      int.profile.plot(currentR, gcdata, trackname=paste(ref.track, ref.barcode, sep="."), tit="quantiles-normalized", track.type="Ref")
      par(mfrow=c(1,1))
      if (mod == "chip") {
        layout(matrix(c(1,1,2,3), 2, 2, byrow=T))
        l2r.profile.plot(current, gcdata, tit=paste(samplename, " quantiles-normalized", sep=""))
        MA.df <- int.plot(10^currentT$int, 10^currentR$int, zval=3, tit="quantiles-normalized")
        MA.plot(MA.df$M, MA.df$A, zval=3, tit="normalized")
        layout(matrix(c(1), 1, 1, byrow=T))
      } else if (mod == "cgh") l2r.profile.plot(current, gcdata, tit=paste(samplename, " quantiles-normalized", sep=""))
    }
  }
  
  #####################
  ## Optional GCfits ##
  #####################
  if ( (gcfit) || (cngfit) ) {
    cat("Performing GC/CnGfit(s) ...\n")
    par(mfrow=c(1,1))
    current <- l2r.gcloop(current, gcdata)
  }
  
  ## ADDING FINAL PROFILE VALUES TO PQC
  #   cat("PQC\n")
  pqc.header.full <- c(pqc.header.full, "L2R.final.MAD", "L2R.final.SSAD")
  pqc.value.full <- c(pqc.value.full, current$l2r.mad, current$rm.mad)
  pqc.header.out <- c(pqc.header.out, "L2R.final.MAD", "L2R.final.SSAD")
  pqc.value.out <- c(pqc.value.out, current$l2r.mad, current$rm.mad)
  
  ####################
  ## Centralization ##
  ####################
  
  autoind <- which(gcdata$Chr < (max(cs$cytob$chr)-2))
  rmC <- as.numeric(runmed(current$l2r[autoind],smo))
  den <- density(rmC)
  fcword = "median-centered"
  if (center) {
    
    cat("Centralization ...\n")
    fcword = "centered"
    
    ## Centralization for CGH and SNP
    fx = ""
    
    if (mod == "cgh") {
#       denf <- data.frame(x=den$x, y=den$y, sign=c(1,sign(diff(den$y))))
#       dcenter <- ( quantile(rmC, 0.90) + quantile(rmC, 0.10) ) / 2
#       denf$dtc <- abs(denf$x - dcenter)
#       denf3 <- denf[which(denf$y >= max(den$y)*mpf),]
#       fx <- fdtc <- 9e+06
#       for (k in 2:nrow(denf3)) {
#         if (denf3$sign[k] < denf3$sign[k-1]) {
#           if (denf3$dtc[k-1] < fdtc) {
#             fdtc <- denf3$dtc[k-1]
#             fx <- denf3$x[k-1]
#           }
#         }
#       }
      
      denf <- data.frame(x=den$x, y=den$y, sign=c(1,sign(diff(den$y))))
      denf$sign2 <- c(diff(denf$sign), 0)
      rrr <- rle(denf$sign2)
      repr <- data.frame(values=rrr$values, start=c(0, (cumsum(rrr$lengths[-c(length(rrr$lengths))])+1)), end=cumsum(rrr$lengths), stringsAsFactors=F)
      npk <- which(repr$values == -2)
#       cat("NPK ini", npk, "\n")
      if (length(npk) == 1) {
        cat("Found a single peak.")
        fx <- den$x[repr$start[npk[1]]]
      }
      else {
        ## Filtering with absolute criterion : peak pop < 10% of the total pop
        parea <- sapply(npk, function(r) {
          sum(den$y[repr$start[r-1]:repr$end[r+1]]) / sum(den$y)
        })
        parea.rel <- parea / max(parea)
        npk <- npk[which((parea > .1) & (parea.rel >= mpf))]
        
#         cat("PAREA", parea, "\n")
#         cat("PAREA.REL", parea.rel, "\n")
#         cat("NPK SELEC", npk, "\n")
        ## Filtering with relative criterion : peak pop < 33% of the most pop one.
        peaklist <- repr$start[npk]
#         cat("PEAKLIST", peaklist, "\n")
        fx <- denf$x[peaklist[which.min(sapply(peaklist, function(p) { abs(sum(denf$y[1:(p-1)]) - sum(denf$y[(p+1):nrow(denf)]))/sum(denf$y) }))]]
      }
    } else if (mod == "chip") fx <- den$x[which.max(den$y)]
    
    ## MODIFYING THE L2R OBJECT
    current$l2r <- current$l2r - fx
    current$rm <- current$rm - fx
    
    den$x <- den$x - fx
  }
  
  #######################
  ## Profile inversion ##
  #######################
  
  # 	if (inv) {
  # 		cat("Inverting profile as requested ...\n")
  # 		current$l2r <- -current$l2r
  # 		current$rm <- as.numeric(runmed(current$l2r, smo))
  # 	}
  
  ############################################
  ## Final pangenomic plot +density +qqnorm ##
  ############################################
  
  cat("Plotting final genomic plot +density +qqnorm ...\n")
  layout(matrix(c(1,1,2,3), 2, 2, byrow=T))
  l2r.profile.plot(current, gcdata, tit=paste(samplename, " final, centered", sep=""))
  plot(den, type="l", main=paste(samplename, ", ", fcword, " distribution of probes' log2(ratio)", sep=""), xlab="log2(ratio)", ylab="fractions")
  abline(v=0, lty=2)
  qqn.l2r <- qqnorm(current$l2r, plot.it=F)
  plot(qqn.l2r, main="QQ-plot of the normalized log2(ratio) distribution", xlab="Theoretical quantiles", ylab="Sample quantiles", pch=20)
  qqline(current$l2r, col=2)
  layout(matrix(c(1), 1, 1, byrow=T))
  
  #######################
  ## CGH-specific part ##
  #######################
  
  if (mod == "cgh") {
    
    scut <<- nrf*current$l2r.mad
    
    pqc.header.full <- c(pqc.header.full, "Aberration.L2R.cutoff")
    pqc.value.full <- c(pqc.value.full, scut)
    pqc.header.out <- c(pqc.header.out, "Aberration.L2R.cutoff")
    pqc.value.out <- c(pqc.value.out, scut)
    
    ## SEGMENTATION
    cat("Segmenting final profile ...\n")
    
    DYE.fed <- cbind(gcdata[,c(1,3,4)], current$l2r)
    colnames(DYE.fed) <- c("Clone", "Chromosome", "Position", "Log2Rat")
    DYE.fed$Position <- DYE.fed$Position
    
    DYE.CNA <- CNA(as.matrix(DYE.fed$Log2Rat), DYE.fed$Chromosome, DYE.fed$Position, data.type = "logratio", sampleid = paste("Sample", samplename, sep="_"))
    if (uSD == 0) { seg.DYE.CNA <- segment(DYE.CNA, min.width=3, nperm=20000, alpha = alpha) } else { seg.DYE.CNA <- segment(DYE.CNA, undo.splits = "sdundo", undo.SD = uSD, alpha = alpha, min.width=3, nperm=20000) }
    
    ## Getting real segments ends
    ## Probes index for start and end
    seg.DYE.CNA$output$mark.start.i <- c(1, sapply(c(1:(nrow(seg.DYE.CNA$output)-1)), function(x) { (sum(seg.DYE.CNA$output$num.mark[1:x])+1) }))
    seg.DYE.CNA$output$mark.end.i <- c(seg.DYE.CNA$output$mark.start.i[2:length(seg.DYE.CNA$output$mark.start.i)]-1, sum(seg.DYE.CNA$output$num.mark))
    
    ## Getting new ends from gcdata, thanks to this index
    seg.DYE.CNA$output$loc.end <- gcdata$End[seg.DYE.CNA$output$mark.end.i]
    seg.DYE.CNA$output$seg.med <- apply(seg.DYE.CNA$output, 1, function(x) { median(current$l2r[x[7]:x[8]])})
    
    cat("Creating DYE.seg.out\n")
    DYE.seg.out <- seg.DYE.CNA$output
    
    DYE.seg.out$genstart <- DYE.seg.out$loc.start + cs$lchrsum[DYE.seg.out$chrom]
    DYE.seg.out$genend <- DYE.seg.out$loc.end + cs$lchrsum[DYE.seg.out$chrom]
    TEMP.normal.index <- which( (DYE.seg.out$seg.med > -scut) & (DYE.seg.out$seg.med < scut) )
    
    ## SECOND CENTRALIZATION ON NORMAL SEGMENTS
    cat("Refining centralization ...\n")
    
    centZ <- median(DYE.seg.out$seg.med[TEMP.normal.index])
    current$l2r <- current$l2r - centZ
    current$runmed <- as.numeric(runmed(current$l2r, smo))
    DYE.seg.out$seg.med <- DYE.seg.out$seg.med - centZ
    
    DYE.gain.index <- which(DYE.seg.out$seg.med >= scut)
    DYE.loss.index <- which(DYE.seg.out$seg.med <= -scut)
    DYE.normal.index <- which( (DYE.seg.out$seg.med > -scut) & (DYE.seg.out$seg.med < scut) )
    
    ## SEGMENTED FINAL PANGENOMIC PLOT
    cat("Plotting segmented final genomic plot ...\n")
    l2r.profile.plot(current, gcdata, tit=paste(samplename, " final, centered", sep=""))
    l2r.seg.add(DYE.seg.out, gcdata, scut)
    ## DUPED AS A SINGLE PNG
    png(file=paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_nrf", nrf, "_", gv, "_pang.png", sep=""), width=1600, height=900)
    par(xaxs="i")
    l2r.profile.plot(current, gcdata, tit=paste(add.title, " ", samplename, " final, centered", sep=""))
    l2r.seg.add(DYE.seg.out, gcdata, scut)
    dev.off()
    
    ## KARYOTYPE-LIKE PLOT
    cat("Plotting karyotypic plot ...\n")
    l2r.karyo.plot(current, gcdata, DYE.seg.out, scut, sp)
    ## DUPED AS A SINGLE PNG
    png(file=paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_nrf", nrf, "_", gv, "_karyo.png", sep=""), width=1600, height=900)
    l2r.karyo.plot(current, gcdata, DYE.seg.out, scut, sp)
    dev.off()
  }
  
  ## CHROMOSOMAL PLOTS
  for (k in sort(unique(cs$cytob$chr))) {
    l2r.profile.chr.plot(current, gcdata, k)
    if (mod == "cgh") l2r.seg.chr.add(DYE.seg.out, gcdata, scut, k)
  }
  
  ## Closing the main graphics output
  dev.off()
  
  ## DUMPING SEGMENTATION RESULTS FOR CGH MODE
  if (mod == "cgh") {
    ## REFORMATING THE SEGMENTATION RESULTS
    DYE.seg.out$ID <- samplename
    DYE.seg.out <- DYE.seg.out[, -c(6:8,10,11)]
    colnames(DYE.seg.out) <- c(samplename, "Chr", "Start", "End", "Probes", "Log2Ratio")
    ## DUMPING TO NOCUT CBS
    write.table(DYE.seg.out, file=paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_NoCut", "_", gv, ".cbs", sep=""), sep="\t", quote=F, row.names=F)
    ## DUMPING TO NOCUT BEDGRAPH FOR IGV
    bedgraph.head <- paste("track type=bedGraph name=\"", samplename, "_GC5_", gv, " l2r\" description=\"", samplename, "_GC5_", gv, " log2(ratio)\" graphtype=bar color=255,0,0", sep="")
    write.table(bedgraph.head, paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_NoCut", "_", gv, ".bedgraph", sep=""), col.names=F, row.names=F, sep="\t", quote=F)
    bedgraph.df <- DYE.seg.out[,-c(1,5)]
#     for (i in 1:nrow(bedgraph.df)) { bedgraph.df[["Chr"]][i] <- chrconv("chr2chrom", sp, bedgraph.df[["Chr"]][i]) }
    bedgraph.df[["Chr"]] <- sapply(1:nrow(bedgraph.df), function(x) { cs$chr2chrom[[bedgraph.df$Chr[x]]] })
    write.table(bedgraph.df, paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_NoCut", "_", gv, ".bedgraph", sep=""), col.names=F, row.names=F, sep="\t", quote=F, append=T)
    if (nrow(DYE.seg.out) == 0) {
      write.table(DYE.seg.out, file=paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_nrf", nrf, "_", gv, ".cbs", sep=""), sep="\t", quote=F, row.names=F)
    } else {
      ## log2(ratio) values of Normal segment are replaced with zero.
      DYE.seg.out[["Log2Ratio"]][DYE.normal.index] <- 0
      ## DUMPING TO CBS
      write.table(DYE.seg.out, file=paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_nrf", nrf, "_", gv, ".cbs", sep=""), sep="\t", quote=F, row.names=F)
      ## DUMPING TO BEDGRAPH FOR IGV
      bedgraph.head <- paste("track type=bedGraph name=\"", samplename, "_GC5_", gv, " l2r\" description=\"", samplename, "_GC5_", gv, " log2(ratio)\" graphtype=bar color=255,0,0", sep="")
      write.table(bedgraph.head, paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_nrf", nrf, "_", gv, ".bedgraph", sep=""), col.names=F, row.names=F, sep="\t", quote=F)
      bedgraph.df <- DYE.seg.out[,-c(1,5)]
#       for (i in 1:nrow(bedgraph.df)) { bedgraph.df[["Chr"]][i] <- chrconv("chr2chrom", sp, bedgraph.df[["Chr"]][i]) }
      bedgraph.df[["Chr"]] <- sapply(1:nrow(bedgraph.df), function(x) { cs$chr2chrom[[bedgraph.df$Chr[x]]] })
      write.table(bedgraph.df, paste(tstamp, "/", samplename, "_GC5", suffix, "_uSD", uSD, "_a", alpha, "_nrf", nrf, "_", gv, ".bedgraph", sep=""), col.names=F, row.names=F, sep="\t", quote=F, append=T)
    }
  }
  
  
  
  #################
  ## Dumping GCX ##
  #################
  GCfex <- data.frame(ProbeName=gcdata$Identifier, Chr=gcdata$Chr, Start=gcdata$Start, End=gcdata$End)
  GCfex$LogRatio <- current$l2r
  write.table(GCfex, file=paste(tstamp, "/", samplename, "_GC5", suffix, "_", gv, ".gcx", sep=""), sep="\t", quote=F, row.names=F)
  system(paste("bzip2 -9 ", tstamp, "/", samplename, "_GC5", suffix, "_", gv, ".gcx", sep=""))
  
  #################
  ## Dumping WIG ##
  #################
  wighead <- paste(paste("track type=wiggle_0 name=\"", samplename, "_GC5_", gv, " l2r\" description=\"", samplename, "_GC5_", gv, " log2(ratio)\" graphtype=line", sep=""), paste("#genome_version=", gv, sep=""), paste("#chrom", "chromStart", "chromEnd", "Log2Ratio", sep="\t"), sep="\n")
  write.table(wighead, paste(tstamp, "/", samplename, "_GC5", suffix, "_", gv, ".wig", sep=""), col.names=F, row.names=F, sep="\t", quote=F)
  wigdf <-data.frame(chrom=gcdata$Chrom, chromStart=gcdata$Start, chromEnd=gcdata$End, l2r=current$l2r)
  write.table(wigdf, paste(tstamp, "/", samplename, "_GC5", suffix, "_", gv, ".wig", sep=""), col.names=F, row.names=F, sep="\t", quote=F, append=T)
  system(paste("bzip2 -9 ", tstamp, "/", samplename, "_GC5", suffix, "_", gv, ".wig", sep=""))

  ## ALSO Z-SCORE WIG IF CHIP MODE
  if (mod == "chip") {
    wighead <- paste(paste("track type=wiggle_0 name=\"", samplename, "_GC5_", gv, " Z\" description=\"", samplename, "_GC5_", gv, " Z-score\" graphtype=line", sep=""), paste("#genome_version=", gv, sep=""), paste("#chrom", "chromStart", "chromEnd", "Z-score", sep="\t"), sep="\n")
    write.table(wighead, paste(tstamp, "/", samplename, "_GC5", suffix, "_", gv, "_Z.wig", sep=""), col.names=F, row.names=F, sep="\t", quote=F)
    wigdf <-data.frame(chrom=gcdata$Chrom, chromStart=gcdata$Start, chromEnd=gcdata$End, l2r=(current$l2r-median(current$l2r))/sd(current$l2r))
    write.table(wigdf, paste(tstamp, "/", samplename, "_GC5", suffix, "_", gv, "_Z.wig", sep=""), col.names=F, row.names=F, sep="\t", quote=F, append=T)
    system(paste("bzip2 -9 ", tstamp, "/", samplename, "_GC5", suffix, "_", gv, "_Z.wig", sep=""))
  }
  #######################
  ## Dumping PQC files ##
  #######################
  
  ## Shortened version
  pqc.df.out <- data.frame(TERM=pqc.header.out, VALUE=pqc.value.out, check.names=F, stringsAsFactors=F)
  write.table(pqc.df.out, file=paste(tstamp, "/", samplename, ".pqc", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  
  ## Full version
  pqc.df.full <- data.frame(TERM=pqc.header.full, VALUE=pqc.value.full, check.names=F, stringsAsFactors=F)
  write.table(pqc.df.full, file=paste(tstamp, "/", samplename, "_full.pqc", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  
  ## IMAGE CONVERSION WHEN NEEDED
  if ((of == "tif") | (of == "tiff")){
    system(paste("gc5psconv ", outname, ".ps", sep=""))
    if (mysystem == "unix") system(paste("rm ", outname, ".ps", sep="")) else if (mysystem == "windows") system(paste("del ", outname, ".ps", sep=""))
  }
  return(pqc.value.out)
}
