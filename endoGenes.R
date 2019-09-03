#!/usr/bin/env Rscript

### Compute relative expression values for realtime quantitative RT-PCRdata

## Compute relative expression values for realtime quantitative RT-PCR data based on Ct or take-off values, respectively. The computations use the PCR efficiency.


#Install and Load Multiple R Packages
if (!require("SLqPCR", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = 'https://cran-r.c3sl.ufpr.br/')
  }
  BiocManager::install("SLqPCR")
  library(SLqPCR)
}


Install_Multiples_Packages <- function(packages) {
  pack <- packages[!(packages %in% installed.packages()[,'Package'])];
  if (length(pack)) {
    install.packages(pack, repos = 'https://cran.rstudio.com/')
  }

  for (package_i in packages) {
    suppressPackageStartupMessages(library(package_i, character.only = TRUE, quietly = TRUE))
    }

}



Install_Multiples_Packages(c('RColorBrewer', 'RankAggreg', 'gplots', 'ctrlGene', 'optparse'))

# suppressPackageStartupMessages(library('optparse'))
# suppressPackageStartupMessages(library('SLqPCR'))
# suppressPackageStartupMessages(library('RColorBrewer'))
# suppressPackageStartupMessages(library("RankAggreg"))
# suppressPackageStartupMessages(library("gplots"))
# suppressPackageStartupMessages(library("ctrlGene"))



option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = "endogenous_ct.txt",
              help = "The filename of the dataset file [default %default]",
              dest = "samplesFile"),
  make_option(c("-e", "--efficiency"), type = "character", default = "efficiencies_list.txt",
              help = "The filename of the gene efficience file [default %default]",
              dest = "efficiencyList"),
  make_option(c("-o", "--output"), type="character", default="01-results",
              help = "output folder [default %default]",
              dest = "output"),
  make_option(c("-m", "--method"), type = "character", default = "CE",
              help = "method to be used to perform rank aggregation: Cross Entropy Monte Carlo (CE) or Genetic Algorithm (GA) in Rankaggreg analysis [default %default]",
              dest = "method"),
  make_option(c("-i", "--iteraction"), type = "numeric", default = 1000,
              help = "the maximum number of iterations allowed in Rankaggreg analysis [default %default]",
              dest = "iteraction"),
  make_option(c("-d", "--distance"), type = "character", default = "Spearman",
              help = "distance which 'measures' the similarity between the ordered lists ('Spearman' or 'Kendall') in Rankaggreg analysis [default %default]",
              dest = "distance"),
  make_option(c("-g", "--group"), type = "logical", default = "TRUE",
              help = "Analysis of Normfinder should be performed in groups? [default %default]",
              dest = "group"),
  make_option(c("-n", "--normfinderFunction"), type = "character", default = "RnormfinderFunction.txt",
              help = "Normfinder Function file [default %default]",
              dest = "normfinderFunction"),
  make_option(c("-c", "--ctvalue"), type = "logical", default = "TRUE",
              help = "Analysis of Normfinder should be performed using raw ct values? [default %default]",
              dest = "ctVal")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, description =  paste('Authors: OLIVEIRA, H.C., sep = "\n", Version: 0.0.3', sep = "\n")))

## Preparing dataset
final_folder <- opt$output
if(!file.exists(file.path(final_folder))) dir.create(file.path(final_folder), recursive = TRUE, showWarnings = FALSE)


if (!file.exists(opt$samplesFile)) {
  write(paste("Sample file", opt$samplesFile, "does not exist\n"), stderr())
  stop()
}

if (!file.exists(opt$efficiencyList)) {
  write(paste("Sample file", opt$efficiencyList, "does not exist\n"), stderr())
  stop()
}

ctvalue <- read.table(opt$samplesFile, header = TRUE, dec = ',', row.names = 1); str(ctvalue); summary(ctvalue); head(ctvalue)


ctvalue_new <- ctvalue[,-length(ctvalue)] #
eficiencia <- read.table(opt$efficiencyList, dec = ',', header = FALSE, col.names = 'EFFICIENCY')  # Efficiency values
rel_values <- function(dados, efi_list){
  res <- list()
  for(i in 1:dim(dados)[2]){
    res[[i]] <- relQuantPCR(dados[,i], E = efi_list[i], na.rm = FALSE)

  }
  res_matrix <- do.call('cbind',res);
}

qvalue <- rel_values(ctvalue_new, eficiencia$EFFICIENCY); nomes_col <- colnames(ctvalue_new[,1:ncol(ctvalue_new)]); colnames(qvalue) <- nomes_col; head(qvalue)

write.csv(qvalue, paste(opt$output, '/', 'results_qvalues.csv', sep = ''), quote = FALSE, row.names = FALSE)

#---------------------------------------REFERENCE GENES----------------------------------


#-------------------------------------------GeNorm---------------------------------------
#Using SLqPCR package
## Selection of reference/housekeeping genes
paste("GeNorm analysis - Using SLqPCR package")

res.gene <- selectHKgenes(qvalue, method = "Vandesompele", minNrHK = 2, geneSymbol = names(ctvalue_new), trace = TRUE, na.rm = FALSE)


## Ranking Stability Gene

rankingeral <- data.frame(c(1, 1, 3:ncol(qvalue)), res.gene$ranking, MvalorG = sort(c(res.gene$meanM[length(res.gene$meanM)], res.gene$meanM))); names(rankingeral) <- c("rank", 'Gene', 'MValue'); rankingeral

write.csv(rankingeral, paste(opt$output, '/', 'general_ranking_genorm.csv', sep = ''), quote = FALSE, row.names = FALSE)


## Plot stability Genes

png(paste(opt$output, '/', 'Rplot_gene_stability_by_genorm.png', sep = ''),
    width = 1280,
    height = 720,
    pointsize = 15)
mypalette <- brewer.pal(3, "Set1")
    matplot(cbind(rankingeral$MValue), type = "b", ylab = "Average expression stability M", xlab = "<==== Most Stable Gene   Least Stable Gene ====>", axes = FALSE, pch = 19, col = mypalette, ylim = c(0, max(rankingeral$MValue)), lty = 1, lwd = 2, main = "Gene stability measure by SLqPCR Package")
    axis(1, at = 1:nrow(rankingeral), labels = as.character(rankingeral$Gene))
    axis(2, at = min(rankingeral$MValue):max(rankingeral$MValue), labels = as.character())
    box()
    abline(h = seq(0, max(rankingeral$MValue), by = 0.2), lty = 2, lwd = 1, col = "grey")
dev.off()


### Plot PairWise variation
###

png(paste(opt$output, '/', 'Rplot_gene_variation_by_genorm.png', sep = ''),
    width = 1280,
    height = 720,
    pointsize = 15)
mypalette <- brewer.pal(8, "Spectral")
    barplot(cbind(res.gene$variation), beside = TRUE, col = mypalette, space = c(0, 0.5), ylim = c(0, 0.25), xlab = c("PairWise Variations"))
    axis(1, at = 1:length(names(res.gene$variation)), labels = as.character(names(res.gene$variation)))
    abline(h = seq(0, max(rankingeral$MValue), by = 0.2), lty = 2, col = "grey")
    abline(h = 0.15, lty = 1, col = "black")
dev.off()





### ---------------------------------------NORMFINDER ---------------------------------------
#Using r.NormOldStab5.txt function file
paste("NormFinder analysis - Using r.NormOldStab5.txt function")
source(opt$normfinderFunction)

#### Analysis using ctvalue
#
#

write.table(t(ctvalue), paste(opt$output, '/', 'normdata_trans_normfinder.txt', sep = ''), quote = FALSE); normdata = read.table(paste(opt$output, '/', 'normdata_trans_normfinder.txt', sep = ''), header = TRUE); head(normdata); tail(normdata)

Resulttotal <- Normfinder(paste(opt$output, '/', 'normdata_trans_normfinder.txt', sep = ''), Groups = opt$group, ctVal = opt$ctVal); Resulttotal$Ordered

write.csv(Resulttotal$Ordered, paste(opt$output, '/', 'ranking_Ordered_normfinder.csv', sep = ''), quote = FALSE, row.names = TRUE)

write.csv(Resulttotal$UnOrdered, paste(opt$output, '/', 'ranking_UnOrdered_normfinder.csv', sep = ''), quote = FALSE, row.names = TRUE)

write.csv(Resulttotal$PairOfGenes, paste(opt$output, '/', 'PairOfGenes_normfinder.csv', sep = ''), quote = FALSE, row.names = TRUE)

png(paste(opt$output, '/', 'Rplot_gene_stability_by_NormFinder.png', sep = ''),
    width = 1280,
    height = 720,
    pointsize = 15)
mypalette <- brewer.pal(3, "Set2")
    matplot(cbind(Resulttotal$Ordered$Stability), type = "b", ylab = "Average expression stability M", xlab = "<==== Most Stable Gene   Least Stable Gene ====>", axes = FALSE, pch = 19, col = mypalette, ylim = c(0, max(Resulttotal$Ordered)), lty = 1, lwd = 2, main = "Gene Stability Measure by NormFinder")
    axis(1, at = 1:nrow(Resulttotal$Ordered), labels = as.character(rownames(Resulttotal$Ordered)))
    axis(2, at = min(Resulttotal$Ordered):max(Resulttotal$Ordered$Stability), labels = as.character())
    box()
    #abline(h = seq(0.2, 1.0, by = 0.2), lty = 1, lwd = 1, col = "grey")
dev.off()

# png(paste(opt$output, '/', 'Rplot_geneGroupDif_GroupSD_stability__by_NormFinder.png', sep = ''),
#     width = 1280,
#     height = 720,
#     pointsize = 15)
# mypalette <- brewer.pal(3, "Set2")
# matplot(cbind(Resulttotal$Ordered), type = "b", ylab = "Average expression stability M", xlab = "<==== Most Stable Gene   Least Stable Gene ====>", axes = FALSE, pch = 19, col = mypalette, ylim = c(0, max(Resulttotal$Ordered)), lty = 1, lwd = 2, main = "Gene Stability Measure by NormFinder")
# axis(1, at = 1:nrow(Resulttotal$Ordered), labels = as.character(rownames(Resulttotal$Ordered)))
# axis(2, at = min(Resulttotal$Ordered):max(Resulttotal$Ordered), labels = as.character(colnames(Resulttotal$Ordered)))
#
# box()
# #abline(h = seq(0.2, 1.0, by = 0.2), lty = 1, lwd = 1, col = "grey")
# dev.off()


#---------------------------------------BESTKEEPER---------------------------------------
#Using ctrlGene package
paste("BestKeeper analysis - Using ctrlGene package")
bestkeeper_results = bestKeeper(ctvalue[,1:ncol(ctvalue) - 1], ctVal = TRUE)
bestkeeper_results

write.csv(bestkeeper_results$pair.Wise.cor, paste(opt$output, '/', 'Bestkeeper_results_pair_wise_correlation.csv', sep = ''), quote = FALSE)

write.csv(bestkeeper_results$HKG.vs.BestKeeper, paste(opt$output, '/', 'Bestkeeper_results_HKG_vs_BestKeeper.csv', sep = ''), quote = FALSE)

write.csv(bestkeeper_results$CP.statistics, paste(opt$output, '/', 'Bestkeeper_results_CP.statistics.csv', sep = ''), quote = FALSE)


bestkeeper_genes <- t(sort(bestkeeper_results$CP.statistics[6,])); best_genes <- t(bestkeeper_genes); colnames(best_genes) <- paste('SD_Value')

write.csv(best_genes, paste(opt$output, '/', 'Bestkeeper_best_genes_ordered.csv', sep = ''), quote = FALSE)

png(paste(opt$output, '/', 'Rplot_gene_stability_by_BestKeeper.png', sep = ''),
    width = 1280,
    height = 720,
    pointsize = 15)
mypalette <- brewer.pal(3, "Set2")
matplot(cbind(best_genes), type = "b", ylab = "Average expression stability M", xlab = "<==== Most Stable Gene   Least Stable Gene ====>", axes = FALSE, pch = 19, col = mypalette, ylim = c(0, max(best_genes)), lty = 1, lwd = 2, main = "Gene Stability Measure by NormFinder")
axis(1, at = 1:nrow(best_genes), labels = as.character(rownames(best_genes)))
axis(2, at = min(best_genes):max(best_genes), labels = as.character())
box()
dev.off()
#abline(h = seq(0.2, 1.0, by = 0.2), lty = 1, lwd = 1, col = "grey")



#---------------------------------------RANKAGGREG---------------------------------------
paste("Final ranking - using RankAggreg package")

#
#

ratools <- t(data.frame(NormFinder = rownames(Resulttotal$Ordered), GeNorm = rankingeral$Gene, Bestkeeper = rownames(best_genes))); colnames(ratools) <- paste(1:ncol(ratools)); write.csv(ratools, paste(opt$output, '/', 'ratools_gene_list.csv', sep = ''), quote = FALSE)


ra_weight_ratools <- t(data.frame(NormFinder = Resulttotal$Ordered$Stability, GeNorm = rankingeral$MValue, Bestkeeper = best_genes[,1])); colnames(ra_weight_ratools) <- paste(1:ncol(ratools)); write.csv(ra_weight_ratools, paste(opt$output, '/', 'ratools_gene_weight.csv', sep = ''), quote = FALSE)


k <- ncol(ratools)

if (k > 10) {
  write(paste('Using RankAggreg function - via the Cross-Entropy Monte Carlo algorithm or the Genetic Algorithm.'), stderr())
    (CESP <- RankAggreg(ratools, k = ncol(ratools), ra_weight_ratools, method = opt$method, distance = opt$distance, weight = .25, rho = .1, maxIter = opt$iteraction, verbose = TRUE))
  pdf(paste(opt$output, '/', "RankAggreg_Iteraction_plot.pdf", sep = ''))
  plot(CESP)
  dev.off()
}else {
  write(paste('Using BruteAggreg function - via the brute force approach.'), stderr())
    (CESP <- BruteAggreg(ratools, k = ncol(ratools), ra_weight_ratools, distance=opt$distance))
  pdf(paste(opt$output, '/', "RankAggreg_Iteraction_plot.pdf", sep = ''))
  plot(CESP)
  dev.off()
}

final_res <- data.frame(Rank = paste(1:length(CESP$top.list)), Gene = CESP$top.list)

write.csv(final_res, paste(opt$output, '/', 'Final_ranking.csv', sep = ''), quote = FALSE, row.names = FALSE)



paste("It's done. See ", opt$output, "folder")
