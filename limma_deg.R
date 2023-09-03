library(tidyverse)
library(edgeR)
source("./functions_bulk.R")
library(docopt)
#' @docopt
#' Usage:
#'   script.R --meta4batching=<meta4batching> --matrix_counts=<matrix_counts>
#'
#' Options:
#'   --meta4batching=<meta4batching>    Path to meta4batching file.
#'   --matrix_counts=<matrix_counts>    Path to matrix counts file.
#'
#' This script performs differential gene expression analysis using edgeR package.
#' It takes two input files: meta4batching and matrix counts. The output is a table of differentially expressed genes.
#'
#' @param meta4batching Path to meta4batching file.
#' @param matrix_counts Path to matrix counts file.
#'
#' @return A table of differentially expressed genes.
#'
#' @examples
#' # Run the script
#' Rscript script.R --meta4batching=meta4batching.txt --matrix_counts=matrix_counts.txt
#'
#' # View the output
#' head output.txt
#'


doc <- "
Usage:
  script.R --meta4batching=<meta4batching> --matrix_counts=<matrix_counts>
  
Options:
  --meta4batching=<meta4batching>    Path to meta4batching file.
  --matrix_counts=<matrix_counts>    Path to matrix counts file.
"

args <- docopt(doc)


meta4batching_file_path <- args$meta4batching
matrix_counts_file_path <- args$matrix_counts

b <- read_file(meta4batching_file_path)
d <- b %>% arrange(REPOSID) %>% select(REPOSID, Group)

counts3 <- read_file_remove_first_col(matrix_counts_file_path)

isexpr = filter_by_expr(counts3)
dgelist <- calc_norm_factors(DGEList(as.matrix(counts3)[isexpr,]))

tt <- cbind(filter_by_expr(counts3), rownames(counts3))
rownames4filter <- as.tibble(tt) %>% filter(V1 == 'TRUE')

log2_cpms <- cpm(dgelist, log = TRUE)
rna <- as.tibble(log2_cpms)

rownames(rna) <- rownames4filter$V2



mds <- cmdscale(dist(t(rna)))



design <- create_design_matrix(d)


cont.matrix <- makeContrasts(con1 = GroupKSplus - GroupKS, levels = design)
efit <- create_efit(dgelist, design, cont.matrix)

png("mod2.png")
plotSA(efit, main = "Final model: Mean-variance trend")
dev.off()

png("mds.png")
plotMDS(log2_cpms, labels = d$Group)
dev.off()

write.table(summary(decideTests(efit)), "summary.txt")


ctpair <- as.data.frame(topTable(efit, coef = 1, number = Inf))

write.table(ctpair, "DEG.txt")




