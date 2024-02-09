library(tidyverse)
library(edgeR)
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

b <- read_tsv(args$meta4batching)
d <- b %>% arrange(REPOSID) %>% select(REPOSID, Group)

counts <- read_tsv(args$matrix_counts)

counts <- as.tibble(counts)
counts2 <- counts %>% select(-1)

rownames(counts2) <- make.names(counts$ID, unique = TRUE)

counts3 <- counts2 %>% mutate_if(is.character, as.numeric)
rownames(counts3) <- make.names(counts$ID, unique = TRUE)

isexpr = filterByExpr(counts3)
dgelist <- calcNormFactors(DGEList(as.matrix(counts3)[isexpr,]), method = "TMM")

tt <- cbind(filterByExpr(counts3), rownames(counts3))
rownames4filter <- as.tibble(tt) %>% filter(V1 == 'TRUE')

log2_cpms <- cpm(dgelist, log = TRUE)
rna <- as.tibble(log2_cpms)

rownames(rna) <- rownames4filter$V2

png("box.png")
boxplot(rna)
dev.off()

mds <- cmdscale(dist(t(rna)))

design <- model.matrix(~ 0 + Group, data = d)

cont.matrix = makeContrasts(
  con1 = GroupKSplus - GroupKS,
  levels = design
)

v <- voom(dgelist, design, plot = TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = cont.matrix)
efit <- eBayes(vfit)


png("mod1.png")
voom(dgelist, design, plot = TRUE)
dev.off()


png("mod2.png")
plotSA(efit, main = "Final model: Mean-variance trend")
dev.off()

png("mds.png")
plotMDS(log2_cpms, labels = d$Group)
dev.off()

  write.table(summary(decideTests(efit)), "summary.txt")


ctpair <- as.data.frame(topTable(efit, coef = 1, number = Inf))

write.table(ctpair, "DEG.txt")


##############################
library(tidyverse)
library(fgsea)
tm<- gmtPathways("./BTM_for_GSEA_20131008.gmt")

ch<-ctpair%>%  arrange(t)
res2<-cbind(rownames(ch), as.numeric(as.character(ch$t)))



res2<-as.tibble(res2)
res2$V2 <- as.numeric(res2$V2)

ranks<-deframe(as.data.frame(res2))

fgseaRes <-fgsea(pathways=tm, stats=ranks, nperm=10000, maxSize = 500) 

Res<-fgseaRes   %>% arrange(padj)
library(data.table)
fwrite(Res, "PATH.txt")
