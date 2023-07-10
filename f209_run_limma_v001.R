library(limma)
library(tidyverse)
library(vroom)
library(edgeR)
library(docopt)
library(janitor)
library(tools)

# Define the command-line interface
doc <- "
Usage:
 f209_run_limma_v001.R [--file=<file>] [--t=<threshold>] [--output-dir=<output_dir>]

Options:
  --file=<file>          Input file [default: ./DATA/merged_20230615_Ramya_freq. Live_singlet.txt]
  --t=<threshold>        Threshold [default: 0.05]
  --output-dir=<output_dir>    Output directory [default: ./OUT2/]
"

# Parse the command-line arguments
args <- docopt(doc)

file <- args$'--file'
threshold <- as.numeric(args$'--t')
output_dir <- args$'--output-dir'

file<-"DATA/merged_20230615_Ramya_freq. main population.txt"

# Read data
aa <- read_tsv(file) %>% filter(!is.na(REPOSID))
aaa <- aa %>% mutate(SampleID = REPOSID)
count <- aaa %>% select(42:139)
rownames(count) <- aaa$SampleID
count2 <- count %>%
  rownames_to_column() %>%
  gather(var, value, -rowname) %>%
  spread(rowname, value)

rownames(count2) <- count2$var
count3 <- count2[, -1] %>% mutate_if(is.character, as.numeric)

count3 <- log(abs(count3)) * sign(count3)
rownames(count3) <- count2$var

# Design matrix
a <- aa %>%
  as.tibble() %>%
  mutate(Group2 = SampleGroup) %>%
  mutate(Group2 = gsub("KS\\+MCD\\/KICS", "KS_MCD_KICS", Group2)) %>%
  select(REPOSID, `KS RESPONSE STATUS`, GROUP, Group2, plate, `CD4 baseline`, `Prior systemic KS therapy`, `HIV VL BL`, Age) %>%
  mutate(SampleID = REPOSID) %>%
  select(SampleID, `KS RESPONSE STATUS`, Group2, plate, `CD4 baseline`, `Prior systemic KS therapy`, `HIV VL BL`, Age) %>%
  arrange(SampleID)

a <- clean_names(a)


design <- model.matrix(~0 + group2 + plate, data = a)

# Perform LIMMA analysis
corfit <- duplicateCorrelation(count3, design)

cont.matrix <- makeContrasts(delta1 = (group2KS_MCD_KICS) - (group2KS), levels = design)

fit <- lmFit(count3, design = design, block = aaa$SampleID, correlation = corfit$consensus)
contrast_fit <- contrasts.fit(fit, contrasts = cont.matrix)
ebays_fit <- eBayes(contrast_fit, robust = TRUE)

baseline <- as.data.frame(topTable(ebays_fit, coef = 1, number = Inf))

p <- baseline %>% filter(adj.P.Val < threshold)

# Write results to file
output_file <- file.path(output_dir, paste0("outp_", tools::file_path_sans_ext(basename(file)), sep = ""))
write.table(p, output_file, sep = "\t", col.names = NA)
