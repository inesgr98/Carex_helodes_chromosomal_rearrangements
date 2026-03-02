### Enrichment analysis ########################################################

setwd("C:/Users/rogel/Desktop/tesis/Tesis/genome_features")

install.packages("BiocManager")
BiocManager::install("regioneR")
library(regioneR)
library(GenomicRanges)
library (dplyr)
library(purrr)
library(ggplot2)

### gene, cds, gc proportion and overall repeat elements #######################

genes <- toGRanges("genes_DTL.txt")
cds <- toGRanges("cds_DTL.txt")
repeats <- toGRanges("repeats_DTL.txt")
df_gc <- read.table("gc_content_DTL.txt", header = FALSE)
colnames(df_gc) <- c("chr", "start", "end", "gc")
gc_gr <- GRanges(df_gc$chr, IRanges(df_gc$start + 1, df_gc$end), gc = df_gc$gc)
breakpoints <- toGRanges("breakpoints_DTL.txt")
genome_df <- read.table("chrm_sizes_DTL.txt", header = FALSE)

head(genes)
head(cds)
head(breakpoints)
head(repeats)
head(genome_df)
head(gc_gr)

res_genes <- permTest(
  A = breakpoints,
  B = genes,
  ntimes = 1000,
  evaluate.function = numOverlaps,
  randomize.function = randomizeRegions,
  genome = genome_df
)

plot(res_genes)

res_cds <- permTest(
  A = breakpoints,
  B = cds,
  ntimes = 1000,
  evaluate.function = numOverlaps,
  randomize.function = randomizeRegions,
  genome = genome_df
)

plot(res_cds)

res_repeats <- permTest(
  A = breakpoints,
  B = repeats,
  ntimes = 1000,
  evaluate.function = numOverlaps,
  randomize.function = randomizeRegions,
  genome = genome_df
)

plot(res_repeats)

res_gc <- permTest(
  A = breakpoints,
  x = gc_gr,
  evaluate.function = meanInRegions,
  ntimes = 1000,
  randomize.function = randomizeRegions,
  genome = genome_df
)

plot(res_gc)

### Repeat family analysis ####################################################

rep_families <- read.table("repeats_family_DTL.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(rep_families) <- c("chr", "start", "end", "family")

rep_families_gr <- GRanges(
  seqnames = rep_families$chr,
  ranges = IRanges(start = rep_families$start + 1, end = rep_families$end),
  family = rep_families$family
)
unique(rep_families$family)

rep_by_family <- split(rep_families_gr, rep_families_gr$family)

results_repeats_by_family <- lapply(names(rep_by_family), function(fam) {
  message("Procesando familia: ", fam)
  res <- permTest(
    A = breakpoints,
    B = rep_by_family[[fam]],
    ntimes = 1000,
    evaluate.function = numOverlaps,
    randomize.function = randomizeRegions,
    genome = genome_df
  )
  return(list(family = fam, result = res))
})

summary_df <- do.call(rbind, lapply(results_repeats_by_family, function(x) {
  data.frame(
    family = x$family,
    pval = x$result$numOverlaps$pval,
    zscore = x$result$numOverlaps$zscore
  )
}))

summary_df <- summary_df[order(summary_df$pval), ]
print(summary_df)

res_LTR_Gypsy <- results_repeats_by_family[[which(sapply(results_repeats_by_family, 
                                                   function(x) x$family == "LTR/Gypsy") )]]$result
res_DNA_hAT_Tag1 <- results_repeats_by_family[[which(sapply(results_repeats_by_family, 
                                                         function(x) x$family == "DNA/hAT-Tag1") )]]$result
res_DNA_hAT_Ac <- results_repeats_by_family[[which(sapply(results_repeats_by_family, 
                                                             function(x) x$family == "DNA/hAT-Ac") )]]$result

plot(res_LTR_Gypsy)
plot(res_DNA_hAT_Tag1)
plot(res_DNA_hAT_Ac)