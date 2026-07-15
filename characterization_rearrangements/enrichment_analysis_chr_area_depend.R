### Enrichment analysis ########################################################
library(regioneR)
library(readxl)
library(GenomicRanges)
library(rtracklayer)
library (dplyr)
library(purrr)
library(ggplot2)
library(Biostrings)


### gene, cds, gc proportion and overall repeat elements #######################
setwd("C:/Users/inesg/Desktop/tesis/helodes")

load("window_analyisis/enrichment/results_permtest_4JMC19C.RData")
load("window_analyisis/enrichment/results_permtest_1JMC18.RData")
load("C:/Users/inesg/Desktop/tesis/helodes/window_analyisis/enrichment/variables_permtest.RData")

breakpoints_1JMC18 <- toGRanges("window_analyisis/enrichment/breakpoints_1JMC18.txt")
breakpoints_4JMC19C <- toGRanges("window_analyisis/enrichment/breakpoints_4JMC19C.txt")
breakpoints_mid_1JMC18 <- resize(breakpoints_1JMC18, width=1, fix="center")
breakpoints_40kb_1JMC18 <- resize(breakpoints_1JMC18, width=40000, fix="center")
breakpoints_mid_4JMC19C <- resize(breakpoints_4JMC19C, width=1, fix="center")
breakpoints_40kb_4JMC19C <- resize(breakpoints_4JMC19C, width=40000, fix="center")

chr_length_1JMC18<-read.table("genome_assembly/ensamblaje_helodes_final/qc_final/1JMC18/1JMC18_contig_length.txt",header=F)
chr_length_4JMC19C<-read.table("genome_assembly/ensamblaje_helodes_final/qc_final/4JMC19C/4JMC19C_contig_length.txt",header=F)

fasta_1JMC18<-readDNAStringSet("genome_assembly/ensamblaje_helodes_final/clean_assemblies/1JMC18/1JMC18.fasta")
fasta_4JMC19C<-readDNAStringSet("repetitive_region_analysis/repeat_explorer/4JMC19C/4JMC19C_revcomp_mask.fasta")

custom_genome_1JMC18 <- data.frame(chr = names(fasta_1JMC18), 
                            size = width(fasta_1JMC18))

custom_genome_4JMC19C <- data.frame(chr = names(fasta_4JMC19C), 
                                  size = width(fasta_4JMC19C))


##Mask telomeric regions

mask_regions_1JMC18 <- GRanges()

mask_regions_4JMC19C <- GRanges()

telomere_size <- 10000

for(i in 1:nrow(custom_genome_1JMC18)) {
  chr <- custom_genome_1JMC18$chr[i]
  len <- custom_genome_1JMC18$size[i]
  
  mask_regions_1JMC18 <- c(mask_regions_1JMC18, 
                    GRanges(chr, IRanges(start=1, end=telomere_size)),
                    GRanges(chr, IRanges(start=len - telomere_size, end=len)))
}

for(i in 1:nrow(custom_genome_4JMC19C)) {
  chr <- custom_genome_4JMC19C$chr[i]
  len <- custom_genome_4JMC19C$size[i]
  
  mask_regions_4JMC19C <- c(mask_regions_4JMC19C, 
                           GRanges(chr, IRanges(start=1, end=telomere_size)),
                           GRanges(chr, IRanges(start=len - telomere_size, end=len)))
}

##Resample distinguishing subtelomeric and interstitial regions

prop <- 0.05 

subtel_regions_1JMC18 <- GRanges()
interst_regions_1JMC18 <- GRanges()


for(i in 1:nrow(custom_genome_1JMC18)) {
  scaff_name <- custom_genome_1JMC18$chr[i]
  scaff_len  <- custom_genome_1JMC18$size[i]
  limit      <- scaff_len * prop
  
  # Construir regiones subteloméricas (P y Q)
  subtel_regions_1JMC18 <- c(subtel_regions_1JMC18, 
                      GRanges(scaff_name, IRanges(start=1, end=limit)),
                      GRanges(scaff_name, IRanges(start=scaff_len - limit, end=scaff_len)))
  
  # Construir región intersticial (el centro del scaffold)
  interst_regions_1JMC18 <- c(interst_regions_1JMC18, 
                       GRanges(scaff_name, IRanges(start=limit + 1, end=scaff_len - limit - 1)))
}

subtel_regions_4JMC19C <- GRanges()
interst_regions_4JMC19C <- GRanges()

for(i in 1:nrow(custom_genome_4JMC19C)) {
  scaff_name <- custom_genome_4JMC19C$chr[i]
  scaff_len  <- custom_genome_4JMC19C$size[i]
  limit      <- scaff_len * prop
  
  # Construir regiones subteloméricas (P y Q)
  subtel_regions_4JMC19C <- c(subtel_regions_4JMC19C, 
                             GRanges(scaff_name, IRanges(start=1, end=limit)),
                             GRanges(scaff_name, IRanges(start=scaff_len - limit, end=scaff_len)))
  
  # Construir región intersticial (el centro del scaffold)
  interst_regions_4JMC19C <- c(interst_regions_4JMC19C, 
                              GRanges(scaff_name, IRanges(start=limit + 1, end=scaff_len - limit - 1)))
}

bp_sub_1JMC18 <- subsetByOverlaps(breakpoints_40kb_1JMC18, subtel_regions_1JMC18)

bp_int_1JMC18 <- subsetByOverlaps(breakpoints_40kb_1JMC18, interst_regions_1JMC18)

 bp_sub_4JMC19C <- subsetByOverlaps(breakpoints_40kb_4JMC19C, subtel_regions_4JMC19C)

bp_int_4JMC19C <- subsetByOverlaps(breakpoints_40kb_4JMC19C, interst_regions_4JMC19C)

perm_test_resuls_1JMC18z<-as.list

##run test per chromosome

features_1JMC18_gr<-list(gene=braker_1JMC18_gr_filt,
                         re=re_1JMC18_gr,
                         te=dante_1JMC18_gr,
                         ltr=ltr_1JMC18_gr,
                         ty1_copia=ltr_1JMC18_gr_ty1,
                         ty3_gypsy=ltr_1JMC18_gr_ty3,
                         dna_trans=dna_trans_1JMC18_gr,
                         line=dante_class1_nonltr_1JMC18_gr,
                         sine=subset(re_class_1JMC18,subclass=="SINE"),
                         unknown=subset(re_class_1JMC18,subclass=="Unknown"),
                         trc=trc_1JMC18_gr,
                         trc1=trc1_1JMC18_gr,
                         gc=gc_1JMC18_gr
)

permtest_perchr_1JMC18_subtel<-as.list(1:length(features_1JMC18_gr))
names(permtest_perchr_1JMC18_subtel)<-names(features_1JMC18_gr)
permtest_perchr_1JMC18_int<-as.list(1:length(features_1JMC18_gr))
names(permtest_perchr_1JMC18_int)<-names(features_1JMC18_gr)
permtest_perchr_4JMC19C_subtel<-as.list(1:length(features_4JMC19C_gr))
names(permtest_perchr_4JMC19C_subtel)<-names(features_4JMC19C_gr)
permtest_perchr_4JMC19C_int<-as.list(1:length(features_4JMC19C_gr))
names(permtest_perchr_4JMC19C_int)<-names(features_4JMC19C_gr)


permtest_perchr_1JMC18_subtel[1:12]<-lapply(features_1JMC18_gr[1:12],FUN=function(feat){
  regioneR::permTest(A = bp_sub_1JMC18, B=feat,
                     ntimes = 1000,
                     mask = mask_regions_1JMC18,
                     universe = subtel_regions_1JMC18, 
                     evaluate.function = numOverlaps,
                     randomize.function = randomizeRegions,
                     per_chromosome=T,
                     genome = chr_length_1JMC18) 
})

permtest_perchr_1JMC18_int[1:12]<-lapply(features_1JMC18_gr[1:12],FUN=function(feat){
  regioneR::permTest(A = bp_int_1JMC18, B=feat,
                     ntimes = 1000,
                     mask = mask_regions_1JMC18,
                     universe = interst_regions_1JMC18, 
                     evaluate.function = numOverlaps,
                     randomize.function = randomizeRegions,
                     per_chromosome=T,
                     genome = chr_length_1JMC18) 
})


permtest_perchr_1JMC18_subtel$gc<-permTest(A = resize(bp_sub_1JMC18, width=20000, fix="center"), 
                                    x =gc_1JMC18_gr, # Replace with your LINEs/GC object
                                    genome = custom_genome_1JMC18,
                                    mask = mask_regions_1JMC18,
                                    universe = subtel_regions_1JMC18, 
                                    randomize.function = randomizeRegions, 
                                    per.chromosome = T, # Corrects for scaffold length
                                    evaluate.function = meanInRegions, 
                                    ntimes = 1000)


permtest_perchr_1JMC18_int$gc<-permTest(A = resize(bp_int_1JMC18, width=20000, fix="center"), 
                                           x =gc_1JMC18_gr, # Replace with your LINEs/GC object
                                           genome = custom_genome_1JMC18,
                                           mask = mask_regions_1JMC18,
                                           universe =  interst_regions_1JMC18, 
                                           randomize.function = randomizeRegions, 
                                           per.chromosome = T, # Corrects for scaffold length
                                           evaluate.function = meanInRegions, 
                                           ntimes = 1000)


features_4JMC19C_gr<-list(gene=braker_4JMC19C_gr_filt,
                          re=re_4JMC19C_gr,
                          te=dante_4JMC19C_gr,
                          ltr=ltr_4JMC19C_gr,
                          ty1_copia=ltr_4JMC19C_gr_ty1,
                          ty3_gypsy=ltr_4JMC19C_gr_ty3,
                          dna_trans=dna_trans_4JMC19C_gr,
                          line=dante_class1_nonltr_4JMC19C_gr,
                          sine=subset(re_class_4JMC19C,subclass=="SINE"),
                          unknown=subset(re_class_4JMC19C,subclass=="Unknown"),
                          trc=trc_4JMC19C_gr,
                          trc1=trc1_4JMC19C_gr,
                          gc=gc_4JMC19C_gr
)



permtest_perchr_4JMC19C_subtel[1:12]<-lapply(features_4JMC19C_gr[1:12],FUN=function(feat){
  regioneR::permTest(A = bp_sub_4JMC19C, B=feat,
                     ntimes = 1000,
                     mask = mask_regions_4JMC19C,
                     universe = subtel_regions_4JMC19C, 
                     evaluate.function = numOverlaps,
                     randomize.function = randomizeRegions,
                     per_chromosome=T,
                     genome = chr_length_4JMC19C) 
})

permtest_perchr_4JMC19C_int[1:12]<-lapply(features_4JMC19C_gr[1:12],FUN=function(feat){
  regioneR::permTest(A = bp_int_4JMC19C, B=feat,
                     ntimes = 1000,
                     mask = mask_regions_4JMC19C,
                     universe = interst_regions_4JMC19C, 
                     evaluate.function = numOverlaps,
                     randomize.function = randomizeRegions,
                     per_chromosome=T,
                     genome = chr_length_4JMC19C) 
})


permtest_perchr_4JMC19C_subtel$gc<-permTest(A = bp_sub_4JMC19C, 
                                           x =gc_4JMC19C_gr, # Replace with your LINEs/GC object
                                           genome = custom_genome_4JMC19C,
                                           mask = mask_regions_4JMC19C,
                                           universe = subtel_regions_4JMC19C, 
                                           randomize.function = randomizeRegions, 
                                           per.chromosome = T, # Corrects for scaffold length
                                           evaluate.function = meanInRegions, 
                                           ntimes = 1000)


permtest_perchr_4JMC19C_int$gc<-permTest(A = bp_int_4JMC19C, 
                                        x =gc_1JMC18_gr, # Replace with your LINEs/GC object
                                        genome = custom_genome_4JMC19C,
                                        mask = mask_regions_4JMC19C,
                                        universe =  interst_regions_4JMC19C, 
                                        randomize.function = randomizeRegions, 
                                        per.chromosome = T, # Corrects for scaffold length
                                        evaluate.function = meanInRegions, 
                                        ntimes = 1000)


lz_score_1JMC18_int_gc<-localZScore(pt = res_inster_1JMC18_gr, A=resize(bp_int_1JMC18[1,], width=100000, fix="center"),
                                  window=100000,x =gc_1JMC18_gr)

res_inster_1JMC18_gene_gr <- permTest(A = bp_int_1JMC18, 
                                 B =features_1JMC18_gr$ty3_gypsy, # Replace with your LINEs/GC object
                                 genome = custom_genome_1JMC18,
                                 mask = mask_regions_1JMC18,
                                 universe = interst_regions_1JMC18, 
                                 randomize.function =randomizeRegions, 
                                 per.chromosome = F, # Corrects for scaffold length
                                 evaluate.function =numOverlaps, 
                                 ntimes = 100)

lz_score_1JMC18_int_gc<-localZScore(pt = res_inster_1JMC18_gr, A=resize(bp_int_1JMC18[1,], width=100000, fix="center"),
                                    window=100000,x =gc_1JMC18_gr)

res_subtel_1JMC18_gr <- permTest(A = bp_sub_1JMC18, 
                              x =gc_1JMC18_gr, # Replace with your LINEs/GC object
                              genome = custom_genome_1JMC18,
                              mask = mask_regions_1JMC18,
                              universe = subtel_regions_1JMC18, 
                              randomize.function = randomizeRegions, 
                              per.chromosome = TRUE, # Corrects for scaffold length
                              evaluate.function = meanInRegions, 
                              ntimes = 1000)

results_permtest_1JMC18<-lapply(features_1JMC18_gr,FUN=function(feat){
  regioneR::permTest(A = breakpoints_1JMC18, B=feat,
  ntimes = 1000,
  evaluate.function = numOverlaps,
  randomize.function = randomizeRegions,
  genome = chr_length_1JMC18) 
})




regioneR::permTest(A = breakpoints_1JMC18, B=trc1_1JMC18_gr,
                   ntimes = 1000,
                   evaluate.function = numOverlaps,
                   randomize.function = randomizeRegions,
                   genome = chr_length_1JMC18) 





features_4JMC19C_gr<-list(gene=braker_4JMC19C_gr_filt,
                         re=re_4JMC19C_gr,
                         te=dante_4JMC19C_gr,
                         ltr=ltr_4JMC19C_gr,
                         ty1_copia=ltr_4JMC19C_gr_ty1,
                         ty3_gypsy=ltr_4JMC19C_gr_ty3,
                         dna_trans=dna_trans_4JMC19C_gr,
                         line=dante_class1_nonltr_4JMC19C_gr,
                         sine=subset(re_class_4JMC19C,subclass=="SINE"),
                         unknown=subset(re_class_4JMC19C,subclass=="Unknown"),
                         trc=trc_4JMC19C_gr,
                         trc1=trc1_4JMC19C_gr
)


results_permtest_4JMC19C<-lapply(features_4JMC19C_gr,FUN=function(feat){
  regioneR::permTest(A = breakpoints_4JMC19C, B=feat,
                     ntimes = 1000,
                     evaluate.function = numOverlaps,
                     randomize.function = randomizeRegions,
                     genome = chr_length_4JMC19C) 
})


###For gc content is a bit different because it is quantitative

#create function for gc content
res_1JMC18_gc <- permTest(
  A = breakpoints_1JMC18,
  x = gc_1JMC18_gr,
  evaluate.function = meanInRegions,
  ntimes = 1000,
  randomize.function = randomizeRegions,
  genome = chr_length_1JMC18
)

res_4JMC19C_gc <- permTest(
  A = breakpoints_4JMC19C,
  x = gc_4JMC19C_gr,
  evaluate.function = meanInRegions,
  ntimes = 1000,
  randomize.function = randomizeRegions,
  genome = chr_length_4JMC19C
)

plot(res_gc)




fold_1JMC18_int<-rep(NA,length(permtest_perchr_1JMC18_int))
fold_1JMC18_subtel<-rep(NA,length(permtest_perchr_1JMC18_subtel))

fold_4JMC19C_int<-rep(NA,length(permtest_perchr_1JMC18_int))
fold_4JMC19C_subtel<-rep(NA,length(permtest_perchr_4JMC19C_subtel))

zscore_1JMC18_int<-rep(NA,length(permtest_perchr_1JMC18_int))
zscore_1JMC18_subtel<-rep(NA,length(permtest_perchr_1JMC18_subtel))

zscore_4JMC19C_int<-rep(NA,length(permtest_perchr_1JMC18_int))
zscore_4JMC19C_subtel<-rep(NA,length(permtest_perchr_4JMC19C_subtel))

pvalue_1JMC18_int<-rep(NA,length(permtest_perchr_1JMC18_int))
pvalue_1JMC18_subtel<-rep(NA,length(permtest_perchr_1JMC18_subtel))

pvalue_4JMC19C_int<-rep(NA,length(permtest_perchr_1JMC18_int))
pvalue_4JMC19C_subtel<-rep(NA,length(permtest_perchr_4JMC19C_subtel))

observed_1JMC18_int<-rep(NA,length(permtest_perchr_1JMC18_int))
observed_1JMC18_subtel<-rep(NA,length(permtest_perchr_1JMC18_subtel))

observed_4JMC19C_int<-rep(NA,length(permtest_perchr_1JMC18_int))
observed_4JMC19C_subtel<-rep(NA,length(permtest_perchr_4JMC19C_subtel))

ci_low_1JMC18_int<-rep(NA,13)
ci_low_1JMC18_subtel<-rep(NA,13)
ci_high_1JMC18_int<-rep(NA,13)
ci_high_1JMC18_subtel<-rep(NA,13)

ci_low_4JMC19C_int<-rep(NA,13)
ci_low_4JMC19C_subtel<-rep(NA,13)
ci_high_4JMC19C_int<-rep(NA,13)
ci_high_4JMC19C_subtel<-rep(NA,13)


for (i in 1:12){
  zscore_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$numOverlaps$zscore
  zscore_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$zscore
  
  zscore_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$numOverlaps$zscore
  zscore_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$zscore
  
  pvalue_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$numOverlaps$pval
  pvalue_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$pval
  
  pvalue_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$numOverlaps$pval
  pvalue_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$pval
  
  observed_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$numOverlaps$observed
  observed_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$observed
  
  observed_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$numOverlaps$observed
  observed_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$observed
  
  fold_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$numOverlaps$observed/median(permtest_perchr_1JMC18_int[[i]]$numOverlaps$permuted)
  fold_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$observed/median(permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$permuted)
  
  fold_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$numOverlaps$observed/median(permtest_perchr_4JMC19C_int[[i]]$numOverlaps$permuted)
  fold_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$observed/median(permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$permuted)
  
  ci_low_1JMC18_int[i]<-quantile(permtest_perchr_1JMC18_int[[i]]$numOverlaps$permuted, 0.025)
  ci_high_1JMC18_int[i]<-quantile(permtest_perchr_1JMC18_int[[i]]$numOverlaps$permuted, 0.975)
  ci_low_4JMC19C_int[i]<-quantile(permtest_perchr_4JMC19C_int[[i]]$numOverlaps$permuted, 0.025)
  ci_high_4JMC19C_int[i]<-quantile(permtest_perchr_4JMC19C_int[[i]]$numOverlaps$permuted, 0.975)
  
  ci_low_1JMC18_subtel[i]<-quantile(permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$permuted, 0.025)
  ci_high_1JMC18_subtel[i]<-quantile(permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$permuted, 0.975)
  ci_low_4JMC19C_subtel[i]<-quantile(permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$permuted, 0.025)
  ci_high_4JMC19C_subtel[i]<-quantile(permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$permuted, 0.975)
}

i<-13

zscore_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$meanInRegions$zscore
zscore_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$meanInRegions$zscore

zscore_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$meanInRegions$zscore
zscore_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$meanInRegions$zscore

pvalue_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$meanInRegions$pval
pvalue_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$meanInRegions$pval

pvalue_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$meanInRegions$pval
pvalue_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$meanInRegions$pval

observed_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$meanInRegions$observed
observed_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$meanInRegions$observed

observed_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$meanInRegions$observed
observed_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$meanInRegions$observed

fold_1JMC18_int[i]<-permtest_perchr_1JMC18_int[[i]]$meanInRegions$observed/median(permtest_perchr_1JMC18_int[[i]]$meanInRegions$permuted)
fold_1JMC18_subtel[i]<-permtest_perchr_1JMC18_subtel[[i]]$meanInRegions$observed/median(permtest_perchr_1JMC18_subtel[[i]]$meanInRegions$permuted)

fold_4JMC19C_int[i]<-permtest_perchr_4JMC19C_int[[i]]$meanInRegions$observed/median(permtest_perchr_4JMC19C_int[[i]]$meanInRegions$permuted)
fold_4JMC19C_subtel[i]<-permtest_perchr_4JMC19C_subtel[[i]]$meanInRegions$observed/median(permtest_perchr_4JMC19C_subtel[[i]]$meanInRegions$permuted)

ci_low_1JMC18_int[i]<-quantile(permtest_perchr_1JMC18_int[[i]]$meanInRegions$permuted, 0.025)
ci_high_1JMC18_int[i]<-quantile(permtest_perchr_1JMC18_int[[i]]$meanInRegions$permuted, 0.975)
ci_low_4JMC19C_int[i]<-quantile(permtest_perchr_4JMC19C_int[[i]]$meanInRegions$permuted, 0.025)
ci_high_4JMC19C_int[i]<-quantile(permtest_perchr_4JMC19C_int[[i]]$meanInRegions$permuted, 0.975)

ci_low_1JMC18_subtel[i]<-quantile(permtest_perchr_1JMC18_subtel[[i]]$meanInRegions$permuted, 0.025)
ci_high_1JMC18_subtel[i]<-quantile(permtest_perchr_1JMC18_subtel[[i]]$meanInRegions$permuted, 0.975)
ci_low_4JMC19C_subtel[i]<-quantile(permtest_perchr_4JMC19C_subtel[[i]]$meanInRegions$permuted, 0.025)
ci_high_4JMC19C_subtel[i]<-quantile(permtest_perchr_4JMC19C_subtel[[i]]$meanInRegions$permuted, 0.975)

#Define zero inflated distributions: 20 % 0s


prop_zeros_1JMC18_int<-rep(NA,12)
prop_zeros_1JMC18_subtel<-rep(NA,12)
prop_zeros_4JMC19C_int<-rep(NA,12)
prop_zeros_4JMC19C_subtel<-rep(NA,12)

for ( i in 1:12){
  prop_zeros_1JMC18_int[i]<-length(which(permtest_perchr_1JMC18_int[[i]]$numOverlaps$permuted==0))/1000
  prop_zeros_1JMC18_subtel[i]<-length(which(permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$permuted==0))/1000
  prop_zeros_4JMC19C_int[i]<-length(which(permtest_perchr_4JMC19C_int[[i]]$numOverlaps$permuted==0))/1000
  prop_zeros_4JMC19C_subtel[i]<-length(which(permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$permuted==0))/1000
  }
#zero inflated variables are the same between the two genomes but not the two zones

for (i in which(prop_zeros_1JMC18_int>0.2)){
  fold_1JMC18_int[i]<-(1+permtest_perchr_1JMC18_int[[i]]$numOverlaps$observed)/(1+mean(permtest_perchr_1JMC18_int[[i]]$numOverlaps$permuted))
  fold_4JMC19C_int[i]<-(1+permtest_perchr_4JMC19C_int[[i]]$numOverlaps$observed)/(1+median(permtest_perchr_4JMC19C_int[[i]]$numOverlaps$permuted))
}
for (i in which(prop_zeros_1JMC18_subtel>0.2)){
  fold_1JMC18_subtel[i]<-(1+permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$observed)/(1+mean(permtest_perchr_1JMC18_subtel[[i]]$numOverlaps$permuted))
  fold_4JMC19C_subtel[i]<-(1+permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$observed)/(1+median(permtest_perchr_4JMC19C_subtel[[i]]$numOverlaps$permuted))
}
