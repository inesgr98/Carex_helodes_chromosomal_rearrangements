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

permtest_perchr_1JMC18<-as.list(1:length(features_1JMC18_gr))
names(permtest_perchr_1JMC18)<-names(features_1JMC18_gr)

permtest_perchr_4JMC19C<-as.list(1:length(features_4JMC19C_gr))
names(permtest_perchr_4JMC19C)<-names(features_4JMC19C_gr)


permtest_perchr_1JMC18[1:12]<-lapply(features_1JMC18_gr[1:12],FUN=function(feat){
  regioneR::permTest(A = breakpoints_40kb_1JMC18, B=feat,
                     ntimes = 1000,
                     mask = mask_regions_1JMC18,
                     evaluate.function = numOverlaps,
                     randomize.function = randomizeRegions,
                     per_chromosome=T,
                     genome = chr_length_1JMC18) 
})


permtest_perchr_1JMC18$gc<-permTest(A = breakpoints_40kb_1JMC18, 
                                    x =gc_1JMC18_gr, # Replace with your LINEs/GC object
                                    genome = custom_genome_1JMC18,
                                    mask = mask_regions_1JMC18,
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



permtest_perchr_4JMC19C[1:12]<-lapply(features_4JMC19C_gr[1:12],FUN=function(feat){
  regioneR::permTest(A =  breakpoints_40kb_4JMC19C, B=feat,
                     ntimes = 1000,
                     mask = mask_regions_4JMC19C,
                     evaluate.function = numOverlaps,
                     randomize.function = randomizeRegions,
                     per_chromosome=T,
                     genome = chr_length_4JMC19C) 
})



permtest_perchr_4JMC19C$gc<-permTest(A = breakpoints_40kb_4JMC19C, 
                                           x =gc_4JMC19C_gr, # Replace with your LINEs/GC object
                                           genome = custom_genome_4JMC19C,
                                           mask = mask_regions_4JMC19C,
                                           randomize.function = randomizeRegions, 
                                           per.chromosome = T, # Corrects for scaffold length
                                           evaluate.function = meanInRegions, 
                                           ntimes = 1000)




fold_1JMC18<-rep(NA,length(permtest_perchr_1JMC18))

fold_4JMC19C<-rep(NA,length(permtest_perchr_1JMC18))

zscore_1JMC18<-rep(NA,length(permtest_perchr_1JMC18))

zscore_4JMC19C<-rep(NA,length(permtest_perchr_1JMC18))

pvalue_1JMC18<-rep(NA,length(permtest_perchr_1JMC18))

pvalue_4JMC19C<-rep(NA,length(permtest_perchr_1JMC18))

observed_1JMC18<-rep(NA,length(permtest_perchr_1JMC18))

observed_4JMC19C<-rep(NA,length(permtest_perchr_1JMC18))

ci_low_1JMC18<-rep(NA,13)
ci_high_1JMC18<-rep(NA,13)

ci_low_4JMC19C<-rep(NA,13)
ci_high_4JMC19C<-rep(NA,13)


for (i in 1:12){
  zscore_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$numOverlaps$zscore

  zscore_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$numOverlaps$zscore

  pvalue_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$numOverlaps$pval

  pvalue_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$numOverlaps$pval

  observed_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$numOverlaps$observed

  observed_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$numOverlaps$observed

  fold_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$numOverlaps$observed/median(permtest_perchr_1JMC18[[i]]$numOverlaps$permuted)

  fold_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$numOverlaps$observed/median(permtest_perchr_4JMC19C[[i]]$numOverlaps$permuted)

  ci_low_1JMC18[i]<-quantile(permtest_perchr_1JMC18[[i]]$numOverlaps$permuted, 0.025)
  ci_high_1JMC18[i]<-quantile(permtest_perchr_1JMC18[[i]]$numOverlaps$permuted, 0.975)
  ci_low_4JMC19C[i]<-quantile(permtest_perchr_4JMC19C[[i]]$numOverlaps$permuted, 0.025)
  ci_high_4JMC19C[i]<-quantile(permtest_perchr_4JMC19C[[i]]$numOverlaps$permuted, 0.975)
  

}

i<-13

zscore_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$meanInRegions$zscore

zscore_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$meanInRegions$zscore

pvalue_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$meanInRegions$pval

pvalue_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$meanInRegions$pval

observed_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$meanInRegions$observed

observed_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$meanInRegions$observed

fold_1JMC18[i]<-permtest_perchr_1JMC18[[i]]$meanInRegions$observed/median(permtest_perchr_1JMC18[[i]]$meanInRegions$permuted)

fold_4JMC19C[i]<-permtest_perchr_4JMC19C[[i]]$meanInRegions$observed/median(permtest_perchr_4JMC19C[[i]]$meanInRegions$permuted)

ci_low_1JMC18[i]<-quantile(permtest_perchr_1JMC18[[i]]$meanInRegions$permuted, 0.025)
ci_high_1JMC18[i]<-quantile(permtest_perchr_1JMC18[[i]]$meanInRegions$permuted, 0.975)
ci_low_4JMC19C[i]<-quantile(permtest_perchr_4JMC19C[[i]]$meanInRegions$permuted, 0.025)
ci_high_4JMC19C[i]<-quantile(permtest_perchr_4JMC19C[[i]]$meanInRegions$permuted, 0.975)


#Define zero inflated distributions: 20 % 0s


prop_zeros_1JMC18<-rep(NA,12)
prop_zeros_4JMC19C<-rep(NA,12)

for ( i in 1:12){
  prop_zeros_1JMC18[i]<-length(which(permtest_perchr_1JMC18[[i]]$numOverlaps$permuted==0))/1000
  prop_zeros_4JMC19C[i]<-length(which(permtest_perchr_4JMC19C[[i]]$numOverlaps$permuted==0))/1000
  }
#zero inflated variables are the same between the two genomes but not the two zones

for (i in which(prop_zeros_1JMC18>0.2)){
  fold_1JMC18[i]<-(1+permtest_perchr_1JMC18[[i]]$numOverlaps$observed)/(1+mean(permtest_perchr_1JMC18[[i]]$numOverlaps$permuted))
  fold_4JMC19C[i]<-(1+permtest_perchr_4JMC19C[[i]]$numOverlaps$observed)/(1+median(permtest_perchr_4JMC19C[[i]]$numOverlaps$permuted))
}
