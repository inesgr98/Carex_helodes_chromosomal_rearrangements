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

breakpoints_1JMC18 <- toGRanges("window_analyisis/enrichment/breakpoints_1JMC18.txt")
breakpoints_4JMC19C <- toGRanges("window_analyisis/enrichment/breakpoints_4JMC19C.txt")

chr_length_1JMC18<-read.table("genome_assembly/ensamblaje_helodes_final/qc_final/1JMC18/1JMC18_contig_length.txt",header=F)
chr_length_4JMC19C<-read.table("genome_assembly/ensamblaje_helodes_final/qc_final/4JMC19C/4JMC19C_contig_length.txt",header=F)

fasta_1JMC18<-readDNAStringSet("genome_assembly/ensamblaje_helodes_final/clean_assemblies/1JMC18/1JMC18.fasta")
fasta_4JMC19C<-readDNAStringSet("repetitive_region_analysis/repeat_explorer/4JMC19C/4JMC19C_revcomp_mask.fasta")



###Guardar los outputs de braker y repeatmasker para GenomicRanges

braker_1JMC18_gr<-import("braker/after_earlgrey/C_helodes_1JMC18_braker.gff3", format = "gff3")
braker_1JMC18_gr_filt<-subset(braker_1JMC18_gr,type=="gene")
braker_4JMC19C_gr<-import("braker/after_earlgrey/C_helodes_4JMC19C_braker.gff3", format = "gff3")
braker_4JMC19C_gr_filt<-subset(braker_4JMC19C_gr,type=="gene")

re_1JMC18_gr<- import("repetitive_region_analysis/after_earlgrey/1JMC18/1JMC18.fasta.out.gff", format = "gff")
re_4JMC19C_gr<-import("repetitive_region_analysis/after_earlgrey/4JMC19C/4JMC19C.fasta.out.gff", format = "gff")

dante_1JMC18_gr<- import("repetitive_region_analysis/repeat_explorer/1JMC18/Galaxy16-[DANTE_on_data_2,_full_output].gff3", format = "gff")
dante_4JMC19C_gr<- import("repetitive_region_analysis/repeat_explorer/4JMC19C/Galaxy16-[DANTE_on_data_14,_full_output].gff3", format = "gff")

ltr_1JMC18_gr<-import("repetitive_region_analysis/repeat_explorer/1JMC18/Galaxy22-[LTR_retrotransposons_annotation_(GFF3)_________based_on_DANTE_annotation_16_and_reference_2].gff3",format="gff3")
ltr_4JMC19C_gr<-import("repetitive_region_analysis/repeat_explorer/4JMC19C/Galaxy27-[LTR_retrotransposons_annotation_(GFF3)_________based_on_DANTE_annotation_16_and_reference_14].gff3",format="gff3")
ltr_1JMC18_gr<-ltr_1JMC18_gr[mcols(ltr_1JMC18_gr)$source == "dante_ltr", ]
ltr_4JMC19C_gr<-ltr_4JMC19C_gr[mcols(ltr_4JMC19C_gr)$source == "dante_ltr", ]

ltr_1JMC18_gr_ty1 <- ltr_1JMC18_gr[ grepl("Ty1/copi", ltr_1JMC18_gr$Name, ignore.case = TRUE) ]
ltr_4JMC19C_gr_ty1 <- ltr_4JMC19C_gr[ grepl("Ty1/copi", ltr_4JMC19C_gr$Name, ignore.case = TRUE) ]
ltr_1JMC18_gr_ty3 <- ltr_1JMC18_gr[ grepl("Ty3/gyps", ltr_1JMC18_gr$Name, ignore.case = TRUE) ]
ltr_4JMC19C_gr_ty3 <- ltr_4JMC19C_gr[ grepl("Ty3/gyps", ltr_4JMC19C_gr$Name, ignore.case = TRUE) ]

dna_trans_1JMC18_gr<-import("repetitive_region_analysis/repeat_explorer/1JMC18/dante_classII_DNA_transposons.gff3",format="gff3")
dna_trans_4JMC19C_gr<-import("repetitive_region_analysis/repeat_explorer/4JMC19C/dante_classII_DNA_transposons.gff3",format="gff3")

dna_trans_1JMC18_gr<-import("repetitive_region_analysis/repeat_explorer/1JMC18/dante_classII_DNA_transposons.gff3",format="gff3")
dna_trans_4JMC19C_gr<-import("repetitive_region_analysis/repeat_explorer/4JMC19C/dante_classII_DNA_transposons.gff3",format="gff3")

dante_class1_nonltr_1JMC18_gr<-import("repetitive_region_analysis/repeat_explorer/1JMC18/dante_classI_nonltr.gff3",format="gff3")
dante_class1_nonltr_4JMC19C_gr<-import("repetitive_region_analysis/repeat_explorer/4JMC19C/dante_classI_nonltr.gff3",format="gff3")


trc_1JMC18_gr<-import("repetitive_region_analysis/repeat_explorer/1JMC18/Galaxy13-[TideCluster_on_data_2__GFF3_TideCluster_Output].gff3",format="gff3")
trc_4JMC19C_gr<-import("repetitive_region_analysis/repeat_explorer/4JMC19C/Galaxy21-[TideCluster_on_data_14__GFF3_TideCluster_Output].gff3",format="gff3")

trc1_1JMC18_gr<-trc_1JMC18_gr[mcols(trc_1JMC18_gr)$Name == "TRC_1", ]
trc1_4JMC19C_gr<-trc_4JMC19C_gr[mcols(trc_4JMC19C_gr)$Name == "TRC_1", ]

re_class_1JMC18<-read.table("repetitive_region_analysis/after_earlgrey/1JMC18/1JMC18.fasta.out",fill=T)
re_class_1JMC18<-re_class_1JMC18[-c(1,2),c(5:7,9,10,11)]
colnames(re_class_1JMC18)<-c("chr","start","end","strand","re","re_type")
re_class_1JMC18<-re_class_1JMC18[-which(re_class_1JMC18$re_type=="ARTEFACT"),]
re_class_1JMC18$subclass<-rep(NA,nrow(re_class_1JMC18))
re_class_1JMC18$start<-as.numeric(re_class_1JMC18$start)
re_class_1JMC18$end<-as.numeric(re_class_1JMC18$end)

re_class_1JMC18[which(re_class_1JMC18$re_type %in% grep("^LTR", unique(re_class_1JMC18$re_type), value = TRUE)),]$subclass<-"LTR"
re_class_1JMC18[which(re_class_1JMC18$re_type %in% grep("^LINE", unique(re_class_1JMC18$re_type), value = TRUE)),]$subclass<-"LINE"
re_class_1JMC18[which(re_class_1JMC18$re_type %in% grep("^DNA", unique(re_class_1JMC18$re_type), value = TRUE)),]$subclass<-"TIR"
re_class_1JMC18[which(re_class_1JMC18$re_type %in% c("tRNA",grep("^SINE", unique(re_class_1JMC18$re_type), value = TRUE))),]$subclass<-"SINE"
re_class_1JMC18[which(re_class_1JMC18$re_type=="RC/Helitron"),]$subclass<-"Helitron"
re_class_1JMC18[which(re_class_1JMC18$re_type=="Unknown"),]$subclass<-"Unknown"
re_class_1JMC18[which(re_class_1JMC18$re_type=="Simple_repeat"),]$subclass<-"Simple_repeat"
re_class_1JMC18[which(re_class_1JMC18$re_type=="Low_complexity"),]$subclass<-"Low_complexity"
re_class_1JMC18[which(re_class_1JMC18$re_type=="Satellite"),]$subclass<-"Satellite"

re_class_1JMC18_1<-re_class_1JMC18[which(is.na(re_class_1JMC18$subclass)),]

unique_types_1JMC18<- unique(re_class_1JMC18$subclass)
re_types_granges_1JMC18<-as.list(unique(re_class_1JMC18$subclass))
for (type in unique_types_1JMC18) {
  re_class_1JMC18_i<-subset(re_class_1JMC18,subclass==type)
  re_types_granges_1JMC18[[type]]<-GRanges(
    seqnames = re_class_1JMC18_i$chr,
    ranges = IRanges(start =  re_class_1JMC18_i$start, end =  re_class_1JMC18_i$end),
    re_type = re_class_1JMC18_i[[type]])
}

re_types_granges_1JMC18<-re_types_granges_1JMC18[10:18]

re_class_4JMC19C<-read.table("repetitive_region_analysis/after_earlgrey/4JMC19C/4JMC19C.fasta.out",fill=T)
re_class_4JMC19C<-re_class_4JMC19C[-c(1,2),c(5:7,9,10,11)]
colnames(re_class_4JMC19C)<-c("chr","start","end","strand","re","re_type")
re_class_4JMC19C<-re_class_4JMC19C[-c(11,which(re_class_4JMC19C$re_type=="ARTEFACT")),]#11 esta vacio
re_class_4JMC19C$subclass<-rep(NA,nrow(re_class_4JMC19C))
re_class_4JMC19C$start<-as.numeric(re_class_4JMC19C$start)
re_class_4JMC19C$end<-as.numeric(re_class_4JMC19C$end)

re_class_4JMC19C[which(re_class_4JMC19C$re_type %in% grep("^LTR", unique(re_class_4JMC19C$re_type), value = TRUE)),]$subclass<-"LTR"
re_class_4JMC19C[which(re_class_4JMC19C$re_type %in% grep("^LINE", unique(re_class_4JMC19C$re_type), value = TRUE)),]$subclass<-"LINE"
re_class_4JMC19C[which(re_class_4JMC19C$re_type %in% grep("^DNA", unique(re_class_4JMC19C$re_type), value = TRUE)),]$subclass<-"TIR"
re_class_4JMC19C[which(re_class_4JMC19C$re_type %in% c("tRNA",grep("^SINE", unique(re_class_4JMC19C$re_type), value = TRUE))),]$subclass<-"SINE"
re_class_4JMC19C[which(re_class_4JMC19C$re_type=="RC/Helitron"),]$subclass<-"Helitron"
re_class_4JMC19C[which(re_class_4JMC19C$re_type=="Unknown"),]$subclass<-"Unknown"
re_class_4JMC19C[which(re_class_4JMC19C$re_type=="Simple_repeat"),]$subclass<-"Simple_repeat"
re_class_4JMC19C[which(re_class_4JMC19C$re_type=="Low_complexity"),]$subclass<-"Low_complexity"
re_class_4JMC19C[which(re_class_4JMC19C$re_type=="Satellite"),]$subclass<-"Satellite"

re_class_4JMC19C_1<-re_class_4JMC19C[which(is.na(re_class_4JMC19C$subclass)),]
#> unique(re_class_4JMC19C_1$re_type)#[1] ""
#> unique(re_class_4JMC19C_1$re)#[1] ""
re_class_4JMC19C<-re_class_4JMC19C[-which(is.na(re_class_4JMC19C$subclass)),]


unique_types_4JMC19C<- unique(re_class_4JMC19C$subclass)
re_types_granges_4JMC19C<-as.list(unique_types_4JMC19C)
for (type in unique_types_4JMC19C) {
  re_class_4JMC19C_i<-subset(re_class_4JMC19C,subclass==type)
  re_types_granges_4JMC19C[[type]]<-GRanges(
    seqnames = re_class_4JMC19C_i$chr,
    ranges = IRanges(start =  re_class_4JMC19C_i$start, end =  re_class_4JMC19C_i$end),
    re_type = re_class_4JMC19C_i[[type]]) 
  
}
re_types_granges_4JMC19C<-re_types_granges_4JMC19C[10:18]

gc_1JMC18<-read.table("window_analyisis/100kb_mod/100kb_window_1JMC18_rearr_points.gc.tsv")[,c(1:4,6)]
gc_4JMC19C<-read.table("window_analyisis/100kb_mod/100kb_window_4JMC19C_rearr_points_fixed_mod.gc.tsv")[,c(1:4,6)]
colnames(gc_1JMC18)<-c("chr","start","end","rearr_points","gc")
colnames(gc_4JMC19C)<-c("chr","start","end","rearr_points","gc")
gc_1JMC18_gr <- GRanges(gc_1JMC18$chr, IRanges(gc_1JMC18$start + 1, gc_1JMC18$end), gc = gc_1JMC18$gc)
gc_4JMC19C_gr <- GRanges(gc_4JMC19C$chr, IRanges(gc_4JMC19C$start + 1, gc_4JMC19C$end), gc = gc_4JMC19C$gc)



###PREPARE LIST FOR PERMTEST###


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
     trc1=trc1_1JMC18_gr
     )


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