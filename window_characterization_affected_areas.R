 library(BiocManager)
library(GenomicRanges)
library(dplyr)
library(xlsx)
library(readxl)
library(rtracklayer)
 #Lo mismo que el otro pero considerando los extremos involucrados en reearrangements


setwd("C:/Users/inesg/Desktop/tesis/helodes")

###Guardar datos de sintenia

syn_window_1JMC18<-read_excel("window_analyisis/100kb_mod/syn_100kb_window_1JMC18_v2.xlsx")
syn_window_4JMC19C<-read_excel("window_analyisis/100kb_mod/syn_100kb_window_4JMC19C_v2.xlsx")

affected_pos_1JMC18<-which(syn_window_1JMC18$affected_areas=="affected")
affected_pos_4JMC19C<-which(syn_window_4JMC19C$affected_areas=="affected")

recomb_window_1JMC18<-read.csv("recombination/1JMC18/recombination_100kbwindow_1JMC18.csv")[,-1]
recomb_window_4JMC19C<-read.csv("recombination/4JMC19C/recombination_100kbwindow_4JMC19C.csv")[,-1]


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

dna_trans_1JMC18_gr<-import("repetitive_region_analysis/repeat_explorer/1JMC18/dante_classII_DNA_transposons.gff3",format="gff3")
dna_trans_4JMC19C_gr<-import("repetitive_region_analysis/repeat_explorer/4JMC19C/dante_classII_DNA_transposons.gff3",format="gff3")

trc_1JMC18_gr<-import("repetitive_region_analysis/repeat_explorer/1JMC18/Galaxy13-[TideCluster_on_data_2__GFF3_TideCluster_Output].gff3",format="gff3")
trc_4JMC19C_gr<-import("repetitive_region_analysis/repeat_explorer/4JMC19C/Galaxy21-[TideCluster_on_data_14__GFF3_TideCluster_Output].gff3",format="gff3")

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

#create the windows for analysis 100kb
chr_length_1JMC18<-as.numeric(read.table("genome_assembly/ensamblaje_helodes_final/qc_final/1JMC18/1JMC18_contig_length.txt")[,2])
#modifocar para incluir el ulttimo tile del scaffold 5
#chr_length_1JMC18[5]<-14412708

chr_mid_1JMC18<-chr_length_1JMC18/2
names(chr_length_1JMC18)<-paste("scaffold",1:35,sep="_")
tiles_1JMC18<-tileGenome(chr_length_1JMC18, tilewidth=100000)

  
  chr_length_4JMC19C<-as.numeric(read.table("genome_assembly/ensamblaje_helodes_final/qc_final/4JMC19C/4JMC19C_contig_length.txt")[,2])
names(chr_length_4JMC19C)<-paste("scaffold",1:36,sep="_")
tiles_4JMC19C<-tileGenome(chr_length_4JMC19C, tilewidth=100000, cut.last.tile.in.chrom=F)



###obtain number of elements from each window
total_re_window_1JMC18 <- countOverlaps(tiles_1JMC18, re_1JMC18_gr)
total_re_window_4JMC19C <- countOverlaps(tiles_4JMC19C, re_4JMC19C_gr)

unique_types_1JMC18<- unique(re_class_1JMC18$subclass)
re_types_windows_1JMC18<-as.list(unique_types_1JMC18)
for (type in unique_types_1JMC18) {
  re_types_windows_1JMC18[[type]]<-countOverlaps(tiles_1JMC18,re_types_granges_1JMC18[[type]])
}
re_types_windows_1JMC18<-re_types_windows_1JMC18[10:18]

unique_types_4JMC19C<- unique(re_class_4JMC19C$subclass)
re_types_windows_4JMC19C<-as.list(unique_types_4JMC19C)
for (type in unique_types_4JMC19C) {
  re_types_windows_4JMC19C[[type]]<-countOverlaps(tiles_4JMC19C,re_types_granges_4JMC19C[[type]])
}
re_types_windows_4JMC19C<-re_types_windows_4JMC19C[10:18]

dante_window_1JMC18<-countOverlaps(tiles_1JMC18, dante_1JMC18_gr)
dante_window_4JMC19C<-countOverlaps(tiles_4JMC19C, dante_4JMC19C_gr)

ltr_window_1JMC18<-countOverlaps(tiles_1JMC18, ltr_1JMC18_gr)
ltr_window_4JMC19C<-countOverlaps(tiles_4JMC19C, ltr_4JMC19C_gr)

dna_trans_window_1JMC18<-countOverlaps(tiles_1JMC18, dna_trans_1JMC18_gr)
dna_trans_window_4JMC19C<-countOverlaps(tiles_4JMC19C, dna_trans_4JMC19C_gr)

trc_window_1JMC18<-countOverlaps(tiles_1JMC18, trc_1JMC18_gr)
trc_window_4JMC19C<-countOverlaps(tiles_4JMC19C, trc_4JMC19C_gr)

total_genes_window_1JMC18<-countOverlaps(tiles_1JMC18, braker_1JMC18_gr_filt)
total_genes_window_4JMC19C<-countOverlaps(tiles_4JMC19C, braker_4JMC19C_gr_filt)



##Graficos exploratorios


par(mgp = c(3, 1.5, 0)) 
boxplot(
  total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")],
  total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "unaffected")], 
  total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "affected")],
  total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "unaffected")], 
  at = c(1, 2, 4, 5), # Posición de los grupos
  names = c(" \naffected\n 1JMC18", " \nunaffected\n 1JMC18", 
            " \naffected\n 4JMC19C", " \nunaffected\n 4JMC19C"), # Etiquetas
  ylab = "Number of Genes", # Etiqueta del eje Y
  col = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de las cajas
  whiskcol = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de los whiskers
  whisklty = 1, # Tipo de línea para whiskers
  whisklwd = 3, # Grosor de los whiskers
  border = "white", # Bordes blancos para las cajas
  outline = TRUE, # Mostrar outliers
  outbg = c("#FF9999", "#9999FF"),
  outcol="white",
  outpch = 21, # Forma de los puntos (16 = círculo sólido)
  cex.lab = 1.5, # Tamaño del texto de las etiquetas
  cex.axis = 1, # Tamaño del texto de los ejes
  las = 1 # Rotar las etiquetas de los ejes para mayor legibilidad
)


#statistical analysis for 4JMC19C
shapiro.test(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# no normal
shapiro.test(total_genes_window_4JMC19C[which(syn_window_1JMC18$affected_areas== "affected")])# normal

t.test(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],total_genes_window_4JMC19C[affected_pos_4JMC19C])#pvalue0.43
mean(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# 8.8
mean(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "affected")])#6.8

wilcox.test(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "affected")])#pvalue0.1214
median(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])#9
median(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "affected")])#9

#> quantile(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "affected")], probs = c(0.25, 0.75))
#25% 75% 
#  7  11 

#> quantile(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "unaffected")], probs = c(0.25, 0.75))
#25% 75% 
#  7  11

#difieren en valores extremos solo, los unaffected contienen tramos con mucha densidad de genes

#statistical analysis for 1JMC18
shapiro.test(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# no normal
shapiro.test(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")])# normal

t.test(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")],total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")])#pvalue0.3
mean(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# 8.8
mean(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")])#7.5

wilcox.test(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")],total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")])#pvalue0.38
median(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])#9
median(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")])#9


#> quantile(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")], probs = c(0.25, 0.75))
#25% 75% 
#  6  12 
#> quantile(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas== "unaffected")], probs = c(0.25, 0.75))
#25% 75% 
#  7  11 

#appear to be non significant difference: only extreme 

#excluding rearr points

wilcox.test(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],total_genes_window_4JMC19C
            [-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])#0.08

summary(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   8.000  10.000   9.718  11.000  23.000 
summary(total_genes_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   7.000   9.000   8.854  11.000  66.000 


wilcox.test(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],total_genes_window_4JMC19C
            [-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])#0.08

summary(total_genes_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   8.000  10.000   9.718  11.000  23.000 
summary(total_genes_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   7.000   9.000   8.854  11.000  66.000 

wilcox.test(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")],total_genes_window_1JMC18
            [-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])#0.08

summary(total_genes_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   8.000  10.000   9.718  11.000  23.000 
summary(total_genes_window_1JMC18[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   7.000   9.000   8.854  11.000  66.000 



###TOTAL RE EARLGREY

par(mgp = c(3, 1.5, 0)) 
boxplot(
  total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")],
  total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")], 
  total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")], 
  total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],
  at = c(1, 2, 4, 5), # Posición de los grupos
  names = c(" \naffected\n 1JMC18", " \nunaffected\n 1JMC18", 
            " \naffected\n 4JMC19C", " \nunaffected\n 4JMC19C"), # Etiquetas
  ylab = "RE content", # Etiqueta del eje Y
  col = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de las cajas
  whiskcol = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de los whiskers
  whisklty = 1, # Tipo de línea para whiskers
  whisklwd = 3, # Grosor de los whiskers
  border = "white", # Bordes blancos para las cajas
  outline = TRUE, # Mostrar outliers
  outbg = c("#FF9999", "#9999FF"),
  outcol="white",
  outpch = 21, # Forma de los puntos (16 = círculo sólido)
  cex.lab = 1.5, # Tamaño del texto de las etiquetas
  cex.axis = 1, # Tamaño del texto de los ejes
  las = 1 # Rotar las etiquetas de los ejes para mayor legibilidad
)

#statistical analysis for 4JMC19C
shapiro.test(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# no normal
shapiro.test(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])# normal

wilcox.test(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#pvalue0.54
median(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# 177
median(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#164
quantile(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "affected")], probs = c(0.25, 0.75))
#25% 75% 
#156.0 193.5 
quantile(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "unaffected")], probs = c(0.25, 0.75))
#25% 75% 
#155 198 


#statistical analysis for 1JMC18
shapiro.test(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# no normal
shapiro.test(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])# normal

wilcox.test(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")],total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])#pvalue0.09
median(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# 177
median(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])#171

quantile(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")], probs = c(0.25, 0.75))
#25% 75% 
#  150 195 
quantile(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas== "unaffected")], probs = c(0.25, 0.75))
#25% 75% 
#  157 199 
#affected tienen lower re content en 1JMC18 que en 4JMC19C


wilcox.test(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],total_re_window_4JMC19C
            [-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])#0.66

summary(total_re_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#60.0   166.5   179.0   176.0   197.5   257.0  
summary(total_re_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.0   155.0   177.0   175.7   198.0   302.0 

wilcox.test(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")],total_re_window_1JMC18
            [-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])#0.99

summary(total_re_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   8.000  10.000   9.718  11.000  23.000 
summary(total_re_window_1JMC18[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   7.000   9.000   8.854  11.000  66.000 


##TOTAL DANTE REPEATEXPLORER

par(mgp = c(3, 1.5, 0)) 
boxplot(
  dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")],
  dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")], 
  dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")], 
  dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],
  at = c(1, 2, 4, 5), # Posición de los grupos
  names = c(" \naffected\n 1JMC18", " \nunaffected\n 1JMC18", 
            " \naffected\n 4JMC19C", " \nunaffected\n 4JMC19C"), # Etiquetas
  ylab = "Tandem repeat content (DANTE)", # Etiqueta del eje Y
  col = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de las cajas
  whiskcol = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de los whiskers
  whisklty = 1, # Tipo de línea para whiskers
  whisklwd = 3, # Grosor de los whiskers
  border = "white", # Bordes blancos para las cajas
  outline = TRUE, # Mostrar outliers
  outbg = c("#FF9999", "#9999FF"),
  outcol="white",
  outpch = 21, # Forma de los puntos (16 = círculo sólido)
  cex.lab = 1.5, # Tamaño del texto de las etiquetas
  cex.axis = 1, # Tamaño del texto de los ejes
  las = 1 # Rotar las etiquetas de los ejes para mayor legibilidad
)

#statistical analysis for 4JMC19C
shapiro.test(dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# no normal
shapiro.test(dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])# normal

wilcox.test(dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#pvalue 0.8

quantile(dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "affected")], probs = c(0.25, 0.5,0.75))
#25% 50% 75% 
# 5   8  12 
quantile(dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "unaffected")], probs = c(0.25, 0.5,0.75))
#25% 50% 75% 
# 5   8  12 

###statistical analysis for 1JMC18
shapiro.test(dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# no normal
shapiro.test(dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])# normal

wilcox.test(dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")],dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])#pvalue0.95
quantile(dante_window_1JMC18[which(syn_window_1JMC18$affected_areas== "affected")], probs = c(0.25,0.5, 0.75))
#25% 50% 75% 
#5   8  13 
quantile(dante_window_1JMC18[which(syn_window_1JMC18$affected_areas== "unaffected")], probs = c(0.25,0.5, 0.75))
#25% 50% 75% 
#5   8  13 
#no hay diferencia significativas


wilcox.test(dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],dante_window_4JMC19C
            [-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])#0.66

summary(dante_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#60.0   166.5   179.0   176.0   197.5   257.0  
summary(dante_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.0   155.0   177.0   175.7   198.0   302.0 

wilcox.test(dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")],dante_window_1JMC18
            [-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])#0.99

summary(dante_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   8.000  10.000   9.718  11.000  23.000 
summary(dante_window_1JMC18[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   7.000   9.000   8.854  11.000  66.000 

###DANTE LTR TRANSPOSONS

par(mgp = c(3, 1.5, 0)) 
boxplot(
  ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")],
  ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")], 
  ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")], 
  ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],
  at = c(1, 2, 4, 5), # Posición de los grupos
  names = c(" \naffected\n 1JMC18", " \nunaffected\n 1JMC18", 
            " \naffected\n 4JMC19C", " \nunaffected\n 4JMC19C"), # Etiquetas
  ylab = "LTR content", # Etiqueta del eje Y
  col = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de las cajas
  whiskcol = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de los whiskers
  whisklty = 1, # Tipo de línea para whiskers
  whisklwd = 3, # Grosor de los whiskers
  border = "white", # Bordes blancos para las cajas
  outline = TRUE, # Mostrar outliers
  outbg = c("#FF9999", "#9999FF"),
  outcol="white",
  outpch = 21, # Forma de los puntos (16 = círculo sólido)
  cex.lab = 1.5, # Tamaño del texto de las etiquetas
  cex.axis = 1, # Tamaño del texto de los ejes
  las = 1 # Rotar las etiquetas de los ejes para mayor legibilidad
)

#statistical analysis for 4JMC19C
shapiro.test(ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# no normal
shapiro.test(ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])# no normal

wilcox.test(ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#pvalue0.12
quantile(ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas== "unaffected")], probs = c(0.25,0.5, 0.75))
#25% 50% 75% 
#0   5  11 
quantile(ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")], probs = c(0.25,0.5, 0.75))
#25%  50%  75% 
#0.00 3.00 8.75 

#statistical analysis for 1JMC18
shapiro.test(ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# no normal
shapiro.test(ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])# no normal

wilcox.test(ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")],ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])#pvalue0.5
quantile(ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")], probs = c(0.25,0.5, 0.75))
#25% 50% 75% 
#0   5  11 
quantile(ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")], probs = c(0.25,0.5, 0.75))
#25% 50% 75% 
#0   4  11 

#lower in affected areas


wilcox.test(ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],ltr_window_4JMC19C
            [-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])#0.04

summary(ltr_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   0.000   3.000   4.231   6.500  18.000 
summary(ltr_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   0.000   5.000   6.886  11.000  55.000 


wilcox.test(ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")],ltr_window_1JMC18
            [-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])#0.44

summary(ltr_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   4.000   5.846  11.000  28.000 
summary(ltr_window_1JMC18[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   0.000   5.000   6.871  11.000  61.000 

###tidecluster trc

# ---- Prepare Data Frame ----
# Combine all values into a long-format dataframe
# Create dataframes for each subset
df1 <- data.frame(
  value = trc_window_1JMC18[syn_window_1JMC18$affected_areas == "affected"],
  group = "affected\n1JMC18"
)

df2 <- data.frame(
  value = trc_window_1JMC18[syn_window_1JMC18$affected_areas == "unaffected"],
  group = "unaffected\n1JMC18"
)

df3 <- data.frame(
  value = trc_window_4JMC19C[syn_window_4JMC19C$affected_areas == "affected"],
  group = "affected\n4JMC19C"
)

df4 <- data.frame(
  value = trc_window_4JMC19C[syn_window_4JMC19C$affected_areas == "unaffected"],
  group = "unaffected\n4JMC19C"
)

# Combine all into one dataframe
df <- rbind(df1, df2, df3, df4)

# Convert group to factor with desired order
df$group <- factor(df$group, levels = c("affected\n1JMC18", "unaffected\n1JMC18", "affected\n4JMC19C", "unaffected\n4JMC19C"))


# ---- Plot ----
ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(aes(color = group),outlier.shape = NA, fill = "white", width = 0.6) +
  geom_jitter(aes(color = group), width = 0.15, size = 2, shape = 21, stroke = 0.5) +
  scale_fill_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
  scale_color_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
  labs(y = "TRC content", x = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )

#statistical analysis
shapiro.test(trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# no normal
shapiro.test(trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])# no normal

trc_window_4JMC19C_df<-as.data.frame(cbind(trc=trc_window_4JMC19C,
                             rearr_points=syn_window_4JMC19C$rearr_points,
                             affected_areas=syn_window_4JMC19C$affected_areas,
                             only_affected=rep(NA,length(trc_window_4JMC19C))))
trc_window_4JMC19C_df[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved"),]$only_affected<-"affected"
trc_window_4JMC19C_df[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved"),]$only_affected<-"unaffected"
# no significativo ninguno

trc_window_1JMC18_df<-as.data.frame(cbind(trc=trc_window_1JMC18,
                                           rearr_points=syn_window_1JMC18$rearr_points,
                                           affected_areas=syn_window_1JMC18$affected_areas,
                                           only_affected=rep(NA,length(trc_window_1JMC18))))
trc_window_1JMC18_df[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved"),]$only_affected<-"affected"
trc_window_1JMC18_df[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved"),]$only_affected<-"unaffected"

summary(zeroinfl(trc ~ rearr_points | 1, data=trc_window_4JMC19C_df ,dist = "poisson"))
summary(zeroinfl(trc ~ rearr_points | 1, data=trc_window_4JMC19C_df ,dist = "poisson"))


wilcox.test(trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#pvalue0.4
summary(trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])
summary(trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])


#statistical analysis for 1JMC18
shapiro.test(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# no normal
shapiro.test(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])# no normal

wilcox.test(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")],trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])#pvalue0.8
summary(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])
summary(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])

###only affected and conserved
#statistical analysis
summary(zeroinfl(trc ~ rearr_points | 1, data=trc_window_4JMC19C_df ,dist = "negbin"))#0.867
summary(zeroinfl(trc ~ affected_areas | 1, data=trc_window_4JMC19C_df ,dist = "negbin"))#0.895
summary(zeroinfl(trc ~ affected_areas | 1, data=trc_window_4JMC19C_df))#0.895



wilcox.test(trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],
            trc_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])#pvalue0.4
summary(trc_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
summary(trc_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])


#statistical analysis for 1JMC18
shapiro.test(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# no normal
shapiro.test(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])# no normal

wilcox.test(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")],trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])#pvalue0.8
summary(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])
summary(trc_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])

###ltr DNA TRANSPOSONS

par(mgp = c(3, 1.5, 0)) 
boxplot(
  dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")],
  dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")], 
  dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")], 
  dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],
  at = c(1, 2, 4, 5), # Posición de los grupos
  names = c(" \naffected\n 1JMC18", " \nunaffected\n 1JMC18", 
            " \naffected\n 4JMC19C", " \nunaffected\n 4JMC19C"), # Etiquetas
  ylab = "DNA transposons content", # Etiqueta del eje Y
  col = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de las cajas
  whiskcol = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de los whiskers
  whisklty = 1, # Tipo de línea para whiskers
  whisklwd = 3, # Grosor de los whiskers
  border = "white", # Bordes blancos para las cajas
  outline = TRUE, # Mostrar outliers
  outbg = c("#FF9999", "#9999FF"),
  outcol="white",
  outpch = 21, # Forma de los puntos (16 = círculo sólido)
  cex.lab = 1.5, # Tamaño del texto de las etiquetas
  cex.axis = 1, # Tamaño del texto de los ejes
  las = 1 # Rotar las etiquetas de los ejes para mayor legibilidad
)

#statistical analysis
shapiro.test(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])# no normal
shapiro.test(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])# no normal

wilcox.test(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#pvalue0.11
summary(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.4277  1.0000  4.0000 
summary(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#0.5
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.5926  1.0000  3.0000 

#statistical analysis for 1JMC18
shapiro.test(dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])# no normal
shapiro.test(dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])# no normal

wilcox.test(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected")],dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected")])#pvalue0.11
summary(dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected")])
median(dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected")])#1




wilcox.test(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],dna_trans_window_4JMC19C
            [-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])#0.19

summary(dna_trans_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#.0000  0.0000  0.0000  0.6154  1.0000  3.0000  
summary(dna_trans_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.4282  1.0000  4.0000 

wilcox.test(dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")],dna_trans_window_1JMC18
            [-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])#0.19

summary(dna_trans_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.0000  0.0000  0.4615  1.0000  3.0000

summary(dna_trans_window_1JMC18[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.4404  1.0000  5.0000 


##RE classes from earlgrey

shapirotest_re_types_1JMC18<-as.list(unique_types_1JMC18)
names(shapirotest_re_types_1JMC18)<-unique_types_1JMC18

wilcoxtest_re_types_1JMC18<-as.list(unique_types_1JMC18)
names(wilcoxtest_re_types_1JMC18)<-unique_types_1JMC18

median_affected_re_types_1JMC18<-as.list(unique_types_1JMC18)
names(median_affected_re_types_1JMC18)<-unique_types_1JMC18

median_unaffected_re_types_1JMC18<-as.list(unique_types_1JMC18)
names(median_unaffected_re_types_1JMC18)<-unique_types_1JMC18


for (i in 1:9){
  shapirotest_re_types_1JMC18[[i]]<-shapiro.test(re_types_windows_1JMC18[[i]][which(syn_window_1JMC18$affected_areas=="affected")])$p.value
  wilcoxtest_re_types_1JMC18[[i]]<-wilcox.test(re_types_windows_1JMC18[[i]][which(syn_window_1JMC18$affected_areas=="affected")],
                                               re_types_windows_1JMC18[[i]][which(syn_window_1JMC18$affected_areas=="unaffected")])$p.value
  median_affected_re_types_1JMC18[i]<-median(re_types_windows_1JMC18[[i]][which(syn_window_1JMC18$affected_areas=="affected")])
  median_unaffected_re_types_1JMC18[i]<-median(re_types_windows_1JMC18[[i]][which(syn_window_1JMC18$affected_areas=="unaffected")])
  }

#diferencia solo en 1JMC18 en TIR con affected lower



shapirotest_re_types_4JMC19C<-as.list(unique_types_4JMC19C)
names(shapirotest_re_types_4JMC19C)<-unique_types_4JMC19C

wilcoxtest_re_types_4JMC19C<-as.list(unique_types_4JMC19C)
names(wilcoxtest_re_types_4JMC19C)<-unique_types_4JMC19C

median_affected_re_types_4JMC19C<-as.list(unique_types_4JMC19C)
names(median_affected_re_types_4JMC19C)<-unique_types_4JMC19C

median_unaffected_re_types_4JMC19C<-as.list(unique_types_4JMC19C)
names(median_unaffected_re_types_4JMC19C)<-unique_types_4JMC19C


for (i in 1:9){
  shapirotest_re_types_4JMC19C[[i]]<-shapiro.test(re_types_windows_4JMC19C[[i]][which(syn_window_4JMC19C$affected_areas=="affected")])$p.value
  wilcoxtest_re_types_4JMC19C[[i]]<-wilcox.test(re_types_windows_4JMC19C[[i]][which(syn_window_4JMC19C$affected_areas=="affected")],
                                                re_types_windows_4JMC19C[[i]][which(syn_window_4JMC19C$affected_areas=="unaffected")])$p.value
  median_affected_re_types_4JMC19C[i]<-median(re_types_windows_4JMC19C[[i]][which(syn_window_4JMC19C$affected_areas=="affected")])
  median_unaffected_re_types_4JMC19C[i]<-median(re_types_windows_4JMC19C[[i]][which(syn_window_4JMC19C$affected_areas=="unaffected")])
}


###only affected and conserved
##1JMC18
wilcoxtest_re_types_1JMC18<-as.list(unique_types_1JMC18)
names(wilcoxtest_re_types_1JMC18)<-unique_types_1JMC18

median_affected_re_types_1JMC18<-as.list(unique_types_1JMC18)
names(median_affected_re_types_1JMC18)<-unique_types_1JMC18

median_unaffected_re_types_1JMC18<-as.list(unique_types_1JMC18)
names(median_unaffected_re_types_1JMC18)<-unique_types_1JMC18


for (i in 1:9){
  wilcoxtest_re_types_1JMC18[[i]]<-wilcox.test(re_types_windows_1JMC18[[i]][which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")],
                                               re_types_windows_1JMC18[[i]][-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])$p.value
  median_affected_re_types_1JMC18[i]<-median(re_types_windows_1JMC18[[i]][which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
  median_unaffected_re_types_1JMC18[i]<-median(re_types_windows_1JMC18[[i]][-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved")])
}

##4JMC19C
wilcoxtest_re_types_4JMC19C<-as.list(unique_types_4JMC19C)
names(wilcoxtest_re_types_4JMC19C)<-unique_types_4JMC19C

median_affected_re_types_4JMC19C<-as.list(unique_types_4JMC19C)
names(median_affected_re_types_4JMC19C)<-unique_types_4JMC19C

median_unaffected_re_types_4JMC19C<-as.list(unique_types_4JMC19C)
names(median_unaffected_re_types_4JMC19C)<-unique_types_4JMC19C


for (i in 1:9){
  wilcoxtest_re_types_4JMC19C[[i]]<-wilcox.test(re_types_windows_4JMC19C[[i]][which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")],
                                               re_types_windows_4JMC19C[[i]][-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])$p.value
  median_affected_re_types_4JMC19C[i]<-median(re_types_windows_4JMC19C[[i]][which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
  median_unaffected_re_types_4JMC19C[i]<-median(re_types_windows_4JMC19C[[i]][-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
}

unlist(wilcoxtest_re_types_4JMC19C)
#Simple_repeat            LTR        Unknown Low_complexity           LINE            TIR      Satellite           SINE 
#0.50454784     0.16459045     0.31735864     0.93476221     0.68764435     0.15598177     0.06714919     0.32554147 
#Helitron 
#0.81627845 

summary(re_types_windows_4JMC19C$Satellite[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.000   3.500   6.000   6.769   9.000  25.000 

summary(re_types_windows_4JMC19C$Satellite[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   5.000   7.000   7.815  10.000  35.000 


##recombination

recomb_window_4JMC19C$recombination1000<-recomb_window_4JMC19C$recombination*1000
recomb_window_1JMC18$recombination1000<-recomb_window_1JMC18$recombination*1000

par(mgp = c(5, 1.5, 0)) 
boxplot(
  recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected"),]$recombination1000,
  recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination1000, 
  recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected"),]$recombination1000,
  recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected"),]$recombination1000, 
  at = c(1, 2, 4, 5), # Posición de los grupos
  names = c(" \nAffected\n 1JMC18", " \nUnaffected\n 1JMC18", 
            " \nAffected\n 4JMC19C", " \nUnaffected\n 4JMC19C"), # Etiquetas
  ylab = "Recombination (cM/kb)", # Etiqueta del eje Y
  col = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de las cajas
  whiskcol = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de los whiskers
  whisklty = 1, # Tipo de línea para whiskers
  whisklwd = 3, # Grosor de los whiskers
  border = "white", # Bordes blancos para las cajas
  outline = TRUE, # Mostrar outliers
  outbg = c("#FF9999", "#9999FF"),
  outcol="white",
  outpch = 21, # Forma de los puntos (16 = círculo sólido)
  cex.lab = 1.5, # Tamaño del texto de las etiquetas
  cex.axis = 0.8, # Tamaño del texto de los ejes
  las = 1 # Rotar las etiquetas de los ejes para mayor legibilidad
)

###afected global

#statistical analysis for 4JMC19C
shapiro.test(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected"),]$recombination)# no normal
shapiro.test(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected"),]$recombination)# no normal

wilcox.test(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected"),]$recombination1000,recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected"),]$recombination1000)#pvalue0.06
summary(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected"),]$recombination1000)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.0000  0.0082  0.0148  0.0216  0.0264  0.3596     408 
summary(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected"),]$recombination1000)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#0.000000 0.003013 0.014024 0.014896 0.017788 0.089578        4 

#statistical analysis for 1JMC18
shapiro.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination)# no normal
shapiro.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected"),]$recombination)# no normal

wilcox.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination1000,recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected"),]$recombination1000)#pvalue4.794e-06
summary(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination1000)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.0000  0.0079  0.0141  0.0214  0.0259  0.3258     360 
summary(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected"),]$recombination1000)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#0.000000 0.003013 0.014024 0.014896 0.017788 0.089578        4

#las rearranged tienen menor recombination

###only affected
#statistical analysis for 4JMC19C
shapiro.test(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="unaffected"),]$recombination)# no normal
shapiro.test(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected"),]$recombination)# no normal

wilcox.test(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved"),]$recombination1000,
            recomb_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved"),]$recombination1000)#pvalue0.06
summary(recomb_window_4JMC19C[which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved"),]$recombination1000)
        #25%          50%          75% 
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#4.200e-07 2.983e-03 1.478e-02 1.470e-02 1.919e-02 8.958e-02 
summary(recomb_window_4JMC19C[-which(syn_window_4JMC19C$affected_areas=="affected" & syn_window_4JMC19C$rearr_points=="conserved"),]$recombination1000)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.0000  0.0082  0.0148  0.0216  0.0264  0.3596     412 
 

#statistical analysis for 1JMC18
shapiro.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination)# no normal
shapiro.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected"),]$recombination)# no normal

wilcox.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved"),]$recombination1000,
            recomb_window_1JMC18[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved"),]$recombination1000)#pvalue0.0004232
summary(recomb_window_1JMC18[-which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved"),]$recombination1000)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.0000  0.0078  0.0141  0.0213  0.0259  0.3258     365 
summary(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected" & syn_window_1JMC18$rearr_points=="conserved"),]$recombination1000)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#2.400e-07 2.122e-03 6.072e-03 1.270e-02 1.776e-02 1.129e-01 

#statistical analysis for 1JMC18
shapiro.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination)# no normal
shapiro.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected"),]$recombination)# no normal

wilcox.test(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination,recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="affected"),]$recombination)#pvalue0.5
quantile(recomb_window_1JMC18[which(syn_window_1JMC18$affected_areas=="unaffected"),]$recombination, probs = c(0.25,0.5, 0.75),na.rm=T)
#25%          50%          75% 
#7.895535e-06 1.410038e-05 2.588239e-05
quantile(recomb_window_1JMC18[which(syn_window_1JMC18$afffected_areas=="affected"),]$recombination, probs = c(0.25,0.5, 0.75),na.rm=T)
#25%  50%  75% 
#7.911307e-07 4.755337e-06 1.733494e-05


##gc content
gc_1JMC18<-read.table("window_analyisis/100kb_mod/100kb_window_1JMC18_affected_areas.gc.tsv")[,c(1:4,6)]
gc_4JMC19C<-read.table("window_analyisis/100kb_mod/100kb_window_4JMC19C_affected_areas.gc.tsv")[,c(1:4,6)]
colnames(gc_1JMC18)<-c("chr","start","end","affected_areas","gc")
colnames(gc_4JMC19C)<-c("chr","start","end","affected_areas","gc")

par(mgp = c(4, 1.5, 0),mar = c(5.1, 6, 4.1, 2.1)) 
boxplot(
  gc_1JMC18[which(gc_1JMC18$affected_areas=="affected"),]$gc,
  gc_1JMC18[which(gc_1JMC18$affected_areas=="unaffected"),]$gc, 
  gc_4JMC19C[which(gc_4JMC19C$affected_areas=="affected"),]$gc,
  gc_4JMC19C[which(gc_4JMC19C$affected_areas=="unaffected"),]$gc, 
  at = c(1, 2, 4, 5), # Posición de los grupos
  names = c(" \nAffected\n 1JMC18", " \nUnaffected\n 1JMC18", 
            " \nAffected\n 4JMC19C", " \nUnaffected\n 4JMC19C"), # Etiquetas
  ylab = "GC proportion", # Etiqueta del eje Y
  col = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de las cajas
  whiskcol = c("#FF9999", "#9999FF", "#FF9999", "#9999FF"), # Colores de los whiskers
  whisklty = 1, # Tipo de línea para whiskers
  whisklwd = 3, # Grosor de los whiskers
  border = "white", # Bordes blancos para las cajas
  outline = TRUE, # Mostrar outliers
  outbg = c("#FF9999", "#9999FF"),
  outcol="white",
  outpch = 21, # Forma de los puntos (16 = círculo sólido)
  cex.lab = 1.5, # Tamaño del texto de las etiquetas
  cex.axis = 0.8, # Tamaño del texto de los ejes
  las = 1 # Rotar las etiquetas de los ejes para mayor legibilidad
)


##statistical analysis for 4JMC19C
shapiro.test(gc_4JMC19C[which(gc_4JMC19C$affected_areas=="affected"),]$gc)# no normal
shapiro.test(gc_4JMC19C[which(gc_4JMC19C$affected_areas=="unaffected"),]$gc)# no normal

wilcox.test(gc_4JMC19C[which(gc_4JMC19C$affected_areas=="affected"),]$gc,gc_4JMC19C[which(gc_4JMC19C$affected_areas=="unaffected"),]$gc)#pvalue0.08
summary(gc_4JMC19C[which(gc_4JMC19C$affected_areas=="affected"),]$gc)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3166  0.3239  0.3316  0.3334  0.3409  0.3607 
summary(gc_4JMC19C[which(gc_4JMC19C$affected_areas=="unaffected"),]$gc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2986  0.3276  0.3354  0.3370  0.3438  0.5467

##statistical analysis for 1JMC18
shapiro.test(gc_1JMC18[which(gc_1JMC18$affected_areas=="affected"),]$gc)#0.25
shapiro.test(gc_1JMC18[which(gc_1JMC18$affected_areas=="unaffected"),]$gc)#no normal

wilcox.test(gc_1JMC18[which(gc_1JMC18$affected_areas=="affected"),]$gc,gc_1JMC18[which(gc_1JMC18$affected_areas=="unaffected"),]$gc)#pvalue0.26
summary(gc_1JMC18[which(gc_1JMC18$affected_areas=="affected"),]$gc)
#      0%      25%      50%      75%     100% 
#0.311289 0.325847 0.333357 0.342450 0.366297
quantile(gc_1JMC18[which(gc_1JMC18$affected_areas=="unaffected"),]$gc)
#0%      25%      50%      75%     100% 
#0.300116 0.327307 0.335227 0.343510 0.546476 


###only affected
##statistical analysis for 4JMC19C

gc_1JMC18_r<-read.table("window_analyisis/100kb_mod/100kb_window_1JMC18_rearr_points.gc.tsv")[,c(1:4,6)]
gc_4JMC19C_r<-read.table("window_analyisis/100kb_mod/100kb_window_4JMC19C_rearr_points.gc.tsv")[,c(1:4,6)]
colnames(gc_1JMC18_r)<-c("chr","start","end","affected_areas","gc")
colnames(gc_4JMC19C_r)<-c("chr","start","end","affected_areas","gc")

wilcox.test(gc_4JMC19C[which(gc_4JMC19C$affected_areas=="affected" & gc_4JMC19C_r$affected_areas=="conserved"),]$gc,
            gc_4JMC19C[-which(gc_4JMC19C$affected_areas=="affected" & gc_4JMC19C_r$affected_areas=="conserved"),]$gc)#pvalue0.08
summary(gc_4JMC19C[which(gc_4JMC19C$affected_areas=="affected" & gc_4JMC19C_r$affected_areas=="conserved"),]$gc)#lower in gc
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3166  0.3240  0.3313  0.3327  0.3398  0.3555 
summary(gc_4JMC19C[-which(gc_4JMC19C$affected_areas=="affected" & gc_4JMC19C_r$affected_areas=="conserved"),]$gc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2986  0.3275  0.3354  0.3369  0.3438  0.5467 

##statistical analysis for 1JMC18
wilcox.test(gc_1JMC18[which(gc_1JMC18$affected_areas=="affected" & gc_1JMC18_r$affected_areas=="conserved"),]$gc,
            gc_1JMC18[-which(gc_1JMC18$affected_areas=="affected" & gc_1JMC18_r$affected_areas=="conserved"),]$gc)#pvalue0.08
summary(gc_1JMC18[which(gc_1JMC18$affected_areas=="affected" & gc_1JMC18_r$affected_areas=="conserved"),]$gc)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3113  0.3254  0.3309  0.3325  0.3390  0.3601 
summary(gc_1JMC18[-which(gc_1JMC18$affected_areas=="affected" & gc_1JMC18_r$affected_areas=="conserved"),]$gc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3001  0.3273  0.3352  0.3366  0.3435  0.5465 


##only affected