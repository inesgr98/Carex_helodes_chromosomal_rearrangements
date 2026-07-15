library(regioneR)
library(GenomicRanges)
library(dplyr)

folders <- c( "inv2", "inv1", "inv3","trans_bp_chr16","trans_bp_chr24","trans_bp_chr5","trans_bp_chr9","fusion_chr17")

custom_genome_4JMC19C<-read.table("4JMC19C_contig_length.txt")
colnames(custom_genome_4JMC19C)<-c("chr","length")

mask_regions_4JMC19C <- GRanges()
telomere_size <- 10000

for(i in 1:nrow(custom_genome_4JMC19C)) {
  chr <- custom_genome_4JMC19C$chr[i]
  len <- custom_genome_4JMC19C$length[i]

  mask_regions_4JMC19C <- c(mask_regions_4JMC19C,
                           GRanges(chr, IRanges(start=1, end=telomere_size)),
                           GRanges(chr, IRanges(start=len - telomere_size, end=len)))
}

 root_dir <- getwd()

for ( i in 2:7){
  zona.gr <- toGRanges(paste(root_dir,"affected_areas_rearr.bed",sep="/"))[i,]

  setwd(file.path(root_dir, folders[i]))

  genes_df<-read.csv("df_input_regioneer.csv")
  genes_gr<-makeGRangesFromDataFrame(genes_df, keep.extra.columns=TRUE)

# 3. Función de Evaluación (Métrica Acumulada)
eval_metrics <- function(A, ...) {
  hits <- subsetByOverlaps(genes_gr, A)
  if(length(hits) == 0) return(c(dn=NA, ds=NA, omega=NA))

  sum_N_sites <- sum(hits$N_sites, na.rm=TRUE)
  sum_S_sites <- sum(hits$S_sites, na.rm=TRUE)
  sum_N_changes <- sum(hits$N_changes, na.rm=TRUE)
  sum_S_changes <- sum(hits$S_changes, na.rm=TRUE)

  if(sum_N_sites == 0 | sum_S_sites == 0) return(c(dn=NA, ds=NA, omega=NA))

  dn_agg <- sum_N_changes / sum_N_sites
  ds_agg <- sum_S_changes / sum_S_sites
  omega_agg <- if(ds_agg > 0) dn_agg / ds_agg else NA

  return(c(dn = dn_agg, ds = ds_agg, omega = omega_agg))
}

eval_dn    <- function(A, ...) eval_metrics(A, ...)[1]
eval_ds    <- function(A, ...) eval_metrics(A, ...)[2]
eval_omega <- function(A, ...) eval_metrics(A, ...)[3]

# 4. Ejecutar Permutación (per chromosome)
pt_dn <- permTest(A=zona.gr,
               evaluate.function=eval_dn,
               randomize.function=randomizeRegions,
               mask=mask_regions_4JMC19C,
               genome=custom_genome_4JMC19C, # Tu archivo .genome
               per.chromosome=TRUE,
               ntimes=100)

pt_ds <- permTest(A=zona.gr,
               evaluate.function=eval_ds,
               randomize.function=randomizeRegions,
               mask=mask_regions_4JMC19C,
               genome=custom_genome_4JMC19C, # Tu archivo .genome
               per.chromosome=TRUE,
               ntimes=100)

pt_omega <- permTest(A=zona.gr,
               evaluate.function=eval_omega,
               randomize.function=randomizeRegions,
               mask=mask_regions_4JMC19C,
               genome=custom_genome_4JMC19C, # Tu archivo .genome
               per.chromosome=TRUE,
               ntimes=100)

# 5. Guardar Resultados
save(pt_dn,pt_ds,pt_omega, file=paste0("permtest_",folders[i], ".RData"))
cat(folders[i])

}
