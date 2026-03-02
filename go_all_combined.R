# --- 2. PREPARE INPUT DATA (REQUIRES: 'eggnog_4JMC19C.csv', 'protein_ids.txt') ---
library(topGO)
library(readxl)
library(xlsx)
library(dplyr)
library(org.At.tair.db)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(viridis)
library(patchwork)

setwd("C:/Users/inesg/Desktop/tesis/helodes/transcriptome/gene_ontology")
aa<-as.list(1:8)

gff_file <- "C:/Users/inesg/Desktop/tesis/helodes/transcriptome/braker/C_helodes_4JMC19C_isoseq.gff3"
gff_4JMC19C <- import(gff_file)
gff_4JMC19C_mrna <- gff_4JMC19C[gff_4JMC19C$type == "mRNA"]  # solo mRNA


eggnog_4JMC19C_file<-"C:/Users/inesg/Desktop/tesis/helodes/gene_ontology/4JMC19C/out.emapper.annotations.xlsx"
eggnog_4JMC19C<-read_excel(eggnog_4JMC19C_file,range = "A3:U23976")


for (i in 1:8){
  bed_file<-"C:/Users/inesg/Desktop/tesis/helodes/synteny_characterization/transform_syntheny_info/all_rearrangements.bed"
  bed <- read.table(bed_file, header=FALSE, stringsAsFactors=FALSE)[i,]
  colnames(bed) <- c("scaffold","start","end")
  
  gr_reg_4JMC19C <- GRanges(
    seqnames = bed$scaffold,
    ranges = IRanges(start = bed$start, end = bed$end),
  )
  ov <- findOverlaps(gff_4JMC19C_mrna, gr_reg_4JMC19C)
  mrna_in_regions <- gff_4JMC19C_mrna[queryHits(ov)]
  mrna_ids <- mcols(mrna_in_regions)$ID

  aa[[i]] <- intersect(mrna_ids, eggnog_4JMC19C$query)#184 out of 269 genes identifie
}

genes_of_interest<-unlist(aa)

# --- 3. CREATE CUSTOM GENE-GO MAPPING (FOR topGO) ---
go_mapping <- list()
for (i in 1:nrow(eggnog_4JMC19C)) {
  gene_id <- eggnog_4JMC19C$query[i]
  go_terms_string <- eggnog_4JMC19C$GOs[i]
  if (!is.na(go_terms_string) && go_terms_string != "" && go_terms_string != "-") {
    go_terms <- strsplit(go_terms_string, ",")[[1]]
    go_mapping[[gene_id]] <- go_terms
  }
}


# --- 3. CREATE CUSTOM GENE-GO MAPPING (FOR topGO) ---
go_mapping <- list()
for (i in 1:nrow(eggnog_4JMC19C)) {
  gene_id <- eggnog_4JMC19C$query[i]
  go_terms_string <- eggnog_4JMC19C$GOs[i]
  if (!is.na(go_terms_string) && go_terms_string != "" && go_terms_string != "-") {
    go_terms <- strsplit(go_terms_string, ",")[[1]]
    go_mapping[[gene_id]] <- go_terms
  }
}

if (length(go_mapping) == 0) {
  stop("ERROR: No GO terms were found in 'eggnog_4JMC19C.csv'. Please check your input file.")
}

# --- 4. PREPARE THE DATA FOR TOPGO ---
gene_list_factor <- as.factor(as.integer(eggnog_4JMC19C$query %in% genes_of_interest))
names(gene_list_factor) <- eggnog_4JMC19C$query

# --- 5. RUN THE GO ENRICHMENT ANALYSIS (Standard topGO) ---
GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = gene_list_factor, annot = annFUN.gene2GO, gene2GO = go_mapping, nodeSize = 5)
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = gene_list_factor, annot = annFUN.gene2GO, gene2GO = go_mapping, nodeSize = 5)
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = gene_list_factor, annot = annFUN.gene2GO, gene2GO = go_mapping, nodeSize = 5)

resultFisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
resultFisher_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
resultFisher_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")

# --- 6. COMBINE AND PROCESS SIGNIFICANT RESULTS ---
GenTable_Filtered <- function(GOdata, resultFisher, ontology_type) {
  GenTable(GOdata, classicFisher = resultFisher, topNodes = length(usedGO(GOdata))) %>%
    mutate(
      function_type = ontology_type,
      p.value = as.numeric(classicFisher) 
    ) %>%
    filter(p.value < 0.05) %>%
    mutate(negLogP = -log10(p.value))
}

combined_go_data <- bind_rows(
  GenTable_Filtered(GOdata_BP, resultFisher_BP, "BP"),
  GenTable_Filtered(GOdata_MF, resultFisher_MF, "MF"),
  GenTable_Filtered(GOdata_CC, resultFisher_CC, "CC")
)

combined_go_data$Term<-Term(GOTERM[combined_go_data$GO.ID])

#write.xlsx(combined_go_data,"go_table_rearr_combined.xlsx")

###Graph combined
bp_terms <- combined_go_data %>% filter(function_type == "BP")
all_go_ids_bp <- unique(bp_terms$GO.ID)
K_CLUSTERS <-  7 #Default number of clusters for BP term families

if (length(all_go_ids_bp) < 2) {
  cat("Less than two significant BP terms found. Cannot perform semantic clustering.\n")
  mds_df_bp <- data.frame(GO.ID = character(), semantic_x = numeric(), semantic_y = numeric(), cluster = factor())
} else {
  organism_db <- "org.At.tair.db" 
  go_data <- GOSemSim::godata(organism_db, ont = "BP", computeIC = TRUE)
  
  cat("1. Calculating BP-ONLY semantic similarity matrix (Resnik measure)...\n")
  sim_matrix_bp <- GOSemSim::goSim(
    GOID1 = all_go_ids_bp,
    GOID2 = all_go_ids_bp,
    measure = 'Resnik',
    semData = go_data
  )
  
  if (is.null(dim(sim_matrix_bp))) {
    n <- length(all_go_ids_bp)
    sim_matrix_bp <- matrix(sim_matrix_bp, nrow = n, ncol = n, dimnames = list(all_go_ids_bp, all_go_ids_bp))
  }
  
  sim_matrix_bp[is.na(sim_matrix_bp)] <- 0
  dist_matrix_bp <- as.dist(1 - sim_matrix_bp)
  
  cat("2. Running Multidimensional Scaling (MDS) on BP terms...\n")
  mds_result_bp <- cmdscale(dist_matrix_bp, k = 2)
  
  mds_df_bp <- data.frame(
    GO.ID = rownames(mds_result_bp),
    semantic_x = mds_result_bp[, 1],
    semantic_y = mds_result_bp[, 2]
  )
  
  # --- 7.5. CLUSTER BP TERMS (FOR FAMILY LABELING) ---
  if (nrow(mds_df_bp) > K_CLUSTERS) {
    cat(paste("3. Running k-means clustering (k=", K_CLUSTERS, ") to define GO families...\n"))
    bp_coords <- mds_df_bp[, c("semantic_x", "semantic_y")]
    kmeans_result <- kmeans(bp_coords, centers = K_CLUSTERS, nstart = 25)
    mds_df_bp$cluster <- as.factor(kmeans_result$cluster)
  } else {
    cat(paste("Not enough BP terms (", nrow(mds_df_bp), ") to perform k-means clustering with k=", K_CLUSTERS, ". Assigning all to one cluster.\n"))
    mds_df_bp$cluster <- as.factor(1) 
  }
}


# --- 8. PREPARE FINAL DATA WITH COORDINATES AND METRICS ---
final_data <- combined_go_data %>%
  left_join(mds_df_bp, by = "GO.ID") %>%
  mutate(
    # Enrichment Factor is used for size across ALL plots
    Enrichment_Factor = Significant / Annotated,  
    Total_Gene_Count = Annotated
  )

# --- 9. PLOT GENERATION FUNCTIONS ---

# 9a. BP Plot: Semantic Space
plot_bp_semantic <- function(data) {
  bp_data <- data %>% filter(function_type == "BP")
  
  ef_range <- range(bp_data$Enrichment_Factor, na.rm = TRUE)
  ef_breaks <- unique(c(round(quantile(ef_range, c(0.1, 0.5, 0.9)), 3), 
                        max(bp_data$Enrichment_Factor, na.rm = TRUE))) %>% 
    round(3)
  
  # Calculate cluster centroids and representative label
  cluster_labels <- bp_data %>%
    filter(!is.na(cluster)) %>%
    group_by(cluster) %>%
    summarise(
      semantic_x = mean(semantic_x),
      semantic_y = mean(semantic_y),
      # Find the term with the highest -log10(p-value) in the cluster to use as the label
      representative_term = Term[which.max(negLogP)],
      .groups = 'drop'
    ) %>%
    # Add a custom family label that includes the cluster number
    mutate(family_label = paste0("Family ", cluster, ": ", representative_term))
  
  
  p <- ggplot(bp_data, aes(x = semantic_x, y = semantic_y)) +
    # Color points by cluster family
    geom_point(aes(size = Enrichment_Factor, fill = negLogP, color = cluster),
               shape = 21, alpha = 0.8) +
    
    # Use different aesthetics for the cluster color boundary
    scale_color_viridis_d(name = "GO Family (Cluster)", 
                          option = "D", begin = 0.2, end = 0.8,
                          guide = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(size = 5, fill="white"))) +
    
    # Label the cluster centroids with the representative family term
    ggrepel::geom_label_repel(data = cluster_labels, 
                              aes(label = family_label),
                              size = 3,
                              box.padding = 0.6,
                              point.padding = 0.5,
                              force = 3,
                              max.overlaps = 100,
                              show.legend = FALSE,
                              fontface = "bold") +
    
    # --- LEGENDS ---
    scale_fill_viridis(option = "plasma", direction = -1,
                       name = "-Log10(p-value)", 
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5, ticks = FALSE)) +
    
    scale_size_continuous(range = c(2, 12),
                          name = "Enrichment Factor (Significant / Annotated)",
                          breaks = ef_breaks,
                          guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
    
    # Scales and Theme
    labs(
      title = "Biological Process (BP) Semantic Clustering and Family Grouping",
      subtitle = paste0("Clustering used k=", K_CLUSTERS, ". Family label is the most significant term in the cluster."),
      x = "Semantic X Coordinate (MDS)",
      y = "Semantic Y Coordinate (MDS)"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", 
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5), # Subtitle size kept default/large as it is short
          legend.title = element_text(face = "bold"))
  
  return(p)
}

# 9b/9c. CC/MF Plot: Annotated vs. Significant
plot_ontology_scatter <- function(data, ontology_name, color_palette = "viridis") {
  ontology_data <- data %>% filter(function_type == ontology_name)
  
  # Select the top 15 most significant terms for labeling (Updated from 10 to 15)
  label_data <- ontology_data %>% 
    arrange(desc(negLogP)) %>% 
    slice_head(n = 15)
  
  ef_range <- range(ontology_data$Enrichment_Factor, na.rm = TRUE)
  ef_breaks <- unique(c(round(quantile(ef_range, c(0.1, 0.5, 0.9)), 3), 
                        max(ontology_data$Enrichment_Factor, na.rm = TRUE))) %>% 
    round(3)
  
  p <- ggplot(ontology_data, aes(x = Annotated, y = Significant)) +
    geom_point(aes(size = Enrichment_Factor, fill = negLogP),
               shape = 21, color = "black", alpha = 0.8) +
    
    # Label the top 15 significant terms
    ggrepel::geom_label_repel(data = label_data, 
                              aes(label = Term),
                              size = 3.5,
                              box.padding = 0.5,
                              point.padding = 0.5,
                              force = 1,
                              max.overlaps = 50,
                              show.legend = FALSE) +
    
    # --- LEGENDS ---
    scale_fill_viridis(option = color_palette, direction = -1,
                       name = "-Log10(p-value)", 
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5, ticks = FALSE)) +
    
    scale_size_continuous(range = c(2, 12),
                          name = "Enrichment Factor (Significant / Annotated)",
                          breaks = ef_breaks,
                          guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
    
    # Scales and Theme
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    
    labs(
      title = paste(ontology_name, "Enrichment Scatter Plot"),
      # Subtitle text updated to reflect 'Top 15'
      subtitle = "Top 15 Terms Labeled. X: Total Annotated Genes. Y: Significant Genes. Size: Enrichment Factor.",
      x = "Total Annotated Genes (Annotated)",
      y = "Significant Genes (Significant)"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", 
          plot.title = element_text(face = "bold", hjust = 0.5),
          # Subtitle font size reduced to prevent overlap (New adjustment)
          plot.subtitle = element_text(hjust = 0.5, size = 11), 
          legend.title = element_text(face = "bold"))
  
  return(p)
}

# --- 10. RENDER AND SAVE THE INDIVIDUAL PLOTS ---

# 10a. BP Plot
p_bp <- plot_bp_semantic(final_data)
pdf("plot_1_BP_semantic_clustering.pdf", width = 12, height = 10)
print(p_bp)
dev.off()
cat("\nSaved plot_1_BP_semantic_clustering.pdf\n")

# 10b. CC Plot
p_cc <- plot_ontology_scatter(final_data, "CC", color_palette = "plasma")
pdf("plot_2_CC_enrichment_scatter.pdf", width = 12, height = 10)
print(p_cc)
dev.off()
cat("Saved plot_2_CC_enrichment_scatter.pdf\n")

# 10c. MF Plot
p_mf <- plot_ontology_scatter(final_data, "MF", color_palette = "inferno")
pdf("plot_3_MF_enrichment_scatter.pdf", width = 12, height = 10)
print(p_mf)
dev.off()
cat("Saved plot_3_MF_enrichment_scatter.pdf\n")


# --- 11. COMBINE ALL THREE PLOTS AND SAVE (FINAL ROBUST PATCHWORK LOGIC with added margins) ---

# 1. Define plots without legends AND with increased internal margins
plot_bp_no_legend <- p_bp + theme(legend.position = "none", plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"))
plot_cc_no_legend <- p_cc + theme(legend.position = "none", plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"))
plot_mf_no_legend <- p_mf + theme(legend.position = "none", plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"))

# 2. Use wrap_plots explicitly to combine CC and MF horizontally (bottom row)
bottom_row <- wrap_plots(plot_cc_no_legend, plot_bp_no_legend, nrow = 1)

# 3. Use wrap_plots explicitly to combine BP (top) and bottom_row (bottom) vertically
combined_plot <- wrap_plots(plot_mf_no_legend, plot_cc_no_legend, plot_bp_no_legend, ncol = 3 )

# 4. Apply global settings (title, collected legend)
#combined_plot <- combined_plot + 
#  plot_layout(guides = "collect") + # Collect all legends
#  plot_annotation(title = "Combined Gene Ontology Enrichment Analysis",
#                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))) & 
#  theme(legend.position = "bottom", legend.box.margin = margin(t=15), text = element_text(size = 14)) # Apply global theme

pdf("plot_4_Combined_GO_Analysis_all_rearr.pdf", width = 16, height = 16)
print(combined_plot)
dev.off()
cat("Saved plot_4_Combined_GO_Analysis.pdf (All three plots combined).\n")
