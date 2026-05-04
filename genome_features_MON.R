library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)

setwd("C:/Users/inesg/Desktop/tesis/helodes/circus_graph")

# === 1. Cargar datos ===
charact <- read.table("feature_table_model_characteristics_4JMC19C.tsv", header = TRUE)
format_table <- read.table("feature_table_characteristics_4JMC19C.tsv", header = TRUE)

df<-cbind(format_table[,1:3],charact[,4:17])

breaks <- read.table("breakpoints_4JMC19C.txt", header = TRUE, sep = "\t")
head(df)


# Convertir posiciones a Mb
df <- df %>%
  mutate(Start = Start / 1e6,
         End = End / 1e6)

breaks <- breaks %>%
  mutate(Start = Start / 1e6,
         End = End / 1e6)

head(df)
head(breaks)

breaks_rearr<-breaks %>% arrange(Chromosome, Start)

# === 2. Calcular posición acumulada ===
# Ordenar por cromosoma
df <- df %>% arrange(Chromosome, Start)
df_rearr<- subset(df,Chromosome %in% c(5,6,8,9,16,17,19,24,32))  %>% arrange(Chromosome, Start)
  
# Calcular longitudes de cada scaffold/cromosoma
chr_lengths <- df %>%
  group_by(Chromosome) %>%
  summarise(chr_len = max(End))
chr_lengths <- chr_lengths %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))
df <- df %>%
  left_join(chr_lengths, by = "Chromosome") %>%
  mutate(Midpoint = (Start + End) / 2,
         GenomePos = Midpoint + offset)
breaks <- breaks %>%
  left_join(chr_lengths, by = "Chromosome") %>%
  mutate(Start_genome = Start + offset,
         End_genome = End + offset)

#For rearr
chr_lengths_rearr <- df_rearr %>%
  group_by(Chromosome) %>%
  summarise(chr_len = max(End))
chr_lengths_rearr <- chr_lengths_rearr %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))
df_rearr <- df_rearr %>%
  left_join(chr_lengths_rearr, by = "Chromosome") %>%
  mutate(Midpoint = (Start + End) / 2,
         GenomePos = Midpoint + offset)
breaks_rearr <- breaks_rearr %>%
  left_join(chr_lengths_rearr, by = "Chromosome") %>%
  mutate(Start_genome = Start + offset,
         End_genome = End + offset)

# === 3. Guardar posiciones de separación entre scaffolds ===
chr_boundaries <- chr_lengths %>%
  mutate(boundary = offset + chr_len)

chr_labels <- chr_lengths %>%
  mutate(midpoint = offset + chr_len / 2)

break_lines <- breaks %>%
  transmute(pos = (Start_genome + End_genome) / 2)

# === 4. Función para graficar (igual que antes) ===
plot_feature <- function(df, yvar, color, fill, ylab, breaks) {
  ggplot(df, aes(x = GenomePos, y = !!sym(yvar), group = Chromosome)) +
    # breakpoints como líneas gruesas
    geom_vline(data = break_lines, aes(xintercept = pos),
               color = "lightblue", linewidth = 0.8, alpha = 0.9) +
    # puntos
    geom_point(alpha = 0.2, size = 0.5) +
    # suavizado LOESS por cromosoma
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                span = 0.3, color = color, fill = fill, linewidth = 0.5) +
    # líneas divisorias entre cromosomas
    geom_vline(data = chr_boundaries, aes(xintercept = boundary),
               linetype = "dashed", color = "black", alpha = 0.6) +
    labs(x = "", y = ylab) +
    coord_cartesian(ylim=(c(min(df[[yvar]],na.rm=T),as.numeric(quantile(df[[yvar]],probs=0.995,na.rm=T))))) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
}
  

# === 5. Crear los paneles ===
p2 <- plot_feature(df, "genes", "purple", "thistle", "Genes", breaks)
p3 <- plot_feature(df, "re", "darkred", "mistyrose", "Repetitive elements", breaks)
p4 <- plot_feature(df, "dante", "darkorange", "moccasin", "Transposable elements", breaks)
p6 <- plot_feature(df, "ltr",  "navy", "lightblue", "LTR", breaks)
p8 <- plot_feature(df, "ty1_copia", "#1565C0", "#CDE4F5", "Ty1/copia", breaks)
p9 <- plot_feature(df, "ty3_gypsy", "#00838F", "#B2EBF2", "Ty3/gypsy", breaks)
p10 <- plot_feature(df, "line", "#7CB342", "#E6F4D7", "LINEs", breaks)
p11 <- plot_feature(df, "sine", "brown4", "bisque", "SINEs", breaks)
p12 <- plot_feature(df, "unknown", "#4E342E", "#D7CCC8", "Unknown repears", breaks)
p7 <- plot_feature(df, "trc",  "deeppink4", "pink", "Tandem.repeats", breaks)


# Panel inferior con etiquetas de cromosomas
p5 <- plot_feature(df, "gc", "darkgreen", "lightgreen", "GC content", breaks) +
  labs(x = "Chromosome") +
  scale_x_continuous(
    breaks = chr_labels$midpoint,
    labels = chr_labels$Chromosome
  ) +
  theme(
    axis.text.x = element_text(angle = 0, size = 10),
    axis.title.x = element_text(size = 10)
  )

# === 6. Combinar ===
combined_plot <- p2 / p4 / p6  / p8/p9 / p10 /p12 / p7 / p5

# === 7. Mostrar ===
print(combined_plot)

ggsave(
  filename = "4JMC19C_genome_features_100kb_genes_te_ltr_withtypes_unk_gc.pdf",
  plot = combined_plot,
  path = "C:/Users/inesg/Desktop/tesis/helodes/circus_graph",
  device = "pdf",
  width = 14,   
  height = 8,   
  units = "in",
  dpi = 300       
)


#####only rearranged chromosomes############

chr_boundaries_rearr <- chr_lengths_rearr %>%
  mutate(boundary = offset + chr_len)

chr_labels_rearr <- chr_lengths_rearr %>%
  mutate(midpoint = offset + chr_len / 2)

break_lines_rearr <- breaks_rearr %>%
  transmute(pos = (Start_genome + End_genome) / 2)

break_lines_rearr<-cbind(break_lines_rearr,chr=breaks_rearr$Chromosome)

# === 4. Función para graficar (igual que antes) ===


plot_feature_rearr <- function(df, yvar, color, fill, ylab, breaks) {
  ggplot(df, aes(x = GenomePos, y = !!sym(yvar), group = Chromosome)) +
    # breakpoints como líneas gruesas
    geom_vline(data = break_lines_rearr, aes(xintercept = pos, color=as.factor(chr)),
               linewidth = 0.8, alpha = 0.9, show.legend = F) +
    scale_color_manual(values=c("17"="lightblue","19"="lightblue",
                                "16"="#CBB69C",
                                "6"="#D8BFD8","8"="#FFA07A","24"="#FFD59A",
                                "5"="#FFD59A","32"="#99FF66","9"="#CBB69C"))+
    # puntos
    geom_point(alpha = 0.2, size = 0.5) +
    # suavizado LOESS por cromosoma
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                span = 0.3, color = color, fill = fill, linewidth = 0.5) +
    # líneas divisorias entre cromosomas
    geom_vline(data = chr_boundaries_rearr, aes(xintercept = boundary),
               linetype = "dashed", color = "black", alpha = 0.6) +
    labs(x = "", y = ylab) +
    coord_cartesian(ylim=(c(min(df[[yvar]],na.rm=T),as.numeric(quantile(df[[yvar]],probs=0.995,na.rm=T))))) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
}


# === 5. Crear los paneles ===
p2 <- plot_feature_rearr(df_rearr, "genes", "purple", "thistle", "Genes", breaks_rearr)
p3 <- plot_feature_rearr(df_rearr, "re", "darkred", "mistyrose", "RE", breaks_rearr)
p4 <- plot_feature_rearr(df_rearr, "dante", "darkorange", "moccasin", "TE", breaks_rearr)
p6 <- plot_feature_rearr(df_rearr, "ltr",  "navy", "lightblue", "LTRs", breaks_rearr)
p8 <- plot_feature_rearr(df_rearr, "ty1_copia", "#1565C0", "#CDE4F5", "Ty1/copia", breaks_rearr)
p9 <- plot_feature_rearr(df_rearr, "ty3_gypsy", "#00838F", "#B2EBF2", "Ty3/gypsy", breaks_rearr)
p10 <- plot_feature_rearr(df_rearr, "line", "#7CB342", "#E6F4D7", "LINEs", breaks_rearr)
p11 <- plot_feature_rearr(df_rearr, "sine", "brown4", "bisque", "SINEs", breaks_rearr)
p12 <- plot_feature_rearr(df_rearr, "unknown", "#4E342E", "#D7CCC8", "Unknown\n repeats", breaks_rearr)
p7 <- plot_feature_rearr(df_rearr, "trc",  "deeppink4", "pink", "Tandem\n repeats", breaks_rearr)


# Panel inferior con etiquetas de cromosomas
p5 <- plot_feature_rearr(df_rearr, "gc", "darkgreen", "lightgreen", "GC", breaks_rearr) +
  labs(x = "Chromosome") +
  scale_x_continuous(
    breaks = chr_labels_rearr$midpoint,
    labels = chr_labels_rearr$Chromosome
  ) +
  theme(
    axis.text.x = element_text(angle = 0, size = 10),
    axis.title.x = element_text(size = 10)
  )

# === 6. Combinar ===
combined_plot <- p2 / p4 / p6 / p8 /p9/p10 /p12 / p7 / p5

# === 7. Mostrar ===
print(combined_plot)

ggsave(
  filename = "4JMC19C_rearr_chr_genome_features_100kb_genes_te_ltr_unk_gc.pdf",
  plot = combined_plot,
  path = "C:/Users/inesg/Desktop/tesis/helodes/circus_graph",
  device = "pdf",
  width = 14,   
  height = 8,   
  units = "in",
  dpi = 300       
)
