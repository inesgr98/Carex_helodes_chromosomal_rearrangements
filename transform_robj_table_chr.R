library(flextable)
library(officer)
library(dplyr)

# 1. Crear el dataframe con tus datos

df_final_1JMC18 <- data.frame(
  Variables = c("Genes", "RE", "TE", "LTR", "Ty1/copia", "Ty3/gypsy", 
                "DNA transp", "LINE", "SINE", "Unknown", "TRC", 
                "Most common TRC", "GC content"),
  zvalue=round(zscore_1JMC18,2),
  Observed=round(observed_1JMC18,2),
  CIs = paste0("[", round(ci_low_1JMC18, 2), ", ", round(ci_high_1JMC18, 2), "]"),
  Fold_enrichment = round(fold_1JMC18, 2),
  pvalue=round(pvalue_1JMC18,3)
)


df_final_4JMC19C <- data.frame(
  Variables = c("Genes", "RE", "TE", "LTR", "Ty1/copia", "Ty3/gypsy", 
                "DNA transp", "LINE", "SINE", "Unknown", "TRC", 
                "Most common TRC", "GC content"),
  zvalue=round(zscore_4JMC19C,2),
  Observed=round(observed_4JMC19C,2),
  CIs = paste0("[", round(ci_low_4JMC19C, 2), ", ", round(ci_high_4JMC19C, 2), "]"),
  Fold_enrichment = round(fold_4JMC19C, 2),
  pvalue= round(pvalue_4JMC19C,3)
)



# Si quieres comparar ambos genotipos en la misma tabla:
# Opción A: Tabla para 1JMC18
df_1JMC18 <- df_final_1JMC18 %>%
  arrange(match(Variables, c("Genes", "RE", "TE", "LTR", "Ty1/copia", "Ty3/gypsy", 
                             "DNA transp", "LINE", "SINE", "Unknown", "TRC", 
                             "Most common TRC", "GC content")))



tabla_jerarquica <- flextable(df_1JMC18) %>%
  # 1. Combinar verticalmente la columna de Variables
  #merge_v(j = 1) %>% 
  
  # 2. Formato de cabeceras
  set_header_labels(
    Variables = "Genomic Feature",
    zvalue = "z-value",
    Observed = "Observed",
    CIs = "95% CIs",
    Fold_enrichment = "Fold Enrichment",
    pvalue = "p-value"
  ) %>%
  
  # 3. Estilo profesional
  theme_booktabs() %>%
  autofit() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  
  # 4. Alineación y estética
  align(j = 1, align = "left", part = "all") %>%
  align(j = 3:6, align = "center", part = "all") %>%
  valign(j = 1, valign = "top") %>%  # El nombre de la variable queda arriba en la celda combinada
  bold(part = "header") %>%

  # 5. Añadir sangría a las sub-familias para mejorar la jerarquía visual
  padding(i = ~ grepl("Ty1|Ty3", Variables), j = 1, padding.left = 15) %>%
  # 7. Línea de separación entre bloques de variables
  hline(i = seq(1, nrow(df_1JMC18), by = 1), 
        border = fp_border(color = "gray80", width = 1)) %>%
  hline_bottom(border = fp_border(color = "black", width = 2), part = "body")

tabla_jerarquica <- set_table_properties(tabla_jerarquica, layout = "autofit")
  
# Creamos un documento Word vacío
doc <- read_docx()

# Añadimos la tabla
doc <- body_add_flextable(doc, value = tabla_jerarquica)

# Guardamos el archivo en tu carpeta de trabajo
print(doc, target = "enrich_perchr_AZN.docx")


##4JMC19C

df_4JMC19C <- df_final_4JMC19C %>%
  arrange(match(Variables, c("Genes", "RE", "TE", "LTR", "Ty1/copia", "Ty3/gypsy", 
                             "DNA transp", "LINE", "SINE", "Unknown", "TRC", 
                             "Most common TRC", "GC content")))

# Reordenamos columnas para que Region esté junto a Variables



tabla_jerarquica <- flextable(df_4JMC19C) %>%
  # 1. Combinar verticalmente la columna de Variables
#  merge_v(j = 1) %>% 
  
  # 2. Formato de cabeceras
  set_header_labels(
    Variables = "Genomic Feature",
    Region = "Region",
    zvalue = "z-value",
    Observed = "Observed",
    CIs = "95% CIs",
    Fold_enrichment = "Fold Enrichment",
    pvalue = "p-value"
  ) %>%
  
  # 3. Estilo profesional
  theme_booktabs() %>%
  autofit() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  
  # 4. Alineación y estética
  align(j = 1:2, align = "left", part = "all") %>%
  align(j = 3:6, align = "center", part = "all") %>%
  valign(j = 1, valign = "top") %>%  # El nombre de la variable queda arriba en la celda combinada
  bold(part = "header") %>%
  
  # 5. Añadir sangría a las sub-familias para mejorar la jerarquía visual
  padding(i = ~ grepl("Ty1|Ty3", Variables), j = 1, padding.left = 15) %>%
  # 7. Línea de separación entre bloques de variables
  hline(i = seq(1, nrow(df_4JMC19C), by = 1), 
        border = fp_border(color = "gray80", width = 1)) %>%
  hline_bottom(border = fp_border(color = "black", width = 2), part = "body")

tabla_jerarquica <- set_table_properties(tabla_jerarquica, layout = "autofit")

# Creamos un documento Word vacío
doc <- read_docx()

# Añadimos la tabla
doc <- body_add_flextable(doc, value = tabla_jerarquica)

# Guardamos el archivo en tu carpeta de trabajo
print(doc, target = "enrich_perchr_MON.docx")
