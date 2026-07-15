library(ggcorrplot)

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

# 2. Preparar los datos (usaremos el dataset 'mtcars' que viene en R)
corr_matrix <- cor(df[,-c(1:3,7,18:21)],method="spearman") # Calculamos la matriz de correlación

# 3. Crear el gráfico de correlación
ggcorrplot(corr_matrix, 
           hc.order = TRUE,           # Reordena las variables para agrupar las similares
           type = "lower",            # Muestra solo la mitad inferior
           lab = TRUE,                # Muestra los coeficientes numéricos
           lab_size = 4,              # Tamaño de los números
           method = "square",         # Forma de los marcadores (puedes usar "circle")
           colors = c("#6D9EC1", "white", "#E46726"), # Colores (Negativo, Neutral, Positivo)
           ggtheme = ggplot2::theme_minimal()
)


#find poiunt of influx:


