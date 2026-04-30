library(ggplot2)
library(patchwork)
library(dplyr)

setwd("C:/Users/inesg/Desktop/tesis/helodes/window_analyisis/100kb_mod")

data_charact_1JMC18<-read.table("feature_table_model_characteristics_1JMC18.tsv")
data_charact_4JMC19C<-read.table("feature_table_model_characteristics_4JMC19C.tsv")

#modify table for graphs

data_graph_1JMC18<- data_charact_1JMC18[,c(2,3,4,5,6,8,9,10,12,13,17)]
data_graph_4JMC19C<- data_charact_4JMC19C[,c(2,3,4,5,6,8,9,10,12,13,17)]

data_graph_1JMC18[which(data_graph_1JMC18$breakpoints=="rearranged"),]$breakpoints<-paste("breakpoints", "1JMC18",sep="\n")
data_graph_1JMC18[which(data_graph_1JMC18$breakpoints=="conserved"),]$breakpoints<-paste("conserved", "1JMC18",sep="\n")

data_graph_1JMC18[which(data_graph_1JMC18$syntheny=="rearranged"),]$syntheny<-paste("rearranged", "1JMC18",sep="\n")
data_graph_1JMC18[which(data_graph_1JMC18$syntheny=="conserved"),]$syntheny<-paste("conserved", "1JMC18",sep="\n")


data_graph_4JMC19C[which(data_graph_4JMC19C$breakpoints=="rearranged"),]$breakpoints<-paste("breakpoints", "4JMC19C",sep="\n")
data_graph_4JMC19C[which(data_graph_4JMC19C$breakpoints=="conserved"),]$breakpoints<-paste("conserved", "4JMC19C",sep="\n")

data_graph_4JMC19C[which(data_graph_4JMC19C$syntheny=="rearranged"),]$syntheny<-paste("rearranged", "4JMC19C",sep="\n")
data_graph_4JMC19C[which(data_graph_4JMC19C$syntheny=="conserved"),]$syntheny<-paste("conserved", "4JMC19C",sep="\n")


data_graph_1JMC18$breakpoints<-as.factor(data_graph_1JMC18$breakpoints)
data_graph_1JMC18$syntheny<-as.factor(data_graph_1JMC18$syntheny)

data_graph_4JMC19C$breakpoints<-as.factor(data_graph_4JMC19C$breakpoints)
data_graph_4JMC19C$syntheny<-as.factor(data_graph_4JMC19C$syntheny)


data_graph<-rbind(data_graph_1JMC18,data_graph_4JMC19C)

colnames(data_graph)<-c("syntheny","breakpoints","Genes","Repetitive elements","TEs","LTRs","Ty1/copia LTRs",
                               "Ty3/gypsy LTRs", "LINEs","Unknown repeats", "GC content")

data_graph$syntheny <- factor(data_graph$syntheny,
                                 levels = c("rearranged\n1JMC18","conserved\n1JMC18","rearranged\n4JMC19C","conserved\n4JMC19C"))

#dont run directly modify step-by-step``like axis and bw in violin

##breakpoints
plots_breakpoints<-as.list(1:9)

for (i in c(1:4,8:9)){
  varname<-colnames(data_graph)[i+2]
  data_graph$y_temp<-data_graph[,i+2]
plots_breakpoints[[i]] <-ggplot(data_graph, aes(x = breakpoints, y = y_temp, fill = breakpoints)) +
  geom_violin(trim = F, alpha = 0.8, color = "white") +
  geom_boxplot(width = 0.15, color = "white", outlier.shape = NA) + 
  coord_cartesian(ylim=(c(min(data_graph$y_temp),as.numeric(quantile(data_graph$y_temp,probs=0.992))))) +
  scale_fill_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
  labs(x = "", y = colnames(data_graph)[i+2]) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    legend.position = "none")
}

for (i in c(5:7)){
  varname<-colnames(data_graph)[i+2]
  data_graph$y_temp<-data_graph[,i+2]
  plots_breakpoints[[i]] <-ggplot(data_graph, aes(x = breakpoints, y = y_temp,color = breakpoints ,fill = breakpoints)) +
    geom_boxplot(aes(color = breakpoints),outlier.shape = NA, fill = "white", width = 0.6) +
    geom_jitter(aes(color = breakpoints), width = 0.15, size = 2, shape = 21, stroke = 0.5) +
    scale_fill_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
    scale_color_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
    coord_cartesian(ylim=(c(min(data_graph$y_temp)-1,as.numeric(quantile(data_graph$y_temp,probs=0.9999))))) +
    labs(y = varname, x = "") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      legend.position = "none"
    )
}


##syntheny rearranged

plots_syntheny<-as.list(1:9)

for (i in c(1:4,8:9)){
  varname<-colnames(data_graph)[i+2]
  data_graph$y_temp<-data_graph[,i+2]
  plots_syntheny[[i]] <-ggplot(data_graph, aes(x = syntheny, y = y_temp, fill = syntheny)) +
    geom_violin(trim = F, color = "white", alpha=0.8) +
    geom_boxplot(width = 0.15, color = "white", outlier.shape = NA) + 
    coord_cartesian(ylim=(c(min(data_graph$y_temp),as.numeric(quantile(data_graph$y_temp,probs=0.995))))) +
    scale_fill_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
    labs(x = "", y = varname) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      legend.position = "none")
}

for (i in c(5:7)){
  varname<-colnames(data_graph)[i+2]
  data_graph$y_temp<-data_graph[,i+2]
  plots_syntheny[[i]] <-ggplot(data_graph, aes(x = syntheny, y = y_temp,color = syntheny ,fill = syntheny)) +
    geom_boxplot(aes(color = syntheny),outlier.shape = NA, fill = "white", width = 0.6) +
    geom_jitter(aes(color = syntheny), width = 0.15, size = 2, shape = 21, stroke = 0.5) +
    scale_fill_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
    scale_color_manual(values = c("#FF9999", "#9999FF", "#FF9999", "#9999FF")) +
    coord_cartesian(ylim=(c(min(data_graph$y_temp)-1,as.numeric(quantile(data_graph$y_temp,probs=0.99))))) +
    labs(y = varname, x = "") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      legend.position = "none"
    )
}

plot_list<-c(plots_breakpoints,plots_syntheny)
final_plot_list<-lapply(plot_list, function(p) {
  p + theme(                     
    axis.title.x = element_text(size = 8)
  )
})

full_plot <- wrap_plots(final_plot_list, ncol = 2, byrow = FALSE)

ggsave("violin_plot_article_big_mod.jpeg", full_plot, width = 33, height = 50, units = "cm", dpi = 300)

