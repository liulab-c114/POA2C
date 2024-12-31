rm(list = ls())
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(limma)
library(tidyverse)
library(xlsx)
library(ggrepel)

#differential analysis function
limma_comparsion <- function(data, meta, group_col, control = 'Control', case = 'Disease'){
  compair.meta <- mutate(meta, cc.group = meta[[group_col]]) %>% 
    filter(cc.group %in% c(control, case))
  compair.data <- select(data, rownames(compair.meta))
  
  group_list <- factor(compair.meta$cc.group,levels = c(case, control))
  design <- model.matrix(~0+group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(compair.data)
  fit <- lmFit(compair.data, design)
  cont.matrix <- makeContrasts(contrasts = paste0(case, '-',control)
                               , levels = design)
  #print(cont.matrix)
  fit2 <- contrasts.fit(fit,cont.matrix) 
  fit2 <- eBayes(fit2) 
  tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC")
  tmpOut$logPval <- -log10(tmpOut$P.Value)
  tmpOut$dir <- ifelse(tmpOut$P.Value > 0.05, 'no'
                       , ifelse(tmpOut$logFC > log(1.5), 'up'
                                , ifelse(tmpOut$logFC < -log(1.5), 
                                         'down', 'no')))
  tmpOut$dir <- factor(tmpOut$dir, levels = c('up', 'down', 'no'))
  tmpOut
  
}
#volcano plot function
plot_volcano <- function(limma_output){
  limma_output$label <- ifelse(limma_output$P.Value <= 0.05 & abs(limma_output$logFC) >= log(1.5)
                               , rownames(limma_output), "")
  limma_output$MD <- ifelse(limma_output$label == "Mdh2", rownames(limma_output), "")
  volcano.p <- ggplot(limma_output, aes(x = logFC, y = logPval, color = dir)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
    geom_vline(xintercept = c(-log(1.5), log(1.5)), linetype = 'dotted') + 
    xlim(c(-max(abs(limma_output$logFC)), max(abs(limma_output$logFC)))) +
    ylim(c(0, max(limma_output$logPval) + 1)) +
    scale_color_manual(values = c('#c0392b', '#2980b9', '#bdc3c7')) + 
    labs(x = "Log Fold Change", y = "-log10(p-value)",
         title = "Volcano Plot") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_text_repel(data = limma_output, aes(label = label), size = 2,
                    box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"),
                    segment.color = "black", show.legend = FALSE)
}

#read expression matrix
cere.data <- read.xlsx('v9.3.xlsx', row.names = 1, sheetIndex = 1)
meta.data <- data.frame(group = sub("_.*", "", colnames(cere.data))
                        , row.names = colnames(cere.data))#metadata construction

#Heatmap for average expression for Lean, Sham and RYGB group
plot.1 <- unique(meta.data$group)#3 groups
heat.1.avg.df <- data.frame(protein_name = rownames(cere.data))#dataframe construction
for(mg in c(plot.1)){
  heat.1.avg.df[[mg]] <- apply(cere.data[rownames(filter(meta.data, group == mg))], 1, mean)
}#mean expression calculation
heat.1.avg.df <- column_to_rownames(heat.1.avg.df, var = 'protein_name')
heat.1.avg.mt <- t(scale(t(heat.1.avg.df)))#data scaling
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#2980b9", "#ecf0f1", "#c0392b"))#Heatmap color
pdf('heatmap_avg_1.pdf', height = 4, width = 3)
Heatmap(heat.1.avg.mt,
        cluster_columns = F,
        col = col_fun,
        show_row_names = F,
        name = 'z-score')#draw heatmap
dev.off()

#Heatmap for individual samples for Lean, Sham and RYGB group
heat.2.df <- select(cere.data, contains(plot.1))#dataframe construction
heat.2.meta <- filter(meta.data, group %in% plot.1) %>% pull(group)
names(heat.2.meta) <- filter(meta.data, group %in% plot.1) %>% rownames()#metadata construction
heat.2.anno <- HeatmapAnnotation(group = heat.2.meta
                                 , col = list(group = c('lean' = '#9b59b6'
                                                        , 'sham' = '#e67e22'
                                                        , 'RYGB' = '#27ae60'
                                                        )))#Heatmap annotation construction
ha = rowAnnotation(foo = anno_mark(at = which(rownames(heat.2.df) %in% c("Mdh2", "Pgk1")), 
                                   labels = c("MDH2", "PGK1")))#Label protein MDH2 and PGK1
heat.2.df <- select(heat.2.df, names(heat.2.meta))
heat.2.mt <- t(scale(t(heat.2.df)))#Data scaling

pdf('heatmap_2.pdf', height = 8, width = 8)
Heatmap(heat.2.mt,
        top_annotation = heat.2.anno, 
        right_annotation = ha,
        cluster_columns = F,
        col = col_fun,
        column_split = c(rep('Lean', 6), rep('Sham', 6), rep('RYGB', 6)
                         ), #Splitting by groups
        show_row_names = F,
        name = 'z-score')
dev.off()

#Heatmap for differentially expressed proteins
sig.dep <- c(RYGB_lean_dep %>% 
               filter(dir %in% c('up', 'down')) %>% 
               rownames(),
             RYGB_sham_dep %>% 
               filter(dir %in% c('up', 'down')) %>% 
               rownames(),
             Lean_sham_dep %>%
               filter(dir %in% c('up', 'down')) %>%
               rownames()
) %>% unique()#Select differentially expressed proteins in 3 groups

sig.dep2.df <- heat.2.df %>% filter(rownames(.) %in% sig.dep) %>% 
  select(names(heat.2.meta))#dataframe construction
heat.3.mt <- t(scale(t(sig.dep2.df)))#Data scaling
heat.11.anno <- HeatmapAnnotation(group = heat.2.meta
                                  , col = list(group = c('lean' = '#9b59b6'
                                                         , 'sham' = '#e67e22'
                                                         , 'RYGB' = '#27ae60'
                                  )))#Heatmap annotation construction
rownames(heat.3.mt) <- toupper(rownames(heat.3.mt))
pdf('heatmap_3.pdf', height = 10, width = 12)
Heatmap(heat.3.mt,
        top_annotation = heat.11.anno, 
        cluster_columns = F,
        col = col_fun,
        column_split = c(rep('Lean', 6), rep('Sham', 6), rep('RYGB', 6)),
        show_row_names = T,
        name = 'z-score',
        row_names_gp = gpar(fontsize = 6))#Heatmap drawing
dev.off()

#Differential expression analysis
RYGB_sham_dep <- limma_comparsion(cere.data, meta.data, group_col = 'group'
                                  , control = 'sham', case = 'RYGB')
RYGB_lean_dep <- limma_comparsion(cere.data, meta.data, group_col = 'group'
                                  , control = 'lean', case = 'RYGB')
Lean_sham_dep <- limma_comparsion(cere.data, meta.data, group_col = 'group'
                                  , control = 'sham', case = 'lean')
write.csv(RYGB_lean_dep, 'RYGB_lean_DEP.csv')
write.csv(RYGB_sham_dep, 'RYGB_sham_DEP.csv')
write.csv(Lean_sham_dep, 'Lean_sham_DEP.csv')#Write DEP results

#Volcano plot drawing
RYGB_lean.vol.p <- plot_volcano(RYGB_lean_dep)
RYGB_sham.vol.p <- plot_volcano(RYGB_sham_dep)
Lean_sham.vol.p <- plot_volcano(Lean_sham_dep)
ggsave('volcano_RYGB_lean.pdf', RYGB_lean.vol.p, width = 5, height = 4)
ggsave('volcano_RYGB_sham.pdf', RYGB_sham.vol.p, width = 5, height = 4)
ggsave('volcano_Lean_sham.pdf', Lean_sham.vol.p, width = 5, height = 4)

#Venn plot for differentially expressed proteins
library(VennDiagram)
RYGB_lean_dep_vector <- filter(RYGB_lean_dep, dir != 'no') %>% rownames()
RYGB_sham_dep_vector <- filter(RYGB_sham_dep, dir != 'no') %>% rownames()
Lean_sham_dep_vector <- filter(Lean_sham_dep, dir != 'no') %>% rownames()
inter_DEP <- get.venn.partitions(x = list(`RYGB vs lean` = RYGB_lean_dep_vector
                                      , `RYGB vs sham` = RYGB_sham_dep_vector
                                      , `Lean vs sham` = Lean_sham_dep_vector))
for (i in 1 : nrow(inter_DEP)){
  inter_DEP[i, 'proteins'] <- paste(inter_DEP[[i, '..values..']], collapse = ', ')
}
write.xlsx(inter_DEP[-c(5, 6)], "vennlist_DEP.xlsx", )#get the detailed information for venn plot
venn.plot <- venn.diagram(
  x = list(`RYGB vs lean` = RYGB_lean_dep_vector
           , `RYGB vs sham` = RYGB_sham_dep_vector
           , `Lean vs sham` = Lean_sham_dep_vector),
  col = "black", 
  alpha = 0.5, 
  fill = c('#9b59b6', '#e67e22', '#27ae60'), 
  margin = 0.2, 
  cat.fontfamily = "Arial", 
  cat.dist = 0.2,
  filename = 'dep_venn.tiff')#Venn plot drawing

#Venn plot for up-regulated proteins
RYGB_lean_up_vector <- filter(RYGB_lean_dep, dir == 'up') %>% rownames()
RYGB_sham_up_vector <- filter(RYGB_sham_dep, dir == 'up') %>% rownames()
Lean_sham_up_vector <- filter(Lean_sham_dep, dir == 'up') %>% rownames()
inter_up <- get.venn.partitions(x = list(`RYGB vs lean` = RYGB_lean_up_vector
                                          , `RYGB vs sham` = RYGB_sham_up_vector
                                          , `Lean vs sham` = Lean_sham_up_vector))
for (i in 1 : nrow(inter_up)){
  inter_up[i, 'proteins'] <- paste(inter_up[[i, '..values..']], collapse = ', ')
}
write.xlsx(inter_up[-c(5, 6)], "vennlist_up.xlsx", )
venn.plot <- venn.diagram(
  x = list(`RYGB vs lean` = RYGB_lean_up_vector
           , `RYGB vs sham` = RYGB_sham_up_vector
           , `Lean vs sham` = Lean_sham_up_vector),
  col = "black", 
  alpha = 0.50, 
  fill = c('#9b59b6', '#e67e22', '#27ae60'), 
  margin = 0.2, 
  cat.fontfamily = "Arial", 
  cat.dist = 0.2,
  filename = 'dep_up_venn.tiff')# Venn plot drawing

#Venn plot for down-regulated proteins
RYGB_lean_down_vector <- filter(RYGB_lean_dep, dir == 'down') %>% rownames()
RYGB_sham_down_vector <- filter(RYGB_sham_dep, dir == 'down') %>% rownames()
Lean_sham_down_vector <- filter(Lean_sham_dep, dir == 'down') %>% rownames()
inter_down <- get.venn.partitions(x = list(`RYGB vs lean` = RYGB_lean_down_vector
                                          , `RYGB vs sham` = RYGB_sham_down_vector
                                          , `Lean vs sham` = Lean_sham_down_vector))
for (i in 1 : nrow(inter_down)){
  inter_down[i, 'proteins'] <- paste(inter_down[[i, '..values..']], collapse = ', ')
}
write.xlsx(inter_down[-c(5, 6)], "vennlist_down.xlsx", )
venn.plot <- venn.diagram(
  x = list(`RYGB vs lean` = RYGB_lean_down_vector
           , `RYGB vs sham` = RYGB_sham_down_vector
           , `Lean vs sham` = Lean_sham_down_vector),
  col = "black", 
  alpha = 0.50, 
  fill = c('#9b59b6', '#e67e22', '#27ae60'), 
  margin = 0.2, 
  cat.fontfamily = "Arial",
  cat.pos = c(-135, -180, -225),
  cat.dist = 0.2,
  filename = 'dep_down_venn.tiff')#Venn plot drawing

