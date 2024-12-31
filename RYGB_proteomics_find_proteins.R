rm(list = ls())
setwd('H:/Bioinformatics/proteomics/RYGB_proteomics_v9.3/')
library(openxlsx)
library(ggplot2)
library(ggrepel)
venn_DEP <- read.xlsx('H:/Bioinformatics/proteomics/RYGB_proteomics_v9.3/vennlist_DEP.xlsx')
venn_meta <- venn_DEP$proteins[2]
venn_meta <- unlist(strsplit(venn_meta, ','))
venn_meta <- sub(" ", "", venn_meta)

RYGB_sham_DEP <- read.csv('H:/Bioinformatics/proteomics/RYGB_proteomics_v9.3/RYGB_sham_DEP.csv')
cp <- RYGB_sham_DEP
RYGB_sham_DEP <- RYGB_sham_DEP[RYGB_sham_DEP$X %in% venn_meta,]
RYGB_sham_DEP <- rbind(RYGB_sham_DEP, cp[cp$dir == 'no',])
RYGB_sham_DEP$label <- ifelse(RYGB_sham_DEP$logPval > 3 & abs(RYGB_sham_DEP$logFC) >= log(1.5),
                              RYGB_sham_DEP$X, "")
RYGB_sham_DEP$change <- ifelse(RYGB_sham_DEP$logPval > 3 & abs(RYGB_sham_DEP$logFC) >= log(1.5),
                             "sig_down", RYGB_sham_DEP$dir)
p1 <- ggplot(RYGB_sham_DEP, aes(x = logFC, y = logPval, color = change)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
  geom_vline(xintercept = c(-log(1.5), log(1.5)), linetype = 'dotted') + 
  xlim(c(-max(abs(RYGB_sham_DEP$logFC)), max(abs(RYGB_sham_DEP$logFC)))) +
  ylim(c(0, max(RYGB_sham_DEP$logPval) + 1)) +
  scale_color_manual(values = c('#2980b9', '#bdc3c7' , 'dodgerblue', '#c0392b')) + 
  labs(x = "Log Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_text_repel(data = RYGB_sham_DEP, aes(label = label), size = 5, color = "black",
                  box.padding = unit(2, "lines"), point.padding = unit(0.8, "lines"),
                  segment.color = "black", show.legend = FALSE)


pdf("vocano plot for DEP.pdf", width = 10, height = 8)
p1
dev.off()


RYGB_sham_DEP <- RYGB_sham_DEP[order(RYGB_sham_DEP$logPval), ]
RYGB_sham_DEP$prank <- seq(from = 1, to = 192, by = 1)
p2 <- ggplot(RYGB_sham_DEP, aes(x = prank, y = logPval, color = logPval)) + 
  geom_point(size = 0.8, alpha = 0.8) + 
  theme_bw() + 
  labs(x = "Rank by pvalue",
       y = "-log10(p)") +
  scale_color_gradient(low = '#EAC0BD', high = '#9C0824') +
  geom_text_repel(data = RYGB_sham_DEP, aes(label = label), size = 3, color = "black",
                  box.padding = unit(1, "lines"), point.padding = unit(0.8, "lines"),
                  segment.color = "black", show.legend = FALSE)


pdf("Find_MDH2/ranked_p.pdf", width = 4, height = 4)
p2
dev.off()

