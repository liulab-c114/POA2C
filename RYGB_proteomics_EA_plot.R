rm(list = ls())
library(ggplot2)
library(openxlsx)
library(tidyverse)
go_enrich = read.xlsx('metascape_result.xlsx', sheet = 2)#read enrichment analysis result
go_enrich <- go_enrich[grepl("Summary", go_enrich$GroupID), ]
go_enrich <- go_enrich[1:10,]
go_enrich$term <- paste(go_enrich$Term, go_enrich$Description, sep = ': ')#enrichment term construction
go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term, ordered = T)
go_enrich$Hits <- as.numeric(sub("/.*", "", go_enrich$InTerm_InList))
p1 <- ggplot(go_enrich,
             aes(x = reorder(term, Hits), y = Hits, fill = Hits)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient(low = "#F4D166", high = "#E8681D", name = "-log10(P)" ) +
  coord_flip() +
  xlab("GO term") +
  ylab("Hits") +
  labs(title = "GO Terms Enrichment") + 
  theme(text = element_text(family = "Times", size = 12, face = "bold"), 
        #legend.position = 'none',
        panel.background = element_rect(fill = NA),
        axis.line = element_line(color = "black")
        )#plot enrichment analysis result


pdf("RYGB_sham_GO_top10_withhit.pdf", width = 15, height = 8)  
p1  
dev.off()
