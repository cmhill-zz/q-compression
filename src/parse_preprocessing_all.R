library(ggplot2)
library(gplots)
library(gtools)
library(xtable)
library(scales)
#install.packages("gplots")

tables = c("~/cbcb2/q-compression/results/rhodo3/results/frag_1.fastq.preprocessing",
           #"~/cbcb2/q-compression/results/staph2/results/frag_1.fastq.preprocessing",
           "~/cbcb2/q-compression/results/ecoli2/results/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.preprocessing",
           "~/cbcb2/q-compression/results/hg14/results/frag_1.fastq.preprocessing",
           "~/cbcb2/q-compression/results/mouse/results//SRR032209.fastq.preprocessing")

dataset_names = c("Rhodobacter sphaeroides",
                  #"S. aureus",
                  "Escherichia coli",
                  "Homo sapiens",
                  "Mus musculus")

dfr_tables <- c()
original_tables <- c()
counter <- 1

for (table in tables) {
  table <- read.table(table, header = T, row.names = 1)
  original <- table[rownames(table) %in% "original", ]$bases_kept
  table <-  table[!rownames(table) %in% c("original","minqual","maxqual"), ]
  table[,5] <- table[,5] - original
  
  dfr <- data.frame(y=seq(1:length(table[,1])), x = table[,5], dataset = dataset_names[[counter]])
  dfr_tables <- rbind(dfr_tables, dfr)
  
  counter <- counter + 1
}

rowLabels = c("2-bin", "regression (0)", "regression (1)", "regression (3)", "regression (5)",
              "regression (7)", "profile (64)", "profile (128)", "profile (256)", "QualComp (6)", "QualComp (10)",
              "QualComp (30)", "QualComp (100)")

library(grid)


pdf("preprocessing_results.pdf", width=2, height=4)
p <- ggplot(dfr_tables, aes(x=jitter(x,factor = 1), y=y, color = dataset)) + geom_point() + scale_y_discrete(labels=rowLabels) + 
  theme_bw() +
  xlab("Proportion of bases kept compared to original") +
  ylab("") + 
  theme_bw() + theme(legend.position="bottom", legend.key = element_blank(), legend.title=element_blank(),
                     legend.text = element_text(face="italic", size=6),
                     axis.text.y = element_text(size=6),
                     axis.title.x = element_text(size=6),
                     axis.text.x = element_text(size=6)) +
  guides(col = guide_legend(ncol = 2, byrow = TRUE)) +  theme(legend.key.height=unit(0.5,"line"))

# Centers the legend at the bottom.
g <- ggplotGrob(p)
id <- which(g$layout$name == "guide-box")
g$layout[id, c("l","r")] <- c(1, ncol(g))
#grid.newpage()
grid.draw(g)

dev.off()
