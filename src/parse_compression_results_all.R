library(ggplot2)

tables = c("~/cbcb2/q-compression/results/rhodo3/results/frag_2.fastq.compression",
           #"~/cbcb2/q-compression/results/staph2/results/frag_1.fastq.compression",
           "~/cbcb2/q-compression//results/ecoli2/results/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.compression",
           "~/cbcb2/q-compression//results/hg14/results/frag_1.fastq.compression",
           "~/cbcb2/q-compression/results/mouse/results//SRR032209.fastq.compression")

dataset_names = c("Rhodobacter sphaeroides",
#"S. aureus",
"Escherichia coli",
"Homo sapiens",
"Mus musculus")

x_breaks <- c(.02,.03,.04,.05,.06,.07,.08,.09,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0,2,3,4,5,6,7,8,9,10)
names = c("2B", "R0", "R1","R3","R5","R7","P64","P128","P256","Q6","Q10","Q30","Q100")
pch_comp = c(1,5,5,5,5,5,2,2,2,0,0,0,0)
pch = c(16,18,18,18,18,18,17,17,17,15,15,15,15)
y_breaks <- c(2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300)

dfr_tables <- c()
original_tables <- c()
counter <- 1

for (table in tables) {
  organism <- read.table(table, row.names = 1, header = T)
  remove = c("original", "maxqual", "minqual")
# remove = c("maxqual", "minqual")
  original_value <- organism[rownames(organism) %in% "original", ]
  organism <- organism[!rownames(organism) %in% remove, ]
  
  y <- organism$mse
  x_bzip2 <- organism$bits_bp
 
#   x_bzip2 <- apply(cbind(organism$comp_bzip2, organism$orig_bzip2), 1, min, na.rm=TRUE)
#   x_bzip2 <- (x_bzip2 * 8) / organism$bases
#   
#   x_comp <- organism$comp_size
#   x_comp <- (x_comp * 8) / organism$bases
#   
#   dfr <- data.frame(x = x_comp, y = y, pch = pch_comp, names=names, dataset=dataset_names[[counter]])
#   dfr_tables <- rbind(dfr_tables,dfr)
  
  dfr <- data.frame(x = x_bzip2, y = y, pch = pch, names=names, dataset=dataset_names[[counter]])
  dfr_tables <- rbind(dfr_tables, dfr)

  original_tables <- rbind(original_tables, data.frame(x = original_value$bits_bp, y = original_value$mse, dataset=dataset_names[[counter]]))

  
  counter <- counter + 1
}

pdf("compression_results.pdf", width=10, height=5)

# Hack to print the asterisks on the axis.
library(grid)
gglabcol <- 
  function(plot1) 
    
  {
    g <- ggplotGrob(plot1)
    
    # legend grobs
    g.b <- g[["grobs"]][[which(g$layout$name=="guide-box")]]
    l <- g.b[["grobs"]][[1]][["grobs"]]
    
    # get grobs for legend symbols (extract colour)
    lg <- l[sapply(l, function(i) grepl("GRID.text", i))]
    
    # get grobs for legend labels 
    lb <- l[sapply(l, function(i) grepl("guide.label", i))]
    
    # get change colour of labels to colour of symbols
    for(i in seq_along(lg)) {
      
      lb[[i]]$gp$col <-  lg[[i]]$gp$col
      
      g.b[["grobs"]][[1]][["grobs"]][sapply(g.b[["grobs"]][[1]][["grobs"]],
                                            function(i) grepl("guide.label", i))][[i]] <- lb[[i]]
    }
    
    # overwrite original legend
    g[["grobs"]][[which(g$layout$name=="guide-box")]] <- g.b
    
    grid.draw(g)
    
    invisible(g)
  }

### End Hack

p <- ggplot(dfr_tables, aes(x = x, y = y, label=names, col=dataset, shape=pch)) + 
  scale_shape_identity() +
  scale_y_log10(breaks = y_breaks, labels = y_breaks, limits = c(5,190)) +
  scale_x_log10(breaks = x_breaks, labels = x_breaks, limits = c(.027,11)) +
  geom_point(size=5, data = original_tables, aes(x = x, y = y, pch=c(8,8,8,8), label=c("","","",""))) +
  geom_point(data = data.frame(x=c(8),y=c(0)), size=5, aes(x=x, y=y, shape=c(8), label=c("")), color="black") +
  xlab("bits/base-call") + 
  ylab("mean squared error") +
  geom_point(size=2) + 
  geom_text(hjust=-.3, vjust=-0.3, size=2, show_guide = FALSE) + 
  theme_bw() + theme(legend.position="bottom", legend.key = element_blank(), legend.title=element_blank(),
                     legend.text = element_text(face="italic", size=8),
                     axis.text.x = element_text(size=6, angle=30), axis.text.y = element_text(size=6),
                     axis.title.x = element_text(size=8), axis.title.y = element_text(size=8))

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)

dev.off()


