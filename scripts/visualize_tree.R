suppressMessages(library(ggplot2));
suppressMessages(library(ggtree));
suppressMessages(library(treeio));

nwk <- read.newick(snakemake@input[[1]]);
labels <- read.csv(snakemake@input[[2]]);
png <- snakemake@output[[1]];
palette <- snakemake@params$palette;
color_by <- snakemake@params$color_by;


t_labelled = rename_taxa(nwk, labels, X, {{color_by}});

p <- ggtree(t_labelled, size=.5, color='grey'); 

tree <- p +
        scale_color_manual(values=palette) +
        geom_tippoint(aes(color=label), size=3) +
        scale_x_reverse() + 
        coord_flip() +
        theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
        theme(legend.position="none");

ggsave(tree, 
       file = png, 
       height = 9,
       width = 22) 
