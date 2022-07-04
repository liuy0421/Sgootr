suppressMessages(library(ggplot2));
suppressMessages(library(ggtree));
suppressMessages(library(treeio));
suppressMessages(library(reticulate));
np <- import('numpy');

nwk <- read.newick(snakemake@input[[1]]);
cells <- np$load(snakemake@input[[2]])$f[["rows"]];
png <- snakemake@output[[1]];
palette <- snakemake@params$palette;

d <- data.frame(label = 0:(length(cells)-1), 
                cells = cells, 
                lesions = substr(cells,1,2), 
                sampling_locations = sapply(strsplit(cells, '_'), "[[", 1));

t_lesions = rename_taxa(nwk, d, label, lesions);
t_sampling_locations = rename_taxa(nwk, d, label, sampling_locations);

p <- ggtree(t_lesions, size=.5, color='grey'); 
p <- p %<+% d;

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
