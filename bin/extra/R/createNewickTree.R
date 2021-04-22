library(ape)
library(reticulate)

args <- commandArgs(trailingOnly = T)

inputFile = as.character(args[1])
outputFile = as.character(args[2])
  
newick = system(paste0("python bin/extra/python/getNewickTree.py --in ", args[1]), intern = TRUE)


mytree <- read.tree(text=newick)


pdf(file=outputFile,
    width=8,height=6, paper="special",bg="white",
    fonts="Helvetica", colormodel="cmyk", pointsize = .5)

plot(mytree,no.margin = T,direction = "downwards")
nodelabels(text=mytree$node.label,node=1:mytree$Nnode+Ntip(mytree),frame="circle")
add.scale.bar(6)
edgelabels(mytree$edge.length, bg="white", col="black")

dev.off()