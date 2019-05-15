dc.lcm <- data.frame(Response = c(rep("Abundance", length(abun) / 2),
                         rep("Richness", length(abun) / 2),
                         rep("Diversity", length(abun) / 2)), 
                     Difference = c(tapply(abun, l.dat[, "Tree.pairs"], diff),
                         tapply(rich, l.dat[, "Tree.pairs"], diff),
                         tapply(shan, l.dat[, "Tree.pairs"], diff)))

pdf("../results/rln_dotchart.pdf")
dotplot(Response ~ Difference, data = dc.lcm)
dev.off()
