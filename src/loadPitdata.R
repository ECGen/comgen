### script for loading data for the
### ACN study
### Pit data collected 2012


pit <- read.csv('../data/acn/arth_cooc_PIT_Lau.csv')
pit[is.na(pit)] <- 0

### separating out "fungal" necrosis
### necrosis <- pit$fungal
pit <- pit[,colnames(pit) != 'fungal']

## Remove 1007
## Has high variability
remove.1007 <- FALSE
if (remove.1007){pit <- pit[pit[, "geno"] != "1007", ]}

## Remove mite
pit <- pit[, colnames(pit) != "mite"]

## combine pb 
pb <- pit$pb.upper + pit$pb.lower + pit$pb.woody 
pb.pred <- pit$pb.pred + pit$pb.woody.pred + pit$pb.hole
pb.abort <- pit$pb.abort
pit <- pit[,!(grepl('pb',colnames(pit)))]
pit <- data.frame(pit,pb.abort,pb.pred,pb)
pit <- data.frame(pit[,1:6],pit[,ncol(pit):7])
tree.info <- paste(pit[, "tree"], pit[, "geno"], pit[, "leaf.type"])
tree.arth <- pit[, 7:ncol(pit)]
tree.arth <- split(tree.arth, tree.info)

## community matrix
com.acn <- do.call(rbind, lapply(tree.arth, function(x) apply(x, 2, sum)))

## tree level networks for arthropods
cn.acn <- lapply(tree.arth, coNets, ci.p = 95, cond = TRUE)

## tree net distances
d.cn.acn <- netDist(cn.acn)

## network stats
l.cn.acn <- do.call(rbind, lapply(cn.acn, enaR:::structure.statistics))[, "L"]
cen.cn.acn <- unlist(lapply(cn.acn, 
                            function(x) 
                                sna::centralization(x, FUN = sna::degree, normalize = FALSE)), 
                     )
nm.cn.acn <- data.frame(L = l.cn.acn, C = cen.cn.acn)

## tree info
acn.dat <- data.frame(do.call(rbind, strsplit(names(cn.acn), split = " ")))
names(acn.dat) <- c("tree", "geno", "leaf.type")

## Network Ordination
set.seed(1234)
nms.cn.acn <- nmds(d.cn.acn, 2, 2)
ord.cn.acn <- nmds.min(nms.cn.acn)
vec.com.acn <- envfit(ord.cn.acn, com.acn[, apply(com.acn, 2, sum) > 10])
vec.nm.acn <- envfit(ord.cn.acn, nm.cn.acn)

## ## Modularity of bipartite networks
## computeModules(com.acn[grepl("live", rownames(com.acn)), ])
## computeModules(com.acn[grepl("sen", rownames(com.acn)), ])
