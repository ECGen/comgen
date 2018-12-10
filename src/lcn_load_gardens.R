###LCN: loading the onc data
###MKLau
###21Mar2014

### Data objects:
## onc.com = "community" occurrences summed across all cells for each tree
## onc.q = occurrence matrices separated out for each tree
## onc.geno = genotypes
## prb.onc = percent rough bark (averaged between the upper and lower)

### Data notes:
## Trees were removed from the analysis genotype RL6 and N1.31
## No physciods
## Lecanoras merged


pkg.list <- c("vegan", "ecodist", "bipartite", "RColorBrewer", "enaR", "devtools")
                                        # Install packages that are not installed
if (any(!(pkg.list %in% installed.packages()[, 1]))){
    sapply(pkg.list[which(!(pkg.list %in% installed.packages()[, 1]))], 
           install.packages, dependencies = TRUE, repos='http://cran.us.r-project.org')
}
                                        # Check for ComGenR
if (!("ComGenR" %in% installed.packages()[, 1])){
    devtools::install_github("CommunityGeneticsAnalyses/ComGenR")
}
                                        # Library packages
sapply(c(pkg.list, "ComGenR"), library, quietly = TRUE, character.only = TRUE)
                                        # Loading some misc helper functions
source('../bin/helpers.R')
                                        # Loading data
xgal.size <- gdata::read.xls("../data/lichen_networks/ONC_Xgal_SizeData_May2011.xlsx")
garden.data <- read.csv('../data/lichen_networks/LCO_data_ONC_PIT.csv')
                                        # remove genotype RL6 and N1.31
garden.data <- garden.data[garden.data$Geno!='RL6',]
garden.data <- garden.data[garden.data$Tree!='N1.31',]
                                        #separate onc
garden.data[,1] <- as.character(garden.data[,1])
g1 <- substr(garden.data[,1],2,2)
g1[g1!='P'] <- 'onc'
onc <- garden.data[g1 == 'onc',]
pit <- garden.data[g1 == 'P',]
					#tree overlap between years
unique(onc$Tree[onc$Year=='2010']) %in% unique(onc$Tree[onc$Year=='2011'])
unique(onc$Tree[onc$Year=='2011']) %in% unique(onc$Tree[onc$Year=='2010'])
                                        # Checking the data
if (!(all(table(onc[,1])==100))){for (i in 1:1000){print('Warning: check input data!!!')}}
                                        # Separate trees
                                        # onc
colnames(onc)[7:ncol(onc)] <- substr(colnames(onc)[7:ncol(onc)],1,2)
onc.q <- split(onc,paste(onc[,1],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])
                                        # pit
colnames(pit)[7:ncol(pit)] <- substr(colnames(pit)[7:ncol(pit)],1,2)
pit.q <- split(pit,paste(pit[,1],pit[,2]))
pit.q <- lapply(pit.q,function(x) x[,7:ncol(x)])
                                        # Get genotype
onc.geno <- unlist(sapply(names(onc.q),function(x) strsplit(x,split=' ')[[1]][2]))
pit.geno <- unlist(sapply(names(pit.q),function(x) strsplit(x,split=' ')[[1]][2]))

                                        # Xgal size data
xgs <- xgal.size[-1:-4, -(ncol(xgal.size) - 1):-ncol(xgal.size)]
xgs.cols <- xgal.size[4, -(ncol(xgal.size) - 1):-ncol(xgal.size)]
colnames(xgs) <- gsub("\\#", "", as.character(unlist(xgs.cols)))
xgs <- xgs[, 1:13]
xgs <- apply(xgs, 2, gsub, pattern = "\\,", replacement = "") 
xgs.dim <- xgs[, "Measurement"]
xgs.geno <- xgs[, "Genotype"]
xgs.tree <- xgs[, "Tree"]
xgs <- xgs[, grep("Thallus", colnames(xgs))] 
                                        # Coercing to numeric
xgs <- apply(xgs, 2, as.numeric) 
                                        # Dealing with NA values
xgs.geno <- xgs.geno[grep("Dimension", xgs.dim)]
xgs.tree <- xgs.tree[grep("Dimension", xgs.dim)]
xgs <- xgs[grep("Dimension", xgs.dim), ]
xgs.dim <- xgs.dim[grep("Dimension", xgs.dim)]
                                        # Convert to cm
xgs <- xgs * 0.1
xgs.ellipse <- pi * xgs[xgs.dim == "Dimension 1", ] * xgs[xgs.dim == "Dimension 2", ] 
xgs.geno <- xgs.geno[xgs.dim == "Dimension 1"]
xgs.tree <- xgs.tree[xgs.dim == "Dimension 1"]
                                        # package all xgs related data
xgs.data <- data.frame(tree = xgs.tree, geno = xgs.geno, 
                       mean.thallus = apply(xgs.ellipse, 1, mean, na.rm = TRUE),
                       median.thallus = apply(xgs.ellipse, 1, median, na.rm = TRUE),
                       xgs.ellipse)
                                        # remove trees not done (i.e. all NA)
xgs.data <- xgs.data[apply(xgs.data[, grep("Thallus", colnames(xgs.data))], 1, function(x) !(all(is.na(x)))),]

                                        # Roughness in the Garden
rough <- read.csv('../data/lichen_networks/ONC_raw_roughness.csv')
                                        # Isolate roughness
rough <- rough[, 1:5]
                                        # Isolate north quadrats
rough <- rough[grepl("North", rough[,3]), ]
                                        # Average roughness
avg.rough <- tapply(rough[,5], rough[,1], mean)
r.tree <- names(avg.rough)
r.tree <- sub('-', '\\.', r.tree)
r.tree <- sub('\\.0', '\\.', r.tree)
names(avg.rough) <- r.tree
                                        # match roughness to to ses values
load('../data/lichen_networks/lcn_onc_ses.rda')
onc.ses <- unlist(os[,1])
onc.ses[is.na(onc.ses)] <- 0
names(onc.ses) <- rownames(os)
if (!(all(names(onc.ses)==names(onc.q)))){print('Holy crap!')}
ses.tree <- as.character(sapply(names(onc.ses),function(x) unlist(strsplit(x,split=' '))[1]))
onc.rough <- avg.rough[match(ses.tree, r.tree)]
if (!(all(ses.tree==names(onc.rough)))){print('Holy Crap!')}
                                        #RFLP distance values from Zink from Martinsen
rflp.d <- readLines('../data/acn/rflp_queller_goodnight.txt')
rflp.d <- lapply(rflp.d,strsplit,split='\t')
rflp.d <- lapply(rflp.d,function(x) x[[1]])
rflp.d[[61]] <- c(rflp.d[[61]],"")
rflp.d <- do.call(rbind,rflp.d)
rflp.n <- rflp.d[1,-1]
rflp.d <- rflp.d[-1,-1]
diag(rflp.d) <- 1
rflp.d <- matrix(as.numeric(rflp.d),nrow=nrow(rflp.d))
rownames(rflp.d) <- colnames(rflp.d) <- rflp.n
rflp.d <- rflp.d[rownames(rflp.d)%in%unique(onc.geno),colnames(rflp.d)%in%unique(onc.geno)]
rflp.d <- rflp.d[match(unique(onc.geno),rownames(rflp.d)),match(unique(onc.geno),rownames(rflp.d))]
if (!(all(rownames(rflp.d)==unique(onc.geno)))){print('Holy crap, rflp.d names match error')}
rflp.d <- as.dist(rflp.d)
                                        # Lichen Network Models
                                        # onc
cn.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), coNets, ci.p = 95)
cn.sign.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), coNets, ci.p = 95, return.sign = TRUE)
cn.d.onc <- netDist(cn.onc, method = "bc")
fn.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), freqNet, zero.diag = TRUE)
fn.d.onc <- netDist(fn.onc, method = "bc")
                                        # pit
cn.pit <- lapply(split(pit[, -1:-6], pit[, "Tree"]), coNets, ci.p = 95)
cn.sign.pit <- lapply(split(pit[, -1:-6], pit[, "Tree"]), coNets, ci.p = 95, return.sign = TRUE)
cn.d.pit <- netDist(cn.pit, method = "bc")
fn.pit <- lapply(split(pit[, -1:-6], pit[, "Tree"]), freqNet, zero.diag = TRUE)
fn.d.pit <- netDist(fn.pit, method = "bc")
                                        # genotype means and mean distances
onc.tree <- do.call(rbind, strsplit(names(onc.geno), " "))[, 1]
cn.mu.onc <- list()
for (i in 1:length(unique(onc.geno))){
    cn.mu.onc[[i]] <- netMean(cn.onc[onc.geno == unique(onc.geno)[i]])
}
cn.mu.d.onc <- netDist(cn.mu.onc, method = "bc")
names(cn.mu.onc) <- unique(onc.geno)
                                        # network statistics
ns.onc <- lapply(cn.onc, pack)
ns.onc <- lapply(ns.onc, enaR::enaStructure)
ns.onc <- do.call(rbind, lapply(ns.onc, function(x) x[[2]]))
                                        # graph level centralization 
dcen.onc <- unlist(lapply(cn.onc, function(x) 
    sna::centralization(x, FUN = sna::degree, normalize = TRUE)))
onc.ns <- cbind(ns.onc, Cen = dcen.onc)
if (!(all(onc.tree == names(cn.onc)))){print("Danger Will Robinson!")}
                                        # species centralities
cen.spp <- lapply(cn.onc, sna::degree, rescale = FALSE)
cen.spp <- do.call(rbind, cen.spp)
colnames(cen.spp) <- colnames(cn.onc[[1]])
                                        # Community data
onc.com <- do.call(rbind,lapply(onc.q,function(x) apply(x,2,sum)))
onc.R <- apply(sign(onc.com),1,sum)
onc.H <- vegan::diversity(onc.com)
onc.com.gm <- apply(onc.com, 2, function(x, g) tapply(x, g, mean), g = onc.geno)
onc.com.gm.rel <- apply(onc.com.gm, 2, function(x) x/max(x))
onc.com.rel <- apply(onc.com, 2, function(x) x/max(x))
onc.com.rel <- cbind(onc.com.rel, ds = rep(min(onc.com.rel[onc.com.rel != 0]) / 1000, nrow(onc.com.rel)))
onc.com <- cbind(onc.com, ds = rep(min(onc.com[onc.com != 0]) / 1000, nrow(onc.com)))
                                        # pit genotype mean community
pit.com <- do.call(rbind,lapply(pit.q,function(x) apply(x,2,sum)))
pit.com.gm <- apply(pit.com, 2, function(x, g) tapply(x, g, mean), g = pit.geno)
pit.com.gm.rel <- apply(pit.com.gm, 2, function(x) x/max(x))
pit.com.gm.rel[is.na(pit.com.gm.rel)] <- 0
                                        # Lichen community metrics
                                        # Percent Total Cover
ptc.onc <- unlist(lapply(onc.q, function(x) sum(apply(x, 1, function(x) sign(sum(x))))))
                                        # Species richness
spr.onc <- apply(onc.com[, colnames(onc.com) != "ds"], 1, function(x) sum(sign(x)))
                                        # Vectors for network similarity
## ns.vec.onc <- envfit(ord, data.frame(onc.ns[, c("L", "Cen")], R = onc.rough, Cov = ptc.onc))
                                        # Bipartite analysis
nperm <- 999
if (!(file.exists("../data/lichen_networks/nest_rel_onc.rda"))){
    nest.onc <- nestedness(onc.com.rel[, colnames(onc.com.rel) != "ds"], n.nulls = 999)
    dput(nest.onc, "../data/lichen_networks/nest_rel_onc.rda")
}else{
    nest.onc <- dget("../data/lichen_networks/nest_rel_onc.rda")
}
if (!(file.exists("../data/lichen_networks/null_mod_onc.csv"))){
    obs.mod.onc <- bipartite::computeModules(onc.com.rel[, colnames(onc.com.rel) != "ds"])
    mods.onc <- tail(apply(slot(obs.mod.onc, "modules"), 2, 
                           function(x) sum(sign(x[2:length(x)]) *
                                               (1:(length(x) - 1)))),
                     sum(dim(onc.com[, colnames(onc.com) != "ds"])))
    mods.onc <- list(sp = tail(mods.onc, ncol(onc.com[, colnames(onc.com) != "ds"])), 
                     tree = head(mods.onc, nrow(onc.com)))
    sim.onc <- lapply(1:nperm, sim.com, x = onc.q)
    sim.onc <- lapply(sim.onc, function(x) x / max(x))
    nul.mod.onc <- lapply(sim.onc, function(x) bipartite::computeModules(x))
    nul.mod.onc <- unlist(lapply(nul.mod.onc, slot, "likelihood"))
    dput(mods.onc, "../data/lichen_networks/mod_list_onc.rda")
    write.csv(slot(obs.mod.onc, "likelihood"), 
              "../data/lichen_networks/obs_mod_onc.csv", 
              row.names = FALSE)
    write.csv(nul.mod.onc, 
              "../data/lichen_networks/null_mod_onc.csv", 
              row.names = FALSE)
}else{
    obs.mod.onc <- read.csv("../data/lichen_networks/obs_mod_onc.csv")[1]
    nul.mod.onc <- read.csv("../data/lichen_networks/null_mod_onc.csv")[,1]
    z.mod.onc <- (obs.mod.onc - mean(nul.mod.onc)) / sd(nul.mod.onc)
    mods.onc <- dget("../data/lichen_networks/mod_list_onc.rda")
}
pval.mod.onc <- length(nul.mod.onc[nul.mod.onc >= obs.mod.onc]) / length(nul.mod.onc)
if (pval.mod.onc == 0){pval.mod.onc <- 1/nperm}
z.mod.onc <- (obs.mod.onc - mean(nul.mod.onc)) / sd(nul.mod.onc)
bp.mod.onc <- c(nperm = nperm, obs = obs.mod.onc, z = z.mod.onc, pval = pval.mod.onc)
                                        # NMDS ordinations
                                        # community ordination
if (!file.exists("../data/lichen_networks/onc_nmds.csv")){
    nms.info.onc <- capture.output(nms.onc <- nmds.min(nmds(
        vegdist(onc.com.rel), 2, 2)))
    write.csv(nms.onc, "../data/lichen_networks/onc_nmds.csv", row.names = FALSE)
    write.table(nms.info.onc, 
                "../data/lichen_networks/onc_nmds_info.txt", 
                col.names = FALSE, row.names = FALSE)
}else{nms.onc <- read.csv("../data/lichen_networks/onc_nmds.csv")}
                                        # Network ordination
if (!(file.exists("../data/lichen_networks/conet_nmds.csv"))){
    cn.nmds.stats.onc <- capture.output(cn.nms.onc <- nmds.min(nmds(cn.d.onc, 2, 2)))
    write.csv(cn.nms.onc, "../data/lichen_networks/conet_nmds.csv", row.names = FALSE)
    write.table(cn.nmds.stats.onc, 
                "../data/lichen_networks/conet_nmds_info.txt", 
                col.names = FALSE, row.names = FALSE)
}else{cn.nms.onc <- read.csv("../data/lichen_networks/conet_nmds.csv")}
                                        # Vector fitting
nv.onc <- envfit(cn.nms.onc, data.frame(onc.com[, colnames(onc.com) != 'ds'], 
                                        R = onc.rough, 
                                        C = onc.ns[, c("C")], 
                                        A = ptc.onc))
cv.onc <- envfit(nms.onc, data.frame(onc.com[, colnames(onc.com) != 'ds'], 
                                        R = onc.rough, 
                                        C = onc.ns[, c("C")], 
                                     A = ptc.onc))
                                        #genotype means
omu <- apply(onc.com[,colnames(onc.com) != 'ds'], 2, 
             function(x,g) tapply(x,g,mean),g=onc.geno)
oms <- tapply(onc.ses, onc.geno, mean)
oms.d <- dist(oms[match(rownames(as.matrix(rflp.d)),names(oms))])
                                        #bark roughness means
oprbmu <- tapply(onc.rough,onc.geno,mean)
oprbmu <- oprbmu[match(rownames(as.matrix(rflp.d)),names(oprbmu))]
                                        #get araujo coordinates
coord <- read.csv('../data/lichen_networks/lcn_coord_onc.csv')
rownames(coord) <- coord[,1]
coord <- coord[,-1]
                                        # packing into a dataframe
tree <- onc.geno
for (i in 1:length(unique(onc.geno))){
    tree[onc.geno == unique(onc.geno)[i]] <- 1:length(tree[onc.geno == unique(onc.geno)[i]])
}
tree <- factor(tree)
onc.dat <- data.frame(ptc.onc, spr.onc, 
                      geno = factor(onc.geno), tree = tree, 
                      onc.rough, onc.ns[, c("L", "Cen")])
                                        # mean bark roughness calculations
prb.mu.onc <- tapply(onc.rough, onc.geno, mean)
prb.mu.d.onc <- dist(prb.mu.onc)
                                        # Plot calculations
pw.onc <- onc.com.rel[, colnames(onc.com.rel) != "ds"]
pw.onc <- pw.onc[order(apply(pw.onc, 1, sum), decreasing = TRUE), 
                 order(apply(pw.onc, 2, sum), decreasing = TRUE)]
rownames(pw.onc) <- onc.geno
col.pal <- RColorBrewer::brewer.pal((max(unlist(mods.onc))), "Paired")

