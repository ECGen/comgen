### Loading data objects
### Data notes:
# Loading data
# separate onc
xgal.size.in = read.csv("../data/lcn/ONC_Xgal_SizeData_May2011.csv")
garden.data.in = read.csv("../data/lcn/LCO_data_ONC_PIT.csv")
rough.in = read.csv("../data/lcn/ONC_raw_roughness.csv")
onc.nc.in = read.csv("../data/lcn/ONC_phytochem_NC.csv")
onc.tan.in = read.csv("../data/lcn/ONC_phytochem_tannin.csv")
onc.ph.in = read.csv("../data/lcn/ONC_Bark_lichen_pH_data.csv")
wild.dat.in = read.csv("../data/lcn/lco_Apr2012.csv")
env.in = read.csv("../data/lcn/Uinta2012_all_data_from_Lamit.csv")
age.in = read.csv(
    "../data/lcn/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv")

proc_garden_data <- function(garden.data) {
                                        # rm genotype RL6 and N1.31
    garden.data <- garden.data[garden.data$Geno != "RL6", ]
    garden.data <- garden.data[garden.data$Tree != "N1.31", ]
    garden.data[, 1] <- as.character(garden.data[, 1])
    return(garden.data)

}


proc_pit <- function(garden.data) {

    g1 <- substr(garden.data[, 1], 2, 2)
    g1[g1 != "P"] <- "onc"
    pit <- garden.data[g1 == "P", ]
    colnames(pit)[7:ncol(pit)] <- substr(
        colnames(pit)[7:ncol(pit)], 1, 2)
    return(pit)

}

proc_onc <- function(garden.data) {

    g1 <- substr(garden.data[, 1], 2, 2)
    g1[g1 != "P"] <- "onc"
    onc <- garden.data[g1 == "onc", ]
    colnames(onc)[which(colnames(onc) == "Ls")] <- "Lh"
    colnames(onc)[7:ncol(onc)] <- substr(colnames(onc)[7:ncol(onc)], 1, 2)
    return(onc)

}

proc_onc_q <- function(onc) {
    onc.q <- split(onc, paste(onc[, 1], onc[, 2]))
    onc.q <- lapply(onc.q, function(x) x[, 7:ncol(x)])
    return(onc.q)
}

proc_onc_ph <- function(garden.data, rough.in, onc, onc.q,
                        onc.nc, onc.tan, onc.ph) {

onc.geno <- unlist(
    sapply(names(onc.q), 
           function(x) strsplit(x, split = " ")[[1]][2])
)
onc.tree <- do.call(rbind, strsplit(names(onc.geno), " "))[, 1]
                                        # roughness
rough <- rough.in[, 1:5]
rough <- rough[grepl("North", rough[, 3]), ]
avg.rough <- tapply(rough[, 5], rough[, 1], mean)
r.tree <- names(avg.rough)
r.tree <- sub("-", "\\.", r.tree)
r.tree <- sub("\\.0", "\\.", r.tree)
names(avg.rough) <- r.tree
onc.rough <- avg.rough[match(onc.tree, r.tree)]
                                        # network models
onc.split <- split(onc[, -1:-6], onc[, "Tree"], drop = TRUE)
cn.onc <- lapply(onc.split, coNet, ci.p = 95)
                                        # modularity
ns.onc <- lapply(lapply(cn.onc, function(x) {
    abs(sign(x))
}), enaR:::structure.statistics)
ns.onc <- do.call(rbind, ns.onc)
cn.mod.onc <- matrix(nrow = length(cn.onc), ncol = 2)
for (i in 1:length(cn.onc)) {
    if (sum(sign(cn.onc[[i]])) >= 3) {
        mod.tmp <- computeModules(cn.onc[[i]])
        cn.mod.onc[i, 1] <- slot(mod.tmp, "likelihood")
        cn.mod.onc[i, 2] <- nrow(slot(mod.tmp, "modules")) - 1
    } else {
        cn.mod.onc[i] <- NA
    }
}
cn.mod.onc[is.na(cn.mod.onc)] <- 0
names(cn.mod.onc) <- c("mod.lik", "mod.n")
                                        # centrality
dcen.onc <- unlist(lapply(cn.onc, function(x) {
    sna::centralization(x, FUN = sna::degree, normalize = FALSE)
}))
onc.ns <- cbind(ns.onc, Cen = dcen.onc, mod.lik = cn.mod.onc[, 1], mod.n = cn.mod.onc[, 
    2])
onc.com <- do.call(rbind, lapply(onc.q, function(x) apply(x, 2, sum)))
onc.com <- cbind(onc.com, ds = rep(min(onc.com[onc.com != 0])/1000, nrow(onc.com)))
ptc.onc <- unlist(lapply(onc.q, function(x) {
    sum(apply(x, 1, function(x) sign(sum(x))))
}))
spr.onc <- apply(onc.com[, colnames(onc.com) != "ds"], 1, function(x) sum(sign(x)))
spd.onc <- diversity(onc.com[, colnames(onc.com) != "ds"])
spe.onc <- spd.onc/log(specnumber(onc.com[, colnames(onc.com) != "ds"]))
spe.onc[is.na(spe.onc)] <- 0
tree <- onc.geno
for (i in 1:length(unique(onc.geno))) {
    tree[onc.geno == unique(onc.geno)[i]] <- 1:length(tree[onc.geno == unique(onc.geno)[i]])
}
tree <- factor(tree)
tree.id <- do.call(rbind, strsplit(names(ptc.onc), split = " "))[, 1]
onc.nc[, 1] <- as.character(paste0("N", gsub("-", "\\.", onc.nc[, 1])))
onc.tan[, 1] <- as.character(paste0("N", gsub("-", "\\.", onc.tan[, 1])))
colnames(onc.nc)[1:4] <- c("tree.id", "sample.mass", "N", "C")
colnames(onc.tan)[1] <- "tree.id"
colnames(onc.tan)[grep("X.CT", colnames(onc.tan))] <- "CT"
onc.nc$rCN <- onc.nc$N/onc.nc$C
onc.ph[, "tree.id"] <- gsub("-", ".", onc.ph[, "tree.id"])
onc.ph[, "tree.id"] <- gsub("\\.0", "\\.", onc.ph[, "tree.id"])
onc.ph[onc.ph[, "tree.id"] == "N7.16", "tree.id"] <- "N7.10"
onc.ph[!is.na(onc.ph[, "pH2"]), "pH"] <- apply(onc.ph[!is.na(onc.ph[, "pH2"]), c("pH", 
    "pH2")], 1, mean)
onc.dat <- data.frame(tree.id, PC = ptc.onc, SR = spr.onc, SD = spd.onc, SE = spe.onc, 
    geno = factor(onc.geno), tree = tree, BR = onc.rough, onc.ns[, c("L", "Cen")])
onc.ph <- onc.ph[onc.ph[, "tree.id"] %in% onc.dat[, "tree.id"], ]
onc.ph <- onc.ph[match(onc.dat[, "tree.id"], onc.ph[, "tree.id"]), ]
return(onc.ph)

}

proc_onc_dat <- function(garden.data, rough, onc, onc.q,
                         onc.nc, onc.tan, onc.ph) {

onc.geno <- unlist(
    sapply(names(onc.q), 
           function(x) strsplit(x, split = " ")[[1]][2])
)
onc.tree <- do.call(rbind, strsplit(names(onc.geno), " "))[, 1]
                                        # roughness
rough <- rough.in[, 1:5]
rough <- rough[grepl("North", rough[, 3]), ]
avg.rough <- tapply(rough[, 5], rough[, 1], mean)
r.tree <- names(avg.rough)
r.tree <- sub("-", "\\.", r.tree)
r.tree <- sub("\\.0", "\\.", r.tree)
names(avg.rough) <- r.tree
onc.rough <- avg.rough[match(onc.tree, r.tree)]
                                        # network models
onc.split <- split(onc[, -1:-6], onc[, "Tree"], drop = TRUE)
cn.onc <- lapply(onc.split, coNet, ci.p = 95)
                                        # modularity
ns.onc <- lapply(lapply(cn.onc, function(x) {
    abs(sign(x))
}), enaR:::structure.statistics)
ns.onc <- do.call(rbind, ns.onc)
cn.mod.onc <- matrix(nrow = length(cn.onc), ncol = 2)
for (i in 1:length(cn.onc)) {
    if (sum(sign(cn.onc[[i]])) >= 3) {
        mod.tmp <- computeModules(cn.onc[[i]])
        cn.mod.onc[i, 1] <- slot(mod.tmp, "likelihood")
        cn.mod.onc[i, 2] <- nrow(slot(mod.tmp, "modules")) - 1
    } else {
        cn.mod.onc[i] <- NA
    }
}
cn.mod.onc[is.na(cn.mod.onc)] <- 0
names(cn.mod.onc) <- c("mod.lik", "mod.n")
                                        # centrality
dcen.onc <- unlist(lapply(cn.onc, function(x) {
    sna::centralization(x, FUN = sna::degree, normalize = FALSE)
}))
onc.ns <- cbind(ns.onc, Cen = dcen.onc, mod.lik = cn.mod.onc[, 1], mod.n = cn.mod.onc[, 
    2])
onc.com <- do.call(rbind, lapply(onc.q, function(x) apply(x, 2, sum)))
onc.com <- cbind(onc.com, 
                 ds = rep(min(onc.com[onc.com != 0])/1000, 
                          nrow(onc.com)))
ptc.onc <- unlist(lapply(onc.q, function(x) {
    sum(apply(x, 1, function(x) sign(sum(x))))
}))
spr.onc <- apply(onc.com[, colnames(onc.com) != "ds"], 1, 
                 function(x) sum(sign(x)))
spd.onc <- diversity(onc.com[, colnames(onc.com) != "ds"])
spe.onc <- spd.onc/log(
                       specnumber(
                           onc.com[, colnames(onc.com) != "ds"]))
spe.onc[is.na(spe.onc)] <- 0
tree <- onc.geno
for (i in 1:length(unique(onc.geno))) {
    tree[onc.geno == 
         unique(onc.geno)[i]] <- 1:length(tree[onc.geno == 
                                               unique(onc.geno)[i]])
}
tree <- factor(tree)
tree.id <- do.call(rbind, strsplit(names(ptc.onc), split = " "))[, 1]
onc.nc[, 1] <- as.character(paste0("N", gsub("-", "\\.", onc.nc[, 1])))
onc.tan[, 1] <- as.character(paste0("N", gsub("-", "\\.", onc.tan[, 1])))
colnames(onc.nc)[1:4] <- c("tree.id", "sample.mass", "N", "C")
colnames(onc.tan)[1] <- "tree.id"
colnames(onc.tan)[grep("X.CT", colnames(onc.tan))] <- "CT"
onc.nc$rCN <- onc.nc$N/onc.nc$C
onc.dat <- data.frame(tree.id, PC = ptc.onc, SR = spr.onc, SD = spd.onc, SE = spe.onc, 
    geno = factor(onc.geno), tree = tree, BR = onc.rough, onc.ns[, c("L", "Cen")])
onc.dat <- data.frame(onc.dat, C = onc.nc[match(onc.dat[, "tree.id"], onc.nc[, "tree.id"]), 
    "C"], N = onc.nc[match(onc.dat[, "tree.id"], onc.nc[, "tree.id"]), "N"], CN = onc.nc[match(onc.dat[, 
    "tree.id"], onc.nc[, "tree.id"]), "rCN"], CT = onc.tan[match(onc.dat[, "tree.id"], 
    onc.tan[, "tree.id"]), "CT"], pH = onc.ph[, "pH"])
return(onc.dat)

}

garden.data <- proc_garden_data(garden.data.in)
pit <- proc_pit(garden.data.in)
onc <- proc_onc(garden.data.in)
onc.q <- proc_onc_q(onc)
onc.ph <- proc_onc_ph(garden.data.in, 
                      rough.in, onc, onc.q, 
                      onc.nc.in, onc.tan.in, onc.ph.in)
onc.dat <- proc_onc_dat(garden.data.in, 
                        rough.in, onc, onc.q, 
                        onc.nc.in, onc.tan.in, onc.ph)



onc <- garden.data[g1 == "onc", ]
colnames(onc)[which(colnames(onc) == "Ls")] <- "Lh"
pit <- garden.data[g1 == "P", ]
# tree overlap between years
unique(onc$Tree[onc$Year == "2010"]) %in%
  unique(onc$Tree[onc$Year == "2011"])
unique(onc$Tree[onc$Year == "2011"]) %in%
  unique(onc$Tree[onc$Year == "2010"])
# Checking the data
if (!(all(table(onc[, 1]) == 100))) {
  for (i in 1:1000) {
    print("Warning: check input data!!!")
  }
}
# Separate trees
# onc
colnames(onc)[7:ncol(onc)] <- substr(colnames(onc)[7:ncol(onc)], 1, 2)
onc.q <- split(onc, paste(onc[, 1], onc[, 2]))
onc.q <- lapply(onc.q, function(wild.dat) x[, 7:ncol(x)])
# pit
colnames(pit)[7:ncol(pit)] <- substr(colnames(pit)[7:ncol(pit)], 1, 2)
pit.q <- split(pit, paste(pit[, 1], pit[, 2]))
pit.q <- lapply(pit.q, function(x) x[, 7:ncol(x)])
# Get genotype
onc.geno <- unlist(sapply(
  names(onc.q),
  function(x) strsplit(x, split = " ")[[1]][2]
))
pit.geno <- unlist(sapply(
  names(pit.q),
  function(x) strsplit(x, split = " ")[[1]][2]
))
# Xgal size data
xgs <- xgal.size[-1:-7, -(ncol(xgal.size) - 1):-ncol(xgal.size)]
xgs.cols <- xgal.size[7, -(ncol(xgal.size) - 1):-ncol(xgal.size)]
colnames(xgs) <- gsub("\\#", "", as.character(unlist(xgs.cols)))
xgs <- xgs[, 1:13]
xgs <- apply(xgs, 2, gsub, pattern = "\\,", replacement = "")
xgs.dim <- xgs[, "Measurement"]
xgs.geno <- xgs[, "Genotype"]
xgs.tree <- xgs[, "Tree"]
xgs <- xgs[, grep("Thallus", colnames(xgs))]
# fix genotypes
# t6
xgs.geno[grep("T6", xgs.geno)] <- "T6"
xgs.geno[grep("H10", xgs.geno)] <- "H-10"
# Coercing to numeric
xgs <- apply(xgs, 2, as.numeric)
# Dealing with NA values
xgs.geno <- xgs.geno[grep("Dimension", xgs.dim)]
xgs.tree <- xgs.tree[grep("Dimension", xgs.dim)]
xgs <- xgs[grep("Dimension", xgs.dim), ]
xgs.dim <- xgs.dim[grep("Dimension", xgs.dim)]
# Convert to cm
xgs <- xgs * 0.1
xgs.ellipse <- pi * xgs[xgs.dim == "Dimension 1", ] *
  xgs[xgs.dim == "Dimension 2", ]
xgs.geno <- xgs.geno[xgs.dim == "Dimension 1"]
xgs.tree <- xgs.tree[xgs.dim == "Dimension 1"]
# package all xgs related data
xgs.data <- data.frame(
  tree = xgs.tree, geno = xgs.geno,
  mean.thallus = apply(xgs.ellipse, 1,
    mean,
    na.rm = TRUE
  ),
  median.thallus = apply(xgs.ellipse, 1,
    median,
    na.rm = TRUE
  ),
  xgs.ellipse
)
# remove trees not done (i.e. all NA)
xgs.data <- xgs.data[apply(
  xgs.data[, grep(
    "Thallus",
    colnames(xgs.data)
  )], 1,
  function(x) !(all(is.na(x)))
), ]
# Isolate roughness
rough <- rough[, 1:5]
# Isolate north quadrats
rough <- rough[grepl("North", rough[, 3]), ]
# Average roughness
avg.rough <- tapply(rough[, 5], rough[, 1], mean)
r.tree <- names(avg.rough)
r.tree <- sub("-", "\\.", r.tree)
r.tree <- sub("\\.0", "\\.", r.tree)
names(avg.rough) <- r.tree
# match roughness to to ses values
load("../data/lcn/lcn_onc_ses.rda")
onc.ses <- unlist(os[, 1])
onc.ses[is.na(onc.ses)] <- 0
names(onc.ses) <- rownames(os)
if (!(all(names(onc.ses) == names(onc.q)))) {
  print("Holy crap!")
}
ses.tree <- as.character(sapply(
  names(onc.ses),
  function(x) unlist(strsplit(x, split = " "))[1]
))
onc.rough <- avg.rough[match(ses.tree, r.tree)]
if (!(all(ses.tree == names(onc.rough)))) {
  print("Holy Crap!")
}

# Lichen Network Models
# onc
cn.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), coNet,
  ci.p = 95
)
cn.sign.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), coNet,
  ci.p = 95
)
cn.d.onc <- distNet(cn.onc, method = "euclidean")
# pit
cn.pit <- lapply(
  split(pit[, -1:-6], pit[, "Tree"]), 
  coNet, ci.p = 95
)
cn.sign.pit <- lapply(
  split(pit[, -1:-6], pit[, "Tree"]), 
  coNet, ci.p = 95
)
cn.d.pit <- distNet(cn.pit, method = "bc")
# genotype means and mean distances
cn.mu.onc <- list()
for (i in 1:length(unique(onc.geno))) {
  cn.mu.onc[[i]] <- meanNet(cn.onc[onc.geno == unique(onc.geno)[i]])
}
names(cn.mu.onc) <- unique(onc.geno)
cn.mu.d.onc <- distNet(cn.mu.onc, method = "bc")
# mean bark roughness calculations
prb.mu.onc <- tapply(onc.rough, onc.geno, mean)
prb.mu.d.onc <- dist(prb.mu.onc)
# network statistics
ns.onc <- lapply(lapply(cn.onc, function(x) {
  abs(sign(x))
}), enaR:::structure.statistics)
ns.onc <- do.call(rbind, ns.onc)
# Ratio P / N
ns.rpn <- unlist(lapply(cn.onc, function(x) {
  mean(x[x > 0]) / mean(x[x < 0])
}))

# modularity
cn.mod.onc <- matrix(nrow = length(cn.onc), ncol = 2)
for (i in 1:length(cn.onc)) {
  if (sum(sign(cn.onc[[i]])) >= 3) {
    mod.tmp <- computeModules(cn.onc[[i]])
    cn.mod.onc[i, 1] <- slot(mod.tmp, "likelihood")
    cn.mod.onc[i, 2] <- nrow(slot(mod.tmp, "modules")) - 1
  } else {
    cn.mod.onc[i] <- NA
  }
}
cn.mod.onc[is.na(cn.mod.onc)] <- 0
names(cn.mod.onc) <- c("mod.lik", "mod.n")
# graph level centralization
dcen.onc <- unlist(lapply(cn.onc, function(x) {
  sna::centralization(x, FUN = sna::degree, normalize = FALSE)
}))
onc.ns <- cbind(
  ns.onc,
  Cen = dcen.onc,
  mod.lik = cn.mod.onc[, 1], mod.n = cn.mod.onc[, 2]
)
if (!(all(onc.tree == names(cn.onc)))) {
  print("Danger Will Robinson!")
}
# Community data
onc.com <- do.call(rbind, lapply(onc.q, function(x) apply(x, 2, sum)))
onc.R <- apply(sign(onc.com), 1, sum)
onc.H <- vegan::diversity(onc.com)
onc.com.gm <- apply(onc.com, 2,
  function(x, g) tapply(x, g, mean),
  g = onc.geno
)
onc.com.gm.rel <- apply(onc.com.gm, 2, function(x) x / max(x))
onc.com.rel <- apply(onc.com, 2, function(x) x / max(x))
onc.com.rel <- cbind(onc.com.rel,
  ds = rep(
    min(onc.com.rel[onc.com.rel != 0]) / 1000,
    nrow(onc.com.rel)
  )
)
onc.com <- cbind(onc.com,
  ds = rep(
    min(onc.com[onc.com != 0]) / 1000,
    nrow(onc.com)
  )
)
# pit genotype mean community
pit.com <- do.call(rbind, lapply(pit.q, function(x) apply(x, 2, sum)))
pit.com.gm <- apply(pit.com, 2,
  function(x, g) tapply(x, g, mean),
  g = pit.geno
)
pit.com.gm.rel <- apply(pit.com.gm, 2, function(x) x / max(x))
pit.com.gm.rel[is.na(pit.com.gm.rel)] <- 0
# Lichen community metrics
# Percent Total Cover
ptc.onc <- unlist(lapply(
  onc.q,
  function(x) {
    sum(apply(
      x, 1,
      function(x) sign(sum(x))
    ))
  }
))
# Species richness
spr.onc <- apply(
  onc.com[, colnames(onc.com) != "ds"], 1,
  function(x) sum(sign(x))
)
# Diversity
spd.onc <- diversity(onc.com[, colnames(onc.com) != "ds"])
# Evenness
spe.onc <- spd.onc / log(specnumber(onc.com[, colnames(onc.com) != "ds"]))
spe.onc[is.na(spe.onc)] <- 0

## "mean" distance matrices
cn.mu.d <- as.matrix(cn.mu.d.onc)
prb.mu.d <- as.matrix(prb.mu.d.onc)
prb.mu.d <- prb.mu.d[
  match(rownames(cn.mu.d), rownames(prb.mu.d)),
  match(rownames(cn.mu.d), rownames(prb.mu.d))
]
prb.mu.d <- as.dist(prb.mu.d)
onc.com.mu <- apply(onc.com[, -ncol(onc.com)], 2,
  function(x, g) tapply(x, g, mean),
  g = onc.geno
)
onc.com.mu <- onc.com.mu[match(rownames(cn.mu.d), 
                               rownames(onc.com.mu)), ]
onc.com.mu.d <- vegdist(onc.com.mu)
if (!(all(rownames(as.matrix(prb.mu.d)) == 
          rownames(as.matrix(cn.mu.d.onc))))) {
  warning("Warning: distance matrices are not aligned!")
} 

# Bipartite analysis
nperm <- 20
if (!(file.exists("../data/lcn/nest_rel_onc.rda"))) {
  nest.onc <- nestedness(onc.com.rel[, colnames(onc.com.rel) != "ds"],
    n.nulls = 999
  )
  dput(nest.onc, "../data/lcn/nest_rel_onc.rda")
} else {
  nest.onc <- dget("../data/lcn/nest_rel_onc.rda")
}
if (!(file.exists("../data/lcn/null_mod_onc.csv"))) {
  obs.mod.onc <- bipartite::computeModules(
    onc.com.rel[, colnames(onc.com.rel) != "ds"]
  )
  mods.onc <- tail(
    apply(
      slot(obs.mod.onc, "modules"), 2,
      function(x) {
        sum(sign(x[2:length(x)]) *
          (1:(length(x) - 1)))
      }
    ),
    sum(dim(onc.com[, colnames(onc.com) != "ds"]))
  )
  mods.onc <- list(
    sp = tail(
      mods.onc,
      ncol(onc.com[, colnames(onc.com) != "ds"])
    ),
    tree = head(mods.onc, nrow(onc.com))
  )
  sim.onc <- lapply(1:nperm, sim.com, x = onc.q)
  sim.onc <- lapply(sim.onc, function(x) x / max(x))
  nul.mod.onc <- lapply(sim.onc, function(x) bipartite::computeModules(x))
  nul.mod.onc <- unlist(lapply(nul.mod.onc, slot, "likelihood"))
  dput(mods.onc, "../data/lcn/mod_list_onc.rda")
  write.csv(slot(obs.mod.onc, "likelihood"),
    "../data/lcn/obs_mod_onc.csv",
    row.names = FALSE
  )
  write.csv(nul.mod.onc,
    "../data/lcn/null_mod_onc.csv",
    row.names = FALSE
  )
} else {
  obs.mod.onc <- read.csv("../data/lcn/obs_mod_onc.csv")[1]
  nul.mod.onc <- read.csv("../data/lcn/null_mod_onc.csv")[, 1]
  z.mod.onc <- (obs.mod.onc - mean(nul.mod.onc)) / sd(nul.mod.onc)
  mods.onc <- dget("../data/lcn/mod_list_onc.rda")
}


## NMDS ordinations
# community ordination
if (!file.exists("../data/lcn/onc_nmds.csv")) {
  nms.info.onc <- capture.output(nms.onc <- nmds.min(nmds(
    vegdist(onc.com.rel), 2, 2
  )))
  write.csv(nms.onc, "../data/lcn/onc_nmds.csv",
    row.names = FALSE
  )
  write.table(nms.info.onc,
    "../data/lcn/onc_nmds_info.txt",
    col.names = FALSE, row.names = FALSE
  )
} else {
  nms.onc <- read.csv("../data/lcn/onc_nmds.csv")
}
# Network ordination
if (!(file.exists("../data/lcn/conet_nmds.csv"))) {
  cn.nmds.stats.onc <- capture.output(
    cn.nms.onc <- nmds.min(nmds(cn.d.onc, 2, 2))
  )
  write.csv(cn.nms.onc, "../data/lcn/conet_nmds.csv",
    row.names = FALSE
  )
  write.table(cn.nmds.stats.onc,
    "../data/lcn/conet_nmds_info.txt",
    col.names = FALSE, row.names = FALSE
  )
} else {
  cn.nms.onc <- read.csv("../data/lcn/conet_nmds.csv")
}
# Vector fitting
nv.onc <- envfit(cn.nms.onc, data.frame(onc.com[, colnames(onc.com) != "ds"],
  R = onc.rough,
  C = onc.ns[, c("C")],
  A = ptc.onc
))
cv.onc <- envfit(nms.onc, data.frame(onc.com[, colnames(onc.com) != "ds"],
  R = onc.rough,
  C = onc.ns[, c("C")],
  A = ptc.onc
))

# packing into a dataframe
tree <- onc.geno
for (i in 1:length(unique(onc.geno))) {
  tree[onc.geno == unique(onc.geno)[i]] <-
    1:length(tree[onc.geno == unique(onc.geno)[i]])
}
tree <- factor(tree)
tree.id <- do.call(rbind, strsplit(names(ptc.onc), split = " "))[, 1]

# add chemistry data
onc.nc[, 1] <- as.character(paste0("N", gsub("-", "\\.", onc.nc[, 1])))
onc.tan[, 1] <- as.character(paste0("N", gsub("-", "\\.", onc.tan[, 1])))
# rename headers
# mass is in mg
colnames(onc.nc)[1:4] <- c("tree.id", "sample.mass", "N", "C")
colnames(onc.tan)[1] <- "tree.id"
colnames(onc.tan)[grep("X.CT", colnames(onc.tan))] <- "CT"
# add C:N ratio
onc.nc$rCN <- onc.nc$N / onc.nc$C
# pH data
onc.ph[, "tree.id"] <- gsub("-", ".", onc.ph[, "tree.id"])
onc.ph[, "tree.id"] <- gsub("\\.0", "\\.", onc.ph[, "tree.id"])
# N7.16 is possibly N7.10
onc.ph[onc.ph[, "tree.id"] == "N7.16", "tree.id"] <- "N7.10"
# updated pH from Lamit
onc.ph[!is.na(onc.ph[, "pH2"]), "pH"] <-
  apply(onc.ph[!is.na(onc.ph[, "pH2"]), c("pH", "pH2")], 1, mean)
# collect into a single df
onc.dat <- data.frame(tree.id,
  PC = ptc.onc, SR = spr.onc,
  SD = spd.onc, SE = spe.onc,
  geno = factor(onc.geno), tree = tree,
  BR = onc.rough, onc.ns[, c("L", "Cen")]
)
# species centralities
cen.spp <- lapply(cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id],
  sna::degree,
  rescale = FALSE
)
cen.spp <- do.call(rbind, cen.spp)
colnames(cen.spp) <- colnames(cn.onc[[1]])
# Get to match onc.dat
onc.ph <- onc.ph[onc.ph[, "tree.id"] %in% onc.dat[, "tree.id"], ]
onc.ph <- onc.ph[match(onc.dat[, "tree.id"], onc.ph[, "tree.id"]), ]
# append chemistry to onc.dat
onc.dat <- data.frame(onc.dat,
  C = onc.nc[match(
    onc.dat[, "tree.id"],
    onc.nc[, "tree.id"]
  ), "C"],
  N = onc.nc[match(
    onc.dat[, "tree.id"],
    onc.nc[, "tree.id"]
  ), "N"],
  CN = onc.nc[match(
    onc.dat[, "tree.id"],
    onc.nc[, "tree.id"]
  ), "rCN"],
  CT = onc.tan[match(
    onc.dat[, "tree.id"],
    onc.tan[, "tree.id"]
  ), "CT"],
  pH = onc.ph[, "pH"]
)
# Plot calculations
pw.onc <- onc.com.rel[, colnames(onc.com.rel) != "ds"]
pw.onc <- pw.onc[
  order(apply(pw.onc, 1, sum), decreasing = TRUE),
  order(apply(pw.onc, 2, sum), decreasing = TRUE)
]
rownames(pw.onc) <- onc.geno
col.pal <- RColorBrewer::brewer.pal((max(unlist(mods.onc))), "Paired")
# Network list NA removed
cn.d.onc.na <- distNet(cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id],
  method = "euclidean"
)
# Figure ordinations
# Communities
if (file.exists("../data/lcn/nms_com_onc.rda")) {
  nms.com <- dget(file = "../data/lcn/nms_com_onc.rda")
} else {
  set.seed(12345)
  nms.com <- nmds(vegdist(onc.com.rel), 2, 3)
  dput(nms.com, file = "../data/lcn/nms_com_onc.rda")
}
# Networks
if (file.exists("../data/lcn/nms_cn_onc.rda")) {
  nms.cn <- dget(file = "../data/lcn/nms_cn_onc.rda")
} else {
  set.seed(12345)
  nms.cn <- nmds(cn.d.onc.na, 1, 2)
  dput(nms.cn, file = "../data/lcn/nms_cn_onc.rda")
}
ord.com <- nmds.min(nms.com, 3)
ord.cn <- nmds.min(nms.cn, 2)
# Vectors for plotting
# Composition
vec.com <- envfit(ord.com,
  env = onc.dat[, c("CT", "SR", "SE", "Cen"), ],
  perm = 10000,
  choices = c(1, 2), na.rm = TRUE
)
# Network similarity
vec.cn <- envfit(ord.cn,
  env = na.omit(onc.dat)[, c("CT", "SR", "SE", "Cen")],
  perm = 10000,
  choices = c(1, 2)
)
# onc
if (!("mod_obsval_onc.csv" %in% dir("../data/lcn"))) {
  mod.onc <- slot(
    bipartite::computeModules(rel(onc.com[, -ncol(onc.com)]),
      deep = FALSE
    ),
    "likelihood"
  )
  write.csv(mod.onc, file = "../data/lcn/mod_obsval_onc.csv", row.names = FALSE)
} else {
  mod.onc <- read.csv(file = "../data/lcn/mod_obsval_onc.csv")[, 1]
}
if (!("mod_simvals_onc.csv" %in% dir("../data/lcn"))) {
  onc.sweb <- simulate(vegan::nullmodel(onc.com[, -ncol(onc.com)],
    method = "swsh_samp_c"
  ), 99)
  for (i in 1:dim(onc.sweb)[3]) {
    onc.sweb[, , i] <- rel(onc.sweb[, , i])
  }
  onc.smod <- apply(onc.sweb, 3, bipartite::computeModules)
  mods.onc.sweb <- unlist(lapply(onc.smod, slot, name = "likelihood"))
  write.csv(mods.onc.sweb,
    file = "../data/lcn/mod_simvals_onc.csv",
    row.names = FALSE
  )
  # nest.onc <- bipartite::nestedness(onc.com.rel)
} else {
  mods.onc.sweb <- read.csv(
    file = "../data/lcn/mod_simvals_onc.csv"
  )[, 1]
}

# pit
if (!("mod_obsval_pit.csv" %in% dir("../data/lcn"))) {
  mod.pit <- slot(
    bipartite::computeModules(rel(pit.com), deep = FALSE),
    "likelihood"
  )
  write.csv(mod.pit,
    file = "../data/lcn/mod_obsval_pit.csv",
    row.names = FALSE
  )
} else {
  mod.pit <- read.csv(
    file = "../data/lcn/mod_obsval_pit.csv"
  )[, 1]
}
if (!("mod_simvals_pit.csv" %in% dir("../data/lcn"))) {
  pit.sweb <- simulate(vegan::nullmodel(pit.com, method = "swsh_samp_c"), 99)
  for (i in 1:dim(pit.sweb)[3]) {
    pit.sweb[, , i] <- rel(pit.sweb[, , i])
  }
  pit.smod <- apply(pit.sweb, 3, bipartite::computeModules)
  mods.pit.sweb <- unlist(lapply(pit.smod, slot, name = "likelihood"))
  write.csv(mods.pit.sweb,
    file = "../data/lcn/mod_simvals_pit.csv",
    row.names = FALSE
  )
  # nest.pit <- bipartite::nestedness(pit.com.rel)
} else {
  mods.pit.sweb <- read.csv(
    file = "../data/lcn/mod_simvals_pit.csv"
  )[, 1]
}

### Wild data

# remove notes
wild.dat <- wild.dat[, colnames(wild.dat) != "NOTES."]
wild.dat <- wild.dat[, colnames(wild.dat) != "dead"]
#
wild.dat <- na.omit(wild.dat)
# remove gnu.44 = FREMONT
wild.dat <- wild.dat[wild.dat$tree != "gnu.44", ]
# rm ll.6, tree with super smooth bark
wild.dat <- wild.dat[wild.dat$tree != "ll.6", ]
wild.dat$tree <- factor(as.character(wild.dat$tree))
# condense species
# lecanora, there can be only one!
lec.sp <- apply(wild.dat[, c(6, 8, 10, 18)], 1, 
                function(x) sign(any(x != 0)))
# no physcioids!
# phy.spp <- apply(x[,c(13,14,15,16)],
# 1,function(x) sign(any(x!=0)))
wild.dat <- cbind(wild.dat, lec = lec.sp)
wild.dat <- wild.dat[, -c(6, 8, 10, 18)]
wild.dat <- wild.dat[, colnames(wild.dat) != "physcioid"]
# break into quadrat list (x.q)
quads <- paste(wild.dat$tree, wild.dat$quadrat)
colnames(wild.dat)[5:ncol(wild.dat)] <- c(
  "Xg", "Cs", "Xm", "fgb", "Rs",
  "Pm", "Pa", "Pu", "Ch", "Ls"
)
wild.dat <- wild.dat[colnames(wild.dat) != "fgb"]
wild.dat.q <- split(wild.dat, quads)
wild.com <- split(wild.dat, wild.dat$tree)
wild.com <- do.call(rbind, lapply(wild.com, function(x) 
    apply(x[, -1:-4], 2, sum)))
wild.com.rel <- apply(wild.com, 2, function(x) x / max(x))
wild.com.rel[is.na(wild.com.rel)] <- 0
wild.q <- lapply(split(wild.dat, wild.dat$tree), 
                 function(x) x[, -1:-4])

# data from lamit

env <- env[is.na(env$Pct.Roughness) == FALSE, ]
env[, 1] <- sub(
  "\\?", "",
  sub(
    "\\.0", "\\.",
    sub(
      "\\_", "\\.",
      sub("\\-", "\\.", tolower(as.character(env[, 1])))
    )
  )
)
env[env[, 1] == "ll.6_(znu.29)", 1] <- "ll.6"
env[env[, 1] == "gnu.85.1ftaway", 1] <- "gnu.85"
env$Quad.Loc <- as.character(sapply(
  as.character(env$Quad.Loc),
  function(x) {
    unlist(strsplit(x, split = "_"))[2]
  }
))
env$Quad.Loc <- sub("\\-", "\\.", env$Quad.Loc)
env$Quad.Loc <- paste("n", env$Quad.Loc, sep = "")
# remove southern aspect
env <- env[env$Aspect != "South", ]
env.tid <- paste(env$Tree.ID, env$Quad.Loc)
# check that the datasets are compatible
all(names(wild.dat.q) %in% env.tid)
# match observations
all(env.tid[match(names(wild.dat.q), env.tid)] == names(x.q))
# delimit co-occurrence and match
env <- env[match(names(wild.dat.q), env.tid), ]
wild.dat.split <- paste(wild.dat$tree, wild.dat$quadrat, 
                        sep = "_")
env.split <- paste(env$Tree.ID, env$Quad.Loc)
wild.dat.split <- as.character(wild.dat$tree)
env.split <- as.character(env$Tree.ID)
# percent rough bark
prb.wild <- tapply(env$Pct.Roughness, env.split, mean)

# age
dbh <- age$DBH.cm_01
age.final <- age$AgeFinal.U
age <- data.frame(tree.id = age[, 1], age.final = age$AgeFinal.U)
age[, 1] <- tolower(age[, 1])
age[, 1] <- sub("_", "\\.", age[, 1])
age[, 1] <- sub("-", "\\.", age[, 1])
age[, 1] <- sub("\\?", "", age[, 1])
age[, 1] <- sub("\\.0", "\\.", age[, 1])
age[age[, 1] == "gnu.85.1ftaway", 1] <- "gnu.85"
# predict age
gnu19.dbh <- dbh[age$tree.id == "gnu.19"]
new <- data.frame(dbh = seq(min(dbh), max(dbh), by = 0.1))
age.final <- na.omit(age.final)
pred.age <- predict(lm(age.final ~ dbh, data = age), new)
gnu19.age <- as.numeric(pred.age[new[, 1] == gnu19.dbh])
#
tree.age <- age[match(names(prb.wild), age[, 1]), 2]
tree.age[is.na(tree.age)] <- gnu19.age
names(tree.age) <- age[match(names(prb.wild), age[, 1]), 1]
age <- tree.age
# percent cover
pc.wild <- unlist(lapply(
  wild.q,
  function(x) {
    sum(apply(
      x, 1,
      function(x) sign(sum(x))
    ))
  }
))
# richness
sr.wild <- unlist(lapply(
  wild.q, function(x) sum(sign(apply(x, 2, sum))))
  )
# networks
cn.wild <- lapply(wild.q, coNet)
cn.mu.wild <- meanNet(cn.wild)
cn.d.wild <- distNet(cn.wild, method = "bc")
# network stats
ns.wild <- do.call(rbind, lapply(lapply(cn.wild, function(x) {
  abs(sign(x))
}), enaR:::structure.statistics))
# centralization
dcen.wild <- unlist(lapply(cn.wild, function(x) {
  sna::centralization(x, FUN = sna::degree, normalize = FALSE)
}))
# wild data frame
wild.dat <- data.frame(
  tree = names(tree.age),
  age = tree.age, BR = prb.wild,
  PC = pc.wild, SR = sr.wild,
  L = ns.wild[, "L"], Cen = dcen.wild
)


