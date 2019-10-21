### Loading data objects
### Data notes:
# Loading data
# separate onc

if (100 < 1){

    xgal.size.in = read.csv("./data/lcn/ONC_Xgal_SizeData_May2011.csv")
    garden.data.in = read.csv("./data/lcn/LCO_data_ONC_PIT.csv")
    rough.in = read.csv("./data/lcn/ONC_raw_roughness.csv")
    onc.nc.in = read.csv("./data/lcn/ONC_phytochem_NC.csv")
    onc.tan.in = read.csv("./data/lcn/ONC_phytochem_tannin.csv")
    onc.ph.in = read.csv("./data/lcn/ONC_Bark_lichen_pH_data.csv")
    wild.dat.in = read.csv("./data/lcn/lco_Apr2012.csv")
    env.in = read.csv("./data/lcn/Uinta2012_all_data_from_Lamit.csv")
    age.in = read.csv(
    "./data/lcn/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv")
    garden.data = proc_garden_data(garden.data.in)
    pit = proc_pit(garden.data.in)
    onc = proc_onc(garden.data)
    onc.q = proc_onc_q(onc)
    onc.ph = proc_onc_ph(garden.data, 
                         rough.in, 
                         onc, onc.q, 
                         onc.nc.in, onc.tan.in, onc.ph.in)
    onc.dat = proc_onc_dat(garden.data, rough.in, 
                           onc, onc.q, 
                           onc.nc.in, onc.tan.in, 
                           onc.ph)
    cn.onc = proc_cn_onc(onc)
    onc.ns = proc_onc_ns(cn.onc)
    cn.d.onc.na = proc_cn_d_onc(cn.onc, onc.dat, rm.na = TRUE)
    onc.com = proc_onc_com(garden.data, onc, onc.q)
    onc.com.rel = proc_onc_com_rel(onc.com)

    h2.tab <- table_h2(onc.dat, cn.d.onc.na, onc.ns, onc.com.rel)
}


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

proc_onc_dat <- function(garden.data, rough.in, onc, onc.q,
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
    onc.ns <- cbind(ns.onc, Cen = dcen.onc, mod.lik = cn.mod.onc[, 1],
                    mod.n = cn.mod.onc[, 2])
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
    onc.dat <- data.frame(tree.id, 
                          PC = ptc.onc, SR = spr.onc, SD = spd.onc, SE = spe.onc, 
                          geno = factor(onc.geno), tree = tree, 
                          BR = onc.rough, onc.ns[, c("L", "Cen")])
    onc.dat <- data.frame(onc.dat, 
                          C = onc.nc[match(onc.dat[, "tree.id"], 
                                           onc.nc[, "tree.id"]), "C"], 
                          N = onc.nc[match(onc.dat[, "tree.id"], 
                                           onc.nc[, "tree.id"]), "N"], 
                          CN = onc.nc[match(onc.dat[,"tree.id"], 
                                            onc.nc[, "tree.id"]), "rCN"], 
                          CT = onc.tan[match(onc.dat[, "tree.id"], 
                                             onc.tan[, "tree.id"]), "CT"], 
                          pH = onc.ph[, "pH"])
    return(onc.dat)
}

proc_onc_com <- function(garden.data, onc, onc.q) {
    g1 <- substr(garden.data[, 1], 2, 2)
    g1[g1 != "P"] <- "onc"
    onc <- garden.data[g1 == "onc", ]
    onc.q <- split(onc, paste(onc[, 1], onc[, 2]))
    onc.q <- lapply(onc.q, function(x) x[, -1:-6])
    onc.com <- do.call(rbind, lapply(onc.q, function(x) apply(x, 2, sum)))
    onc.com <- cbind(onc.com, 
                     ds = rep(min(onc.com[onc.com != 0]) / 1000, nrow(onc.com)))
    return(onc.com)
}


proc_onc_com_rel <- function(onc.com) {
    if (any(tolower(colnames(onc.com)) == "ds")){
        onc.com <- onc.com[, tolower(colnames(onc.com)) != "ds"]
    }
    onc.com.rel <- apply(onc.com, 2, function(x) x / max(x))
    onc.com.rel <- cbind(onc.com.rel, 
                         ds = rep(min(onc.com.rel[onc.com.rel != 0]) / 1000, 
                                  nrow(onc.com.rel)))
    return(onc.com.rel)
}


proc_cn_onc <- function(onc, ci.p = 95) {
    lapply(split(onc[, -1:-6], onc[, "Tree"]), 
           coNet, ci.p = ci.p)
}

proc_cn_d_onc <- function(cn.onc, onc.dat, rm.na = TRUE){
    if (rm.na){
        cn.d.onc.na <- distNet(cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id], 
                               method = "euclidean")
    }else{
        cn.d.onc.na <- distNet(cn.onc, method = "euclidean")
    }
    return(cn.d.onc.na)
}

proc_onc_ns <- function(cn.onc) {
    ns.onc <- lapply(
        lapply(cn.onc, 
               function(x) {
                   abs(sign(x))
               }
               ), 
        enaR:::structure.statistics
    )
    ns.onc <- do.call(rbind, ns.onc)
    cn.mod.onc <- matrix(nrow = length(cn.onc), ncol = 2)
    for (i in 1:length(cn.onc)) {
        if (sum(sign(cn.onc[[i]])) >= 3) {
            mod.tmp <- computeModules(cn.onc[[i]])
            cn.mod.onc[i, 1] <- slot(mod.tmp, "likelihood")
            cn.mod.onc[i, 2] <- nrow(slot(mod.tmp, "modules")) - 1
        }
        else {
            cn.mod.onc[i] <- NA
        }
    }
    cn.mod.onc[is.na(cn.mod.onc)] <- 0
    names(cn.mod.onc) <- c("mod.lik", "mod.n")
    dcen.onc <- unlist(lapply(cn.onc, function(x) {
        sna::centralization(x, FUN = sna::degree, normalize = FALSE)
    }))
    onc.ns <- cbind(ns.onc, 
                    Cen = dcen.onc, 
                    mod.lik = cn.mod.onc[, 1], 
                    mod.n = cn.mod.onc[, 2])
return(onc.ns)

}

table_h2 <- function(onc.dat, cn.d.onc.na, onc.ns, onc.com.rel){
    h2.tab <- matrix("", 1, 4)
    colnames(h2.tab) <- c("Response", "H2", "R2", "p-value")
    ptc.reml <- lme4::lmer(I(PC^(1 / 2)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    ptc.reml.pval <- RLRsim::exactRLRT(ptc.reml)
    ptc.reml.result <- c("Percent Lichen Cover", H2(ptc.reml, g = onc.dat$geno), R2(ptc.reml), ptc.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, ptc.reml.result)
    spr.reml <- lme4::lmer(I(SR^(1 / 2)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    spr.reml.pval <- RLRsim::exactRLRT(spr.reml)
    spr.reml.result <- c("Lichen Species Richness", H2(spr.reml, g = onc.dat$geno), R2(spr.reml), spr.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, spr.reml.result)
    prb.reml <- lme4::lmer(I(BR^(1 / 2)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    prb.reml.pval <- RLRsim::exactRLRT(prb.reml)
    prb.reml.result <- c("Percent Rough Bark", H2(prb.reml, g = onc.dat$geno), R2(prb.reml), prb.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, prb.reml.result)
    ph.reml <- lme4::lmer(I(pH^(1 / 2)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    ph.reml.pval <- RLRsim::exactRLRT(ph.reml)
    ph.reml.result <- c("pH", H2(ph.reml, g = onc.dat$geno), R2(ph.reml), ph.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, ph.reml.result)
    ct.reml <- lme4::lmer(I(CT^(1 / 4)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    ct.reml.pval <- RLRsim::exactRLRT(ct.reml)
    ct.reml.result <- c("Condensed Tannins (CT)", H2(ct.reml, g = onc.dat$geno), R2(ct.reml), ct.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, ct.reml.result)
    cnr.reml <- lme4::lmer(I(CN^(1)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    cnr.reml.pval <- RLRsim::exactRLRT(cnr.reml)
    cnr.reml.result <- c("Carbon-Nitrogen (CN) Ratio", H2(cnr.reml, g = onc.dat$geno), R2(cnr.reml), cnr.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, cnr.reml.result)
    rcom.perm <- vegan::adonis2(onc.com.rel^(1 / 1) ~ geno + BR + PC + SR, data = onc.dat, perm = 10000, mrank = TRUE)
    h2.tab[4, "p-value"] <- unlist(rcom.perm)["Pr(>F)1"]
    h2.tab[4, "H2"] <- H2(rcom.perm, g = onc.dat$geno)
    h2.tab[4, "R2"] <- R2(rcom.perm)
    h2.tab[4, "Response"] <- "Lichen Community Composition"
    cn.perm <- vegan::adonis2(
                          cn.d.onc.na ~ geno + 
                              BR + 
                              pH + CN + CT + 
                              PC + SR + SE, by = "term", 
                          data = na.omit(onc.dat), permutations = 10000)
    h2.tab[5, "p-value"] <- as.matrix(cn.perm)[1, "Pr(>F)"]
    h2.tab[5, "H2"] <- H2(cn.perm, g = onc.dat[, "geno"], perm = 10000)
    h2.tab[5, "R2"] <- R2(cn.perm)
    h2.tab[5, "Response"] <- "Lichen Network"
    link.reml <- lme4::lmer(I(log(L + 1e-08)) ~ (1 | geno), data = onc.dat, REML = TRUE)
    link.reml.pval <- RLRsim::exactRLRT(link.reml, nsim = 50000)
    link.reml.result <- c("Number of Network Links", H2(link.reml, g = onc.dat$geno), R2(link.reml), link.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, link.reml.result)
    cen.reml <- lme4::lmer(I(Cen^(1 / 2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
    cen.reml.pval <- RLRsim::exactRLRT(cen.reml, nsim = 50000)
    cen.reml.result <- c("Network Centrality", H2(cen.reml, g = onc.dat$geno), R2(cen.reml), cen.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, cen.reml.result)
    mod.reml <- lme4::lmer(I(onc.ns[, "mod.lik"]^(1 / 4)) ~ (1 | geno), data = onc.dat, REML = TRUE)
    mod.reml.pval <- RLRsim::exactRLRT(mod.reml)
    mod.reml.result <- c("Network Modularity", H2(mod.reml, g = onc.dat$geno), R2(mod.reml), mod.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, mod.reml.result)
    spd.reml <- lme4::lmer(I(SD^(1 / 2)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    spd.reml.pval <- RLRsim::exactRLRT(spd.reml)
    spd.reml.result <- c("Lichen Species Diversity", H2(spd.reml, g = onc.dat$geno), R2(spd.reml), spd.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, spd.reml.result)
    spe.reml <- lme4::lmer(I(SE^(1 / 4)) ~ (1 | geno), data = na.omit(onc.dat), REML = TRUE)
    spe.reml.pval <- RLRsim::exactRLRT(spe.reml)
    spe.reml.result <- c("Lichen Species Evenness", H2(spe.reml, g = onc.dat$geno), R2(spe.reml), spe.reml.pval$p.value)
    h2.tab <- rbind(h2.tab, spe.reml.result)
    h2.tab[, "H2"] <- round(as.numeric(h2.tab[, "H2"]), digits = 5)
    h2.tab[, "R2"] <- round(as.numeric(h2.tab[, "R2"]), digits = 5)
    h2.tab[, "p-value"] <- round(as.numeric(h2.tab[, "p-value"]), digits = 5)
    h2.tab <- h2.tab[order(h2.tab[, "H2"], decreasing = TRUE), ]
    return(h2.tab)
}

