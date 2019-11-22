### Loading data objects
### Data notes:
# Loading data
# separate onc

proc_garden_data <- function(garden.data, rm.geno, rm.tree) {
                                        # rm genotype RL6 and N1.31
    garden.data <- garden.data[!(garden.data$Geno %in% rm.geno), ]
    garden.data <- garden.data[!(garden.data$Geno %in% rm.tree), ]
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


proc_onc_dat <- function(garden.data, rough.in, onc, onc.q,
                         onc.nc.in, onc.tan.in, onc.ph.in, rm.na = TRUE) {
    ## Genotype and Tree vectors
    onc.geno <- unlist(sapply(
        names(onc.q), 
        function(x) 
            strsplit(x, split = " ")[[1]][2]))
    onc.tree <- do.call(rbind, 
                        strsplit(names(onc.geno), " "))[, 1]
    ## Roughness data
    rough <- rough.in[, 1:5]
    rough <- rough[grepl("North", rough[, 3]), ]
    avg.rough <- tapply(rough[, 5], rough[, 1], mean)
    r.tree <- names(avg.rough)
    r.tree <- sub("-", "\\.", r.tree)
    r.tree <- sub("\\.0", "\\.", r.tree)
    names(avg.rough) <- r.tree
    onc.rough <- avg.rough[match(onc.tree, r.tree)]
    onc.split <- split(onc[, -1:-6], onc[, "Tree"], drop = TRUE)
    ## Network data and metrics
    cn.onc <- lapply(onc.split, coNet, ci.p = 95)
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
    onc.com <- do.call(rbind, lapply(onc.q, function(x) apply(x, 2, sum)))
    onc.com <- cbind(onc.com, 
                     ds = rep(min(onc.com[onc.com != 0]) / 1000, 
                              nrow(onc.com)))
    ptc.onc <- unlist(lapply(onc.q, function(x) {
        sum(apply(x, 1, function(x) sign(sum(x))))
    }))
    spr.onc <- apply(onc.com[, colnames(onc.com) != "ds"], 
                     1, 
                     function(x) sum(sign(x)))
    spd.onc <- diversity(onc.com[, colnames(onc.com) != "ds"])
    spe.onc <- spd.onc / 
        log(specnumber(onc.com[, colnames(onc.com) != "ds"]))
    spe.onc[is.na(spe.onc)] <- 0
    tree <- onc.geno
    for (i in seq_along(unique(onc.geno))) {
        tree[onc.geno == unique(onc.geno)[i]] <- seq_along(
            tree[onc.geno == unique(onc.geno)[i]])
    }
    tree <- factor(tree)
    tree.id <- do.call(rbind, strsplit(names(ptc.onc), split = " "))[, 1]
    ## Chemistry Data: C:N, tannins and pH
    onc.nc <- onc.nc.in
    onc.tan <- onc.tan.in
    onc.ph <- onc.ph.in
    onc.nc[, 1] <- as.character(
        paste0("N", gsub("-", "\\.", onc.nc[, 1])))
    onc.tan[, 1] <- as.character(
        paste0("N", gsub("-", "\\.", onc.tan[, 1])))
    colnames(onc.nc)[1:4] <- c("tree.id", "sample.mass", "N", "C")
    colnames(onc.tan)[1] <- "tree.id"
    colnames(onc.tan)[grep("X.CT", colnames(onc.tan))] <- "CT"
    onc.nc$rCN <- onc.nc$N / onc.nc$C
    onc.ph[, "tree.id"] <- gsub("-", ".", onc.ph[, "tree.id"])
    onc.ph[, "tree.id"] <- gsub("\\.0", "\\.", onc.ph[, "tree.id"])
    onc.ph[onc.ph[, "tree.id"] == "N7.16", "tree.id"] <- "N7.10"
    onc.ph[!is.na(onc.ph[, "pH2"]), "pH"] <- apply(
        onc.ph[!is.na(onc.ph[, "pH2"]), c("pH", "pH2")], 1, mean)
    ## Create a data frame for export
    onc.dat <- data.frame(
        tree.id, 
        PC = ptc.onc, 
        SR = spr.onc, 
        SD = spd.onc, 
        SE = spe.onc, 
        geno = factor(onc.geno), 
        tree = tree, 
        BR = onc.rough, 
        onc.ns[, c("L", "Cen")])
    onc.ph <- onc.ph[onc.ph[, "tree.id"] %in% onc.dat[, "tree.id"], ]
    onc.ph <- onc.ph[match(onc.dat[, "tree.id"], onc.ph[, "tree.id"]), ]
    onc.geno <- unlist(sapply(
        names(onc.q), function(x) 
            strsplit(x, split = " ")[[1]][2]))
    onc.tree <- do.call(rbind, strsplit(names(onc.geno), " "))[, 1]
    rough <- rough.in[, 1:5]
    rough <- rough[grepl("North", rough[, 3]), ]
    avg.rough <- tapply(rough[, 5], rough[, 1], mean)
    r.tree <- names(avg.rough)
    r.tree <- sub("-", "\\.", r.tree)
    r.tree <- sub("\\.0", "\\.", r.tree)
    names(avg.rough) <- r.tree
    onc.rough <- avg.rough[match(onc.tree, r.tree)]
    onc.split <- split(onc[, -1:-6], onc[, "Tree"], drop = TRUE)
    cn.onc <- lapply(onc.split, coNet, ci.p = 95)
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
    onc.com <- do.call(rbind, lapply(onc.q, function(x) apply(x, 2, sum)))
    onc.com <- cbind(onc.com, 
                     ds = rep(min(onc.com[onc.com != 0]) / 1000, 
                              nrow(onc.com)))
    ptc.onc <- unlist(lapply(onc.q, function(x) {
        sum(apply(x, 1, function(x) sign(sum(x))))
    }))
    spr.onc <- apply(onc.com[, colnames(onc.com) != "ds"], 
                     1, 
                     function(x) sum(sign(x)))
    spd.onc <- diversity(onc.com[, colnames(onc.com) != "ds"])
    spe.onc <- spd.onc / log(specnumber(onc.com[, colnames(onc.com) != "ds"]))
    spe.onc[is.na(spe.onc)] <- 0
    tree <- onc.geno
    for (i in seq_along(unique(onc.geno))) {
        tree[onc.geno == unique(onc.geno)[i]] <- seq_along(
            tree[onc.geno == unique(onc.geno)[i]])
    }
    tree <- factor(tree)
    tree.id <- do.call(rbind, strsplit(names(ptc.onc), split = " "))[, 1]
    onc.dat <- data.frame(
        tree.id, 
        PC = ptc.onc, 
        SR = spr.onc, 
        SD = spd.onc, 
        SE = spe.onc, 
        geno = factor(onc.geno), 
        tree = tree, 
        BR = onc.rough, 
        onc.ns[, c("L", "Cen", "mod.lik")],
        C = onc.nc[match(onc.dat[, "tree.id"], 
                         onc.nc[, "tree.id"]), "C"], 
        N = onc.nc[match(onc.dat[, "tree.id"], 
                         onc.nc[, "tree.id"]), "N"], 
        CN = onc.nc[match(onc.dat[, "tree.id"], 
                          onc.nc[, "tree.id"]), "rCN"], 
        CT = onc.tan[match(onc.dat[, "tree.id"], 
                           onc.tan[, "tree.id"]), "CT"], 
        pH = onc.ph[, "pH"]
    )
    if (rm.na){onc.dat <- na.omit(onc.dat)}
    return(onc.dat)
}

proc_onc_com <- function(garden.data, onc, onc.q, onc.dat, rm.na = TRUE) {
    g1 <- substr(garden.data[, 1], 2, 2)
    g1[g1 != "P"] <- "onc"
    onc <- garden.data[g1 == "onc", ]
    onc.q <- split(onc, paste(onc[, 1], onc[, 2]))
    onc.q <- lapply(onc.q, function(x) x[, -1:-6])
    onc.com <- do.call(
        rbind, 
        lapply(onc.q, 
               function(x) apply(x, 2, sum))
    )
    onc.com <- cbind(
        onc.com, 
        ds = rep(
            min(onc.com[onc.com != 0]) / 1000, nrow(onc.com))
    )
    if (rm.na){
        onc.com <- onc.com[match(
            rownames(onc.dat), rownames(onc.com)), ]
    }
    return(onc.com)
}


proc_onc_com_rel <- function(onc.com, onc.dat, rm.na = TRUE) {
    if (rm.na){
        onc.com <- onc.com[match(
            rownames(onc.dat), rownames(onc.com)), ]
    }
    if (any(tolower(colnames(onc.com)) == "ds")){
        onc.com <- onc.com[, tolower(colnames(onc.com)) != "ds"]
    }
    onc.com.rel <- apply(onc.com, 2, function(x) x / max(x))
    onc.com.rel <- cbind(onc.com.rel, 
                         ds = rep(min(onc.com.rel[onc.com.rel != 0]) / 1000, 
                                  nrow(onc.com.rel)))
    return(onc.com.rel)
}


proc_cn_onc <- function(onc, onc.dat, ci.p = 95) {
    lapply(split(onc[, -1:-6], onc[, "Tree"]), 
           coNet, ci.p = ci.p)
}

proc_cn_d_onc <- function(cn.onc, onc.dat, rm.na = TRUE){
    if (rm.na){
        cn.d.onc <- distNet(
            cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id], 
            method = "euclidean"
        )
    }else{
        cn.d.onc <- distNet(cn.onc, method = "euclidean")
    }
    return(cn.d.onc)
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

check_fligner <- function(onc.dat){
    y <- onc.dat[, c(-1, -6, -7)]
    flig.y <- do.call(rbind, apply(y, 2, 
                         fligner.test, 
                         g = onc.dat[, "geno"]))
    flig.y2 <- do.call(rbind, apply(y^(2), 2, 
                         fligner.test, 
                         g = onc.dat[, "geno"]))
    flig.r2y <- do.call(rbind, apply(y^(1/2), 2, 
                         fligner.test, 
                         g = onc.dat[, "geno"]))
    flig.r4y <- do.call(rbind, apply(y^(1/4), 2, 
                         fligner.test, 
                         g = onc.dat[, "geno"]))
    flig.logy <- do.call(rbind, apply(log(y + 0.0001), 2, 
                         fligner.test, 
                         g = onc.dat[, "geno"]))
    out <- list(y = flig.y, 
                y2 = flig.y2,
                r2y = flig.r2y, 
                r4y = flig.r4y, 
                logy = flig.logy)
    out <- melt.list(out)
    out <- out[, -7]
    out <- cbind(out[, "L1"], out[, colnames(out) != "L1"])
    colnames(out)[1] <- "transformation"
    return(out)
}


check_shapiro <- function(reml.reml){
    reml.resids <- lapply(reml.reml, residuals)
    reml.shapiros <- lapply(reml.resids, shapiro.test)
    reml.forms <- lapply(reml.reml, formula)
    reml.forms <- lapply(reml.forms, as.character)
    reml.forms <- lapply(reml.forms, paste0, collapse = "")
    formula <- unlist(reml.forms)
    out <- do.call(rbind, reml.shapiros)
    out <- out[, -ncol(out)]
    out <- cbind(formula, out)
    out <- as.data.frame(out)
    return(out)
}

run_spp_centrality <- function(cn.onc, onc.dat){
    cen.spp <- lapply(cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id],
                      sna::degree, 
                      rescale = TRUE
                      )
    cen.spp <- do.call(rbind, cen.spp)
    cen.spp[is.na(cen.spp)] <- 0
    colnames(cen.spp) <- colnames(cn.onc[[1]])
    if (any(rownames(cen.spp) != onc.dat[, "tree.id"])){stop()}
    colnames(cen.spp) <- paste0("cen_", colnames(cen.spp))
    cen.spp <- cen.spp[, apply(sign(cen.spp), 2, sum) > 3]
    cen.spp.reml <- list()
    for (i in seq_along(colnames(cen.spp))){
        cen.reml <- lme4::lmer(I(cen.spp[ ,i]^(1 / 4)) ~ (1 | geno), 
                               data = onc.dat, REML = TRUE)
        reml.pval <- RLRsim::exactRLRT(cen.reml)
        cen.spp.reml[[i]] <- c(
            "spp centrality", 
            H2 = H2(cen.reml, g = onc.dat[, "geno"]), 
            R2 = R2(cen.reml), 
            p.value = reml.pval[["p.value"]]
            )
    }
    cen.spp.table <- do.call(rbind, cen.spp.reml)[, -1]
    cen.spp.table <- apply(cen.spp.table, 2, as.numeric)
    rownames(cen.spp.table) <- colnames(cen.spp)
    out <- list(cen.spp = cen.spp, cen.spp.reml = cen.spp.table)
}

run_reml <- function(onc.dat, rm.na = TRUE, raw.reml = FALSE){
    if (rm.na){onc.dat <- na.omit(onc.dat)}
    ## tree traits
    prb.reml <- lme4::lmer(I(BR^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    prb.reml.pval <- RLRsim::exactRLRT(prb.reml)
    prb.reml.result <- c("Percent Rough Bark", 
                         H2(prb.reml, g = onc.dat$geno), 
                         R2(prb.reml), prb.reml.pval$p.value)
    ph.reml <- lme4::lmer(I(pH^(1 / 4)) ~ (1 | geno), 
                          data = onc.dat, REML = TRUE)
    ph.reml.pval <- RLRsim::exactRLRT(ph.reml)
    ph.reml.result <- c("pH", H2(ph.reml, g = onc.dat$geno), 
                        R2(ph.reml), ph.reml.pval$p.value)
    ct.reml <- lme4::lmer(I(CT^(1 / 4)) ~ (1 | geno), 
                          data = onc.dat, REML = TRUE)
    ct.reml.pval <- RLRsim::exactRLRT(ct.reml)
    ct.reml.result <- c("Condensed Tannins (CT)", 
                        H2(ct.reml, g = onc.dat$geno), 
                        R2(ct.reml), ct.reml.pval$p.value)
    cnr.reml <- lme4::lmer(I(log(CN^(1 / 1) + 0.001)) ~ (1 | geno), 
                            data = onc.dat, REML = TRUE)
    cnr.reml.pval <- RLRsim::exactRLRT(cnr.reml)
    cnr.reml.result <- c("Carbon-Nitrogen (CN) Ratio", 
                         H2(cnr.reml, g = onc.dat$geno), R2(cnr.reml), 
                         cnr.reml.pval$p.value)
    ## lichen community metrics
    ptc.reml <- lme4::lmer(I(PC^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    ptc.reml.pval <- RLRsim::exactRLRT(ptc.reml)
    ptc.reml.result <- c("Percent Lichen Cover", 
                         H2(ptc.reml, g = onc.dat$geno), 
                         R2(ptc.reml), ptc.reml.pval$p.value)
    spr.reml <- lme4::lmer(I(log(SR^(1 / 1) + 0.001)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    spr.reml.pval <- RLRsim::exactRLRT(spr.reml)
    spr.reml.result <- c("Lichen Species Richness", 
                         H2(spr.reml, g = onc.dat$geno), 
                         R2(spr.reml), spr.reml.pval$p.value)
    spd.reml <- lme4::lmer(I(SD^(1 / 4)) ~  (1 | geno), 
                           data = onc.dat, REML = TRUE)
    spd.reml.pval <- RLRsim::exactRLRT(spd.reml)
    spd.reml.result <- c("Lichen Species Diversity", 
                         H2(spd.reml, g = onc.dat$geno), 
                         R2(spd.reml), spd.reml.pval$p.value)
    spe.reml <- lme4::lmer(I(SE^(1 / 4)) ~  (1 | geno), 
                           data = onc.dat, REML = TRUE)
    spe.reml.pval <- RLRsim::exactRLRT(spe.reml)
    spe.reml.result <- c("Lichen Species Evenness", 
                         H2(spe.reml, g = onc.dat$geno), 
                         R2(spe.reml), spe.reml.pval$p.value)
                                        # Lichen network metrics
    link.reml <- lme4::lmer(I(L^(1 / 4)) ~ (1 | geno), 
                            data = onc.dat, REML = TRUE)
    link.reml.pval <- RLRsim::exactRLRT(link.reml, nsim = 50000)
    link.reml.result <- c("Number of Network Links", 
                          H2(link.reml, g = onc.dat$geno), 
                          R2(link.reml), link.reml.pval$p.value)
    cen.reml <- lme4::lmer(I(Cen^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.reml.pval <- RLRsim::exactRLRT(cen.reml, nsim = 50000)
    cen.reml.result <- c("Network Centrality", 
                         H2(cen.reml, g = onc.dat$geno), 
                         R2(cen.reml), cen.reml.pval$p.value)
    mod.reml <- lme4::lmer(I(mod.lik^(2 / 1)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    mod.reml.pval <- RLRsim::exactRLRT(mod.reml)
    mod.reml.result <- c("Network Modularity", 
                         H2(mod.reml, g = onc.dat$geno), 
                         R2(mod.reml), mod.reml.pval$p.value)
    if (raw.reml){
        out <- list(prb.reml, 
                     ph.reml,
                     ct.reml,    
                     cnr.reml,
                     ptc.reml,
                     spr.reml,
                     spe.reml,
                     spd.reml,
                     link.reml,
                     mod.reml,
                     cen.reml)
    }else{
        out <- rbind(prb.reml.result, 
                     ph.reml.result,
                     ct.reml.result,    
                     cnr.reml.result,
                     ptc.reml.result,
                     spr.reml.result,
                     spe.reml.result,
                     spd.reml.result,
                     link.reml.result,
                     mod.reml.result,
                     cen.reml.result)

    }
    return(out)
}

g.mat <- model.matrix(~ g - 1)
gm.d <- dist(g.mat)
t.mat <- onc.dat[, c("BR", "CT")]
t.mat <- apply(t.mat, 2, function(x) x/max(x))
tm.d <- dist(t.mat)

MRM(tm.d ~ gm.d)
MRM(cn.d.onc ~ gm.d + tm.d)
MRM(cn.d.onc ~  tm.d + gm.d)
MRM(cn.d.onc ~  gm.d)


run_perm <- function(onc.dat, onc.com, cn.d.onc){
    com.perm <- vegan::adonis2((onc.com^(1/4)) ~ geno + 
                                   BR + pH + CN + CT + 
                                       PC + SR + SE,
                               data = onc.dat, 
                               by = "term",
                               mrank = TRUE,
                               perm = 10000)
    cn.perm <- vegan::adonis2(cn.d.onc ~ geno + 
                                  BR + pH + CN + CT +
                                      PC + SR + SE + 
                                         L + mod.lik + Cen,
                              by = "term", 
                              data = onc.dat, 
                              mrank = TRUE,
                              permutations = 10000)
    out <- list(com = com.perm, 
                cn = cn.perm)
    return(out)
}

make_tables <- function(onc.dat, reml.results, perm.results, digits = 4){
    ## Heritability table
    h2.tab <- reml.results
    colnames(h2.tab) <- c("Response", "H2", "R2", "p-value")
    ## Add PERMANOVA results
    com.perm.h2 <- c("Community Composition", 
                     H2(perm.results[["com"]], g = onc.dat[["geno"]]),
                     R2(perm.results[["com"]]),
                     unlist(perm.results[["com"]])["Pr(>F)1"]
                     )
    cn.perm.h2 <- c("Lichen Network Similarity", 
                    H2(perm.results[["cn"]], g = onc.dat[, "geno"], perm = 10000),
                    R2(perm.results[["cn"]]),
                    unlist(perm.results[["cn"]])["Pr(>F)1"]
                    )
    h2.tab <- rbind(h2.tab, 
                    cn.perm.h2,
                    com.perm.h2)
    ## Format Heritability Table
    h2.tab[, c("H2", "R2", "p-value")] <- apply(
        h2.tab[, c("H2", "R2", "p-value")], 
        2, 
        function(x, digits) round(as.numeric(x), digits = digits), 
        digits = digits
    )
    h2.tab <- na.omit(h2.tab)
    ## Organize H2 table
    h2.tab <- h2.tab[c(12, 9, 11, 10, 5:8, 13, 1, 3, 4, 2), ]
                                        # Lichen Networks
                                        # Lichen Network Metrics
                                        # Lichen Community
                                        # Tree Traits
    ## Remove R2 from H2 table
    h2.tab <- h2.tab[, colnames(h2.tab) != "R2"]
    ## Format lichen network permanova table
    cn.perm <- as.data.frame(perm.results[["cn"]])
    rownames(cn.perm) <- c("Genotype", "Bark Roughness", "pH", 
                           "C:N Ratio", "Condensed Tannins", 
                           "Percent Cover", "Species Richness",
                           "Species Evenness", "Number of Links", 
                           "Network Modularity", "Network Centrality", 
                           "Residual", "Total")
    colnames(cn.perm) <- c("df", "SS", "R2", "F", "p-value")
    ## Create the latex
    tab.h2 <- xtable::xtable(
       h2.tab,
       caption = "Genotypic effects on the associated lichen community.",
       label = "tab:h2_table",
       type = "latex",
       include.rownames = FALSE,
       include.colnames = TRUE, 
        digits = digits
       )
    tab.cn.perm <- xtable::xtable(
        cn.perm,
        caption = "PERMANOVA Pseudo-F Table of lichen network similarity.",
        label = "tab:cn_perm",
        type = "latex", 
        include.rownames = TRUE,
        include.colnames = TRUE, 
        digits = digits
        )
    tab.com.perm <- xtable::xtable(
        perm.results[["com"]],
        caption = "Pseudo-F Table of lichen community similarity PERMANOVA.",
        label = "tab:com_perm",
        type = "latex", 
        include.rownames = TRUE,
        include.colnames = TRUE, 
        digits = digits
        )
    out <- list(h2_reml = tab.h2, 
                cn = tab.cn.perm, 
                com = tab.com.perm)
    return(out)
}

run_nms <- function(d, vec.data, dim = 2, seed = 12345){
    set.seed(seed)
    nms <- nmds(d, dim, dim)
    nms.out <- capture.output(nms <- nmds.min(nms))
    vec <- vf(nms, vec.data)
    out <- list(nms = nms, vec = vec, report = nms.out)
    return(out)
}

run_sppcen_aov <- function(spp.cen){
    sppcen.aov <- aov(value ~ X2, data = melt(spp.cen[["cen.spp"]]))
    sppcen.mct <- TukeyHSD(sppcen.aov)
    out <- list(aov = sppcen.aov, mct = sppcen.mct)
    return(out)
}

plot_sppcen <- function(spp.cen, file = "results/spp_cen.pdf"){
    dat <- melt(spp.cen[["cen.spp"]])
    spp <- do.call(rbind, strsplit(as.character(dat[, "X2"]), split = "_"))[, 2]
    cen <- dat[, "value"]
    mu <- tapply(cen, spp, mean)
    se <- tapply(cen, spp, sd) / sqrt(nrow(dat))
    se <- se[order(mu, decreasing = TRUE)]
    mu <- mu[order(mu, decreasing = TRUE)]
    pdf(file = file)
    barplot2(mu, plot.ci = TRUE, ci.u = mu + se, ci.l = mu - se,
             ylab = "Centrality", xlab = "Lichen Species")
    dev.off()
}

plot_netsim <- function(ord, onc.dat, sig.alpha = 1, plot.vectors = FALSE,
                        file = "./cn_chplot.pdf"){
    par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
    pdf(file)
    chp.coord <- ch.plot(ord[["nms"]], onc.dat[, "geno"],
                         cex = 2.5, lwd = 2.5, mu.pch = 15,
                         pt.col = "white",
                         bar.col = grey(0.5)
                         )
    text(chp.coord, labels = rownames(chp.coord), cex = 0.55)
    if (plot.vectors){
        plot(ord[["vec"]], pval = sig.alpha, col = grey(0.01), 
             lwd = 1.0, ascale = 1.25, cex = 1)
    }
    dev.off()
}

plot_mdc <- function(onc.dat, file = "./cn_metrics.pdf"){
    ## Significant Genotype and Network Effects
    pdf(file)
    mdc.plot(onc.dat[, "geno"], onc.dat[, "CT"],
             ylim = c(-1.25, 3),
             xlab = "Tree Genotype", ylab = "Standardized Metric",
             xlas = 2, 
             ord = order(tapply(onc.dat[, "CT"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE)
             )
    mdc.plot(onc.dat[, "geno"], onc.dat[, "BR"],
             add = TRUE, pch = 1,
             ord = order(tapply(onc.dat[, "CT"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE), xjit = 0.005, xlas = 2
             )
    legend("topright", 
           legend = c("Condensed Tannins", "Bark Roughness"), 
           pch = c(19, 1), bty = "none")
    dev.off()
}

plot_h2 <- function(ord, onc.dat, sig.alpha = 1, plot.vectors = FALSE,
                    file = "./cn_trait_h2.pdf"){
    ## Significant Genotype and Network Effects
    pdf(file, width = 9, height = 4.5)
    par(mfrow = c(1, 2), mar = c(5.1, 4.1, 4.1, 2.1))
    chp.coord <- ch.plot(ord[["nms"]], onc.dat[, "geno"],
                         cex = 2.5, lwd = 2.5, mu.pch = 15,
                         pt.col = "white",
                         bar.col = grey(0.5)
                         )
    text(chp.coord, labels = rownames(chp.coord), cex = 0.55)
    legend("topleft", "A", bty = "n", text.font = 2)
    if (plot.vectors){
        plot(ord[["vec"]], pval = sig.alpha, 
             col = grey(0.01), 
             lwd = 1.0, ascale = 2.0, cex = 1)
    }
    ## MDC Plot
    mdc.plot(onc.dat[, "geno"], onc.dat[, "CT"],
             ylim = c(-1.25, 3),
             xlab = "Tree Genotype", ylab = "Standardized Metric",
             xlas = 2, 
             ord = order(tapply(onc.dat[, "CT"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE)
             )
    mdc.plot(onc.dat[, "geno"], onc.dat[, "BR"],
             add = TRUE, pch = 1,
             ord = order(tapply(onc.dat[, "CT"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE), xjit = 0.005, xlas = 2
             )
    legend("topright", 
           legend = c("Condensed Tannins", "Bark Roughness"), 
           pch = c(19, 1), bty = "none")
    legend("topleft", "B", bty = "n", text.font = 2)
    dev.off()
}

## 
plot_nets <- function(cn.onc, onc.dat, file = "./cn_onc.pdf"){
    pdf(file)
    cn.onc <- cn.onc[match(onc.dat[, "tree.id"], names(cn.onc))]
    cn.mu.onc <- tapply(cn.onc, onc.dat[, "geno"], meanNet)
    par(mfrow = c(2, 2), mar = c(0, 0.1, 1.0, 0.1))
    set.seed(123)
    net.col <- sign(meanNet(cn.onc))
    net.col[net.col == -1] <- 2
    net.col[net.col == 1] <- 1
    net.elwd <- (abs(meanNet(cn.onc)) * 10)^2
    coord <- gplot(abs(meanNet(cn.onc)),
                   gmode = "digraph",
                   displaylabels = TRUE,
                   edge.lwd = net.elwd,
                   edge.col = net.col,
                   vertex.col = "black",
                   vertex.cex = 0.5,
                   arrowhead.cex = 0.5,
                   label.cex = 1,
                   main = "All Genotypes"
                   )
    cn.mu.plot <- cn.mu.onc[names(cn.mu.onc) %in%
                            c("996", "WC5", "1008")]
    cn.mu.plot <- cn.mu.plot[order(unlist(lapply(
        cn.mu.plot, function(x) sum(abs(sign(x)))
        )))]
    for (i in 1:length(cn.mu.plot)) {
        net.col <- sign(cn.mu.plot[[i]])
        net.col[net.col == -1] <- 2
        net.col[net.col == 1] <- 1
        net.elwd <- (abs(cn.mu.plot[[i]]) * 10)^2
        set.seed(123)
        gplot(abs(cn.mu.plot[[i]]),
              gmode = "digraph",
              displaylabels = TRUE,
              coord = coord,
              edge.lwd = net.elwd,
              edge.col = net.col,
              vertex.col = "black",
              vertex.cex = 0.5,
              arrowhead.cex = 0.5,
              label.cex = 1,
              main = names(cn.mu.plot)[i]
              )
    }
    dev.off()
}

## Supplementary Figures
proc_size <- function(xgal.size.in){
    xgal.size <- xgal.size.in
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
    return(xgs.data)
}
## X. galericulata size analysis
run_xgsize <- function(xgs.data){
    xgs.median.reml <- lme4::lmer(I(median.thallus^(1 / 4)) ~ (1 | geno),
                                  data = xgs.data[xgs.data$geno %in%
                                      names(which(table(xgs.data$geno) > 2)), ],
                                  REML = TRUE
                                  )
    reml.results <- RLRsim::exactRLRT(xgs.median.reml)
    check.norm <- shapiro.test(residuals(xgs.median.reml))
    check.var <- fligner.test(xgs.data$median.thallus, xgs.data$geno)
    out <- list(reml = reml.results, 
                shapiro = check.norm, 
                fligner = check.var)
    return(out)
}

## Run trait regressions on network metrics
run_trait_nm <- function(onc.dat){
    out <- list()
    out[[1]] <- summary(lm(I(L^(1/4)) ~ BR + CT, data = onc.dat))
    out[[2]] <- summary(lm(I(Cen^(1/4)) ~ BR + CT, data = onc.dat))
    out[[3]] <- summary(lm(I(mod.lik^(1/4)) ~ BR + CT, data = onc.dat))
    names(out) <- c("L", "Cen", "mod.lik")
    return(out)
}

## X. galericulata size plot
plot_xg_size <- function(xgs.data, file = "./xg_size.pdf"){
    pdf(file)
    plot(density(xgs.data$median.thallus),
         xlab = "Median Lichen Thallus Area (cm^2)",
         main = ""
         )
    abline(v = median(xgs.data$median.thallus, 
               na.rm = TRUE), lty = 2)
    dev.off()
}

## Updates the manuscript
update_manuscript <- function(files, dir, file.tex = "main.tex"){
    files = paste0("results/", names(files))
    if (dir.exists(dir)){
        file.copy(
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE,
            from = files,
            to = dir
            )
        texi2pdf(paste(dir, file.tex, sep = "/"), clean = TRUE)
        file.pdf <- gsub(".tex", ".pdf", file.tex, ignore.case = TRUE)
        file.rename(from = file.pdf, to = paste(dir, file.pdf, sep = "/"))
    }else{
        print("Manuscript directory not present.")
    }
}
