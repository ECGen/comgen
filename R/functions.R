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

proc_onc_q <- function(garden.data, rm.zeros = FALSE, rm.nless = 2) {
    g1 <- substr(garden.data[, 1], 2, 2)
    g1[g1 != "P"] <- "onc"
    onc <- garden.data[g1 == "onc", ]
    colnames(onc)[which(colnames(onc) == "Ls")] <- "Lh"
    colnames(onc)[7:ncol(onc)] <- substr(colnames(onc)[7:ncol(onc)], 1, 2)
    onc.q <- split(onc, paste(onc[, 1], onc[, 2]))
    onc.q <- lapply(onc.q, function(x) x[, 7:ncol(x)])
    if (rm.zeros){
        onc.q <- onc.q[unlist(lapply(onc.q, sum)) != 0]
    }
    sr <- unlist(lapply(onc.q, 
                        function(x) sum(sign(apply(x, 2, sum)))))
    onc.q <- onc.q[sr > rm.nless]
    return(onc.q)
}


proc_onc_dat <- function(garden.data, rough.in, onc.q,
                         onc.nc.in, onc.tan.in, onc.ph.in, 
                         rm.na = TRUE) {
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
    ## Network data and metrics
    cn.onc <- lapply(onc.q, coNet, ci.p = 95)
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
        sna::centralization(x, FUN = sna::degree, normalize = TRUE)
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
    spd.onc <- vegan::diversity(onc.com[, colnames(onc.com) != "ds"])
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
    cn.onc <- lapply(onc.q, coNet, ci.p = 95)
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
        sna::centralization(x, FUN = sna::degree, normalize = TRUE)
    }))
    dcen.in.onc <- unlist(lapply(cn.onc, function(x) {
        sna::centralization(x, FUN = sna::degree, cmode = "indegree", normalize = TRUE)
    }))
    dcen.out.onc <- unlist(lapply(cn.onc, function(x) {
        sna::centralization(x, FUN = sna::degree, cmode = "outdegree", normalize = TRUE)
    }))
    dcen.inp.onc <- unlist(lapply(cn.onc, function(x) {
        centralization_signed(x, mode = "in", type = "pos")
    }))
    dcen.outp.onc <- unlist(lapply(cn.onc, function(x) {
        centralization_signed(x, mode = "out", type = "pos")
    }))
    dcen.inn.onc <- unlist(lapply(cn.onc, function(x) {
        centralization_signed(x, mode = "in", type = "neg")
    }))
    dcen.outn.onc <- unlist(lapply(cn.onc, function(x) {
        centralization_signed(x, mode = "out", type = "neg")
    }))
    onc.ns <- cbind(ns.onc, 
                    Cen = dcen.onc, 
                    Cen.in = dcen.in.onc,
                    Cen.out = dcen.out.onc,
                    Cen.in.pos = dcen.inp.onc,
                    Cen.out.pos = dcen.outp.onc,
                    Cen.in.neg = dcen.inn.onc,
                    Cen.out.neg = dcen.outn.onc,
                    mod.lik = cn.mod.onc[, 1], 
                    mod.n = cn.mod.onc[, 2])
    ## adds ascendency
    cn.onc.net <- lapply(cn.onc, as.network)
    asc.l <- lapply(cn.onc.net, enaAscendency)
    asc.df <- do.call(rbind, asc.l)
    onc.ns <- cbind(onc.ns, asc.df)
    ## zero NA network statistics
    onc.ns[is.na(onc.ns)] <- 0
    ## com
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
    spd.onc <- vegan::diversity(onc.com[, colnames(onc.com) != "ds"])
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
        onc.ns[, c("L", 
                   "Cen", "Cen.in", "Cen.out",
                   "Cen.in.pos", "Cen.in.neg", 
                   "Cen.out.pos", "Cen.out.neg", 
                   "mod.lik", 
                   "AMI", "ASC")],
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

proc_onc_com <- function(garden.data, onc.q, onc.dat, rm.na = TRUE) {
    g1 <- substr(garden.data[, 1], 2, 2)
    g1[g1 != "P"] <- "onc"
    onc <- garden.data[g1 == "onc", ]
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


proc_cn_onc <- function(onc.q, ci.p = 95) {
    out <- lapply(onc.q, coNet, ci.p = ci.p)
    names(out) <- do.call(rbind, 
                          strsplit(names(out), split = " "))[, 1]
    return(out)
}

proc_cn_d_onc <- function(cn.onc, onc.dat, method = "euclidean", rm.na = TRUE){
    if (rm.na){
        cn.d.onc <- distNet(
            cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id], 
            method = method
        )
    }else{
        cn.d.onc <- distNet(cn.onc, method = method)
    }
    return(cn.d.onc)
}

trans_cn_d <- function(cn.d.onc, root = 1){
    out <- as.matrix(cn.d.onc)^(1/root)
    out <- as.dist(out)
    return(out)
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

run_spp_centrality <- function(cn.onc, onc.dat, rescale = FALSE,
                               cmode = c("freeman", "in", "out"), 
                               type = c("pos", "neg"), nsim = 100000, 
                               rlrt.seed = 2623, digits = 15){
    if (cmode == "freeman"){
        cen.spp <- lapply(cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id],
                          sna::degree, 
                          rescale = rescale,
                          cmode = cmode
                          )
    }else if (cmode == "in" | cmode == "out"){
        cen.spp <- lapply(cn.onc[names(cn.onc) %in% na.omit(onc.dat)$tree.id],
                          centrality_signed,
                          mode = cmode,
                          type = type
                          )
    }else {
        warning("Error: unknown mode or type.")
    }
    cen.spp <- do.call(rbind, cen.spp)
    cen.spp[is.na(cen.spp)] <- 0
    colnames(cen.spp) <- colnames(cn.onc[[1]])
    if (any(rownames(cen.spp) != onc.dat[, "tree.id"])){stop()}
    colnames(cen.spp) <- paste0("cen_", colnames(cen.spp))
    cen.spp.reml <- list()
    for (i in seq_along(colnames(cen.spp))){
        if (sum(cen.spp[ ,i]) > 0){
            onc.dat[, "cs"] <- cen.spp[ ,i]
            cen.reml <- lme4::lmer(cs^(1 / 4) ~ (1 | geno), 
                                   data = onc.dat, REML = TRUE)
            reml.pval <- RLRsim::exactRLRT(cen.reml, nsim = nsim, seed = rlrt.seed)
            cen.spp.reml[[i]] <- c(
                "spp centrality", 
                mean(cen.spp[, i]),
                reml.pval["statistic"],
                H2 = H2(cen.reml, g = onc.dat[, "geno"]), 
                p.value = reml.pval[["p.value"]]
            )
        }else{cen.spp.reml[[i]] <- c("spp centrality", mean(cen.spp[, i]), rep(NA, 3))}
    }
    spp <- c("X. galericulata", "C. subdeflexa", "L. spp.", "C. holocarpa", "X. montana", 
             "P. melanchra", "P. adscendens", "P. undulata", "R. sp.")
    cen.spp.table <- do.call(rbind, cen.spp.reml)[, -1]
    cen.spp.table <- apply(cen.spp.table, 2, as.numeric)
    cen.spp.table <- round(cen.spp.table, digits)
    cen.spp.table <- cbind(spp, cen.spp.table)
    colnames(cen.spp.table) <- c("lichen species", "mean", "statistic", "H2", "p-value")
    rownames(cen.spp.table) <- colnames(cen.spp)
    out <- list(cen.spp = cen.spp, cen.spp.reml = cen.spp.table)
}

run_reml <- function(onc.dat, trait.results, rm.na = TRUE, raw.reml = FALSE, nsim = 100000, rlrt.seed = 2623){
    if (rm.na){onc.dat <- na.omit(onc.dat)}
    ## tree traits
    prb.reml <- lme4::lmer(I(BR^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    prb.reml.pval <- RLRsim::exactRLRT(prb.reml, nsim = nsim, seed = rlrt.seed)
    prb.reml.result <- c("Percent Rough Bark", 
                         prb.reml.pval["statistic"],
                         H2(prb.reml, g = onc.dat$geno), 
                         R2(prb.reml), prb.reml.pval$p.value)
    ph.reml <- lme4::lmer(I(pH^(1 / 4)) ~ (1 | geno), 
                          data = onc.dat, REML = TRUE)
    ph.reml.pval <- RLRsim::exactRLRT(ph.reml, nsim = nsim, seed = rlrt.seed)
    ph.reml.result <- c("pH", 
                        ph.reml.pval["statistic"],
                        H2(ph.reml, g = onc.dat$geno), 
                        R2(ph.reml), ph.reml.pval$p.value)
    ct.reml <- lme4::lmer(I(CT^(1 / 4)) ~ (1 | geno), 
                          data = onc.dat, REML = TRUE)
    ct.reml.pval <- RLRsim::exactRLRT(ct.reml, nsim = nsim, seed = rlrt.seed)
    ct.reml.result <- c("Condensed Tannins (CT)", 
                        ct.reml.pval["statistic"],
                        H2(ct.reml, g = onc.dat$geno), 
                        R2(ct.reml), ct.reml.pval$p.value)
    cnr.reml <- lme4::lmer(I(log(CN^(1 / 1) + 0.001)) ~ (1 | geno), 
                            data = onc.dat, REML = TRUE)
    cnr.reml.pval <- RLRsim::exactRLRT(cnr.reml, nsim = nsim, seed = rlrt.seed)
    cnr.reml.result <- c("Carbon-Nitrogen (CN) Ratio", 
                         cnr.reml.pval["statistic"],
                         H2(cnr.reml, g = onc.dat$geno), 
                         R2(cnr.reml), 
                         cnr.reml.pval$p.value)
    ## lichen community metrics
    ptc.reml <- lme4::lmer(I(PC^(1 / 2)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    ptc.reml.pval <- RLRsim::exactRLRT(ptc.reml, nsim = nsim, seed = rlrt.seed)
    ptc.reml.result <- c("Percent Lichen Cover", 
                         ptc.reml.pval["statistic"],
                         H2(ptc.reml, g = onc.dat$geno), 
                         R2(ptc.reml), ptc.reml.pval$p.value)
    spr.reml <- lme4::lmer(I(SR^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    spr.reml.pval <- RLRsim::exactRLRT(spr.reml, nsim = nsim, seed = rlrt.seed)
    spr.reml.result <- c("Lichen Species Richness", 
                         spr.reml.pval["statistic"],
                         H2(spr.reml, g = onc.dat$geno), 
                         R2(spr.reml), spr.reml.pval$p.value)
    spd.reml <- lme4::lmer(I(SD^(1 / 4)) ~  (1 | geno), 
                           data = onc.dat, REML = TRUE)
    spd.reml.pval <- RLRsim::exactRLRT(spd.reml, nsim = nsim, seed = rlrt.seed)
    spd.reml.result <- c("Lichen Species Diversity", 
                         spd.reml.pval["statistic"],
                         H2(spd.reml, g = onc.dat$geno), 
                         R2(spd.reml), spd.reml.pval$p.value)
    spe.reml <- lme4::lmer(I(SE^(1 / 4)) ~  (1 | geno), 
                           data = onc.dat, REML = TRUE)
    spe.reml.pval <- RLRsim::exactRLRT(spe.reml, nsim = nsim, seed = rlrt.seed)
    spe.reml.result <- c("Lichen Species Evenness", 
                         spe.reml.pval["statistic"],
                         H2(spe.reml, g = onc.dat$geno), 
                         R2(spe.reml), spe.reml.pval$p.value)
                                        # Lichen network metrics
    link.reml <- lme4::lmer(I(L^(1 / 4)) ~ (1 | geno), 
                            data = onc.dat, REML = TRUE)
    link.reml.pval <- RLRsim::exactRLRT(link.reml, nsim = nsim, seed = rlrt.seed)
    link.reml.result <- c("Number of Network Links", 
                          link.reml.pval["statistic"],
                          H2(link.reml, g = onc.dat$geno), 
                          R2(link.reml), link.reml.pval$p.value)
    cen.reml <- lme4::lmer(I(Cen^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.reml.pval <- RLRsim::exactRLRT(cen.reml, nsim = nsim, seed = rlrt.seed)
    cen.reml.result <- c("Degree Centralization", 
                         cen.reml.pval["statistic"],
                         H2(cen.reml, g = onc.dat$geno), 
                         R2(cen.reml), cen.reml.pval$p.value)
    cen.in.reml <- lme4::lmer(I(Cen.in^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.in.reml.pval <- RLRsim::exactRLRT(cen.in.reml, nsim = nsim, seed = rlrt.seed)
    cen.in.reml.result <- c("In-degree Centralization", 
                            cen.in.reml.pval["statistic"],
                         H2(cen.in.reml, g = onc.dat$geno), 
                         R2(cen.in.reml), cen.in.reml.pval$p.value)
    cen.out.reml <- lme4::lmer(I(Cen.out^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.out.reml.pval <- RLRsim::exactRLRT(cen.out.reml, nsim = nsim, seed = rlrt.seed)
    cen.out.reml.result <- c("Out-degree Centralization", 
                             cen.out.reml.pval["statistic"],
                         H2(cen.out.reml, g = onc.dat$geno), 
                         R2(cen.out.reml), cen.out.reml.pval$p.value)
    cen.inp.reml <- lme4::lmer(I(Cen.in.pos^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.inp.reml.pval <- RLRsim::exactRLRT(cen.inp.reml, nsim = nsim, seed = rlrt.seed)
    cen.inp.reml.result <- c("In-Positive Centralization", 
                             cen.inp.reml.pval["statistic"],
                             H2(cen.inp.reml, g = onc.dat$geno), 
                             R2(cen.inp.reml), cen.inp.reml.pval$p.value)
    cen.inn.reml <- lme4::lmer(I(Cen.in.neg^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.inn.reml.pval <- RLRsim::exactRLRT(cen.inn.reml, nsim = nsim, seed = rlrt.seed)
    cen.inn.reml.result <- c("In-Negative Centralization", 
                             cen.inn.reml.pval["statistic"],
                             H2(cen.inn.reml, g = onc.dat$geno), 
                             R2(cen.inn.reml), cen.inn.reml.pval$p.value)
    cen.outp.reml <- lme4::lmer(I(Cen.out.pos^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.outp.reml.pval <- RLRsim::exactRLRT(cen.outp.reml, nsim = nsim, seed = rlrt.seed)
    cen.outp.reml.result <- c("Out-Positive Centralization", 
                              cen.outp.reml.pval["statistic"],
                              H2(cen.outp.reml, g = onc.dat$geno), 
                              R2(cen.outp.reml), 
                              cen.outp.reml.pval$p.value)
    cen.outn.reml <- lme4::lmer(I(Cen.out.neg^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    cen.outn.reml.pval <- RLRsim::exactRLRT(cen.outn.reml, nsim = nsim, seed = rlrt.seed)
    cen.outn.reml.result <- c("Out-Negative Centralization", 
                              cen.outn.reml.pval["statistic"],
                              H2(cen.outn.reml, g = onc.dat$geno), 
                              R2(cen.outn.reml), cen.outn.reml.pval$p.value)
    ami.reml <- lme4::lmer(I(AMI^(1 / 4)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    ami.reml.pval <- RLRsim::exactRLRT(ami.reml, nsim = nsim, seed = rlrt.seed)
    ami.reml.result <- c("Average Mutual Information", 
                         ami.reml.pval["statistic"],
                         H2(ami.reml, g = onc.dat$geno), 
                         R2(ami.reml), 
                         ami.reml.pval$p.value)
    asc.reml <- lme4::lmer(I(ASC^(1 / 2)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    asc.reml.pval <- RLRsim::exactRLRT(asc.reml, nsim = nsim, seed = rlrt.seed)
    asc.reml.result <- c("Network Ascendency", 
                         asc.reml.pval["statistic"],
                         H2(asc.reml, g = onc.dat$geno), 
                         R2(asc.reml), 
                         asc.reml.pval$p.value)
    mod.reml <- lme4::lmer(I(mod.lik^(2 / 1)) ~ (1 | geno), 
                           data = onc.dat, REML = TRUE)
    mod.reml.pval <- RLRsim::exactRLRT(mod.reml, nsim = nsim, seed = rlrt.seed)
    mod.reml.result <- c("Network Modularity", 
                         mod.reml.pval["statistic"],
                         H2(mod.reml, g = onc.dat$geno), 
                         R2(mod.reml), 
                         mod.reml.pval$p.value)
    ## Residual bark roughness effect on net
    res.df <- data.frame(geno = onc.dat[, "geno"],
                    resL = residuals(trait.results[["br_L"]]),
                    resCen = residuals(trait.results[["br_Cen"]]),
                    resAMI = residuals(trait.results[["br_AMI"]]))
    resL.reml <- lme4::lmer(resL ~ (1 | geno), 
                            data = res.df, REML = TRUE)
    resL.reml.pval <- RLRsim::exactRLRT(resL.reml, nsim = nsim, seed = rlrt.seed)
    resL.reml.result <- c("BR-L Residuals", 
                          resL.reml.pval["statistic"],
                          H2(resL.reml, g = res.df$geno), 
                          R2(resL.reml), 
                          resL.reml.pval$p.value)
    resCen.reml <- lme4::lmer(resCen ~ (1 | geno), 
                           data = res.df, REML = TRUE)
    resCen.reml.pval <- RLRsim::exactRLRT(resCen.reml, nsim = nsim, seed = rlrt.seed)
    resCen.reml.result <- c("BR-Cen Residuals", 
                            resCen.reml.pval["statistic"],
                          H2(resCen.reml, g = res.df$geno), 
                          R2(resCen.reml), 
                          resCen.reml.pval$p.value)
    resAMI.reml <- lme4::lmer(resAMI ~ (1 | geno), 
                           data = res.df, REML = TRUE)
    resAMI.reml.pval <- RLRsim::exactRLRT(resAMI.reml, nsim = nsim, seed = rlrt.seed)
    resAMI.reml.result <- c("BR-AMI Residuals", 
                            resAMI.reml.pval["statistic"],
                          H2(resAMI.reml, g = res.df$geno), 
                          R2(resAMI.reml), 
                          resAMI.reml.pval$p.value)
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
                    cen.reml,
                    cen.in.reml,
                    cen.out.reml,
                    cen.inp.reml,
                    cen.inn.reml,
                    cen.outp.reml,
                    cen.outn.reml,
                    ami.reml,
                    asc.reml,
                    resL.reml, 
                    resCen.reml,
                    resAMI.reml)
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
                     cen.reml.result,
                     cen.in.reml.result,
                     cen.out.reml.result,
                     cen.inp.reml.result,
                     cen.inn.reml.result,
                     cen.outp.reml.result,
                     cen.outn.reml.result,
                     ami.reml.result,
                     asc.reml.result,
                     resL.reml.result,
                     resCen.reml.result,
                     resAMI.reml.result)
        colnames(out) <- c("response", "statistic", "H2", "R2", "p-value")
    }
    return(out)
}

std <- function(x){
    (x - mean(x)) / sqrt(length(x))
}

run_trait_path <- function(onc.dat){
    out <- list()
    ## Direct effect of BR on network
    out[["br_L"]] <- lm(L ~ BR, data = onc.dat)
    out[["br_Cen"]] <- lm(Cen ~ BR, data = onc.dat)
    out[["br_AMI"]] <- lm(AMI ~ BR, data = onc.dat)
    ## Direct effect of CT on network
    out[["ct_L"]] <- lm(L ~ CT, data = onc.dat)
    out[["ct_Cen"]] <- lm(Cen ~ CT, data = onc.dat)
    out[["ct_AMI"]] <- lm(AMI ~ CT, data = onc.dat)
    ## Direct effect of pH on network
    out[["ph_L"]] <- lm(L ~ pH, data = onc.dat)
    out[["ph_Cen"]] <- lm(Cen ~ pH, data = onc.dat)
    out[["ph_AMI"]] <- lm(AMI ~ pH, data = onc.dat)
    ## Direct effect of CN on network
    out[["cn_L"]] <- lm(L ~ CN, data = onc.dat)
    out[["cn_Cen"]] <- lm(Cen ~ CN, data = onc.dat)
    out[["cn_AMI"]] <- lm(AMI ~ CN, data = onc.dat)
    return(out)
}


run_SEM <- function(onc.dat, cn.d.onc, np = 100000){
    
    CT.d <- dist(std(onc.dat[, "CT"]))
    BR.d <- dist(std(onc.dat[, "BR"]))
    PC.d <- dist(std(onc.dat[, "PC"]))
    SR.d <- dist(std(onc.dat[, "SR"]))
    SE.d <- dist(std(onc.dat[, "SE"]))
    pH.d <- dist(std(onc.dat[, "pH"]))
    CN.d <- dist(std(onc.dat[, "CN"]))
    g.m <- model.matrix(~ geno - 1, data = onc.dat)
    geno.d <- dist(g.m)

    adonis2(cn.d.onc^(1/4) ~ geno, 
                          data = onc.dat, mrank = TRUE, nperm = np)
    
    ## geno -> net
    mantel(cn.d.onc^(1/4) ~ geno.d, nperm = np)
    mantel(cn.d.onc^(1/4) ~ geno.d + SR.d + SE.d + 
               PC.d + CT.d + BR.d + pH.d + CN.d, 
           nperm = np)
    ## geno -> lichen stats
    mantel(SR.d^(1/4) ~ geno.d, nperm = np)
    mantel(SE.d^(1/4) ~ geno.d, nperm = np)
    mantel(PC.d^(1/4) ~ geno.d, nperm = np)
    ## geno -> trait
    mantel(CT.d^(1/1) ~ geno.d, nperm = np)
    mantel(BR.d^(1/1) ~ geno.d, nperm = np)
    mantel(pH.d^(1/4) ~ geno.d, nperm = np)
    mantel(CN.d^(1/4) ~ geno.d, nperm = np)
    ## Trait -> SR
    set.seed(70); mantel(SR.d ~ CT.d + pH.d + CN.d + BR.d, nperm = np)
    set.seed(70); mantel(SR.d ~ CT.d + pH.d + CN.d + BR.d + geno.d, nperm = np)
    set.seed(70); mantel(SR.d ~ BR.d + CT.d + pH.d + CN.d + geno.d, nperm = np)
    set.seed(70); mantel(SR.d ~ CT.d + pH.d + CN.d + BR.d + geno.d, nperm = np)
    set.seed(70); mantel(SR.d ~ pH.d + CN.d + BR.d + CT.d + geno.d, nperm = np)
    set.seed(70); mantel(SR.d ~ CN.d + BR.d + CT.d + pH.d + geno.d, nperm = np)
    ## SR -> net
    set.seed(70); mantel(cn.d.onc^(1/4) ~ SR.d + PC.d + SE.d + 
                             BR.d + pH.d + CN.d, nperm = np)
    set.seed(70); mantel(cn.d.onc^(1/4) ~ SR.d + PC.d + SE.d + 
                             BR.d + pH.d + CN.d + CT.d + geno.d, nperm = np)
    set.seed(70); mantel(cn.d.onc^(1/4) ~ PC.d + SR.d + SE.d, nperm = np)
    set.seed(70); mantel(cn.d.onc^(1/4) ~ SE.d + PC.d + SR.d, nperm = np)
    ## Trait -> net
    set.seed(70); mantel(cn.d.onc^(1/4) ~ BR.d + CT.d + pH.d + CN.d, nperm = np)
    set.seed(70); mantel(cn.d.onc^(1/4) ~ CT.d + pH.d + CN.d + BR.d, nperm = np)
    set.seed(70); mantel(cn.d.onc^(1/4) ~ pH.d + CN.d + BR.d + CT.d, nperm = np)
    set.seed(70); mantel(cn.d.onc^(1/4) ~ CN.d + BR.d + CT.d + pH.d, nperm = np)

}

run_perm <- function(onc.dat, onc.com, cn.d.onc){
    com.perm <- vegan::adonis2((onc.com^(2)) ~ geno,
                               data = onc.dat, 
                               by = "term",
                               mrank = TRUE,
                               perm = 100000)
    cn.perm <- vegan::adonis2(cn.d.onc^(2) ~ geno,
                              by = "term", 
                              data = onc.dat, 
                              mrank = TRUE,
                              permutations = 100000)
    cn.trait.perm <- vegan::adonis2(
                                cn.d.onc^(2)~BR+CT+pH+CN,
                                by = "term", 
                              data = onc.dat, 
                              mrank = TRUE,
                              permutations = 100000)
    out <- list(com = com.perm, 
                cn = cn.perm,
                cn.trait = cn.trait.perm)
    return(out)
}

make_table_sppcen <- function(spp.cen.pos.in, spp.cen.pos.out, spp.cen.neg.in, spp.cen.neg.out, 
                              xtab = TRUE, digits = 4){
    sc <- list(c("Positive", rep(NA, (ncol(spp.cen.pos.in[["cen.spp.reml"]]) - 1))),
               c("In-Degree", rep(NA, (ncol(spp.cen.pos.in[["cen.spp.reml"]]) - 1))),
               spp.cen.pos.in[["cen.spp.reml"]], 
               c("Out-Degree", rep(NA, (ncol(spp.cen.pos.in[["cen.spp.reml"]]) - 1))),
               spp.cen.pos.out[["cen.spp.reml"]], 
               c("Negative", rep(NA, (ncol(spp.cen.pos.in[["cen.spp.reml"]]) - 1))),
               c("In-Degree", rep(NA, (ncol(spp.cen.pos.in[["cen.spp.reml"]]) - 1))),
               spp.cen.neg.in[["cen.spp.reml"]], 
               c("Out-Degree", rep(NA, (ncol(spp.cen.pos.in[["cen.spp.reml"]]) - 1))),
               spp.cen.neg.out[["cen.spp.reml"]]
               )
    out <- do.call(rbind, sc)
    if (xtab){
        out  <-  xtable(
            out, 
            caption = "REML tests of the effect of tree genotype on lichen species centrality.",
            label = "tab:sppcen",
            type = "latex",
            digits = digits
        )
    }
    return(out)
}

make_table_vectors <- function(vec, xtab = TRUE, digits = 3){
    out <- vec[, c("r", "pval")]
    rownames(out) <- c("Bark Roughness",
                       "Number of Links",
                       "Centralization",
                       "AMI")
    colnames(out) <- c("r", "p-value")
    if (xtab){
        out  <-  xtable(
            out, 
            caption = "Correlation tests for vectors displayed in NMDS ordination of network similarity.",
            label = "tab:vec",
            type = "latex",
            digits = digits
        )
    }
    return(out)
}

make_table_path<- function(trait.results, onc.dat, digits = 7, xtab = TRUE){
    ## genotype -> bark roughness
    ## genotype (-> br) -> net(L, Cen, AMI)
    ## correlations
    cm <- cor.mat(onc.dat[, c("BR", "CT", "pH", "CN", "L", "Cen", "AMI")], 
                  digits = digits, sig.only = FALSE, p.val = TRUE)[["r"]]
    r <- as.vector(cm[c("L", "Cen", "AMI"), c("BR", "CT", "pH", "CN")])
    ## trait -> net(L, Cen, AMI)
    out <- lapply(trait.results, function(x) as.matrix(summary(x)[["coefficients"]])[2, ])
    r2 <- sapply(trait.results, function(x) summary(x)[["r.squared"]])
    ## output binding
    out <- cbind(r, r2, do.call(rbind, out))
    colnames(out) <- c("r", "R2", "estimate", "SE", "t", "p-value")
    if (xtab){
        out  <-  xtable(
            out, 
            caption = "Tests of the correlation between tree bark traits and lichen network structure",
            label = "tab:trait_path",
            type = "latex",
            include.rownames = TRUE,
            include.colnames = TRUE
        )
    }
    return(out)
}

make_tables <- function(onc.dat, reml.results, perm.results, digits = 4){
    ## Heritability table
    h2.tab <- reml.results
    ## Add PERMANOVA results
    com.perm.h2 <- c("Community Composition", 
                     as.data.frame(perm.results[["com"]])["geno", "F"],
                     H2(perm.results[["com"]], 
                        g = onc.dat[["geno"]]),
                     R2(perm.results[["com"]]),
                     as.data.frame(perm.results[["com"]])["geno", "Pr(>F)"]
                     )
    cn.perm.h2 <- c("Lichen Network Similarity", 
                    as.data.frame(perm.results[["cn"]])["geno", "F"],
                    H2(perm.results[["cn"]], 
                       g = onc.dat[, "geno"], 
                       perm = 10000),
                    R2(perm.results[["cn"]]),
                    as.data.frame(perm.results[["cn"]])["geno", "Pr(>F)"]
                    )
    h2.tab <- rbind(h2.tab, 
                    cn.perm.h2,
                    com.perm.h2)
    ## Format Heritability Table
    h2.tab[, c("statistic", "H2", "R2", "p-value")] <- apply(
        h2.tab[, c("statistic", "H2", "R2", "p-value")], 
        2, 
        function(x, digits) round(as.numeric(x), digits = digits), 
        digits = digits
    )
    h2.tab <- na.omit(h2.tab)
    ## Oragnize H2 table
    h2.tab <- h2.tab[c("cn.perm.h2",
                       "ami.reml.result",
                       "cen.reml.result",
                       "cen.in.reml.result",
                       "cen.out.reml.result",
                       "cen.inp.reml.result",
                       "cen.inn.reml.result",
                       "cen.outp.reml.result",
                       "cen.outn.reml.result",
                       "link.reml.result",
                       "ptc.reml.result",
                       "spd.reml.result",
                       "spr.reml.result",
                       "spe.reml.result",
                       "prb.reml.result",
                       "ph.reml.result",
                       "cnr.reml.result",
                       "ct.reml.result",
                       "resL.reml.result",
                       "resCen.reml.result",
                       "resAMI.reml.result"), ]
                                        # Lichen Networks
                                        # Lichen Network Metrics
                                        # Lichen Community
                                        # Tree Traits
    ## Remove R2 from H2 table
    h2.tab <- h2.tab[, colnames(h2.tab) != "R2"]
    ## Format lichen network permanova table
    cn.perm <- as.data.frame(perm.results[["cn"]])
    cn.trait.perm <- as.data.frame(perm.results[["cn.trait"]])
    ## rownames(cn.perm) <- c("Genotype", "Bark Roughness", "pH", 
    ##                        "C:N Ratio", "Condensed Tannins", 
    ##                        "Percent Cover", "Species Richness",
    ##                        "Species Evenness", "Number of Links", 
    ##                        "Network Modularity", "Network Centrality", 
    ##                        "Residual", "Total")
    colnames(cn.perm) <- c("df", "SS", "R2", "F", "p-value")
    ## Create the latex
    tab.h2 <- xtable::xtable(
       h2.tab,
       caption = "Genotypic effects on tree traits and bark lichen.",
       label = "tab:h2_table",
       type = "latex",
       include.rownames = FALSE,
       include.colnames = TRUE, 
        digits = digits
       )
    tab.h2.net <- xtable::xtable(
       h2.tab[c("cn.perm.h2",
                "ami.reml.result",
                "cen.reml.result",
                "cen.in.reml.result",
                "cen.inp.reml.result",
                "cen.inn.reml.result",
                "cen.out.reml.result",
                "cen.outp.reml.result",
                "cen.outn.reml.result",
                "link.reml.result"), ],
       caption = "Genotypic effects on the associated lichen network structure.",
       label = "tab:h2_net",
       type = "latex",
       include.rownames = FALSE,
       include.colnames = TRUE, 
        digits = digits
       )
    tab.h2.trait <- xtable::xtable(
       h2.tab[c("prb.reml.result",
                "ph.reml.result",
                "cnr.reml.result",
                "ct.reml.result",
                "resL.reml.result",
                "resCen.reml.result",
                "resAMI.reml.result"), ],
       caption = "Genotypic effects on tree traits and residuals from trait regressions of lichen network structure.",
       label = "tab:h2_trait",
       type = "latex",
       include.rownames = FALSE,
       include.colnames = TRUE, 
        digits = digits
       )
     tab.cn.perm <- xtable::xtable(
        cn.perm,
        caption = "PERMANOVA Pseudo-F Table of lichen network similarity to genotype.",
        label = "tab:cn_perm",
        type = "latex", 
        include.rownames = TRUE,
        include.colnames = TRUE, 
        digits = digits
        )
    tab.cn.trait.perm <- xtable::xtable(
        cn.trait.perm,
        caption = "PERMANOVA Pseudo-F Table of lichen network similarity response to bark traits.",
        label = "tab:cn_trait_perm",
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
                h2_net = tab.h2.net, 
                h2_trait = tab.h2.trait, 
                cn = tab.cn.perm, 
                cn_trait = tab.cn.trait.perm,
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
    spp.cen.dat <- melt(spp.cen[["cen.spp"]])
    colnames(spp.cen.dat) <- c("tree.id", "species", "centrality")
    sppcen.aov <- lm(centrality ~ species, data = spp.cen.dat)
    sppcen.mct <- TukeyHSD(aov(centrality ~ species, data = spp.cen.dat))
    out <- list(aov = sppcen.aov, mct = sppcen.mct)
    return(out)
}

plot_sppcen <- function(spp.cen, file = "results/spp_cen.pdf", ylab = "Centrality"){
    dat <- melt(spp.cen[["cen.spp"]])
    spp <- do.call(rbind, strsplit(as.character(dat[, "X2"]), split = "_"))[, 2]
    cen <- dat[, "value"]
    mu <- tapply(cen, spp, mean)
    se <- tapply(cen, spp, sd) / sqrt(nrow(dat))
    se <- se[order(mu, decreasing = TRUE)]
    mu <- mu[order(mu, decreasing = TRUE)]
    pdf(file = file)
    barplot2(mu, plot.ci = TRUE, ci.u = mu + se, ci.l = mu - se,
             ylab = ylab, xlab = "Lichen Species")
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
             lwd = 1.0, ascale = 1.00, cex = 1)
    }
    dev.off()
}

plot_mdc <- function(onc.dat, file = "./cn_metrics.pdf"){
    ## Significant Genotype and Network Effects
    pdf(file)
    mdc.plot(onc.dat[, "geno"], onc.dat[, "AMI"],
             ylim = c(-1.25, 3),
             xlab = "Tree Genotype", ylab = "Standardized Metric",
             xlas = 2, 
             ord = order(tapply(onc.dat[, "AMI"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE)
             )
    mdc.plot(onc.dat[, "geno"], onc.dat[, "BR"],
             add = TRUE, pch = 1,
             ord = order(tapply(onc.dat[, "AMI"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE), xjit = 0.005, xlas = 2
             )
    legend("topright", 
           legend = c("AMI", "Bark Roughness"), 
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
             lwd = 1.0, cex = 0.75)
    }
    ## MDC Plot
    mdc.plot(onc.dat[, "geno"], onc.dat[, "AMI"],
             ylim = c(-1.5, 2.5),
             xlab = "Tree Genotype", ylab = "Standardized Metric",
             xlas = 2, 
             ord = order(tapply(onc.dat[, "AMI"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE),
             std = TRUE
             )
    mdc.plot(onc.dat[, "geno"], onc.dat[, "Cen"],
             add = TRUE, pch = 1,
             ord = order(tapply(onc.dat[, "AMI"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE), 
             xjit = 0.01, xlas = 2,
             std = TRUE
             )
    mdc.plot(onc.dat[, "geno"], onc.dat[, "L"],
             add = TRUE, pch = 2,
             ord = order(tapply(onc.dat[, "AMI"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE), 
             xjit = 0.02, xlas = 2,
             std = TRUE
             )
    legend("topright", 
           legend = c("AMI", "Centralization", "Number of Links"), 
           pch = c(19, 1, 2), bty = "none")
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

plot_br_net <- function(onc.dat, file = "./results/br_net.pdf", cex = 2.5, lwd = 1.5, lab.cex = 1.5, no.geno.line = FALSE){
    gmu <- data.frame(apply(onc.dat[, c("BR", "L", "Cen", "AMI")], 2, function(x, y) tapply(x, y, mean), y = onc.dat[, "geno"]))
    pdf(file, width = 15, height = 5)
    par(mfrow = c(1, 3), mar = c(5.1, 4.1, 4.1, 2.1), cex.lab = 1.5, cex.axis = 1.25)
    chp.coord <- ch.plot(onc.dat[, c("BR", "L")], onc.dat[, "geno"],
                         cex = cex, lwd = lwd, mu.pch = 15,
                         pt.col = "white",
                         bar.col = "black",
                         xlab = "Bark Roughness", ylab = "Number of Links (L)"
                         )
    text(chp.coord, labels = rownames(chp.coord), cex = lab.cex)
    abline(lm(L ~ BR, data = gmu), lty = 1)
    if (no.geno.line){abline(lm(L ~ BR, data = onc.dat), lty = 2)}
    legend("topleft", "A", bty = "n", text.font = 2, cex = 2)
    chp.coord <- ch.plot(onc.dat[, c("BR", "Cen")], onc.dat[, "geno"],
                         cex = cex, lwd = lwd, mu.pch = 15,
                         pt.col = "white",
                         bar.col = "black",
                         xlab = "Bark Roughness", ylab = "Centralization"
                         )
    text(chp.coord, labels = rownames(chp.coord), cex = lab.cex)
    abline(lm(Cen ~ BR, data = gmu), lty = 1)
    if (no.geno.line){abline(lm(Cen ~ BR, data = onc.dat), lty = 2)}
    legend("topleft", "B", bty = "n", text.font = 2, cex = 2)
    chp.coord <- ch.plot(onc.dat[, c("BR", "AMI")], onc.dat[, "geno"],
                         cex = cex, lwd = lwd, mu.pch = 15,
                         pt.col = "white",
                         bar.col = "black",
                         xlab = "Bark Roughness", ylab = "Average Mutual Information (AMI)"
                         )
    text(chp.coord, labels = rownames(chp.coord), cex = lab.cex)
    abline(lm(AMI ~ BR, data = gmu), lty = 1)
    if (no.geno.line){abline(lm(AMI ~ BR, data = onc.dat), lty = 2)}
    legend("topleft", "C", bty = "n", text.font = 2, cex = 2)
    dev.off()
}

## Supplementary Figures
plot_geno_sppcen <- function(onc.dat, spp.cen.pos.in, spp.cen.pos.out, file = "./results/geno_sppcen.pdf"){
    ylim.min <- 0
    ylim.max <- max(cbind(spp.cen.pos.in[["cen.spp"]][, c("cen_Ch", "cen_Xm")], 
                          spp.cen.pos.out[["cen.spp"]][, c("cen_Ch", "cen_Xm")]))
    pdf(file, width = 9, height = 4.5)
    par(mfrow = c(1, 2), mar = c(5.1, 4.1, 4.1, 2.1), cex.lab = 1.0, cex.axis = 1.0)
    mdc.plot(onc.dat[, "geno"], spp.cen.pos.in[["cen.spp"]][, "cen_Ch"],
             ylim = c(ylim.min, ylim.max),
             xlab = "Tree Genotype", ylab = "Centraliity (In)",
             ord = order(tapply(spp.cen.pos.in[["cen.spp"]][, "cen_Ch"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE),
             std = FALSE,
             xjit = -0.005,
             xlas = 2
             )
    mdc.plot(onc.dat[, "geno"], spp.cen.pos.in[["cen.spp"]][, "cen_Xm"],
             add = TRUE, pch = 1,
             ord = order(tapply(spp.cen.pos.in[["cen.spp"]][, "cen_Ch"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE), 
             std = FALSE,
             xjit = 0.005, 
             xlas = 2
             )
    legend("topright", 
           legend = c("C. holocarpa", "X. montana"), 
           pch = c(19, 1), bty = "none")
    legend("topleft", "A", bty = "n", text.font = 2)
    mdc.plot(onc.dat[, "geno"], spp.cen.pos.out[["cen.spp"]][, "cen_Ch"],
             ylim = c(ylim.min, ylim.max),
             xlab = "Tree Genotype", ylab = "Centraliity (Out)",
             ord = order(tapply(spp.cen.pos.out[["cen.spp"]][, "cen_Ch"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE),
             std = FALSE, 
             xjit = -0.005,
             xlas = 2
             )
    mdc.plot(onc.dat[, "geno"], spp.cen.pos.out[["cen.spp"]][, "cen_Xm"],
             add = TRUE, pch = 1,
             ord = order(tapply(spp.cen.pos.out[["cen.spp"]][, "cen_Ch"], 
                 onc.dat[, "geno"], mean), 
                 decreasing = TRUE), 
             std = FALSE,
             xjit = 0.005, 
             xlas = 2
             )
    legend("topright", 
           legend = c("C. holocarpa", "X. montana"), 
           pch = c(19, 1), bty = "none")
    legend("topleft", "B", bty = "n", text.font = 2)
    dev.off()
    
}

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
run_xgsize <- function(xgs.data, nsim = 100000, rlrt.seed = 2623){
    xgs.median.reml <- lme4::lmer(I(median.thallus^(1/4)) ~ (1 | geno),
                                  data = xgs.data[xgs.data$geno %in%
                                      names(which(table(xgs.data$geno) > 2)), ],
                                  REML = TRUE
                                  )
    reml.results <- RLRsim::exactRLRT(xgs.median.reml, nsim = nsim, seed = rlrt.seed)
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

## Species accumulation curves by genotype
emp_spac <- function(x){
    spac <- lapply(x, specaccum)
    rich <- lapply(spac, function(x) x[["richness"]])
    mu.rich <- apply(do.call(rbind, rich), 2, mean)
    sd.rich <- apply(do.call(rbind, rich), 2, sd)
    freq <- ldply(lapply(spac, function(x) x[["freq"]]), rbind)
    freq[is.na(freq)] <- 0
    mu.freq <- apply(freq[, -1], 2, mean)
    out <- spac[[1]]
    out[["richness"]] <- mu.rich
    out[["sd"]] <- sd.rich
    out[["freq"]] <- mu.freq
    return(out)
}

spac_geno <- function(onc.q, onc.dat){
    tapply(onc.q, onc.dat[, "geno"], emp_spac)
}

plot_spag <- function(spac.geno, file = "./spac_geno.pdf", y.max = 10){
    if (!(exists("y.max"))){
        mu.max <- max(do.call(rbind, lapply(spac.geno, function(x) x[["richness"]])))
        sd.max <- max(do.call(rbind, lapply(spac.geno, function(x) x[["sd"]])))
        y.max <- mu.max + (1.96 * (sd.max / sqrt(length(spac.geno))))
    }
    pdf(file)
    ## plot(spac.geno[[1]], ci.type = "polygon", ci.col = "lightgrey", ylim = c(0, y.max))
    ## lapply(spac.geno[-1], plot, add = TRUE, ci.type = "polygon", ci.col = "lightgrey")
    plot(spac.geno[[1]], ci.type = "bar", 
         ylim = c(0, y.max), xlab = expression("Cumulative Area Sampled cm"^2), ylab = "Lichen Species")
    lapply(seq(2, length(spac.geno)), 
               function(i, x, col) 
                   plot(x[[i]], add = TRUE, ci.type = "bar", col = col[i]),
               x = spac.geno, col = rainbow(length(spac.geno) - 1)
           )
    dev.off()
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

## cortest matrix
cor.mat <- function(x, digits = 2, sig.only = TRUE, alpha = 0.05, p.val = FALSE){
    cm <- array(NA, dim = rep(ncol(x), 2))
    rownames(cm) <- colnames(cm) <- colnames(x)
    cm.p <- array(NA, dim = rep(ncol(x), 2))
    for (i in seq(1, ncol(x))){
        for (j in seq(1, ncol(x))){
        ct <- cor.test(x[, i], x[, j])
        cm[i, j] <- unlist(ct["estimate"])
        cm.p[i, j] <- unlist(ct["p.value"])
        }
    }
    cm <- round(cm, digits)
    if (sig.only){
        out <- cm
        out[cm.p >= alpha] <- NA
    }else if (p.val){
        out <- list(r = cm, p = cm.p)
    }else{
        out <- cm
    }
    return(out)
}

## cormat table
cormat_tab <- function(onc.dat, upper = TRUE, xtab = TRUE, digits = 2){
    out <- cor.mat(onc.dat[,c("BR", "CT", "pH", "CN", "PC","SR","SE","SD","L","Cen","AMI")])
    if (upper){
        out[lower.tri(out)] <- NA
        diag(out) <- NA
    }
    if (xtab){
        out <- xtable(out, 
               caption = "Matrix of correlations among tree traits, lichen community metrics and network metrics",
               label = "tab:cormat",
               type = "latex",               
               digits = digits)
    }
   return(out)
}

### Signed network metrics
as_graph_signed <- function(x){
    g = graph_from_adjacency_matrix(x, mode = "directed", weighted = TRUE)
    sv  <- numeric()
    for (i in seq(1, nrow(x))){
        for (j in seq(1, ncol(x))){
            if (x[i, j] != 0){sv  <- append(sv, sign(x[i, j]))}
        }
    }
    E(g)$sign <- sv
    return(g)
}

freeman <- function(x){
    sum(max(x) - x) 
}

centrality_signed <- function(x, mode = c("in", "out"), type = c("pos", "neg", "ratio")){
    if (sum(abs(sign(x))) == 0){
        out <- rep(0, nrow(x))
    }else{
        g <- as_graph_signed(x)    
        out <- signnet::degree_signed(g, mode = mode, type = type)
    }
    return(out)
}

centralization_signed <- function(x, mode = c("in", "out"), type = c("pos", "neg", "ratio")){
    if (sum(abs(sign(x))) == 0){
        out <- 0
    }else{
        g <- as_graph_signed(x)    
        out <- freeman(signnet::degree_signed(g, mode = mode, type = type))
    }
    return(out)
}

## Updates the manuscript
update_manuscript <- function(files, dir, file.tex = "main.tex", render = FALSE){
    files = paste0("results/", names(files))
    if (dir.exists(dir)){
        file.copy(
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE,
            from = files,
            to = dir
            )
        if (render){
            texi2pdf(paste(dir, file.tex, sep = "/"), clean = TRUE)
            file.pdf <- gsub(".tex", ".pdf", file.tex, ignore.case = TRUE)
            file.rename(from = file.pdf, to = paste(dir, file.pdf, sep = "/"))
        }
    }else{
        print("Manuscript directory not present.")
        print("Setting up submodule:")
        system("git add submodule https://github.com/ECGen/lcn_manuscript docs/lcn_manuscript")
        print("Re-run make.R")
    }
}

