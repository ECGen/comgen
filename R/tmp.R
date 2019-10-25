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
    onc.ph[!is.na(onc.ph[, "pH2"]), "pH"] <- apply(
        onc.ph[!is.na(onc.ph[, "pH2"]), c("pH", "pH2")], 1, mean)

    onc.dat <- data.frame(tree.id, 
                          PC = ptc.onc, 
                          SR = spr.onc, 
                          SD = spd.onc, 
                          SE = spe.onc, 
                          geno = factor(onc.geno), 
                          tree = tree, 
                          BR = onc.rough, 
                          onc.ns[, c("L", "Cen")])
    onc.ph <- onc.ph[onc.ph[, "tree.id"] %in% 
                     onc.dat[, "tree.id"], ]
    onc.ph <- onc.ph[match(onc.dat[, "tree.id"], 
                           onc.ph[, "tree.id"]), ]

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

