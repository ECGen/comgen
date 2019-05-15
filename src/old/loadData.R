
### Garden
gard <- read.csv('../data/lichen_networks/LCO_data_ONC_PIT.csv')
g.tree <- gard[, 1:4]
g.net <- gard[, 7:ncol(gard)]
g.net <- split(g.net, gard$Tree)
g.net <- lapply(g.net, as.matrix)
g.net <- lapply(g.net, coNets)
g.nd <- netDist(g.net)
g.mn <- Reduce("+", g.net) / sum(unlist(g.net))
g.info <- split(g.tree, gard$Tree)
g.info <- do.call(rbind, lapply(g.info, head, n = 1))
garden <- grepl("NP", g.info[, "Tree"])
garden[garden == FALSE] <- "onc"
garden[garden == "TRUE"] <- "pit"
g.info <- cbind(g.info, Garden = garden)

### Wild
                                        #remove notes and NA
                                        #remove gnu.44  =  FREMONT and ll.6  =  weird
                                        #condense species
                                        #use only the 45.55 quadrats
wild <- read.csv('../data/lichen_networks/lco_Apr2012.csv')
wild <- wild[, colnames(wild) != 'NOTES.']
wild <- wild[, colnames(wild) !=  'dead']
wild <- na.omit(wild)
wild <- wild[wild$tree != 'gnu.44', ]
wild <- wild[wild$tree != 'll.6', ]
lec.sp <- apply(wild[, c(6, 8, 10, 18)], 1, function(x) sign(any(x != 0)))
wild <- cbind(wild, lec = lec.sp)
wild <- wild[, -c(6, 8, 10, 18)]
wild <- wild[wild$quadrat  !=  'n80.90', ]
w.tree <- wild[, 1:4]
w.net <- wild[, 5:ncol(wild)]
w.net <- split(w.net, wild$tree)
w.net <- lapply(w.net, as.matrix)
w.net <- lapply(w.net, coNets)
w.nd <- netDist(w.net)
w.mn <- Reduce("+", w.net) / sum(unlist(w.net))
w.info <- split(w.tree, wild$tree)
w.info <- do.call(rbind, lapply(w.info, head, n = 1))
