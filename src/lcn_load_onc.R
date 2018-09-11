###LCN: loading the onc data
###MKLau
###21Mar2014

library(vegan)
library(bipartite)
library(ComGenR)
source('../bin/helpers.R')
source('../bin/cooc/src/cooc.R')

## source('~/projects/packages/ComGenR/R/CoCo.R')
cs <- function(x){nestedchecker(x)[[1]][1]}
mm <- function(x){slot(computeModules(x),'likelihood')}
garden.data <- read.csv('../data/lichen_networks/LCO_data_ONC_PIT.csv')
                                        #remove genotype RL6 and N1.31
garden.data <- garden.data[garden.data$Geno!='RL6',]
garden.data <- garden.data[garden.data$Tree!='N1.31',]
                                        #separate onc
garden.data[,1] <- as.character(garden.data[,1])
g1 <- substr(garden.data[,1],2,2)
g1[g1!='P'] <- 'onc'
onc <- garden.data[g1=='onc',]
					#tree overlap between years
unique(onc$Tree[onc$Year=='2010']) %in% unique(onc$Tree[onc$Year=='2011'])
unique(onc$Tree[onc$Year=='2011']) %in% unique(onc$Tree[onc$Year=='2010'])
                                        #
if (all(table(onc[,1])==100)){print('Good to go!')}else{for (i in 1:1000){print('Warning: check input data!!!')}}
                                        #separate trees
colnames(onc)[7:ncol(onc)] <- substr(colnames(onc)[7:ncol(onc)],1,2)
onc.q <- split(onc,paste(onc[,1],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])

                                        #get genotype
onc.geno <- unlist(sapply(names(onc.q),function(x) strsplit(x,split=' ')[[1]][2]))
                                        #Roughness in the Garden
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
                                        #match roughness to to ses values
load('../data/lichen_networks/lcn_onc_ses.rda')
onc.ses <- unlist(os[,1])
onc.ses[is.na(onc.ses)] <- 0
names(onc.ses) <- rownames(os)
if (all(names(onc.ses)==names(onc.q))){print('Good to go!')}else{print('Holy crap!')}
ses.tree <- as.character(sapply(names(onc.ses),function(x) unlist(strsplit(x,split=' '))[1]))
onc.rough <- avg.rough[match(ses.tree, r.tree)]
if (all(ses.tree==names(onc.rough))){print('Good to go!')}else{print('Holy Crap!')}

                                        #Microsat data from Nash
## gen.d <- read.csv(file='../../lcn/data/ONC_MSAT_datafromnash.csv')[,-1]
## gen.d[is.na(gen.d)] <- 0
## gen.d <- as.dist((gen.d))
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
if (all(rownames(rflp.d)==unique(onc.geno))){print('Good to go!')}else{print('Holy crap, rflp.d names match error')}
rflp.d <- as.dist(rflp.d)

                                        # Lichen Network Models
cn.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), coNets, ci.p = 95)
cn.sign.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), coNets, ci.p = 95, return.sign = TRUE)
cn.d.onc <- netDist(cn.onc, method = "bc")
fn.onc <- lapply(split(onc[, -1:-6], onc[, "Tree"]), freqNet, zero.diag = TRUE)
fn.d.onc <- netDist(fn.onc, method = "bc")

onc.tree <- do.call(rbind, strsplit(names(onc.geno), " "))[, 1]
cn.mu.onc <- list()
for (i in 1:length(unique(onc.geno))){cn.mu.onc[[i]] <- netMean(cn.onc[onc.geno == unique(onc.geno)[i]])}
cn.mu.d.onc <- netDist(cn.mu.onc, method = "bc")
names(cn.mu.onc) <- unique(onc.geno)
if (all(onc.tree == names(cn.onc))){print("Good to go!")}else{
    print("Danger Will Robinson!")
}

                                        #community data
onc.com <- do.call(rbind,lapply(onc.q,function(x) apply(x,2,sum)))
onc.R <- apply(sign(onc.com),1,sum)
onc.H <- vegan::diversity(onc.com)
onc.com.rel <- apply(onc.com,2,function(x) x/max(x))
onc.com.rel <- cbind(onc.com.rel,ds=rep(min(onc.com.rel[onc.com.rel!=0]),nrow(onc.com.rel)))
onc.com <- cbind(onc.com,ds=rep(min(onc.com[onc.com!=0]),nrow(onc.com)))
                                        #nmds and procrustes rotation
## onc.nms <- nmds.min(nmds(vegdist(onc.com.rel),2,2))
## onc.rot <- procrustes(dist(onc.rough),onc.nms,scale=TRUE)
## onc.rot <- t(onc.rot$rotation)
## onc.rot.sem <- onc.rot
## onc.rot <- onc.rot[,(1:2)[abs(cor(cbind(onc.rough,onc.rot))[1,2:3])==max(abs(cor(cbind(onc.rough,onc.rot))[1,2:3]))]]

                                        #genotype means
omu <- apply(onc.com[,colnames(onc.com)!='ds'], 2, 
             function(x,g) tapply(x,g,mean),g=onc.geno)
oms <- tapply(onc.ses, onc.geno, mean)
oms.d <- dist(oms[match(rownames(as.matrix(rflp.d)),names(oms))])
                                        #bark roughness means
oprbmu <- tapply(onc.rough,onc.geno,mean)
oprbmu <- oprbmu[match(rownames(as.matrix(rflp.d)),names(oprbmu))]
                                        #Co-occurrence counts
                                        #co-occurrence patterns
oco <- do.call(rbind,lapply(onc.q,function(x,t) apply(CoCo(x,type=t),2,sum),t='pos'))
och <- do.call(rbind,lapply(onc.q,function(x,t) apply(CoCo(x,type=t),2,sum),t='neg'))
oco <- oco[order(apply(oco, 1, sum), decreasing = TRUE), order(apply(oco, 2, sum), decreasing = TRUE)]
och <- och[order(apply(och, 1, sum), decreasing = TRUE), order(apply(och, 2, sum), decreasing = TRUE)]
                                        #get araujo coordinates
coord <- read.csv('../data/lichen_networks/lcn_coord_onc.csv')
rownames(coord) <- coord[,1]
coord <- coord[,-1]
### Renaming
oq <- onc.q
oc <- onc.com[,colnames(onc.com)!='ds']
os <- onc.ses
og <- onc.geno
osgmu <- tapply(os ,og,mean)
osgse <- tapply(os ,og,function(x)sd(x)/sqrt(length(x)))
prb.onc <- onc.rough
prb.mu.onc <- tapply(prb.onc, og, mean)
prb.mu.d.onc <- dist(prb.mu.onc)
