### Creates the plot of the species accumulation curves
library(vegan)

### Load the data files for ONC
garden.data <- read.csv("../data/lichen_networks/LCO_data_ONC_PIT.csv")
garden.data <- garden.data[garden.data$Geno != "RL6", ]
garden.data <- garden.data[garden.data$Tree != "N1.31", ]
garden.data[, 1] <- as.character(garden.data[, 1])
g1 <- substr(garden.data[, 1], 2, 2)
g1[g1 != "P"] <- "onc"
onc <- garden.data[g1 == "onc", ]
com <- split(onc[, 7:ncol(onc)], paste(onc[, 1], onc[, 3], onc[, 4]))
com <- do.call(rbind, lapply(com, function(x) apply(x, 2, sum)))
com.onc <- com

## Load the data for the wild stand
## The community matrix
x <- read.csv("../data/lichen_networks/lco_Apr2012.csv")
x <- x[, colnames(x) != "physcioid"]
x <- x[,colnames(x) != "dead"]
x <- x[, colnames(x) != "NOTES."]
x <- x[,!(apply(x,2,function(x) all(is.na(x))))]
x <- x[!(apply(x,1,function(x) as.character(x[1]) == "")),]
x <- x[x$tree != "gnu.44", ]
x <- x[x$tree != "ll.6", ]
lec.spp <- apply(x[, c(6, 8, 10, 18)], 1, function(x) sign(any(x != 0)))
x <- cbind(x, lec = lec.spp)
x <- x[, -c(6, 8, 10, 18)]
com <- split(x, paste(x$tree, x$quadrat, sep = "_"))
com <- lapply(com, function(x) apply(x[, -1:-4], 2, sum))
com <- do.call(rbind, com)
com.wild <- com

### The site vector
site <- c(rep("onc", nrow(com.onc)), rep("wild", nrow(com.wild)))

### The species accumulation curve plot
sac.onc <- specaccum(com.onc, method = "exact")
sac.wild <- specaccum(com.wild, method = "exact")

### Plotting
png("../results/spp_accum.png", width = 1000, height = 1000, pointsize = 30)
plot(sac.wild, xlim = c(0,nrow(com.onc)), xlab = "Number of Quadrats", ylab = "Number of Species", lwd = 3)
plot(sac.onc, add = TRUE, col = "darkgrey", lty = 2, lwd = 3)
dev.off()
