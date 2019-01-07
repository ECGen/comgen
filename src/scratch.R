apply(cbind(l.dat[l.dat[, "Moth"] == 0, "Litter.."], l.dat[l.dat[, "Moth"] == 1, "Litter.."], 
      tapply(l.dat[, "Litter.."] , l.dat[, "Tree.pairs"], diff)), 2, mean)
cbind(l.dat[l.dat[, "Moth"] == 0, "Light...average"], l.dat[l.dat[, "Moth"] == 1, "Light...average"], 
      tapply(l.dat[, "Light...average"] , l.dat[, "Tree.pairs"], diff))
cbind(abun[l.dat[, "Moth"] == 0], abun[l.dat[, "Moth"] == 1], 
      tapply(abun , l.dat[, "Tree.pairs"], diff))


