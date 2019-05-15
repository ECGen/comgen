
## Equal sample sizes
## Uses Euclidean distance
g <- c(rep("A", 20), rep("B", 20))
y <- c(rnorm(20, 5), rnorm(20, 10))
d <- dist(y)
anova(aov(y ~ g)) # One-Way ANOVA 
adonis(d ~ g) # perMANOVA Anderson 2001
adonis2(d ~ g) # perMANOVA McArdle and Anderson 2001
anova(dbrda(d ~ g)) # db-RDA McArdle and Anderson 2001

# Uses Bray Curtis Distance
d <- vegan::vegdist(as.matrix(y))
r.aov <- aov(y ~ g)
r.perm <- adonis(d ~ g)
r.perm2 <- adonis2(d ~ g)
r.dbrd <- anova(dbrda(d ~ g))

## Un-equal sample size
g <- c(rep("A", 30), rep("B", 10), rep("C", 15))
y <- c(rnorm(30, 5), rnorm(10, 10), rnorm(15, 5))
d <- dist(y)
r.aov <- aov(y ~ g)
r.reml <- lme4::lmer(y ~ (1 | g), REML = TRUE)
r.perm <- adonis(d ~ g)
r.perm2 <- adonis2(d ~ g)
r.dbrd <- dbrda(d ~ g)

H2(r.aov, g)
H2(r.reml, g)
H2(r.perm2, g)
H2(r.dbrd, g)

## Un-equal sample size and using Bray Curtis
d <- vegan::vegdist(as.matrix(y))
anova(aov(y ~ g))
adonis(d ~ g)
adonis2(d ~ g)
anova(dbrda(d ~ g))
