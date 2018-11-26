### script for loading data for the
### ACN study
### Pit data collected 2012

source('../bin/helpers.R')
library(sna)
library(vegan)
source("loadPitdata.R")

### Genotype effect
vegan::adonis2(nd.acn.pit ~ geno * leaf.type, data = tree.info, sqrt.dist = TRUE)
vegan::adonis2(nd.acn.liv.pit ~ geno, data = tree.info[tree.info[, "leaf.type"] == "live", ], 
               sqrt.dist = TRUE)
vegan::adonis2(nd.acn.sen.pit ~ geno, data = tree.info[tree.info[, "leaf.type"] == "sen", ], 
               sqrt.dist = TRUE)


reml.dc <- lme4::lmer(I(all.dc^(1/1)) ~ (1 | geno) * leaf.type, 
                       data = tree.info,
                       REML = TRUE)
reml.pval.dc <- RLRsim::exactRLRT(reml.dc)
reml.pval.dc

reml.dc.liv <- lme4::lmer(I(liv.dc^(1/1)) ~ (1 | geno), 
                       data = tree.info[tree.info[, "leaf.type"] == "live", ],
                       REML = TRUE)
reml.pval.dc.liv <- RLRsim::exactRLRT(reml.dc.liv)
reml.pval.dc.liv

reml.dc.sen <- lme4::lmer(I(sen.dc^(1/1)) ~ (1 | geno), 
                       data = tree.info[tree.info[, "leaf.type"] == "sen", ],
                       REML = TRUE)
reml.pval.dc.sen <- RLRsim::exactRLRT(reml.dc.sen)
reml.pval.dc.sen
