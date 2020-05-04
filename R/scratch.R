## Results

## Genotype affects networks via bark roughness
### Genotype predicted bark roughness
### Genotype predicted lichen network metrics
### Genotype did not predict network metrics after controlling for bark roughness
### Genotype did no predict any other bark traits
### No other bark traits predicted lichen networks
### Bark roughness was negatively correlated with L, Cen and AMI

## Lichen community variation, which did not respond to genotype or
## bark roughness, correlated with lichen network structure
## Bark roughness not correlated with lichen PC, SR, SE, SD
### PC negatively correlated with Cen
### SR positively correlated with L
### SE positively correlated with L
### SD positively correlated with L, Cen


## Discussion

### Network structure genetically based 

### Genetically based trait, bark roughness, drove pattern, trees with
### rougher bark tended to have networks with fewer, more diffuse
### interactions

### Network structure also varied with other community dynamics that
### were not genetically based, driven by the dominance of one lichen
### species, X. galericulata

xg.reml <- lmer(onc.com[,"Xg"]^(1/4) ~ (1 | geno), 
                      data = onc.dat, REML = TRUE)
ls.reml <- lmer(onc.com[,"Ch"]^(1/4) ~ (1 | geno), 
                      data = onc.dat, REML = TRUE)
ls.reml <- lmer(onc.com[,"Ls"]^(1/4) ~ (1 | geno), 
                      data = onc.dat, REML = TRUE)

RLRsim::exactRLRT(ch.reml)
RLRsim::exactRLRT(xg.reml)
RLRsim::exactRLRT(ls.reml)
