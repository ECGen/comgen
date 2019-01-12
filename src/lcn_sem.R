### Structural Equation Modeling
### For the following LCN model

## G -> PRB -> PTC -> COM-NMS1
##              \     ^
##               V   /     
##               SPR -> NET-NMS2 

## G needs to be made into a dummy matrix

geno.dm <- model.matrix(~ geno, data = onc.dat)
lcn.sem.dat <- cbind(geno.dm, onc.dat, 
                     com.nms1 = ord.com[, 1],
                     cn.nms1 = ord.cn[, 1], cn.nms2 = ord.cn[, 2])
lcn.sem.dat <- lcn.sem.dat[, !(colnames(lcn.sem.dat) %in% 
                                 c("tree", "geno", "L", "Cen", "(Intercept)"))]
lcn.sem.dat <- apply(lcn.sem.dat, 2, as.numeric)
lcn.cov <- cov(lcn.sem.dat)

### Using SEM
lcn.sem.results <- sem(specifyModel("./lcn_sem_model.txt"), S = lcn.cov, N = 500)

lcn.sem.results <- sem(specifyModel("./lcn_sem_model_composite.txt"), S = lcn.cov, N = 500)

summary(lcn.sem.results)
semPaths(lcn.sem.results, 
         "equality", "estimates", 
         style = "lisrel", 
         nCharNodes = 0, edge.label.cex = 0.35)
modIndices(lcn.sem.results)

### Using Lavann

## https://www.r-bloggers.com/ecological-sems-and-composite-variables-what-why-and-how/
## https://github.com/jebyrnes
## https://github.com/kelpecosystems/observational_data
