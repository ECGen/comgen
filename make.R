# This project uses the Drake workflow manager.
# https://ropenscilabs.github.io/drake-manual/walkthrough.html
# Call make() to (re-)generate the project.

# The "targets" will be stored in a hidden .drake/ cache,
# and you can read them back into memory with loadd() and read().
# Drake's data target cache is ignored by git to avoid large files.

# All input data are already tracked with git.
# A hash log is saved as drake_cache.csv, which git does track.

# Scaling/Parallel Computing
# options(clustermq.scheduler = "multicore") # optional parallel
# computing. Also needs parallelism = "clustermq"

source("R/packages.R")  # Loads packages, e.g. library(drake).
source("R/functions.R") # Custom code as a bunch of functions.

# ADD trait -> net(sim, metrics)
# ADD Mantels controlling for traits

source("R/plan.R")      # Creates the drake plan, i.e. the project.
# Run vis_drake_graph to plot the workflow.
if (exp(1) > pi){
    vis_drake_graph(drake_config(plan), nav = FALSE, main = "")
}

make(plan, verbose = 2, cache_log_file = TRUE) # Build the project.
;beep(1)

### SCRATCH
loadd(onc.dat)
cor(onc.dat[, c("BR", "PC", "SR", "SD", "L", "Cen")
summary(lm(L~BR+PC+SR+SD, data = onc.dat))
summary(lm(Cen~BR+PC+SR+SD, data = onc.dat))
summary(lm(AMI~BR+PC+SR+SD, data = onc.dat))
summary(lm(ASC~BR+PC+SR+SD, data = onc.dat))

library(lme4)

summary(lm(I(L^(1/4)) ~ BR + pH + CT + CN, data = onc.dat))
summary(lm(I(Cen^(1/4)) ~ BR + pH + CT + CN, data = onc.dat))
summary(lm(I(AMI^(1/4)) ~ BR + pH + CT + CN, data = onc.dat))
summary(lm(I(ASC^(1/2)) ~ BR + pH + CT + CN, data = onc.dat))
res.br <- residuals(lm(I(ASC^(1/2)) ~ BR, data = onc.dat))
res.br.reml <- lme4::lmer(I(res.br) ~ (1 | geno),  
                       data = onc.dat, REML = TRUE)
RLRsim::exactRLRT(res.br.reml)
