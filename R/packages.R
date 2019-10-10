## Check for supporting packages
pkg.list <- c("drake", "magrittr", "devtools", "xtable", "reshape",
              "MuMIn", "lme4", "RLRsim", "vegan", "ecodist", 
              "bipartite", "RColorBrewer", "enaR"
              )
## install packages that are not installed
if (any(!(pkg.list %in% installed.packages()[, 1]))){
    sapply(pkg.list[which(!(pkg.list %in% installed.packages()[, 1]))], 
           install.packages, dependencies = TRUE, 
           repos = 'http://cran.us.r-project.org')
}
## Check non-CRAN packages
if (!("ComGenR" %in% installed.packages()[, 1])){
    devtools::install_github("ECGen/ComGenR")
}
if (!("conetto" %in% installed.packages()[, 1])){
    devtools::install_github("ECGen/conetto")
}
## Load libraries
sapply(c(pkg.list, "ComGenR", "conetto"), 
       library, quietly = TRUE, character.only = TRUE)
