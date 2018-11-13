
#Test of different heritability approaches using lichen data



#The methods compared are heritability estimated from aov, REML, adonis, adonis2, dbRDA, permanova in primer. 


Data<-read.csv("../data/HeritTest_Nov2018.csv")
attach(Data)
attributes(Data)

# [1] "Tree.ONC"  "Geno.ONC"  "T.Cov.ONC" "Rough.ONC" "CT.ONC" 

summary(Geno.ONC) #To get genotype sample sizes
#Output
#  10 1005 1007 1008 1012 1017 1023   11  996  999  EC1  H10 HE10  RL6  RM2  T15   T6  WC5 
#   4    4    3    8    4    8    4    5    9    3    3    3    2    3    3    3    3    4 










############################################################################

######Heritability and R2 Estimates for Total Lichen Cover#####



###########################
####one-way anova method####
summary(h2.aov <- aov(T.Cov.ONC~ Geno.ONC))

#output anova table
## > summary(aov(T.Cov.ONC~ Geno.ONC))
##             Df Sum Sq Mean Sq F value  Pr(>F)   
## Geno.ONC    17   2206  129.78    2.74 0.00224 **
## Residuals   58   2747   47.37  


#R2
2206 / (2206+2747)  # R2 = 0.4453866

#Heritability (from Shuster's calculator)
    #H2 = 0.294460929
H2(h2.aov, Geno.ONC)

###########################
####REML method####

summary(h2.reml <- lme4::lmer(T.Cov.ONC~ 1|Geno.ONC, REML=TRUE))

#Output
## Linear mixed model fit by REML ['lmerMod']
## Formula: T.Cov.ONC ~ 1 | Geno.ONC

## REML criterion at convergence: 523

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.5559 -0.6230 -0.2626  0.4320  2.6764 

## Random effects:
##  Groups   Name        Variance Std.Dev.
##  Geno.ONC (Intercept) 19.07    4.366   
##  Residual             47.60    6.899   
## Number of obs: 76, groups:  Geno.ONC, 18

## Fixed effects:
##             Estimate Std. Error t value
## (Intercept)    8.389      1.325   6.332


#Heritability
19.07/(19.07+47.60)   #H2 = 0.2860357

H2(h2.reml)

H2(h2.reml)

###########################
####dbRDA method####

anova(h2.dbrd <- vegan::dbrda(T.Cov.ONC~ Geno.ONC, dist="euclidean"))

#Output
## Model: vegan::dbrda(formula = T.Cov.ONC ~ Geno.ONC, distance = "euclidean")
##          Df Variance      F Pr(>F)   
## Model    17   29.417 2.7398  0.006 **
## Residual 58   36.632  


#R2 or H2(???) using provided varaince values. 
29.417/(29.417+36.632)  # R or H2 = 0.4453815

H2(h2.dbrd)

H2(h2.dbrd)










###########################
####PERMANOVA method using adonis####


vegan::adonis((T.Cov.ONC~ Geno.ONC), method="euclidean")

#Notice how similar the anova table here is to that of the aov function, above. 
#Output
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Geno.ONC  17    2206.3 129.781  2.7398 0.44538  0.004 **
Residuals 58    2747.4  47.369         0.55462          
Total     75    4953.7                 1.00000   




#R2 
2206.3/(2206.3 + 2747.4)  # R2 =  0.4453843. Same as the model output


#Heritability (from Shuster's calculator). This is the same as that using the aov function. The very slight difference is only due to there being more decimal places in the adonis anova table
    #H2 = 0.294457782

















###########################
####PERMANOVA method using adonis####


vegan::adonis2((T.Cov.ONC~ Geno.ONC), method="euclidean")

#Notice how similar the anova table here is to that of the aov function, above. 
#Output
         Df SumOfSqs      R2      F Pr(>F)   
Geno.ONC 17   2206.3 0.44538 2.7398  0.003 **
Residual 58   2747.4 0.55462                 
Total    75   4953.7 1.00000  



#R2 
2206.3/(2206.3 + 2747.4)  # R2 =  0.4453843. Same as the model output


#Heritability (from Shuster's calculator). Same this as regular adonis
    #H2 = 0.294457782











###########################
####PERMANOVA from Primer####

#Ofcourse I didn't run the permanova in R. But, I calculated teh R2 from the anova table,and the H2 from the variance components output from a one-way random effect permanova and H2 using shuster's calculator using the provided anova table. And, below is the anova table copied from Primer. 

          Df     SS        MSE       F        P  
Geno.ONC  17    2206.3    129.78    2.7398  0.004
Residuals 58    2747.4    47.36                       
Total     75    4953.7                 




#R2 calculated from anova table
2206.3/(2206.3 + 2747.4)  # R2 =  0.4453843


#Heritability from Shuster's calculator, using the anova table from primer
    #H2 = 0.294457782



#Heritability from variance components output from Primer. The genotype factor was set as a random effect and primer adjusts for unequal sample sizes. Note how the heritability estimate is equivalent to that from Shuster's calculator.  Also look how similar these variance estimates are to the REML outpout, and REML is the best approach for estimating varaince components with unbalanced designs.
((19.769)/(19.769+47.369))	  #H2 = 0.2944532
