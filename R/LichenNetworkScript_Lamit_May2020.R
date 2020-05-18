#Analyses for Lau's lichen network stuff.  May 15, 2020




###################################################################################
##########Packages
library (lmerTest) #For overall models
library(RLRsim) #For testing random effects
#library (emmeans) #For fixed effects and covariates

# https://www.ssc.wisc.edu/sscc/pubs/MM/MM_TestEffects.html.  #For GLM





###################################################################################
##########Data
load("data/onc_cn_dist.Rdata")
load("data/onc_lcn_data.Rdata")
ls()
attributes(cn.d.onc)
attributes(onc.dat)






###################################################################################
##########PC

PC_Mod1<-lmer(log(PC)~(1|geno), data=onc.dat,  REML=TRUE)

plot(PC_Mod1)

hist(onc.dat$PC)
hist(resid(PC_Mod1)) #Not great but probably fine.
qqnorm(resid(PC_Mod1))
qqline(resid(PC_Mod1))
fligner.test(onc.dat$PC ~ onc.dat$geno)

plot(onc.dat$PC~ onc.dat$geno)

exactRLRT(PC_Mod1) #Test of random effect






JunkFact<-as.factor(onc.dat$Junk)

###################################################################################
##########SR

SR_Mod1<-lmer(SR~(1|geno), data=onc.dat, REML=TRUE)
#SR_Mod1<-glmer(SR~(1|geno), data=onc.dat, family=poisson(link = "log")) #In the form of a GLM

summary(SR_Mod1)
plot(SR_Mod1)

hist(onc.dat$SR)
hist(resid(SR_Mod1)) #Not great but probably fine.
qqnorm(resid(SR_Mod1))
qqline(resid(SR_Mod1))
fligner.test(onc.dat$SR ~ onc.dat$geno)

exactRLRT(SR_Mod1) #Note, doesn't work on glmer, just on lmer









###################################################################################
##########SR

SD_Mod1<-lmer(SD~(1|geno), data=onc.dat, REML=TRUE)
#SR_Mod1<-glmer(SR~(1|geno), data=onc.dat, family=poisson(link = "log")) #In the form of a GLM. The variance associated with genotype is zero 

summary(SD_Mod1)
plot(SD_Mod1)

hist(onc.dat$SD)
hist(resid(SD_Mod1)) #Not great but probably fine.
qqnorm(resid(SD_Mod1))
qqline(resid(SD_Mod1))
fligner.test(onc.dat$SD ~ onc.dat$geno)

exactRLRT(SD_Mod1)












###################################################################################
##########SE

SE_Mod1<-lmer(SE~(1|geno), data=onc.dat, REML=TRUE)
#SE_Mod1<-glmer(SE~(1|geno), data=onc.dat, family=poisson(link = "log")) #In the form of a GLM. Warning message

summary(SE_Mod1)
plot(SE_Mod1)

hist(onc.dat$SE)
hist(resid(SE_Mod1)) 
qqnorm(resid(SE_Mod1))
qqline(resid(SE_Mod1))
fligner.test(onc.dat$SE ~ onc.dat$geno)

exactRLRT(SE_Mod1)












###################################################################################
##########BR

BR_Mod1<-lmer((BR)~(1|geno), data=onc.dat, REML=TRUE)
#BR_Mod1<-glmer(SE~(1|geno), data=onc.dat, family=Gamma(link = "inverse")) #But, not really needed

summary(BR_Mod1)
plot(BR_Mod1)

hist(onc.dat$BR)
hist(resid(BR_Mod1)) 
qqnorm(resid(BR_Mod1))
qqline(resid(BR_Mod1))
fligner.test(onc.dat$BR ~ onc.dat$geno)
plot(onc.dat$BR ~ onc.dat$geno)

exactRLRT(BR_Mod1) #Doesn't work with GLM!  p for lmer = RLRT = 4.3895, p-value = 0.0172. I used this function in the past, so...let's go with it.
lmerTest::ranova(BR_Mod1, )   #Another version of a random effect test. Doesn't work with GLM. P for lmer = 0.05963 ...so close!

102.4/(102.4+269.6) #Herit from regular lm = 0.2752688
0.03621/(0.09144+0.03621) #Herit = 0.2836663 #Heritability from GLM








###################################################################################
##########L.      

L_Mod1<-lmer(log(L+1)~(1|geno), data=onc.dat, REML=TRUE)
#L_Mod1<-glmer(SE~(1|geno), data=onc.dat, family=Gamma(link = "inverse")) 
#L_Mod1<-glmer(SE~(1|geno), data=onc.dat, family=poisson(link = "log")) #Gets errors but seems like poisson distribution makes the most sense

summary(L_Mod1)
plot(L_Mod1)

hist(onc.dat$L)
hist(resid(L_Mod1)) 
qqnorm(resid(L_Mod1))
qqline(resid(L_Mod1))
fligner.test(onc.dat$L ~ onc.dat$geno)
plot((onc.dat$L) ~ onc.dat$geno) 

plot(log(onc.dat$L+1) ~ onc.dat$geno)

exactRLRT(L_Mod1)

#Try the negative binomial model? Or just rank transform?  The distribution is all screwwy. Log+1 transformation maybe fine.









###################################################################################
##########L.      #Try the negative binomial model?

Cen_Mod1<-lmer(log(L+0.1)~(1|geno), data=onc.dat, REML=TRUE)
Cen_Mod1<-glmer(L~(1|geno), data=onc.dat, family=poisson(link = "log")) #In the form of a GLM

summary(Cen_Mod1)
plot(Cen_Mod1)

hist(onc.dat$Cen)
hist(resid(Cen_Mod1)) 
qqnorm(resid(Cen_Mod1))
qqline(resid(Cen_Mod1))
fligner.test(onc.dat$Cen ~ onc.dat$geno)
plot((onc.dat$Cen) ~ onc.dat$geno) 

exactRLRT(Cen_Mod1)


plot(log(onc.dat$Cen +1) ~ onc.dat$geno)
















###################################################################################
##########mod.lik   

mod.lik_Mod1<-lmer(log(mod.lik+1)~(1|geno), data=onc.dat, REML=TRUE)
#mod.lik_Mod1<-glmer(mod.lik ~(1|geno)+(1|as.factor(Junk)), data=onc.dat, family=poisson(link = "log")) #In the form of a GLM

summary(mod.lik_Mod1)
plot(mod.lik_Mod1)

hist(onc.dat$mod.lik)
hist(resid(mod.lik_Mod1)) 
qqnorm(resid(mod.lik_Mod1))
qqline(resid(mod.lik_Mod1))
fligner.test(onc.dat$mod.lik ~ onc.dat$geno)
plot((onc.dat$mod.lik) ~ onc.dat$geno) 

exactRLRT(mod.lik_Mod1)


#very Zero inflated. Hard to test.



















###################################################################################
##########C.     

C_Mod1<-lmer(sqrt(C)~(1|geno), data=onc.dat, REML=TRUE)

summary(C_Mod1)
plot(C_Mod1)

hist(onc.dat$C)
hist(resid(C_Mod1)) 
qqnorm(resid(C_Mod1))
qqline(resid(C_Mod1))
fligner.test(onc.dat$C ~ onc.dat$geno)
plot((onc.dat$C) ~ onc.dat$geno) 

exactRLRT(C_Mod1)













###################################################################################
##########CN     

CN_Mod1<-lmer(log(CN)~(1|geno), data=onc.dat, REML=TRUE)

summary(CN_Mod1)
plot(CN_Mod1)

hist(onc.dat$CN)
hist(resid(CN_Mod1)) 
qqnorm(resid(CN_Mod1))
qqline(resid(CN_Mod1))
fligner.test(onc.dat$CN ~ onc.dat$geno)
plot((onc.dat$CN) ~ onc.dat$geno) 

exactRLRT(CN_Mod1)









###################################################################################
##########N

N_Mod1<-lmer(sqrt(N)~(1|geno), data=onc.dat, REML=TRUE)

summary(N_Mod1)
plot(N_Mod1)

hist(onc.dat$N)
hist(resid(N_Mod1)) 
qqnorm(resid(N_Mod1))
qqline(resid(N_Mod1))
fligner.test(onc.dat$N ~ onc.dat$geno)
plot((onc.dat$N) ~ onc.dat$geno) 

exactRLRT(N_Mod1)

















###################################################################################
##########CT

CT_Mod1<-lmer(rank(CT)~(1|geno), data=onc.dat, REML=TRUE)  #Need to rank transform or something
#CT_Mod1<-glmer(CT~(1|geno), data=onc.dat, family=Gamma(link = "inverse")) #But, not really needed

summary(CT_Mod1)
plot(CT_Mod1)

hist(onc.dat$CT)
hist(resid(CT_Mod1)) 
qqnorm(resid(CT_Mod1))
qqline(resid(CT_Mod1))
fligner.test(onc.dat$CT ~ onc.dat$geno)
plot((onc.dat$CT) ~ onc.dat$geno) 

exactRLRT(CT_Mod1)

#The irritation outlier. But, reducing the dataset more isn't probably good for the lichen variables.










###################################################################################
##########pH
pH_Mod1<-lmer((pH)~(1|geno), data=onc.dat, REML=TRUE)

summary(pH_Mod1)
plot(pH_Mod1)

hist(onc.dat$pH)
hist(resid(pH_Mod1)) 
qqnorm(resid(pH_Mod1))
qqline(resid(pH_Mod1))
fligner.test(onc.dat$pH ~ onc.dat$geno)
plot((onc.dat$pH) ~ onc.dat$geno) 

exactRLRT(pH_Mod1)













###################################################################################
##########Composition      #Try the negative binomial model?

library(vegan)
## adonis((cn.d.onc^.25)~onc.dat$geno, family="bray")
adonis((cn.d.onc^2)~onc.dat$geno, perm = 100000)

txtplot::txtplot(cn.d.onc, cn.d.onc^2)
