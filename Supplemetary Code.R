#Authors: Mads F. Schou
#Date: 30/01/2022

pacman::p_load(MCMCglmm,tidyverse,parallel)

# All MCMCglmm models were run 3 times like this
# m1_3 <- mclapply(1:3, function(i) {
#   MCMCglmm()
# }, mc.cores=3)

# see definations of variables in supplementary materials.

Myburn = 100000
Mythin = 4000
Mynitt = 5100000 + Myburn
Nsamples = (Mynitt-Myburn)/Mythin

#For reproductive success:
# Myburn = 1500000
# Mythin = 10000
# Mynitt = 31500000 + Myburn
# Nsamples = (Mynitt-Myburn)/Mythin


#############################################################
# 1. Thermal plasticity models
#############################################################

### Table S1: Results from model of the change in head surface temperature with increasing/decreasing temperatures (random regression). 

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(1),   
                 nu        = 0.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002)))

MCMCglmm(HeadAve ~ TempCon:TempDir*sex + TempCon:TempDir*Pop + poly(daytime_z,2),
         random = ~ us(1 + TempCon:TempDir):ID + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = T)

### Table S2: Results from model of the change in neck surface temperature with increasing/decreasing temperatures (random regression). 

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(1),   
                 nu        = 0.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002)))

MCMCglmm(NeckAve ~ TempCon:TempDir*sex + TempCon:TempDir*Pop + poly(daytime_z,2),
         random = ~ us(1 + TempCon:TempDir):ID + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = T)

### Table S3: Results from model of the change in surface temperature (neck and head) with increasing/decreasing temperatures (random regression).  

thermalL <- thermal %>%
  gather(key = BodyPart, "SurfaceTemp", NeckAve:HeadAve, factor_key=TRUE)

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(1),   
                 nu        = 0.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002),
         G5=list(V        = diag(1),   
                 nu        = 0.002)))

MCMCglmm(SurfaceTemp ~ BodyPart-1 + TempCon:TempDir*sex*BodyPart + TempCon:TempDir*Pop*BodyPart + poly(daytime_z,2)*BodyPart,
         random = ~ us(1 + TempCon:TempDir):ID + Image + enclosure + date + year,
         data   = thermalL,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S4: Results from model of the change in head surface temperature with increasing/decreasing temperatures (random regression) while accounting for body mass.

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(1),   
                 nu        = 0.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002)))

MCMCglmm(HeadAve ~ TempCon:TempDir*sex*BodyMass_z + TempCon:TempDir*Pop*BodyMass_z + poly(daytime_z,2)*BodyMass_z,
         random = ~ us(1 + TempCon:TempDir):ID + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S5: Results from model of the change in neck surface temperature with increasing/decreasing temperatures (random regression) while accounting for body mass.

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(1),   
                 nu        = 0.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002)))

MCMCglmm(NeckAve ~ TempCon:TempDir*sex*BodyMass_z + TempCon:TempDir*Pop*BodyMass_z + poly(daytime_z,2)*BodyMass_z,
         random = ~ us(1 + TempCon:TempDir):ID + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S6: Results from model of the difference between neck and head surface temperature across cold, benign or hot temperatures (character-state).  

prior.m1 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 n        = 2.002),
         G2=list(V        = diag(1),   
                 n        = 0.002),
         G3=list(V        = diag(1),   
                 n        = 0.002),
         G4=list(V        = diag(1),   
                 n        = 0.002)))

MCMCglmm(Diff ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*sex + TempCat*Pop,
         random = ~ us(TempCat):ID + year + enclosure + date,
         rcov = ~idh(TempCat):units,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)


#############################################################
# 2. Thermal plasticity and reproductive success
#############################################################

### Table S7: Results from analyses of impact of temperature differences between head and neck on egg-laying rate (character-state).

prior.m4 <- list(
  R=list(V = diag(1), nu=0.002, fix = 1), #No fix as I do this in the rcov term
  G=list(G1=list(V        = diag(1),
                 nu        = 0.002),
         G2=list(V        = diag(1),
                 nu        = 0.002),
         G3=list(V        = diag(1),
                 nu        = 0.002),
         G4=list(V        = diag(1),
                 nu        = 0.002)))

MCMCglmm(EggOrNot ~ MaxTempCat+ MaxTempCat:DiffAve_z-1 + Pop + Age2orNot,
         random = ~ ID + year + enclosure + date,
         data   = thermal,
         prior  = prior.m4,
         family = "threshold",
         burnin = Myburn, thin = Mythin,nitt = Mynitt,
         verbose = TRUE,
         pr = FALSE)

#############################################################
# 3. Incubation vs Standing
#############################################################

#----------- Analysis head 

t.test(sitting$Head_sit[sitting$TimeGroup %in% "Afternoon"], sitting$Head_befaf[sitting$TimeGroup %in% "Afternoon"], paired = TRUE, alternative = "two.sided")

t.test(sitting$Head_sit[sitting$TimeGroup %in% "Morning"], sitting$Head_befaf[sitting$TimeGroup %in% "Morning"], paired = TRUE, alternative = "two.sided")

#----------- Analysis neck 

t.test(sitting$Neck_sit[sitting$TimeGroup %in% "Afternoon"], sitting$Neck_befaf[sitting$TimeGroup %in% "Afternoon"], paired = TRUE, alternative = "two.sided")

t.test(sitting$Neck_sit[sitting$TimeGroup %in% "Morning"], sitting$Neck_befaf[sitting$TimeGroup %in% "Morning"], paired = TRUE, alternative = "two.sided")

#----------- Analysis Diff 

t.test(sitting$NeckTempDiff[sitting$TimeGroup %in% "Afternoon"], sitting$HeadTempDiff[sitting$TimeGroup %in% "Afternoon"], paired = TRUE, alternative = "two.sided")

t.test(sitting$NeckTempDiff[sitting$TimeGroup %in% "Morning"], sitting$HeadTempDiff[sitting$TimeGroup %in% "Morning"], paired = TRUE, alternative = "two.sided")


#############################################################
# 4. Evolutionary potential of surface temperatures
#############################################################

### Table S8: Results from animal model of neck-head differences in surface temperature (character-state). 

prior.m1 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 n        = 2.002),
         G2=list(V        = diag(3),
                 n        = 2.002),
         G3=list(V        = diag(1),   
                 n        = 0.002),
         G4=list(V        = diag(1),   
                 n        = 0.002),
         G5=list(V        = diag(1),   
                 n        = 0.002)))

MCMCglmm(Diff ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*sex + TempCat*Pop,
         random = ~ us(TempCat):ID + us(TempCat):animal + year + enclosure + date,
         rcov = ~idh(TempCat):units,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S9: Results from animal model of head surface temperature (character-state). 
prior.m1 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 n        = 2.002),
         G2=list(V        = diag(3),
                 n        = 2.002),
         G3=list(V        = diag(1),   
                 n        = 0.002),
         G4=list(V        = diag(1),   
                 n        = 0.002),
         G5=list(V        = diag(1),   
                 n        = 0.002)))

MCMCglmm(HeadAve ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*sex + TempCat*Pop,
         random = ~ us(TempCat):ID + us(TempCat):animal + year + enclosure + date,
         rcov = ~idh(TempCat):units,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S10: Results from animal model of neck surface temperature (character-state). 

prior.m1 <- list(
  R=list(V = diag(3), nu=2.002, fix = FALSE), 
  G=list(G1=list(V        = diag(3),
                 n        = 2.002),
         G2=list(V        = diag(3),
                 n        = 2.002),
         G3=list(V        = diag(1),   
                 n        = 0.002),
         G4=list(V        = diag(1),   
                 n        = 0.002),
         G5=list(V        = diag(1),   
                 n        = 0.002)))

MCMCglmm(NeckAve ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*sex + TempCat*Pop,
         random = ~ us(TempCat):ID + us(TempCat):animal + year + enclosure + date,
         rcov = ~idh(TempCat):units,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)


### Table S11: Results from animal model of neck-head differences in surface temperature (random regression). 

thermal <- thermal %>%
  mutate(Diff = NeckAve-HeadAve)

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),  #Inverse-Wishart prior with low beleif. 2.002 because 3 in diag!
  G=list(G1=list(V        = diag(3),
                 nu        = 2.002),
         G2=list(V        = diag(3),
                 nu        = 2.002),
         G3=list(V        = diag(1),   
                 nu        = 0.002),
         G4=list(V        = diag(1),   
                 nu        = 0.002),
         G5=list(V        = diag(1),   
                 nu        = 0.002)))

MCMCglmm(Diff ~ TempCon:TempDir*sex + TempCon:TempDir*Pop + poly(daytime_z,2),
         random = ~ us(1 + TempCon:TempDir):ID + us(1+TempCon:TempDir):animal + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = T,
         pl = T)

#############################################################
# 5. Population differences in morphology and origin
#############################################################

### Table S12: Results from analyses of population differences in neck length.

MCMCglmm(necklenbody ~ Pop-1 + sex +height_z:Pop + height_z:sex , data = morphan,nitt=MyItt, thin=MyThin, burnin=MyBurn)

### Table S13: Results from analyses of population differences in neck length to height ratio.

MCMCglmm(neckbody.height ~ Pop-1 + Pop:sex, data = morphan,nitt=MyItt, thin=MyThin, burnin=MyBurn)


#############################################################
# Estimating evolutionary parameters
#############################################################

#Examples from table S10

#----------  evolvability  ----------#

evo_High<-100*m1$VCV[,'TempCatHot:TempCatHot.animal']/
  (posterior.mode(m1$Sol[,"(Intercept)"] + m1$Sol[,"TempCatHot"])^2)
evo_High_est <- data.frame(t(c(posterior.mode(evo_High),HPDinterval(evo_High))))

#----------  repeatability ----------#

rep_high<-(m1$VCV[,'TempCatHot:TempCatHot.animal']+m1$VCV[,'TempCatHot:TempCatHot.ID'])/
  (
    m1$VCV[,'TempCatHot:TempCatHot.animal']+
      m1$VCV[,'TempCatHot:TempCatHot.ID']+
      m1$VCV[,'TempCatHot.units']+
      m1$VCV[,'date']+
      m1$VCV[,'year']+
      m1$VCV[,'enclosure'])
rep_high_est <- data.frame(t(c(posterior.mode(rep_high),HPDinterval(rep_high))))

#----------  heritability ----------#

h2_high<-m1$VCV[,'TempCatHot:TempCatHot.animal']/
  (
    m1$VCV[,'TempCatHot:TempCatHot.animal']+
      m1$VCV[,'TempCatHot:TempCatHot.ID']+
      m1$VCV[,'TempCatHot.units']+
      m1$VCV[,'enclosure']+
      m1$VCV[,'year']+
      m1$VCV[,'date'])

h2_high_est <- data.frame(t(c(posterior.mode(h2_high),HPDinterval(h2_high))))

