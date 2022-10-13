#Authors: Mads F. Schou
#Date: 13/10/2022

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

MCMCglmm(HeadAve ~ TempCon:TempDir*Pop + poly(daytime_z,2),
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

MCMCglmm(NeckAve ~ TempCon:TempDir*Pop + poly(daytime_z,2),
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

MCMCglmm(SurfaceTemp ~ BodyPart-1 + TempCon:TempDir*Pop*BodyPart + poly(daytime_z,2)*BodyPart,
         random = ~ us(1 + TempCon:TempDir):ID + Image + enclosure + date + year,
         data   = thermalL,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)



### Table S4: Results from model of the difference between head and neck surface temperature across cold, benign or hot temperatures (character-state).  

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

MCMCglmm(Diff ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*Pop,
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

### Table S5: Results from analyses of impact of temperature differences between head and neck on egg-laying rate (character-state).

prior.m4 <- list(
  R=list(V = diag(1), nu=0.002, fix = 1),
  G=list(G1=list(V        = diag(1),
                 nu        = 0.002),
         G2=list(V        = diag(1),
                 nu        = 0.002),
         G3=list(V        = diag(1),
                 nu        = 0.002),
         G4=list(V        = diag(1),
                 nu        = 0.002)))

thermal$DiffAve_z <- scale(thermal$DiffAve)[,1]

MCMCglmm(EggOrNot ~ MaxTempCat+ MaxTempCat:DiffAve_z-1 + Pop + Age2orNot,
         random = ~ ID + year + enclosure + date,
         data   = thermal,
         prior  = prior.m4,
         family = "threshold",
         burnin = Myburn, thin = Mythin,nitt = Mynitt,
         verbose = TRUE,
         pr = FALSE)


#############################################################
# 3. Evolutionary potential of surface temperatures
#############################################################

### Table S6: Results from animal model of head-neck differences in surface temperature (character-state). 

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

MCMCglmm(Diff ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*Pop,
         random = ~ us(TempCat):ID + us(TempCat):animal + year + enclosure + date,
         rcov = ~idh(TempCat):units,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S7: Results from animal model of head surface temperature (character-state). 
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

MCMCglmm(HeadAve ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*Pop,
         random = ~ us(TempCat):ID + us(TempCat):animal + year + enclosure + date,
         rcov = ~idh(TempCat):units,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S8: Results from animal model of neck surface temperature (character-state). 

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

MCMCglmm(NeckAve ~ 1 + Temp_zcat*TempCat + poly(daytime_z,2) + TempCat*Pop,
         random = ~ us(TempCat):ID + us(TempCat):animal + year + enclosure + date,
         rcov = ~idh(TempCat):units,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)


### Table S9: Results from animal model of head-neck differences in surface temperature (random regression). 

thermal <- thermal %>%
  mutate(Diff = HeadAve-NeckAve)

prior.m1 <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE), 
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

MCMCglmm(Diff ~ TempCon:TempDir*Pop + poly(daytime_z,2),
         random = ~ us(1 + TempCon:TempDir):ID + us(1+TempCon:TempDir):animal + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         pedigree = MyPed,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = T,
         pl = T)


#############################################################
# 4. Population differences in morphology
#############################################################

### Table S11:  Results from model of the difference between head and neck surface temperature with increasing/decreasing temperatures (random regression)

thermal <- thermal %>%
  mutate(Diff = HeadAve-NeckAve)

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

MCMCglmm(Diff ~ TempCon:TempDir*Pop + poly(daytime_z,2),
         random = ~ us(1 + TempCon:TempDir):ID + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE)


### Table S12: Results from model of the change in head surface temperature with increasing/decreasing temperatures (random regression) while accounting for body mass.

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

MCMCglmm(HeadAve ~ TempCon:TempDir*Pop*BodyMass_z + poly(daytime_z,2)*BodyMass_z,
         random = ~ us(1 + TempCon:TempDir):ID + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)

### Table S13: Results from model of the change in neck surface temperature with increasing/decreasing temperatures (random regression) while accounting for body mass.

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

MCMCglmm(NeckAve ~ TempCon:TempDir*Pop*BodyMass_z + poly(daytime_z,2)*BodyMass_z,
         random = ~ us(1 + TempCon:TempDir):ID + enclosure + date + year,
         data   = thermal,
         prior  = prior.m1,
         thin   = Mythin, burnin = Myburn, nitt   = Mynitt,
         verbose = TRUE,
         pr = FALSE)



#############################################################
# Estimating evolutionary parameters
#############################################################

#Examples from table S6

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

