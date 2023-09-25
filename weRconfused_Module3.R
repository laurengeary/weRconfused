#Georgia, Lauren, Elizabeth, & Weihang
#load libraries
library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

#Question 1: Establish anole.log data tibble
anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")
anole.log <- anole%>%
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%
  na.omit()%>%
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)

#Question 2:Construct two simple linear models that assess the effect of perch diameter and height 
lm.ArbPD <- lm(HTotal~SVL+ArbPD, anole.log)
lm.PH <- lm(HTotal~SVL+PH, anole.log)

#Question 3:Plot residuals of simple linear model to explore perch diameter and height effect hindlimb-SVL relationship
anole.log <- anole.log %>%
  mutate(res.ArbPD=residuals(lm.ArbPD), res.PH=residuals(lm.PH))
anole.log%>%
  dplyr::select(Ecomorph2,res.ArbPD,res.PH)%>%
  pivot_longer(cols=c("res.ArbPD","res.PH"))%>%
  ggplot(aes(x=Ecomorph2,y=value))+geom_boxplot()+facet_grid(name~.,scales = "free_y")+ylab("residual")

#Question 4:Construct phylogenetic least squares models of the hindlimb-SVL relationships
anole.tree <- read.tree("anole.tre")
#PGLS model with the hindlimb-SVL relationship + perch height
pgls.BM1 <- gls(HTotal~SVL+PH, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")
#PGLS model with the hindlimb-SVL relationship + perch diameter
pgls.BM2 <- gls(HTotal~SVL+ArbPD, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")
#PGLS model with the hindlimb-SVL relationship + perch height + perch diameter
pgls.BM3 <- gls(HTotal~SVL+ArbPD+PH, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#Question 5: Assess the fit of each of these three models using AICc and AICw 
anole.phylo.aic <- AICc(pgls.BM1,pgls.BM2,pgls.BM3)
aicw(anole.phylo.aic$AICc)
#comment whether one or both of covariates is a significant predictor of hindlimb length 
anova(pgls.BM3)
#Based on its p-value of 0.0005, perch diameter (ArbPD) is a better predictor of hind-limb length (HTotal)

#Question 6: Produce a plot that visualizes effect of your covatiate(s) and factors on hindlimb residuals of best fitting PGLS model
anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM3))
anole.log%>%
  ggplot(aes(Ecomorph,phylo.res, color="pink"))+geom_boxplot()