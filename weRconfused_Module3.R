#going through project/adding notes and running the code/this is not report
#load libraries
library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

#load data
anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")

#merge anole data tibble with anole.eco tibble
anole2 <- anole%>%
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%
  na.omit()%>%
  print()

anole.log <- anole2%>%
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)

#Assessing hind limb length vs size with ggplot
anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_smooth(method="lm")

#linear relationship = Htotal=α⋅SVL+β
#a is slope and B is the y-intercept

#evaluating fit with simple linear model
anole.lm <- lm(HTotal~SVL,anole2)

coef(anole.lm)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_abline(slope=coef(anole.lm)[2],intercept=coef(anole.lm)[1],col="blue")

#Using linear model to predict HTotal...creating tibble with predictions
SVL2 <- seq(min(anole2$SVL),max(anole2$SVL),0.1)

pred.lm <-tibble(
  SVL=SVL2,
  H.pred=predict(anole.lm,newdata = data.frame(SVL=SVL2))
)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_point(data=pred.lm,aes(SVL,H.pred),col="blue")

summary(anole.lm)

#allometric model/exponenetial equation: Htotal=α⋅SVLβ
# a is the intercept and B is the scaling exponent
#fitting allometric model with nls
anole.allo <- nls(HTotal~a*SVL^b, start=list(b=1, a=1),data = anole2)

summary(anole.allo)

#calculate relative support of each model using AIC weights
#AICc from the MuMIn package
anole.aic <- AICc(anole.lm,anole.allo)

#aicw from the geiger package
anole.aicw <- aicw(anole.aic$AICc)

print(anole.aicw)

#alometric model has highest likelihood and because models have the same number 
#of parameters, is the model that fits best
logLik(anole.lm)
logLik(anole.allo)

#visualizing hindlimb-SVL relationships for each ecomorph in ggplot
anole.log%>%
  ggplot(aes(HTotal,SVL,col=Ecomorph2))+geom_point()+geom_smooth(method="lm")

#considering ecomorph in explaining hindlimb-SVL relationship
anole.log.eco.lm <- lm(HTotal~SVL*Ecomorph2,anole.log)
summary(anole.log.eco.lm)

#two way analysis of covariance...
#assessing effect of categorcal varibale (Ecomorph2) in context of how HTotal covaries with SVL
anova(anole.log.eco.lm)

#establish simple model with log-transformed data and compare it with AIC and AICw to this more complicated model
anole.log.lm  <- lm(HTotal~SVL,anole.log)
anova(anole.log.lm)

anole.log.aic <- AICc(anole.log.lm,anole.log.eco.lm)
aicw(anole.log.aic$AICc)

#compute the residuals based on our global anole.log.lm model
anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm))

#plot the residuals against Ecomorph2
anole.log%>%
  ggplot(aes(Ecomorph2,res))+geom_point()

#plot summary statistics including the median, hinges that correspond to the first and third quartiles, and whiskers that correspond to 1.5 times the interquartile range
p.eco <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=res)) +geom_boxplot()
print(p.eco)

#add a geometry that includes a representation of mean residual for each ecomorph
p.eco+ geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)

#phylogenetic comparative methods (PCMs)...Phylogenetic generalized least squares (PGLS)
anole.tree <- read.tree("anole.tre")
plot(anole.tree,cex=0.4)

#PGLS using simple regression models that don't include ecomorph and then models that do
#PGLS under BM, no ecomorph
pgls.BM1 <- gls(HTotal ~SVL, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under BM, w ecomorph
pgls.BM2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")


#PGLS under OU, no ecomorph
pgls.OU1 <- gls(HTotal ~SVL, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under OU, w, ecomorph
pgls.OU2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#perform AIC operations to see which models fit data best
anole.phylo.aic <- AICc(pgls.BM1,pgls.BM2,pgls.OU1,pgls.OU2)
aicw(anole.phylo.aic$AICc)

#consider if Ecomorph is significant in hindlimb-SVL relationship
anova(pgls.BM2)

#mutate and redefine our anole.log data to include a column for phylogenetically corrected residuals and then plot them against ecomorph just as we did in the case of our simple lm() model that didn’t consider phylogeny
anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM2))

p.eco.phylo <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=phylo.res)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)

print(p.eco.phylo)

#facets are used to break the plot up into a grid according values in the data. Let’s massage the data to a longer format and then use a faceting function from ggplot to plot the phlyognetically corrected and uncorrected residuals against ecomorph with a boxplot
anole.log%>%
  dplyr::select(Ecomorph2,res,phylo.res)%>%
  pivot_longer(cols=c("res","phylo.res"))%>%
  print%>%
  ggplot(aes(x=Ecomorph2,y=value)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)+facet_grid(name~.,scales = "free_y")+ylab("residual")

