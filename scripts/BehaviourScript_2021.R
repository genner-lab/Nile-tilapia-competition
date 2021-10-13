#Install packages

library(ggplot2)
library(ggpubr)
library(MASS)
library(RColorBrewer)
library(dplyr)  
library(ggeffects)
library(effects)
library(lmerTest)
library(AICcmodavg)
library(MuMIn)
library(performance)
library(jtools)
library(sjPlot)

#Split MC and NT behaviour, into two NT densities

BoxPlotData <- read.table("assets/Behaviour_Data.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
BoxPlotData_Fives <- filter(BoxPlotData, NTDensity == "Five")
BoxPlotData_Tens <- filter(BoxPlotData, NTDensity == "Ten")

#Generate the summary figures

brewer.pal(8,"Set1")
RedPalette <- c("orange", "#3787b8", "#red", "#1e1c57")

Part1 <- ggplot(BoxPlotData_Fives, aes(x=Species, y=Movement, fill = Species)) + scale_fill_manual(values=RedPalette) + 
  theme_classic() + geom_boxplot() + ylim(0, 150) + labs(x =" ", y = "Movement (# movement events)")
Part1

Part2 <- ggplot(BoxPlotData_Tens, aes(x=Species, y=Movement, fill = Species)) + scale_fill_manual(values=RedPalette) + 
  theme_classic() + geom_boxplot() + ylim(0, 150) + labs(x =" ", y = " ")
Part2

Part3 <- ggplot(BoxPlotData_Fives, aes(x=Species, y=Chase, fill = Species)) + scale_fill_manual(values=RedPalette) + 
  theme_classic() + geom_boxplot() + ylim(0, 6) + labs(x =" ", y = "Heterospecific chases (# events)")
Part3

Part4 <- ggplot(BoxPlotData_Tens, aes(x=Species, y=Chase, fill = Species)) + scale_fill_manual(values=RedPalette) + 
  theme_classic() + geom_boxplot() + ylim(0, 6) + labs(x =" ", y = " ")
Part4

Part5 <- ggplot(BoxPlotData_Fives, aes(x=Species, y=Shelter, fill = Species)) + scale_fill_manual(values=RedPalette) + 
  theme_classic() + geom_boxplot() + ylim(0, 10) + labs(x ="Species", y = "Shelter use (mean #individuals in shelter)")
Part5

Part6 <- ggplot(BoxPlotData_Tens, aes(x=Species, y=Shelter, fill = Species)) + scale_fill_manual(values=RedPalette) + 
  theme_classic() + geom_boxplot() + ylim(0, 10) + labs(x ="Species", y = " ")
Part6

#Save at 7x9
Figure3 <- ggarrange(Part1, Part2, Part3, Part4, Part5, Part6, labels = c("a", " ", "b"," ", "c", " "), ncol = 2, nrow = 3, common.legend =TRUE, legend = "none",align = "hv",  hjust = -2.5, label.x = 0.12)
Figure3

#Spearman's tests

MC_Fives <- filter(BoxPlotData_Fives, Species == "MC")
cor.test( ~ Movement + Shelter, data=MC_Fives, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Movement + Chase, data=MC_Fives, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Chase + Shelter, data=MC_Fives, method = "spearman", conf.level = 0.95, exact=FALSE)

MC_Tens <- filter(BoxPlotData_Tens, Species == "MC")
cor.test( ~ Movement + Shelter, data=MC_Tens, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Movement + Chase, data=MC_Tens, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Chase + Shelter, data=MC_Tens, method = "spearman", conf.level = 0.95, exact=FALSE)

NT_Fives <- filter(BoxPlotData_Fives, Species == "NT")
cor.test( ~ Movement + Shelter, data=NT_Fives, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Movement + Chase, data=NT_Fives, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Chase + Shelter, data=NT_Fives, method = "spearman", conf.level = 0.95, exact=FALSE)

NT_Tens <- filter(BoxPlotData_Tens, Species == "NT")
cor.test( ~ Movement + Shelter, data=NT_Tens, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Movement + Chase, data=NT_Tens, method = "spearman", conf.level = 0.95, exact=FALSE)
cor.test( ~ Chase + Shelter, data=NT_Tens, method = "spearman", conf.level = 0.95, exact=FALSE)

#Quantifying correlations between environmental variables across all trials

Expt_All<- read.table("assets/Full_Data.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

cor(Expt_All$DO,Expt_All$pH, use="complete.obs", method="spearman")
cor(Expt_All$DO,Expt_All$WT, use="complete.obs", method="spearman")
cor(Expt_All$DO,Expt_All$PC1, use="complete.obs", method="spearman")
cor(Expt_All$pH,Expt_All$WT, use="complete.obs", method="spearman")
cor(Expt_All$pH,Expt_All$PC1, use="complete.obs", method="spearman")
cor(Expt_All$WT,Expt_All$PC1, use="complete.obs", method="spearman")

#MC_Movement models

ModA1 <- glmer(MC_Movement ~ NT_Density + sqrt(NT_Movement) + NT_Chase + PC1 + Delta_Size +(1|POOL), data=Expt_All, family=poisson(link="log"), na.action = "na.fail")
ModA2 <- glmer(MC_Movement ~ sqrt(NT_Movement) + NT_Chase + PC1 + Delta_Size +(1|POOL), data=Expt_All, family=poisson(link="log"))
ModA3 <- glmer(MC_Movement ~ NT_Density + sqrt(NT_Movement) + NT_Chase + Delta_Size +(1|POOL), data=Expt_All, family=poisson(link="log"))
ModA4 <- glmer(MC_Movement ~ NT_Density + sqrt(NT_Movement) + NT_Chase + PC1 + (1|POOL), data=Expt_All, family=poisson(link="log"))
ModA5 <- glmer(MC_Movement ~ NT_Density + PC1 + Delta_Size + NT_Chase + (1|POOL), data=Expt_All, family=poisson(link="log"))
ModA6 <- glmer(MC_Movement ~ NT_Density + PC1 + Delta_Size + sqrt(NT_Movement) +(1|POOL), data=Expt_All, family=poisson(link="log"))

#Compare AIC of models
models <- list(ModA1, ModA2, ModA3, ModA4, ModA5, ModA6)
model.names <- c('ModA1', 'ModA2', 'ModA3', 'ModA4', 'ModA5', 'ModA6')
aictab(cand.set = models, modnames = model.names)

#Key detail on best MC_Movement model
summary(ModA2)
plot(ModA2)
qqnorm(residuals(ModA2))
r.squaredGLMM(ModA2)

#MC_Shelter models

ModB1 <- lmer(MC_Shelter ~ NT_Density + sqrt(NT_Movement) + NT_Chase + PC1 + Delta_Size +(1|POOL), data=Expt_All, REML=F, na.action = "na.fail")
ModB2 <- lmer(MC_Shelter ~ sqrt(NT_Movement) + NT_Chase + PC1 + Delta_Size +(1|POOL), data=Expt_All, REML=F)
ModB3 <- lmer(MC_Shelter ~ NT_Density + sqrt(NT_Movement) + NT_Chase + Delta_Size +(1|POOL), data=Expt_All, REML=F)
ModB4 <- lmer(MC_Shelter ~ NT_Density + sqrt(NT_Movement) + NT_Chase + PC1 + (1|POOL), data=Expt_All, REML=F)
ModB5 <- lmer(MC_Shelter ~ NT_Density + PC1 + Delta_Size + NT_Chase +(1|POOL), data=Expt_All, REML=F)
ModB6 <- lmer(MC_Shelter ~ NT_Density + PC1 + sqrt(NT_Movement) + Delta_Size +(1|POOL), data=Expt_All, REML=F)

#Compare AIC of models
models <- list(ModB1, ModB2, ModB3, ModB4, ModB5, ModB6)
model.names <- c('ModB1', 'ModB2', 'ModB3', 'ModB4', 'ModB5', 'ModB6')
aictab(cand.set = models, modnames = model.names)

#Key detail on best MC_Shelter model
summary(ModB4)
plot(ModB4)
qqnorm(residuals(ModB4))
r.squaredGLMM(ModB4)

#NT_Movement Models

ModC1 <- glmer(NT_Movement ~ NT_Density + sqrt(MC_Movement) + MC_Chase + PC1 + Delta_Size +(1|POOL), data=Expt_All, family=poisson(link="log"),na.action = "na.fail")
ModC2 <- glmer(NT_Movement ~ sqrt(MC_Movement) + MC_Chase + PC1 + Delta_Size +(1|POOL), data=Expt_All,  family=poisson(link="log"))
ModC3 <- glmer(NT_Movement ~ NT_Density + sqrt(MC_Movement) + MC_Chase + Delta_Size +(1|POOL), data=Expt_All, family=poisson(link="log"))
ModC4 <- glmer(NT_Movement ~ NT_Density + sqrt(MC_Movement) + MC_Chase + PC1 +(1|POOL), data=Expt_All, family=poisson(link="log"))
ModC5 <- glmer(NT_Movement ~ NT_Density + PC1 + MC_Chase + Delta_Size +(1|POOL), data=Expt_All, family=poisson(link="log"))
ModC6 <- glmer(NT_Movement ~ NT_Density + sqrt(MC_Movement) + PC1 + Delta_Size +(1|POOL), data=Expt_All, family=poisson(link="log"))

#Compare AIC of models
models <- list(ModC1, ModC2, ModC3, ModC4, ModC5, ModC6)
model.names <- c('ModC1', 'ModC2', 'ModC3', 'ModC4', 'ModC5', 'ModC6')
aictab(cand.set = models, modnames = model.names)

#Key detail on best NT_Movement model
summary(ModC4)
plot(ModC4)
qqnorm(residuals(ModC4))
r.squaredGLMM(ModC4)

#NT_Shelter models

ModD1 <- lmer(NT_Shelter ~ NT_Density + sqrt(MC_Movement) + MC_Chase +  PC1 + Delta_Size +(1|POOL), data=Expt_All,REML=F,na.action = "na.fail")
ModD2 <- lmer(NT_Shelter ~ sqrt(MC_Movement) + MC_Chase +  PC1 + Delta_Size +(1|POOL), data=Expt_All,REML=F)
ModD3 <- lmer(NT_Shelter ~ NT_Density + sqrt(MC_Movement) + MC_Chase +  Delta_Size +(1|POOL), data=Expt_All,REML=F)
ModD4 <- lmer(NT_Shelter ~ NT_Density + sqrt(MC_Movement) + MC_Chase +  PC1 +(1|POOL), data=Expt_All,REML=F)
ModD5 <- lmer(NT_Shelter ~ NT_Density + PC1 + MC_Chase + Delta_Size +(1|POOL), data=Expt_All,REML=F)
ModD6 <- lmer(NT_Shelter ~ NT_Density + PC1 + sqrt(MC_Movement) + Delta_Size +(1|POOL), data=Expt_All,REML=F)

#Compare AIC of models

models <- list(ModD1, ModD2, ModD3, ModD4, ModD5, ModD6)
model.names <- c('ModD1', 'ModD2', 'ModD3', 'ModD4', 'ModD5', 'ModD6')
aictab(cand.set = models, modnames = model.names)

#Key detail on best NT_Shelter model
summary(ModD6)
plot(ModD6)
qqnorm(residuals(ModD6))
r.squaredGLMM(ModD6)

###Plotting Figure 4 (Export as 8x8)

PlotModA2_Mo <- effect_plot(ModA2, pred = NT_Movement, interval = TRUE, plot.points = TRUE, colors = "orange")
PlotModA2_Mo <- PlotModA2_Mo + theme_classic() + labs(x ="NT movement (# events)", y = "MC movement (# events)")
PlotModA2_Mo <- PlotModA2_Mo + xlim(0,150) + ylim(0,60)
PlotModA2_Mo

PlotModA2_Ch <- effect_plot(ModA2, pred = NT_Chase, interval = TRUE, plot.points = TRUE, colors = "orange")
PlotModA2_Ch <- PlotModA2_Ch + theme_classic() + labs(x ="NT chases initiated (# events)", y = "MC movement (# events)")
PlotModA2_Ch <- PlotModA2_Ch + xlim(0,6) + ylim(0,60)
PlotModA2_Ch

PlotModC4_Mo <- effect_plot(ModC4, pred = MC_Movement, interval = TRUE, plot.points = TRUE, colors = "darkblue")
PlotModC4_Mo <- PlotModC4_Mo + theme_classic() + labs(x ="MC movement (# events)", y = "NT movement (# events)")
PlotModC4_Mo <- PlotModC4_Mo + xlim(0,60) + ylim(0,150)
PlotModC4_Mo

PlotModC4_Ch <- effect_plot(ModC4, pred = MC_Chase, interval = TRUE, plot.points = TRUE, colors = "darkblue")
PlotModC4_Ch <- PlotModC4_Ch + theme_classic() + labs(x ="MC chases initiated (# events)", y = "NT movement (# events)")
PlotModC4_Ch <- PlotModC4_Ch + xlim(0,2) + ylim(0,150)
PlotModC4_Ch

Figure4 <- ggarrange(PlotModA2_Mo, PlotModA2_Ch, PlotModC4_Mo, PlotModC4_Ch, labels = c("a", "b", "c","d"), ncol = 2, nrow = 2, common.legend =TRUE, legend = "none",align = "hv",  hjust = -2.5, label.x = 0.12)
Figure4

###Plotting Figure 5 (Export as 8x8)

#MCActivity
PlotModA2_PC <- effect_plot(ModA2, pred =  PC1, interval = TRUE, plot.points = TRUE, colors = "orange")
PlotModA2_PC <- PlotModA2_PC + theme_classic() + labs(x ="PC1", y = "MC movement (# events)")
PlotModA2_PC <- PlotModA2_PC + ylim(0,60)
PlotModA2_PC
#MCZone
PlotModB5_PC <- effect_plot(ModB5, pred = PC1 , interval = TRUE, plot.points = TRUE, colors = "orange")
PlotModB5_PC <- PlotModB5_PC + theme_classic() + labs(x ="PC1", y = "MC shelter use (mean # individuals in shelter)")
PlotModB5_PC
#NTActivity
PlotModC4_PC <- effect_plot(ModC4, pred =  PC1, interval = TRUE, plot.points = TRUE, colors = "darkblue")
PlotModC4_PC <- PlotModC4_PC + theme_classic() + labs(x ="PC1", y = "NT movement (# events)")
PlotModC4_PC
#NTZone
PlotModD6_PC <- effect_plot(ModD6, pred =  PC1, interval = TRUE, plot.points = TRUE, colors = "darkblue")
PlotModD6_PC <- PlotModD6_PC + theme_classic() + labs(x ="PC1", y = "NT shelter use (mean # individuals in shelter)")
PlotModD6_PC

Figure5 <- ggarrange(PlotModA2_PC, PlotModB5_PC, PlotModC4_PC, PlotModD6_PC, labels = c("a", "b", "c","d"), ncol = 2, nrow = 2, common.legend =TRUE, legend = "none",align = "hv",  hjust = -2.5, label.x = 0.12)
Figure5

#Plotting association of environmental variables (export as 8x8)

Env1<- ggplot(Expt_All, aes(x=PC1, y=WT)) +
  theme_classic() + geom_smooth(method="lm", level=0.95, colour="black") + 
  geom_point() + labs(x ="PC1", y = "Temperature °C")
Env2 <- ggplot(Expt_All, aes(x=PC1, y=pH)) +
  theme_classic() + geom_smooth(method="lm", level=0.95, colour="black") + 
  geom_point() + labs(x ="PC1", y = "pH")
Env3 <- ggplot(Expt_All, aes(x=PC1, y=DO)) +
  theme_classic() + geom_smooth(method="lm", level=0.95, colour="black") + 
  geom_point() + labs(x ="PC1", y = "Dissolved oxygen (mg/l)")

ESM3 <- ggarrange(Env1, Env2, Env3, labels = c("a", "b", "c"), ncol = 2, nrow = 2, common.legend =TRUE, legend = "none",align = "hv", label.x = 0.12)
ESM3

#Plotting associations of behaviour variables, export as 12x9

MC_Fives1 <- ggplot(MC_Fives, aes(x=Movement, y=Shelter)) +
  theme_classic() + geom_point() + labs(x ="MC movement (#events)", y = "MC shelter use (mean # individuals in shelter)")
MC_Fives2 <- ggplot(MC_Fives, aes(x=Movement, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="MC movement (#events)", y = "MC chases initiated (# events)")
MC_Fives3 <- ggplot(MC_Fives, aes(x=Shelter, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="MC shelter use (mean # individuals in shelter)", y = "MC chases initiated (# events)")
MC_Tens1 <- ggplot(MC_Tens, aes(x=Movement, y=Shelter)) +
  theme_classic() + geom_point() + labs(x ="MC movement (#events)", y = "MC shelter use (mean # individuals in shelter)")
MC_Tens2 <- ggplot(MC_Tens, aes(x=Movement, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="MC movement (#events)", y = "MC chases initiated (# events)")
MC_Tens3 <- ggplot(MC_Tens, aes(x=Shelter, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="MC shelter use (mean # individuals in shelter)", y = "MC chases initiated (# events)")

ESM4 <- ggarrange(MC_Fives1, MC_Fives2, MC_Fives3, MC_Tens1, MC_Tens2, MC_Tens3, labels = c("a", "b", "c","d","e","f"), ncol = 3, nrow = 2, common.legend =TRUE, legend = "none",align = "hv", vjust =1.5, hjust=3.5, label.x = 0.12)
ESM4

NT_Fives1 <- ggplot(NT_Fives, aes(x=Movement, y=Shelter)) +
  theme_classic() + geom_point() + labs(x ="NT movement (#events)", y = "NT shelter use (mean # individuals in shelter)")
NT_Fives2 <- ggplot(NT_Fives, aes(x=Movement, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="NT movement (#events)", y = "NT chases initiated (# events)")
NT_Fives3 <- ggplot(NT_Fives, aes(x=Shelter, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="NT shelter use (mean # individuals in shelter)", y = "NT chases initiated (# events)")
NT_Tens1 <- ggplot(NT_Tens, aes(x=Movement, y=Shelter)) +
  theme_classic() + geom_point() + labs(x ="NT movement (#events)", y = "NT shelter use (mean # individuals in shelter)")
NT_Tens2 <- ggplot(NT_Tens, aes(x=Movement, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="NT movement (#events)", y = "NT chases initiated (# events)")
NT_Tens3 <- ggplot(NT_Tens, aes(x=Shelter, y=Chase)) +
  theme_classic() + geom_point() + labs(x ="NT shelter use (mean # individuals in shelter)", y = "NT chases initiated (# events)")

ESM5 <- ggarrange(NT_Fives1, NT_Fives2, NT_Fives3, NT_Tens1, NT_Tens2, NT_Tens3, labels = c("a", "b", "c","d","e","f"), ncol = 3, nrow = 2, common.legend =TRUE, legend = "none",align = "hv", vjust =1.5, hjust=3.5, label.x = 0.12)
ESM5

##end of code
