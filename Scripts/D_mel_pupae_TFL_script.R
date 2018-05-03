#### Loading packages and data in ####
rm(list = ls()) #clear R

#Load packages
library(ggplot2)
library('Cairo')
library(dplyr)
library(binom)
library(tidyr)
library(survival)

#Load data
setwd("C:/Users/bwalsh/Documents/R docs")
pupaeTFL<-read.csv("D.mel_pupae_TFL_DATA.csv", header=TRUE)

#### Rename vectors and convert from numeric to factors #####
#pupaeTFL$temp<-as.factor(pupaeTFL$temp)
pupaeTFL$block<-as.factor(pupaeTFL$block)
pupaeTFL$PCR<-as.factor(pupaeTFL$PCR)
pupaeTFL$plate_ID<-as.factor(pupaeTFL$plate_ID)
pupaeTFL$vial_ID<-as.factor(pupaeTFL$vial_ID)
pupaeTFL$emerged<-as.factor(pupaeTFL$emerged)
pupaeTFL$sex<-as.factor(pupaeTFL$sex)
pupaeTFL$survived_to_mating<-as.factor(pupaeTFL$survived_to_mating)
pupaeTFL$M_male_survived_mating<-as.factor(pupaeTFL$M_male_survived_mating)
pupaeTFL$M_female_survived_mating<-as.factor(pupaeTFL$M_female_survived_mating)
pupaeTFL$M_female_survived_laying<-as.factor(pupaeTFL$M_female_survived_laying)
pupaeTFL$M_larval<-as.factor(pupaeTFL$M_larval)
#pupaeTFL$F_female_survived_mating<-as.factor(pupaeTFL$F_male_survived_mating)
pupaeTFL$F_male_survived_mating<-as.factor(pupaeTFL$F_male_survived_mating)
#pupaeTFL$F_female_survived_laying1<-as.factor(pupaeTFL$F_female_survived_laying1)
pupaeTFL$F_larval1<-as.factor(pupaeTFL$F_larval1)
#pupaeTFL$F_female_survived_laying2<-as.factor(pupaeTFL$F_female_survived_laying2)
pupaeTFL$F_larval2<-as.factor(pupaeTFL$F_larval2)
#pupaeTFL$F_female_survived_laying3<-as.factor(pupaeTFL$F_female_survived_laying3)
pupaeTFL$F_larval3<-as.factor(pupaeTFL$F_larval3)

# STATS ####
#### 1. Kruskal Wallis ####
#Kruskal Wallis with larval vs temp
#kruskal.test(M_larval ~ temp, data=pupaeTFL)

#### 2. Binomail GLM ####
#Basic GLM with larval vs temp in males
#M1<-glm(M_larval~temp,data=pupaeTFL, family='binomial')
#anova(M1,test="Chisq")
#summary(M1)

# DATA MANAGEMENT ####
#### 1. Creating subsets of data for plotting ####

#Subset of individuals with emergence data
subemerged<- subset(pupaeTFL,emerged==1 |emerged==0)
subemerged$emerged<-as.numeric(subemerged$emerged)
subemerged$emerged<-subemerged$emerged-1

#Create subset with male larval
subpupaeTFL<- subset(pupaeTFL,M_larval==1 |M_larval==0)
subpupaeTFL$M_larval<-as.numeric(subpupaeTFL$M_larval)
subpupaeTFL$M_larval<-subpupaeTFL$M_larval-1 #changing back to normal values

#Create subset with female larval where data for larval 3 is available
FsubpupaeTFL<- subset(pupaeTFL, F_larval3==1 |F_larval3==0)
FsubpupaeTFL$F_larval1<-as.numeric(FsubpupaeTFL$F_larval1)
FsubpupaeTFL$F_larval2<-as.numeric(FsubpupaeTFL$F_larval2)
FsubpupaeTFL$F_larval3<-as.numeric(FsubpupaeTFL$F_larval3)
FsubpupaeTFL$F_larval1<-FsubpupaeTFL$F_larval1-1
FsubpupaeTFL$F_larval2<-FsubpupaeTFL$F_larval2-1
FsubpupaeTFL$F_larval3<-FsubpupaeTFL$F_larval3-1

Fsublarval1<-subset(pupaeTFL, F_larval1==1 | F_larval1==0)
Fsublarval1$F_larval1<-as.numeric(Fsublarval1$F_larval1)
Fsublarval1$F_larval1<-Fsublarval1$F_larval1-1



#### 2. Calculating proportions for binary data ####
prop_emerged <- summarise(group_by(subemerged, temp), mean(emerged), N=length(emerged))
prop_M_larval <- summarise(group_by(subpupaeTFL, temp), mean(M_larval), N=length(M_larval))
prop_F_larvaltot <- summarise(group_by(FsubpupaeTFL, temp), mean(F_larval1)+mean(F_larval2)+mean(F_larval3), N=length(F_larval1))
prop_F_larval1<-summarise(group_by(Fsublarval1, temp), mean(F_larval1), N=length(F_larval1))

# PLOTS ####

#### 1. Emerged with error bars ####
prop_emerged<-cbind((binom.confint(prop_emerged$`mean(emerged)`*prop_emerged$N, prop_emerged$N, conf.level = 0.95, methods = "logit")),prop_emerged)

pd<-position_dodge(1)

emerged_error<- ggplot(prop_emerged, aes(x=temp, y=`mean(emerged)`, ymin=lower, ymax=upper)) + 
                geom_point(position=pd, size=0.5, colour="black") +
                geom_line(position=pd, size=0.5, colour="black")+
                geom_errorbar(colour="black", position=pd, size=0.5, width=.2)+
                scale_x_continuous('Pupal heat shock temperature', breaks= round(seq(min(28), max(35), by = 1),1)) +
                scale_y_continuous('Proportion of emerging individuals')+
                coord_cartesian(ylim=c(0, 1)) +
                theme(panel.background = element_rect(fill = "white",colour = "gray90"),
                      axis.line=element_line(colour='gray90'),
                      axis.text=element_text(colour='black'),
                      panel.grid.major = element_line("gray90"),
                      axis.title.y = element_text(vjust=1.5),
                      plot.title = element_text(face="bold"),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      legend.position='none')

emerged_error

#### 2. Male fertility vs emergence ####

TFLvsCTL<- ggplot() +
           geom_point(data = prop_emerged, aes(temp, `mean(emerged)`, colour="chartreuse4"))+ 
           geom_point(data = prop_M_larval,aes(temp, `mean(M_larval)`, colour= "red")) + 
           geom_line(data = prop_emerged,  size = 1, aes(temp, `mean(emerged)`, colour="chartreuse4")) + 
           geom_line(data = prop_M_larval, size=1,  aes(temp, `mean(M_larval)`, colour="red")) + 
           xlab("Pupal heat-shock temperature") + 
           ylab ("Proportion") + 
           scale_x_continuous(breaks = round(seq(min(28), max(35), by = 1),1)) +
           #labs(colour='') +
           scale_color_manual("", values = c("red", "chartreuse4"), labels = c("Pupal emergence", "Surviving male fertility")) + 
           coord_cartesian(ylim=c(0, 1)) +
           theme(axis.title.x = element_text(size = 12, vjust=-.2), 
                 axis.title.y = element_text(size = 12, vjust=0.3), 
                 panel.background = element_rect(fill="white", colour="gray90") , 
                 panel.grid.major = element_line("gray90"),
                 legend.position = "top")

         
TFLvsCTL

#Saving high res plot
png(filename="C:/Users/bwalsh/Documents/R docs/TFLvsCTL.png", type="cairo", units="in", width=8, height=6, pointsize=10, res=1000)
print(TFLvsCTL)
dev.off()

#### 3. Male fertility with error bars ####

#Calculate ymin and ymax
prop_M_larval<-cbind((binom.confint(prop_M_larval$`mean(M_larval)`*prop_M_larval$N, prop_M_larval$N, conf.level = 0.95, methods = "logit")),prop_M_larval)

#Set position dodge
pd<-position_dodge(1)

TFL_error<- ggplot(prop_M_larval, aes(x=temp, y=`mean(M_larval)`, ymin=lower, ymax=upper)) + 
            geom_point(position=pd, size=0.5, colour="black") +
            geom_line(position=pd, size=0.5, colour="black")+
            geom_errorbar(colour="black", position=pd, size=0.5, width=.2)+
            scale_x_continuous('Pupal heat shock temperature', breaks= round(seq(min(28), max(35), by = 1),1)) +
            scale_y_continuous('Proportion of surviving males fertile')+
            coord_cartesian(ylim=c(0, 1)) +
            theme(panel.background = element_rect(fill = "white",colour = "gray90"),
                  axis.line=element_line(colour='gray90'),
                  axis.text=element_text(colour='black'),
                  panel.grid.major = element_line("gray90"),
                  axis.title.y = element_text(vjust=1.5),
                  plot.title = element_text(face="bold"),
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  legend.position='none')

TFL_error

png(filename="C:/Users/bwalsh/Documents/R docs/TFL_error.png", type="cairo", units="in", width=8, height=6, pointsize=10, res=1000)
print(TFL_error)
dev.off()

#### 4. Female fertility (larval_1 only) vs emergence ####



TFLvsCTL_female<- ggplot() +
                  geom_point(data = prop_emerged, aes(temp, `mean(emerged)`, colour="chartreuse4"))+ 
                  geom_point(data = prop_F_larval1,aes(temp, `mean(F_larval1)`, colour= "red")) + 
                  geom_line(data = prop_emerged,  size = 1, aes(temp, `mean(emerged)`, colour="chartreuse4")) + 
                  geom_line(data = prop_F_larval1, size=1,  aes(temp, `mean(F_larval1)`, colour="red")) + 
                  xlab("Pupal heat-shock temperature") + 
                  ylab ("Proportion") + 
                  scale_x_continuous(breaks = round(seq(min(28), max(35), by = 1),1)) +
                  #labs(colour='') +
                  scale_color_manual("", values = c("red", "chartreuse4"), labels = c("Pupal emergence", "Surviving female fertility")) + 
                  coord_cartesian(ylim=c(0, 1)) +
                  theme(axis.title.x = element_text(size = 12, vjust=-.2), 
                        axis.title.y = element_text(size = 12, vjust=0.3), 
                        panel.background = element_rect(fill="white", colour="gray90") , 
                        panel.grid.major = element_line("gray90"),
                        legend.position = "top")

TFLvsCTL_female

#High resolution TFL female
png(filename="C:/Users/bwalsh/Documents/R docs/TFL_female.png", type="cairo", units="in", width=8, height=6, pointsize=10, res=1000)
print(TFL_female)
dev.off()


#### 5. Female survival throughout TFL experiment ####

#Create subset with femlaes that have the full data
subF_survival<-subset(pupaeTFL, F_female_survived_laying3==1 |F_female_survived_laying3==0)
#subF_survival$F_female_survived_laying3<-as.numeric(subF_survival$F_female_survived_mating)

#Gather useful columns and assing values based on days
x<- gather(subF_survival, Days, Survival, emerged, survived_to_mating, F_female_survived_mating, F_female_survived_laying1, F_female_survived_laying2, F_female_survived_laying3)
x$Days<-sapply(x$Days,switch,'emerged'=0, 'survived_to_mating'=5,'F_female_survived_mating'=6,'F_female_survived_laying1'=8,'F_female_survived_laying2'=10, 'F_female_survived_laying3'=12)
x$temp<-as.factor(x$temp)

#Work out proportions alive at each time point by temperature
x$Survival<-as.numeric(x$Survival)
prop_F_survival <- summarise(group_by(x, temp, Days), mean(Survival),N=length(Survival))

#Normal ggplot of survival over time
F_survival<- ggplot(prop_F_survival, aes(x=Days, y=`mean(Survival)`, group=temp)) +
             geom_point(aes(colour=temp)) +
             geom_line(aes(colour=temp)) +
             coord_cartesian(ylim=c(0, 1), xlim=c(0,13)) +
             scale_x_continuous('Days after emergence', breaks= round(seq(min(0), max(13), by = 1),1)) +
             scale_y_continuous('Proportion of females alive', breaks=round(seq(min(0), max(1), by =0.1),1)) +
             labs(colour = "Temperature") +
             theme(panel.background = element_rect(fill = "white",colour = "gray90"),
                   axis.line=element_line(colour='gray90'),
                   axis.text=element_text(colour='black'),
                   panel.grid.major = element_line("gray90"),
                   axis.title.y = element_text(vjust=1.5),
                   plot.title = element_text(face="bold"),
                   plot.background = element_rect(fill = "transparent",colour = NA))

F_survival

#High resolution survival
png(filename="C:/Users/bwalsh/Documents/R docs/F_survival.png", type="cairo", units="in", width=8, height=6, pointsize=10, res=1000)
print(F_survival)
dev.off()
