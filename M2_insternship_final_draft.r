##############
## STAGE M2 ##
##############

#write.table(mydata, "c:/mydata.txt",row.names = F) 

### REQUIRED PACKAGE ----
library(multcomp)
library(nlme)
library(vegan)
library(IDPmisc)
library(reshape2)
library(ggplot2)
library(plyr)
library(psych)
library(lme4)
library(car)
library(multcompView)
library(ppcor)
library(gridExtra)
library(stargazer)


#### BIOLOGICAL ######################################################################

# VEGETATION COVER ----
# required dataset ----
cov <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/cover.csv");
cov0=cov;
str(cov);
cov$plot=as.factor(cov$plot);
cov$inner_replicate=as.factor(cov$inner_replicate);
cov$age <- factor(cov$age, levels=c("y","o"));
cov$state <- factor(cov$state, levels=c("N","R1","R2"));
cov$Plantation=paste(cov$state,cov$age); 
cov$Plantation <- factor(cov$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
colnames(cov)[13] <- "sp_8";
colnames(cov)[10] <- "asteraceae_2";
colnames(cov)[14] <- "fabaceae_1";
colnames(cov)[15] <- "fabaceae_2";
colnames(cov)[16] <- "fabaceae_3";
cov=cov[,-c(6,7,20:22)];
cov$total=apply(cov[,6:29],1,sum);
cov$specnumber=specnumber(cov[,6:29]);
cov$shannon=diversity(cov[,6:29]);
cov$pielou = cov$shannon/log(cov$specnumber);
which(cov$shannon==0);cov$pielou[c(8,35)]=0;
cov2=aggregate(cov,by=list(state=cov$state,age=cov$age,plot=cov$plot),mean);
cov2=cov2[,-c(5:8,33)];
cov2$Plantation=paste(cov2$state,cov2$age);
cov2$Plantation <- factor(cov2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
# div indices
# div ----
# total relative cover 
glmm0 = glmer(total ~ 1+(1|plot/inner_replicate),data=cov,family="poisson");
glmm1 = glmer(total ~ Plantation+(1|plot/inner_replicate),data=cov,family="poisson");
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11(); qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(glmm1), cov$Plantation);#x11(); plot(lme1); homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(cov$Plantation),letter=grpplot);
x11();ggplot(cov, aes(x=Plantation, y=total))  + geom_boxplot(aes(fill=age),notch=F) +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=level,y=250,label = letter)) +
  xlab("Plantation") +
  ylab("Total relative cover (%)")+guides(fill=FALSE)+
  theme(axis.title=element_text(size=20),axis.text=(element_text(size=15)));

# richness
lme1 = lme(fixed= specnumber ~ Plantation,random= ~ 1|plot/inner_replicate,data=cov);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), cov$Plantation);#x11();plot(lme1); hetero
lme1 = lme(fixed= log(specnumber) ~ Plantation,random= ~ 1|plot/inner_replicate,data=cov);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));#x11(); qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), cov$Plantation);#x11(); plot(lme1); homo
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters;m=aggregate(log(cov2$specnumber),by=list(cov2$Plantation),mean);
xx=data.frame(level=levels(cov$Plantation),letter=grpplot)
x11();ggplot(cov, aes(x=Plantation, y=specnumber))  + geom_boxplot(aes(fill=age),notch=F) +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=level,y=10,label = letter)) +
  xlab("Plantation") +
  ylab("specnumber")+guides(fill=FALSE);

# shannon
lme1 = lme(fixed= shannon ~ Plantation,random= ~ 1|plot/inner_replicate,data=cov);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), cov$Plantation);#x11();plot(lme1); hetero
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters;m=aggregate(log(cov2$shannon),by=list(cov2$Plantation),mean);
xx=data.frame(level=levels(cov$Plantation),letter=grpplot)
x11();ggplot(cov, aes(x=Plantation, y=shannon))  + geom_boxplot(aes(fill=age),notch=F) +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=level,y=2,label = letter)) +
  xlab("Plantation") +
  ylab("shannon")+guides(fill=FALSE);

# pielou
lme1 = lme(fixed= pielou ~ Plantation,random= ~ 1|plot/inner_replicate,data=cov);
summary(lme1);
anova(lme1);
x11();ggplot(cov, aes(x=Plantation, y=pielou))  + geom_boxplot(aes(fill=age),notch=F) +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("pielou weight (g/m²)")+guides(fill=FALSE);
# multivariate analysis ----
library(ade4);
# coa
coacover<-dudi.coa(cov[,6:29],scannf=F,nf=8); # 9 not good
# Kaiser-Guttman
barplot(coacover$eig, main="Eigenvalues", col="grey")
abline(h=mean(coacover$eig), col="red")
pvp=100*coacover$eig/sum(coacover$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(cov$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);#s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bcacover<-bca(coacover,cov$Plantation,scannf=FALSE);
pvp=100*bcacover$eig/sum(bcacover$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));s.class(bcacover$ls,cov$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);
title(xlab="Dim 1 (46%) ",ylab="Dim 2 (28%)",main="p-value: 0.001",font.main=1);
s.corcircle(bcacover$co,clabel=0.5,box=T,full=T);
# Permutation test
ptcover=rtest(bcacover,999)
ptcover # difference
x11();plot(ptcover,main="Between class inertia")


### LITTER COMPOSITION ----
# required dataset ----
lit <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/litter_0.csv")
lit0=lit
summary(lit)
str(lit)
#sapply(lit, class)
lit$plot=as.factor(lit$plot);
#levels(lit$plot) <- c("N", "R1","R2")
lit$inner_replicate=as.factor(lit$inner_replicate);
lit$total=apply(lit[,6:12],1,sum);
lit$Plantation=paste(lit$state,lit$age); # create a variable with gathered state and age
lit$Plantation=as.factor(lit$Plantation); # define it as factor
levels(lit$Plantation); # give the current level of $Plantation
lit$Plantation <- factor(lit$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
lit$age <- factor(lit$age, levels=c("y","o"));
lit[,6:18]=lit[,6:18]*16; # m²
#x11();plot(lit)
lit2=aggregate(lit,by=list(state=lit$state,age=lit$age,plot=lit$plot),mean); # create a tab with the mean of inner-replicate
lit2=lit2[,-c(4:8)]; # delete columns with is no more needed
lit2$Plantation=paste(lit2$state,lit2$age); # create a variable with gathered state and age
lit2$Plantation=as.factor(lit2$Plantation); # define it as factor
levels(lit2$Plantation); # give the current level of $Plantation
lit2$Plantation <- factor(lit2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
lit2$total=apply(lit2[,4:10],1,sum); # create a column with the total weight of each type (excluding intra compact weight)
#x11();par(mfrow = c(4, 2));plot(lit2$fresh~lit2$Plantation);plot(lit2$entire~lit2$Plantation);plot(lit2$fragment~lit2$Plantation);plot(lit2$skeleton~lit2$Plantation);plot(lit2$wood~lit2$Plantation);plot(lit2$cast~lit2$Plantation);plot(lit2$compact_total~lit2$Plantation);plot(lit2$total~lit2$Plantation)
#par(mfrow = c(1, 1));x11();plot(lit2$total-lit2$fresh-lit2$cast~lit2$Plantation)
## cor ----
pcor(lit[,6:12])$estimate
## total ----
# with Plantation
lme1 = lme(fixed= total ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), lit$Plantation);#x11();plot(lme1); hetero
lme1 = lme(fixed= log(total) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));#x11(); qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), lit$Plantation);#x11(); plot(lme1); homo
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters;m=aggregate(log(lit2$total),by=list(lit2$Plantation),mean);
xx=data.frame(level=levels(lit$Plantation),letter=grpplot)
x11();ggplot(lit, aes(x=Plantation, y=total))  + geom_boxplot(aes(fill=age),notch=F) +
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=level,y=2100,label = letter)) +
  xlab("Plantation") +
  ylab("Total weight (g/m²)")+guides(fill=FALSE);
## cast ----
lme2 = lme(fixed= cast ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit);
summary(lme2);
anova(lme2);
shapiro.test(residuals(lme2));#x11();qqnorm(residuals(lme2));qqline(residuals(lme2)) non normal
bartlett.test(residuals(lme2), lit$Plantation);#x11();plot(lme2); # hetero
lme2 = lme(fixed= log(cast+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit);
summary(lme2);
anova(lme2);
shapiro.test(residuals(lme2));#x11();qqnorm(residuals(lme2));qqline(residuals(lme2)); # normal
bartlett.test(residuals(lme2), lit$Plantation);#x11();plot(lme2); # hetero
lm0=lm(cast~1,data=lit2)
lm1=lm(cast~Plantation,data=lit2)
anova(lm0,lm1)
x11();ggplot(lit, aes(x=Plantation, y=cast))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Cast weight (g/m²)")+guides(fill=F)
## wood ----
# Plantation
lme4 = lme(fixed= wood ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme4)
anova(lme4) # Plantationment effect
shapiro.test(residuals(lme4)); #x11();qqnorm(residuals(lme4));qqline(residuals(lme4)) # non normal
bartlett.test(residuals(lme4), lit$Plantation); #x11();plot(lme4)  # hetero
lme4 = lme(fixed= log(wood) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme4)
anova(lme4) # Plantationment effect
shapiro.test(residuals(lme4));#x11(); qqnorm(residuals(lme4));qqline(residuals(lme4)) #  normal
bartlett.test(residuals(lme4), lit$Plantation); #x11();plot(lme4) # homo
summary(glht(lme4,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(lme4,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1 #  NP different de ctrl
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(lit$Plantation),letter=grpplot)
x11();ggplot(lit, aes(x=Plantation, y=wood))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=650,label = letter)) +
  xlab("Plantation") +
  ylab("Wood weight (g/m²)")+guides(fill=F)
# maybe because of wood
## Total without wood ----
lme1 = lme(fixed= total-wood ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), lit$Plantation);#x11();plot(lme1); hetero
lme1 = lme(fixed= log(total-wood) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));#x11(); qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), lit$Plantation);#x11(); plot(lme1); homo
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters;m=aggregate(log(lit2$total),by=list(lit2$Plantation),mean);
xx=data.frame(level=levels(lit$Plantation),letter=grpplot)
x11();ggplot(lit, aes(x=Plantation, y=total-wood))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=2000,label = letter)) +
  xlab("Plantation") +
  ylab("Total weight-wood (g/m²)")
# fresh ----
lme5 = lme(fixed= fresh ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme5)
anova(lme5) # Plantationment effect
shapiro.test(residuals(lme5)); #x11();qqnorm(residuals(lme5));qqline(residuals(lme5)) # non normal
bartlett.test(residuals(lme5), lit$Plantation); #x11();plot(lme5)  # hetero
lme5 = lme(fixed= log(fresh+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme5)
anova(lme5) # Plantationment effect
shapiro.test(residuals(lme5)); #x11();qqnorm(residuals(lme5));qqline(residuals(lme5)) #  non normal
bartlett.test(residuals(lme5), lit$Plantation); #x11();plot(lme5) # homo
summary(glht(lme5,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme5,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(lit$Plantation),letter=grpplot)
x11();ggplot(lit, aes(x=Plantation, y=fresh))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=350,label = letter)) +
  xlab("Plantation") +
  ylab("Fresh weight (g/m²)")+guides(fill=F)
## entire ----
lme7=lme(fixed= entire ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme7)
anova(lme7) 
shapiro.test(residuals(lme7)); #x11();qqnorm(residuals(lme7));qqline(residuals(lme7)) # non normal
bartlett.test(residuals(lme7), lit$Plantation); #x11();plot(lme7)  # hetero
lme7=lme(fixed= log(entire+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme7)
anova(lme7)
shapiro.test(residuals(lme7)); #x11();qqnorm(residuals(lme7));qqline(residuals(lme7)) # non normal
bartlett.test(residuals(lme7), lit$Plantation); #x11();plot(lme7) # hetero
lm0=lm(entire~1,data=lit2)
lm1=lm(entire~Plantation,data=lit2)
anova(lm0,lm1)
## fragment ----
lme7 = lme(fixed= fragment ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme7)
anova(lme7) # Plantationment effect
shapiro.test(residuals(lme7)); #x11();qqnorm(residuals(lme7));qqline(residuals(lme7)) # non normal
bartlett.test(residuals(lme7), lit$Plantation); #x11();plot(lme7)  # hetero
lme7 = lme(fixed= log(fragment+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme7)
anova(lme7) # Plantationment effect
shapiro.test(residuals(lme7)); x11();qqnorm(residuals(lme7));qqline(residuals(lme7)) #  normal
bartlett.test(residuals(lme7), lit$Plantation); #x11();plot(lme7) # homo
summary(glht(lme7,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme7,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(lit$Plantation),letter=grpplot)
x11();ggplot(lit, aes(x=Plantation, y=fragment))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=60,label = letter)) +
  xlab("Plantation") +
  ylab("Fragment weight (g/m²)")+guides(fill=F)
# skeleton ----
lme7=lme(fixed= skeleton ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme7)
anova(lme7) 
shapiro.test(residuals(lme7)); #x11();qqnorm(residuals(lme7));qqline(residuals(lme7)) # non normal
bartlett.test(residuals(lme7), lit$Plantation); #x11();plot(lme7) # hetero
lme7=lme(fixed= log(skeleton+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme7)
anova(lme7) 
shapiro.test(residuals(lme7)); #x11();qqnorm(residuals(lme7));qqline(residuals(lme7)) # non normal
bartlett.test(residuals(lme7), lit$Plantation); #x11();plot(lme7) # hetero
lm0=lm(log(skeleton+1)~1,data=lit2)
lm1=lm(log(skeleton+1)~Plantation,data=lit2)
anova(lm0,lm1)
shapiro.test(residuals(lm1)); #x11();qqnorm(residuals(lme7));qqline(residuals(lme7)) # non normal
bartlett.test(residuals(lm1), lit2$Plantation); #x11();plot(lme7) # hetero
kruskal.test(lit2$skeleton,lit2$Plantation) # significant 
test=pairwise.wilcox.test(lit2$skeleton, lit2$Plantation, p.adj = "bonf");test
test$p.value # non significant
## compact total ----
lme7=lme(fixed= compact_total ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit)
summary(lme7)
anova(lme7) # no difference
## adjusted p.val  ----
#bonferroni simple sequentialy rejective test
pval=data.frame(test=c("total-cast","wood","total"),p.val=c(0.0002,0.0001,0.0001),factor=c(1:3))
pval$p.val.adjusted=pval$p.val*pval$factor;pval
## GRAPH ----
# STAR PLOT FOR LITTER COMPOSITION BY Plantation
lit3=aggregate(lit2,by=list(state=lit2$state,age=lit2$age),mean); # tab with mean by Plantation
lit3=lit3[,-c(3:5,19)]; # delete NA
lit3$Plantation=paste(lit3$state,lit3$age);
row.names(lit3) <- c("NO","R1O","R2O","R1Y","R2Y"); # define row names
lit3=lit3[c(1,4,2,5,3),]; # put the row in the right order
palette(rainbow(12, s = 0.6, v = 0.75));
x11();stars(lit3[, 3:9], len = 0.8, key.loc = c(4.7, 2.3),main = "Litter structure", draw.segments = TRUE);

library("fmsb")
x11();radarchart(lit3[, 3:9],pfcol=c(1:5),plty=1,pdensity=0,axistype=4,seg=4,plwd=2,maxmin=F)
legend("topleft", legend=unique(lit3$Plantation), seg.len=1.5, title="state", pch=1,bty="n" ,lwd=3, y.intersp=0.5, horiz=FALSE, col=c(1:5))

# STACKED BAR CHART
# without the compact type (should do the same for all the compact type)
library(reshape2)
lit5= melt(lit3[,c(3:9,16)], id.var="Plantation"); # melt to get the long format
colnames(lit5)[3] <- "weight";
lit5$Plantation <- factor(lit5$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
library(ggplot2)
#ggplot(lit5, aes(x = Plantation, y = weight, fill = variable)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))
# with plotted weight
lit5$round=round(lit5$weight,1); # round the weight for light plotting
library(plyr);
lit5=ddply(lit5, "Plantation", transform, label_y=cumsum(weight)); # add the cumsum in order to give the text position
#x11();ggplot(lit5, aes(x = Plantation, y = weight, fill = variable)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
#
lme1 = lme(fixed= log(total) ~ Plantation,random= ~ 1|plot/inner_replicate,data=lit);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni')); grpplot=grp1$mcletters$Letters;m=aggregate(log(lit2$total),by=list(lit2$Plantation),mean); xx=data.frame(level=levels(lit$Plantation),letter=grpplot);
#
ce=arrange(lit5,Plantation,variable); # to have number in the middle of each bar
  ce <- ddply(ce, "Plantation", transform, label_y=cumsum(weight)-0.5*weight);
  x11();ggplot(ce, aes(x = Plantation, y = weight, fill = variable)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)+   geom_text(data=xx, aes(x=level,y=1100,label = letter,fill=NULL));
# by relative composition
lit55=ddply(lit5, "Plantation", transform,percent_weight = weight / sum(weight) * 100); # add % in the tab
lit55$Plantation <- factor(lit55$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
lit55$roundp=round(lit55$percent_weight,1); # round the weight for light plotting
lit55=ddply(lit55, "Plantation", transform, label_y=cumsum(percent_weight)); # add the cumsum in order to give the text position
x11();ggplot(lit55, aes(x = Plantation, y = percent_weight, fill = variable)) + geom_bar(stat = "identity")+guides(fill=guide_legend(reverse=TRUE))+  geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3);
# graph by litter compound
lit6= melt(lit3[,c(3:9,15,16)], id.var="Plantation") # add the total in this tab
colnames(lit6)[3] <- "weight"
lit6$Plantation <- factor(lit6$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#x11();c <- ggplot(lit6, aes(x = Plantation, y = weight));c + facet_wrap(~ variable) + geom_bar(stat = "identity")
# total minus cast
lit4=lit3
lit4$total_minus_cast=lit4$total-lit4$cast
lit7= melt(lit4[,c(3:9,15,16,17)], id.var="Plantation") # add the total in this tab
colnames(lit7)[3] <- "weight"
lit7$Plantation <- factor(lit7$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(lit7, aes(x = Plantation, y = weight));c + facet_wrap(~ variable) + geom_bar(stat = "identity")#----
## Multivariate analysis ----
palette("default")
library(ade4)
# PCA
pcalitter<-dudi.pca(lit[,6:12],scale=T,scannf=F,nf=3);
# Kaiser-Guttman
barplot(pcalitter$eig, main="Eigenvalues", col="grey")
abline(h=mean(pcalitter$eig), col="red")
pvp=100*pcalitter$eig/sum(pcalitter$eig);pvp;
cumsum(pvp)
#x11();s.class(pca2$li,as.factor(lit$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);#s.arrow(pca2$co,boxes=F,add.plot=T)
# BCA
bcalitter<-bca(pcalitter,lit$Plantation,scannf=FALSE);
pvp=100*bcalitter$eig/sum(bcalitter$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));s.class(bcalitter$ls,lit$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);
title(xlab="Dim 1 (42%) ",ylab="Dim 2 (31%)",main="p-value: 0.001",font.main=1);
s.corcircle(bcalitter$co,clabel=0.5,box=T,full=T);
# Permutation test
ptlitter=rtest(bcalitter,999)
ptlitter # difference
x11();plot(ptlitter,main="Between class inertia")

# without cast
pca2<-dudi.pca(lit[,c(6:10,12)],scale=T,scannf=F,nf=3)
pvp=100*pca2$eig/sum(pca2$eig);pvp;
cumsum(pvp)
#x11();s.class(pca2$li,as.factor(lit$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);#s.arrow(pca2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(pca2,lit$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,lit$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5)#;s.arrow(bca2$co,boxes=F,add.plot=T)
title(xlab="Dim 1 (45%) ",ylab="Dim 2 (28%)")
s.corcircle(bca2$co,box=T)
# Permutation test
pt2=rtest(bca2,999)
pt2 # difference
x11();plot(pt2,main="Between class inertia")

### ROOTS ----
roo<- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/roots.csv", dec=",");
str(roo)
roo$plot=as.factor(roo$plot);
roo$inner_replicate=as.factor(roo$inner_replicate);
roo$Plantation=paste(roo$state,roo$age);
roo$Plantation=as.factor(roo$Plantation);
levels(roo$Plantation)
roo$Plantation <- factor(roo$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));

lme1 = lme(fixed= dw_.of.root ~ Plantation,random= ~ 1|plot/inner_replicate,data=roo);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch05$Plantation);#x11();plot(lme1); hetero
lme1 = lme(fixed= log(dw_.of.root) ~ Plantation,random= ~ 1|plot/inner_replicate,data=roo);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), ch05$Plantation);#x11();plot(lme1); homo
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(roo$Plantation),letter=grpplot);
x11();ggplot(roo, aes(x=Plantation, y=dw_.of.root))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, size=9,aes(x=c(1:5),y=7,label = letter)) +
  xlab("Plantation") +
  ylab("Root weight (g)")+guides(fill=F);
# less root in y

## MACROFAUNA ----
# Required dataset ----
md <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_density_0.csv");
md[md$strata=="Litter",c(6:84)]=md[md$strata=="Litter",c(6:84)]*16;
md[md$strata=="Soil",c(6:84)]=md[md$strata=="Soil",c(6:84)]*160;
md$oligochaeta=md$oligochaeta+md$oligochaeta_diapause; # sum oligochaeta
names(md[c(6,12,14,22,24,80:84)]);
md=md[,-c(6,12,14,22,24,80:84)]; # remove mesofauna & egg & oligochaeta_diapause
str(md);
md$plot=as.factor(md$plot);
md$inner_replicate=as.factor(md$inner_replicate);
md$age <- factor(md$age, levels=c("y","o"));
md$Plantation=paste(md$state,md$age);
md$Plantation=as.factor(md$Plantation);
md$Plantation <- factor(md$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
md$sp=apply(md[,20:74],1,sum); # sum sp
md$total=apply(md[,6:74],1,sum);
md$specnumber=specnumber(md[,6:74]);
md$shannon=diversity(md[,6:74]);
md[md$specnumber==0,]$shannon=NA; # non sense
md$pielou = md$shannon/log(md$specnumber); # NaN problem:
which(md$shannon==0);md$pielou[c(1,5,9,31,41,55,61,77,83)]=0;
# functionnal group
# engineer: hymenoptera, isoptera, oligochaeta
# carnivore: araneae, chilopoda
# omnivore: blattodea
# herbivore: gastropoda, hemiptera, lepidoptera, orthoptera 12131719
# detritivore: coleoptera, coleoptera_larva, diplopoda, isopoda 9101115
# + sp
fg=data.frame(md[,c(1:5)],engineer=md[,14]+md[,16]+md[,18],carnivore=md[,6]+md[,8],omnivore=md[,7],herbivore=md[,12]+md[,13]+md[,17]+md[,19],detritivore=md[,9]+md[,10]+md[,11]+md[,15]);
md=data.frame(md,fg[,c(6:10)]);
mdl=subset(md,strata=="Litter");
mds=subset(md,strata=="Soil");
mdtot=data.frame(mds[,c(1:4,75)],mds[,c(6:74,76:78,81:85)]+mdl[,c(6:74,76:78,81:85)]);
mdtot$shannon=diversity(mdtot[,6:74]);
mdtot$pielou= mdtot$shannon/log(mdtot$specnumber);
# sum inner replicate
mdtot2=aggregate(mdtot[,c(6:19,75:84)],by=list(state=mdtot$state,age=mdtot$age,plot=mdtot$plot),median); # create a tab with the mean of inner-replicate
mdtot2$Plantation=paste(mdtot2$state,mdtot2$age); # create a variable with gathered state and age
mdtot2$Plantation=as.factor(mdtot2$Plantation); # define it as factor
levels(mdtot2$Plantation); # give the current level of $Plantation
mdtot2$Plantation <- factor(mdtot2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mdl2=aggregate(mdl[,c(6:19,74,76:85)],by=list(state=mdl$state,age=mdl$age,plot=mdl$plot),median); # create a tab with the mean of inner-replicate
mdl2$Plantation=paste(mdl2$state,mdl2$age); # create a variable with gathered state and age
mdl2$Plantation=as.factor(mdl2$Plantation); # define it as factor
levels(mdl2$Plantation); # give the current level of $Plantation
mdl2$Plantation <- factor(mdl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mds2=aggregate(mds[,c(6:19,74,76:85)],by=list(state=mds$state,age=mds$age,plot=mds$plot),median); # create a tab with the mean of inner-replicate
mds2$Plantation=paste(mds2$state,mds2$age); # create a variable with gathered state and age
mds2$Plantation=as.factor(mds2$Plantation); # define it as factor
levels(mds2$Plantation); # give the current level of $Plantation
mds2$Plantation <- factor(mds2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mdol <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_density_0.csv");
mdol$plot=as.factor(mdol$plot);
mdol$inner_replicate=as.factor(mdol$inner_replicate);
mdol$age <- factor(mdol$age, levels=c("y","o"));
mdol$Plantation=paste(mdol$state,mdol$age);
mdol$Plantation=as.factor(mdol$Plantation);
mdol$Plantation <- factor(mdol$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mdoll=subset(mdol,strata=="Litter");
mdols=subset(mdol,strata=="Soil");
mdtol=data.frame(mdols[,c(1:4,85)],mdols[,c(21,22)]+mdoll[,c(21,22)]);
#
mdtol2=aggregate(mdtol[,6:7],by=list(state=mdtot$state,age=mdtot$age,plot=mdtot$plot),median); # create a tab with the mean of inner-replicate
mdtol2$Plantation=paste(mdtol2$state,mdtol2$age); # create a variable with gathered state and age
mdtol2$Plantation=as.factor(mdtol2$Plantation); # define it as factor
levels(mdtol2$Plantation); # give the current level of $Plantation
mdtol2$Plantation <- factor(mdtol2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mds2=aggregate(mds,by=list(state=mds$state,age=mds$age,plot=mds$plot),mean);
mds2=mds2[,-c(4:8)];
mds2$Plantation=paste(mds2$state,mds2$age); # create a variable with gathered state and age
mds2$Plantation=as.factor(mds2$Plantation); # define it as factor
levels(mds2$Plantation); # give the current level of $Plantation
mds2$Plantation <- factor(mds2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
#
mb <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_biomass_0.csv");
mb[mb$strata=="Litter",c(6:84)]=mb[mb$strata=="Litter",c(6:84)]*16;
mb[mb$strata=="Soil",c(6:84)]=mb[mb$strata=="Soil",c(6:84)]*160;

str(mb);
mb$plot=as.factor(mb$plot);
mb$inner_replicate=as.factor(mb$inner_replicate);

mb$oligochaeta=mb$oligochaeta+mb$oligochaeta_diapause; # sum oligochaeta
names(mb[c(6,12,14,22,24,80:84)]);
mb=mb[,-c(6,12,14,22,24,80:84)]; # remove mesofauna & egg & oligochaeta_diapause
mb$age <- factor(mb$age, levels=c("y","o"));
mb$Plantation=paste(mb$state,mb$age);
mb$Plantation=as.factor(mb$Plantation);
mb$Plantation <- factor(mb$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mb$sp=apply(mb[,20:74],1,sum); # sum sp
mb$biomass=apply(mb[,6:74],1,sum);
#
## functionnal group
# engineer: hymenoptera, isoptera, oligochaeta
# carnivore: araneae, chilopoda
# omnivore: blattodea
# herbivore: gastropoda, hemiptera, lepidoptera, orthoptera 12131719
# detritivore: coleoptera, coleoptera_larva, diplopoda, isopoda 9101115
fg=data.frame(mb[,c(1:5)],engineer=mb[,14]+mb[,16]+mb[,18],carnivore=mb[,6]+mb[,8],omnivore=mb[,7],herbivore=mb[,12]+mb[,13]+mb[,17]+mb[,19],detritivore=mb[,9]+mb[,10]+mb[,11]+mb[,15],sp=mb[,76]);
mb=data.frame(mb,fg[,c(6:10)]);
mbl=subset(mb,strata=="Litter");
mbs=subset(mb,strata=="Soil");
mbtot=data.frame(mbs[,c(1:4,75)],mbs[,c(6:74,76:82)]+mbl[,c(6:74,76:82)]);
#
mbl2=aggregate(mbl,by=list(state=mbl$state,age=mbl$age,plot=mbl$plot),mean);
mbl2=mbl2[,-c(4:8)];
mbl2$Plantation=paste(mbl2$state,mbl2$age); # create a variable with gathered state and age
mbl2$Plantation=as.factor(mbl2$Plantation); # define it as factor
levels(mbl2$Plantation); # give the current level of $Plantation
mbl2$Plantation <- factor(mbl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mdl2=aggregate(mdl,by=list(state=mdl$state,age=mdl$age,plot=mdl$plot),mean);
mdl2=mdl2[,-c(4:8)];
mdl2$Plantation=paste(mdl2$state,mdl2$age); # create a variable with gathered state and age
mdl2$Plantation=as.factor(mdl2$Plantation); # define it as factor
levels(mdl2$Plantation); # give the current level of $Plantation
mdl2$Plantation <- factor(mdl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
#
mbs2=aggregate(mbs,by=list(state=mbs$state,age=mbs$age,plot=mbs$plot),mean); # create a tab with the mean of inner-replicate
mbs2=mbs2[,-c(4:8)]; # delete columns with is no more needed
mbs2$Plantation=paste(mbs2$state,mbs2$age); # create a variable with gathered state and age
mbs2$Plantation=as.factor(mbs2$Plantation); # define it as factor
levels(mbs2$Plantation); # give the current level of $Plantation
mbs2$Plantation <- factor(mbs2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mbol <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_biomass_0.csv");
mbol$plot=as.factor(mbol$plot);
mbol$inner_replicate=as.factor(mbol$inner_replicate);
mbol$age <- factor(mbol$age, levels=c("y","o"));
mbol$Plantation=paste(mbol$state,mbol$age);
mbol$Plantation=as.factor(mbol$Plantation);
mbol$Plantation <- factor(mbol$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mboll=subset(mbol,strata=="Litter");
mbols=subset(mbol,strata=="Soil");
mbtol=data.frame(mbols[,c(1:4,85)],mbols[,c(21,22)]+mboll[,c(21,22)]);
#
mbtol2=aggregate(mbtol,by=list(state=mbtol$state,age=mbtol$age,plot=mbtol$plot),mean); # create a tab with the mean of inner-replicate
mbtol2=mbtol2[,-c(4:8)]; # delete columns with is no more needed
mbtol2$Plantation=paste(mbtol2$state,mbtol2$age); # create a variable with gathered state and age
mbtol2$Plantation=as.factor(mbtol2$Plantation); # define it as factor
levels(mbtol2$Plantation); # give the current level of $Plantation
mbtol2$Plantation <- factor(mbtol2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#

## LITTER DENSITY----
#### litter -----
# density -----
## cor ----
pcor(mdl[,c(6:19,76)])$estimate
## total ----
glmm0 = glmer(total ~ 1+(1|plot/inner_replicate),data=mdl,family="poisson");
glmm1 = glmer(total ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson");
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mdl$Plantation),letter=grpplot);
x11();ggplot(mdl, aes(x=Plantation, y=total))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=600,label = letter)) +
  xlab("Plantation") +
  ylab("Total density (m²)")+guides(fill=F);
# more in R1o less in R2o
# specnumber ----
glmm0 = glmer(specnumber ~ 1+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(specnumber ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # no difference
# no difference between specie total number; idem for litter (mdl) or soil (mds) separated
x11();ggplot(mdl, aes(x=Plantation, y=specnumber))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Taxon richness")
# shannon ----
mdl1=na.exclude(mdl)
lme1 = lme(fixed= shannon ~ Plantation,random= ~ 1|plot/inner_replicate,data=mdl1)
summary(lme1)
anova(lme1) # no difference
# pielou ----
lme1 = lme(fixed= pielou ~ Plantation,random= ~ 1|plot/inner_replicate,data=mdl1)
summary(lme1)
anova(lme1) # no difference
# araneae ----
glmm0 = glmer(araneae ~ 0+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(araneae ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
lm0=lm(mdl2$araneae~1)
lm1=lm(mdl2$araneae~mdl2$Plantation)
anova(lm0,lm1)
# coleoptera ----
glmm0 = glmer(coleoptera ~ 0+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(coleoptera ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mdl$Plantation),letter=grpplot)
x11();ggplot(mdl, aes(x=Plantation, y=coleoptera))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=150,label = letter)) +
  xlab("Plantation") +
  ylab("Coleoptera density (m²)")
# hemiptera ----
glmm0 = glmer(hemiptera ~ 0+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(hemiptera ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # diff
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # non normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # hetero
lm0=lm(mdl2$hemiptera~1)
lm1=lm(mdl2$hemiptera~mdl2$Plantation)
anova(lm0,lm1) # not enough
# hymenoptera ----
glmm0 = glmer(hymenoptera ~ 0+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(hymenoptera ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mdl$Plantation),letter=grpplot)
x11();ggplot(mdl, aes(x=Plantation, y=hymenoptera))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=300,label = letter)) +
  xlab("Plantation") +
  ylab("Hymenoptera density (m²)")
# isopoda ----
glmm0 = glmer(isopoda ~ 0+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(isopoda ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) 
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
lm0=lm(mdl2$isopoda~1)
lm1=lm(mdl2$isopoda~mdl2$Plantation)
anova(lm0,lm1)
# isoptera ----
glmm0 = glmer(isoptera ~ 0+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(isoptera ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
lm0=lm(mdl2$isoptera~1)
lm1=lm(mdl2$isoptera~mdl2$Plantation)
anova(lm0,lm1)
# oligochaeta (sum) ----
glmm0 = glmer(oligochaeta ~ 0+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(oligochaeta ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1) # not enouth sp
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # non normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # hetero
lm0=lm(mdl2$oligochaeta~1)
lm1=lm(mdl2$oligochaeta~mdl2$Plantation)
anova(lm0,lm1)
# Engineer ----
glmm0 = glmer(engineer ~ 1+(1|plot/inner_replicate),data=mdl,family="poisson");
glmm1 = glmer(engineer ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson");
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mdl$Plantation),letter=grpplot);
x11();ggplot(mdl, aes(x=Plantation, y=engineer))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=410,label = letter)) +
  xlab("Plantation") +
  ylab("Engineer density (m²)")+guides(fill=F);
# less in  R1o more in young (hymenoptera)
# carnivore ----
glmm0 = glmer(carnivore ~ 1+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(carnivore ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # non normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
lm0=lm(mdl2$carnivore~1)
lm1=lm(mdl2$carnivore~mdl2$Plantation)
anova(lm0,lm1)
# omnivore ----
glmm0 = glmer(omnivore ~ 1+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(omnivore ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) #significative
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # non normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
lm0=lm(mdl2$omnivore~1)
lm1=lm(mdl2$omnivore~mdl2$Plantation)
anova(lm0,lm1)
# herbivore ----
glmm0 = glmer(herbivore ~ 1+(1|plot/inner_replicate),data=mdl,family="poisson")
glmm1 = glmer(herbivore ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # non normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # hetero
lm0=lm(mdl2$herbivore~1)
lm1=lm(mdl2$herbivore~mdl2$Plantation)
anova(lm0,lm1)
# detritivore ----
glmm0 = glmer(detritivore ~ 1+(1|plot/inner_replicate),data=mdl,family="poisson");
glmm1 = glmer(detritivore ~ Plantation+(1|plot/inner_replicate),data=mdl,family="poisson");
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mdl$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mdl$Plantation),letter=grpplot);
x11();ggplot(mdl, aes(x=Plantation, y=detritivore))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=150,label = letter)) +
  xlab("Plantation") +
  ylab("Detritivore density (m²)")+guides(fille=F);
# graph ----
mdl2=aggregate(mdl,by=list(state=mdl$state,age=mdl$age),mean)
mdl2$Plantation=paste(mdl2$state,mdl2$age)
mdl2$Plantation=as.factor(mdl2$Plantation)
mdl2$Plantation <- factor(mdl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"))
mdl3= melt(mdl2[,c(8:21,77,78)], id.var="Plantation") # melt to get the long format
colnames(mdl3)[3] <- "Density"
colnames(mdl3)[2] <- "Taxon"
mdl3$Plantation <- factor(mdl3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#library(ggplot2)
# with plotted weight
mdl3$round=round(mdl3$Density,1) # round the Density for light plotting
#library(plyr)
mdl3=ddply(mdl3, "Plantation", transform, label_y=cumsum(Density)) # add the cumsum in order to give the text position
ce=arrange(mdl3,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Density)-0.5*Density)
x11();ggplot(ce, aes(x = Plantation, y = Density, fill = Taxon)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mdl33=ddply(mdl3, "Plantation", transform,percent_Density = Density / sum(Density) * 100) # add % in the tab
mdl33$Plantation <- factor(mdl33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mdl33$roundp=round(mdl33$percent_Density,1) # round the Density for light plotting
mdl33=ddply(mdl33, "Plantation", transform, label_y=cumsum(percent_Density)) # add the cumsum in order to give the text position
ce=arrange(mdl33,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Density)-0.5*percent_Density)
x11();ggplot(ce, aes(x = Plantation, y = percent_Density, fill = Taxon)) + geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)
# graph by compound
mdl4= melt(mdl2[,c(8:21,77,78,79)], id.var="Plantation") # melt to get the long format
colnames(mdl4)[3] <- "Density"
colnames(mdl4)[2] <- "Taxon"
mdl4$Plantation <- factor(mdl4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(mdl4, aes(x = Plantation, y = Density));c + facet_wrap(~ Taxon,scales="free_y") + geom_bar(stat = "identity")
#
## same for functional group
mdl6=aggregate(mdl[,c(81:85)],by=list(state=mdl$state,age=mdl$age),mean)
mdl6$Plantation=paste(mdl2$state,mdl2$age)
mdl6$Plantation=as.factor(mdl2$Plantation)
mdl6$Plantation <- factor(mdl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"))
#GRAPH
mdl3= melt(mdl6[,c(3:8)], id.var="Plantation") # melt to get the long format
colnames(mdl3)[3] <- "Density"
colnames(mdl3)[2] <- "Functional_group"
mdl3$Plantation <- factor(mdl3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#library(ggplot2)
# with plotted weight
mdl3$round=round(mdl3$Density,1) # round the Density for light plotting
#library(plyr)
mdl3=ddply(mdl3, "Plantation", transform, label_y=cumsum(Density)) # add the cumsum in order to give the text position
ce=arrange(mdl3,Plantation,Functional_group) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Density)-0.5*Density)
x11();ggplot(ce, aes(x = Plantation, y = Density, fill = Functional_group)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mdl33=ddply(mdl3, "Plantation", transform,percent_Density = Density / sum(Density) * 100) # add % in the tab
mdl33$Plantation <- factor(mdl33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mdl33$roundp=round(mdl33$percent_Density,1) # round the Density for light plotting
mdl33=ddply(mdl33, "Plantation", transform, label_y=cumsum(percent_Density)) # add the cumsum in order to give the text position
ce=arrange(mdl33,Plantation,Functional_group) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Density)-0.5*percent_Density)
x11();ggplot(ce, aes(x = Plantation, y = percent_Density, fill = Functional_group)) + geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)
# graph by compound
mdl5= melt(mdl2[,c(83:87,77,79)], id.var="Plantation") # melt to get the long format
colnames(mdl5)[3] <- "Density"
colnames(mdl5)[2] <- "Taxon"
mdl5$Plantation <- factor(mdl5$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(mdl5, aes(x = Plantation, y = Density));c + facet_wrap(~ Taxon) + geom_bar(stat = "identity")
## LITTER BIOMASS -----
## total ----
lme1 = lme(fixed= biomass ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # no difference
## coleoptera ----
lme1 = lme(fixed= coleoptera ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal QQplot
bartlett.test(residuals(lme1), mbl$Plantation);#x11();plot(lme1); 
lme1 = lme(fixed= log(coleoptera+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), mbl$Plantation);#x11();plot(lme1); hetero
lm0=lm(mbl2$coleoptera~1)
lm1=lm(mbl2$coleoptera~mbl2$Plantation)
anova(lm0,lm1)
## hymenoptera ----
lme1 = lme(fixed= hymenoptera ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # no difference
## isoptera ----
lme1 = lme(fixed= isoptera ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # no difference
## carnivore ----
lme1 = lme(fixed= carnivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # no difference
## herbivore ----
lme1 = lme(fixed= herbivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # no difference
# Engineer ----
lme1 = lme(fixed= engineer ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # no diff
# omnivore ----
lme1 = lme(fixed= omnivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # no diff
# detritivore ----
lme1 = lme(fixed= detritivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # diff
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), mbl$Plantation); # hetero
lme1 = lme(fixed= log(detritivore+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbl);
summary(lme1);
anova(lme1); # diff
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), mbl$Plantation); # hetero
lm0=lm(detritivore~1,data=mbl2)
lm1=lm(detritivore~Plantation,data=mbl2)
anova(lm0,lm1)
shapiro.test(residuals(lm1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lm1), mbl2$Plantation); # hetero
lsmeans(lm1,pairwise~Plantation,data=mbl2,adjust='bonferroni')
lettre=cld(summary(glht(lm1,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")))
lettre$mcletters
#let=scan(what="")
let=c("a","b","a","a","a")
xx=data.frame(level=levels(mbl2$Plantation),letter=let)
x11();ggplot(mbl, aes(x=Plantation, y=detritivore))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=30,label = letter)) +
  xlab("Plantation") +
  ylab("Detritivore biomass (m³)")
## graph -----
mbl2=aggregate(mbl,by=list(state=mbl$state,age=mbl$age),mean)
mbl2$Plantation=paste(mbl2$state,mbl2$age)
mbl2$Plantation=as.factor(mbl2$Plantation)
mbl2$Plantation <- factor(mbl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"))
mbl3= melt(mbl2[,c(8:21,77,78)], id.var="Plantation") # melt to get the long format
colnames(mbl3)[3] <- "Biomass"
colnames(mbl3)[2] <- "Taxon"
mbl3$Plantation <- factor(mbl3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#library(ggplot2)
# with plotted weight
mbl3$round=round(mbl3$Biomass,2) # round the Biomass for light plotting
#library(plyr)
mbl3=ddply(mbl3, "Plantation", transform, label_y=cumsum(Biomass)) # add the cumsum in order to give the text position
ce=arrange(mbl3,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Biomass)-0.5*Biomass)
x11();ggplot(ce, aes(x = Plantation, y = Biomass, fill = Taxon)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mbl33=ddply(mbl3, "Plantation", transform,percent_Biomass = Biomass / sum(Biomass) * 100) # add % in the tab
mbl33$Plantation <- factor(mbl33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mbl33$roundp=round(mbl33$percent_Biomass,1) # round the Biomass for light plotting
mbl33=ddply(mbl33, "Plantation", transform, label_y=cumsum(percent_Biomass)) # add the cumsum in order to give the text position
ce=arrange(mbl33,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Biomass)-0.5*percent_Biomass)
x11();ggplot(ce, aes(x = Plantation, y = percent_Biomass, fill = Taxon)) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)
# graph by compound
mbl22=melt(mbl2[,c(7:21,77,78)], id.var="Plantation")
colnames(mbl22)[3] <- "weight"
mbl22$Plantation <- factor(mbl22$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(mbl22, aes(x = Plantation, y = weight));c + facet_wrap(~ variable) + geom_bar(stat = "identity")
# a lot of diplopoda in R1Y and not in R10
# coleoptera = oligochaeta ?
# isoptera = hymenoptera
#
## same for functional group
mblfg2=aggregate(mbl[,c(78:82,75)],by=list(state=mbl$state,age=mbl$age),mean)
mblfg2$Plantation=paste(mblfg2$state,mblfg2$age)
mblfg2$Plantation=as.factor(mblfg2$Plantation)
mblfg2$Plantation <- factor(mblfg2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"))
#GRAPH
mblfg3= melt(mblfg2[,c(3:8)], id.var="Plantation") # melt to get the long format
colnames(mblfg3)[3] <- "Biomass"
colnames(mblfg3)[2] <- "Functional_group"
mblfg3$Plantation <- factor(mblfg3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#library(ggplot2)
# with plotted weight
mblfg3$round=round(mblfg3$Biomass,3) # round the Biomass for light plotting
#library(plyr)
mblfg3=ddply(mblfg3, "Plantation", transform, label_y=cumsum(Biomass)) # add the cumsum in order to give the text position
ce=arrange(mblfg3,Plantation,Functional_group) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Biomass)-0.5*Biomass)
x11();ggplot(ce, aes(x = Plantation, y = Biomass, fill = Functional_group)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mblfg33=ddply(mblfg3, "Plantation", transform,percent_Biomass = Biomass / sum(Biomass) * 100) # add % in the tab
mblfg33$Plantation <- factor(mblfg33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mblfg33$roundp=round(mblfg33$percent_Biomass,1) # round the Biomass for light plotting
mblfg33=ddply(mblfg33, "Plantation", transform, label_y=cumsum(percent_Biomass)) # add the cumsum in order to give the text position
ce=arrange(mblfg33,Plantation,Functional_group) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Biomass)-0.5*percent_Biomass)
x11();ggplot(ce, aes(x = Plantation, y = percent_Biomass, fill = Functional_group)) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)

# SOIL DENSITY  ----
## total ----
glmm0 = glmer(total ~ 1+(1|plot/inner_replicate),data=mds,family="poisson");
glmm1 = glmer(total ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson");
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mds$Plantation),letter=grpplot);
x11();ggplot(mds, aes(x=Plantation, y=total))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=100000,label = letter)) +
  xlab("Plantation") +
  ylab("Total density (m3)")+guides(fill=F);
# total - hymenoptera - isoptera
glmm0 = glmer((total-hymenoptera-isoptera) ~ 1+(1|plot/inner_replicate),data=mds,family="poisson");
glmm1 = glmer((total-hymenoptera-isoptera) ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson");
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mds$Plantation),letter=grpplot);
x11();ggplot(mds, aes(x=Plantation, y=(total-hymenoptera-isoptera)))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=15000,label = letter)) +
  xlab("Plantation") +
  ylab("Total density without hymenoptera and isoptera (m3)")+guides(fill=F);
# density increase
# specnumber ----
glmm0 = glmer(specnumber ~ 1+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(specnumber ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # no difference
# shannon ----
mds1=na.exclude(mds)
lme1 = lme(fixed= shannon ~ Plantation,random= ~ 1|plot/inner_replicate,data=mds1)
summary(lme1)
anova(lme1) # difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)) # normal
bartlett.test(residuals(lme1), mds$Plantation);#x11();plot(lme1) # homo
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1 #  NP different de ctrl
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mds$Plantation),letter=grpplot)
x11();ggplot(mds, aes(x=Plantation, y=shannon))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=2,label = letter)) +
  xlab("Plantation") +
  ylab("Shannon index (H')")
# pielou ----
lme1 = lme(fixed= pielou ~ Plantation,random= ~ 1|plot/inner_replicate,data=mds1)
summary(lme1)
anova(lme1)
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)) # normal
bartlett.test(residuals(lme1), mds$Plantation);#x11();plot(lme1) # homo
x11();ggplot(mds, aes(x=Plantation, y=pielou))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=1.1,label = letter)) +
  xlab("Plantation") +
  ylab("Pielou index (J')")
# araneae ----
glmm0 = glmer(araneae ~ 1+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(araneae ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1))
bartlett.test(residuals(glmm1), mds$Plantation);
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1 #  NP different de ctrl
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mds$Plantation),letter=grpplot)
x11();ggplot(mds, aes(x=Plantation, y=araneae))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=1200,label = letter)) +
  xlab("Plantation") +
  ylab("Craneae density (m³)")
# chilopoda ----
glmm0 = glmer(chilopoda ~ 1+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(chilopoda ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1 #  NP different de ctrl
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mds$Plantation),letter=grpplot)
x11();ggplot(mds, aes(x=Plantation, y=chilopoda))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=4100,label = letter)) +
  xlab("Plantation") +
  ylab("Chilopoda density (m³)")
# coleoptera ----
glmm0 = glmer(coleoptera ~ 1+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(coleoptera ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1 #  NP different de ctrl
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mds$Plantation),letter=grpplot)
x11();ggplot(mds, aes(x=Plantation, y=coleoptera))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=2500,label = letter)) +
  xlab("Plantation") +
  ylab("Coleoptera density (m³)")
# hymenoptera ----
glmm0 = glmer(hymenoptera ~ 0+(1|plot/inner_replicate),data=mds,family="poisson");
glmm1 = glmer(hymenoptera ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson");
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mds$Plantation),letter=grpplot);
x11();ggplot(mds, aes(x=Plantation, y=hymenoptera))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=20000,label = letter)) +
  xlab("Plantation") +
  ylab("Hymenoptera density (m3)")+guides(fill=F)+ylim(0,20000)
# isopoda ----
glmm0 = glmer(isopoda ~ 0+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(isopoda ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
lm0=lm(mds2$isopoda~1)
lm1=lm(mds2$isopoda~mds2$Plantation)
anova(lm0,lm1)
shapiro.test(residuals(lm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(lm1), mds2$Plantation);#x11();plot(glmm1) # homo
kruskal.test(mds2$isopoda,mds2$Plantation)
test=pairwise.wilcox.test(mds2$isopoda,mds2$Plantation, p.adj = "bonf");test
test$p.value
# isoptera ----
glmm0 = glmer(isoptera ~ 0+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(isoptera ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
lm0=lm(isoptera~1,data=mds2)
lm1=lm(isoptera~Plantation,data=mds2)
anova(lm0,lm1)
# oligochaeta (sum) ----
glmm0 = glmer(oligochaeta ~ 0+(1|plot/inner_replicate),data=mds,family="poisson");
glmm1 = glmer(oligochaeta ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson");
summary(glmm1)
anova(glmm0,glmm1) # diff
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mds$Plantation),letter=grpplot);
x11();ggplot(mds, aes(x=Plantation, y=oligochaeta))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=5000,label = letter)) +
  xlab("Plantation") +
  ylab("Oligochaeta density (m3)")+guides(fill=F);
# Engineer ----
glmm0 = glmer(engineer ~ 1+(1|plot/inner_replicate),data=mds,family="poisson");
glmm1 = glmer(engineer ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson");
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mds$Plantation),letter=grpplot);
x11();ggplot(mds, aes(x=Plantation, y=engineer))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=90000,label = letter)) +
  xlab("Plantation") +
  ylab("Engineer density (m³)")+guides(fill=F);
# carnivore ----
glmm0 = glmer(carnivore ~ 1+(1|plot/inner_replicate),data=mds,family="poisson");
glmm1 = glmer(carnivore ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson");
summary(glmm1)
anova(glmm0,glmm1) # no difference
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mds$Plantation),letter=grpplot);
x11();ggplot(mds, aes(x=Plantation, y=carnivore))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=5000,label = letter)) +
  xlab("Plantation") +
  ylab("Carnivore density (m³)")+guides(fill=F);
# omnivore ----
glmm0 = glmer(omnivore ~ 1+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(omnivore ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1) # not significative
# herbivore ----
glmm0 = glmer(herbivore ~ 1+(1|plot/inner_replicate),data=mds,family="poisson")
glmm1 = glmer(herbivore ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson")
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # hetero
lm0=lm(mds2$araneae~1)
lm1=lm(mds2$araneae~mds2$Plantation)
anova(lm0,lm1)
# detritivore ----
glmm0 = glmer(detritivore ~ 1+(1|plot/inner_replicate),data=mds,family="poisson");
glmm1 = glmer(detritivore ~ Plantation+(1|plot/inner_replicate),data=mds,family="poisson");
summary(glmm1)
anova(glmm0,glmm1)
shapiro.test(residuals(glmm1));#x11();qqnorm(residuals(glmm1));qqline(residuals(glmm1)) # normal
bartlett.test(residuals(glmm1), mds$Plantation);#x11();plot(glmm1) # homo
summary(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(glmm1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mds$Plantation),letter=grpplot);
x11();ggplot(mds, aes(x=Plantation, y=detritivore))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=5200,label = letter)) +
  xlab("Plantation") +
  ylab("Detritivore density (m³)")+guides(fill=F);
# graph ----
mds2=aggregate(mds,by=list(state=mds$state,age=mds$age),mean);
mds2$Plantation=paste(mds2$state,mds2$age);
mds2$Plantation=as.factor(mds2$Plantation);
mds2$Plantation <- factor(mds2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));

mds3= melt(mds2[,c(8:21,77,78)], id.var="Plantation") # melt to get the long format
colnames(mds3)[3] <- "Density"
colnames(mds3)[2] <- "Taxon"
mds3$Plantation <- factor(mds3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#library(ggplot2)
# with plotted weight
mds3$round=round(mds3$Density,1) # round the Density for light plotting
#library(plyr)
mds3=ddply(mds3, "Plantation", transform, label_y=cumsum(Density)) # add the cumsum in order to give the text position
ce=arrange(mds3,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Density)-0.5*Density)
x11();ggplot(ce, aes(x = Plantation, y = Density, fill = Taxon)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mds33=ddply(mds3, "Plantation", transform,percent_Density = Density / sum(Density) * 100) # add % in the tab
mds33$Plantation <- factor(mds33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mds33$roundp=round(mds33$percent_Density,1) # round the Density for light plotting
mds33=ddply(mds33, "Plantation", transform, label_y=cumsum(percent_Density)) # add the cumsum in order to give the text position
ce=arrange(mds33,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Density)-0.5*percent_Density)
x11();ggplot(ce, aes(x = Plantation, y = percent_Density, fill = Taxon)) + geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)
# graph by compound
mds4= melt(mds2[,c(8:21,77,78,79)], id.var="Plantation") # melt to get the long format
colnames(mds4)[3] <- "Density"
colnames(mds4)[2] <- "Taxon"
mds4$Plantation <- factor(mds4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(mds4, aes(x = Plantation, y = Density));c + facet_wrap(~ Taxon,scales="free_y") + geom_bar(stat = "identity")
#
## same for functional group
#GRAPH
mds3= melt(mds2[,c(83:87,77)], id.var="Plantation");
colnames(mds3)[3] <- "Density";
colnames(mds3)[2] <- "Functional_group";
mds3$Plantation <- factor(mds3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
#library(ggplot2)
# with plotted weight
mds3$round=round(mds3$Density,1) # round the Density for light plotting
#library(plyr)
mds3=ddply(mds3, "Plantation", transform, label_y=cumsum(Density)) # add the cumsum in order to give the text position
ce=arrange(mds3,Plantation,Functional_group) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Density)-0.5*Density)
x11();ggplot(ce, aes(x = Plantation, y = Density, fill = Functional_group)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mds33=ddply(mds3, "Plantation", transform,percent_Density = Density / sum(Density) * 100);
mds33$Plantation <- factor(mds33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
mds33$roundp=round(mds33$percent_Density,1);
mds33=ddply(mds33, "Plantation", transform, label_y=cumsum(percent_Density));
ce=arrange(mds33,Plantation,Functional_group);
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Density)-0.5*percent_Density);
x11();ggplot(ce, aes(x = Plantation, y = percent_Density, fill = Functional_group)) + geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3);
# by compound
mds5= melt(mds2[,c(83:87,79,77)], id.var="Plantation") # melt to get the long format
colnames(mds5)[3] <- "Density"
colnames(mds5)[2] <- "Taxon"
mds5$Plantation <- factor(mds5$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(mds5, aes(x = Plantation, y = Density));c + facet_wrap(~ Taxon,scales="free_y") + geom_bar(stat = "identity")
## SOIL BIOMASS ----
## total ----
lme1 = lme(fixed= biomass ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal QQplot
bartlett.test(residuals(lme1), mbs$Plantation);#x11();plot(lme1); 
lme1 = lme(fixed= log(biomass) ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal QQplot
bartlett.test(residuals(lme1), mbs$Plantation);#x11();plot(lme1); 
lm0=lm(biomass~1,data=mbs2)
lm1=lm(biomass~Plantation,data=mbs2)
anova(lm0,lm1)
## coleoptera ----
lme1 = lme(fixed= coleoptera ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
## hymenoptera ----
lme1 = lme(fixed= hymenoptera ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
## isoptera ----
lme1 = lme(fixed= isoptera ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
## oligochaeta (sum)
lme1 = lme(fixed= oligochaeta ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mbs$Plantation);#x11();plot(lme1); 
lme1 = lme(fixed= log(oligochaeta+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal QQplot
bartlett.test(residuals(lme1), mbs$Plantation);#x11();plot(lme1); 
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1 #  NP different de ctrl
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mds$Plantation),letter=grpplot)
x11();ggplot(mds, aes(x=Plantation, y=oligochaeta))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=5000,label = letter)) +
  xlab("Plantation") +
  ylab("Soil oligochaeta biomass (m³)")
## carnivore ----
lme1 = lme(fixed= carnivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
## herbivore ----
lme1 = lme(fixed= herbivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
## engineer ----
lme1 = lme(fixed= engineer ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal QQplot
bartlett.test(residuals(lme1), mbs$Plantation);#x11();plot(lme1); 
lme1 = lme(fixed= log(engineer+1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal QQplot
bartlett.test(residuals(lme1), mbs$Plantation);#x11();plot(lme1); 
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'))
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1 #  NP different de ctrl
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(mds$Plantation),letter=grpplot)
x11();ggplot(mds, aes(x=Plantation, y=engineer))  + 
  geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=90000,label = letter)) +
  xlab("Plantation") +
  ylab("Engineer biomass (m³)")
## omnivore ----
lme1 = lme(fixed= omnivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
## detritivore ----
lme1 = lme(fixed= detritivore ~ Plantation,random= ~ 1|plot/inner_replicate,data=mbs);
summary(lme1);
anova(lme1); # no difference
## graph ----
mbs2=aggregate(mbs,by=list(state=mbs$state,age=mbs$age),mean)
mbs2$Plantation=paste(mbs2$state,mbs2$age)
mbs2$Plantation=as.factor(mbs2$Plantation)
mbs2$Plantation <- factor(mbs2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"))
mbs3= melt(mbs2[,c(8:21,77,78)], id.var="Plantation") # melt to get the long format
colnames(mbs3)[3] <- "Biomass"
colnames(mbs3)[2] <- "Taxon"
mbs3$Plantation <- factor(mbs3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#library(ggplot2)
# with plotted weight
mbs3$round=round(mbs3$Biomass,2) # round the Biomass for light plotting
#library(plyr)
mbs3=ddply(mbs3, "Plantation", transform, label_y=cumsum(Biomass)) # add the cumsum in order to give the text position
ce=arrange(mbs3,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Biomass)-0.5*Biomass)
x11();ggplot(ce, aes(x = Plantation, y = Biomass, fill = Taxon)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mbs33=ddply(mbs3, "Plantation", transform,percent_Biomass = Biomass / sum(Biomass) * 100) # add % in the tab
mbs33$Plantation <- factor(mbs33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mbs33$roundp=round(mbs33$percent_Biomass,1) # round the Biomass for light plotting
mbs33=ddply(mbs33, "Plantation", transform, label_y=cumsum(percent_Biomass)) # add the cumsum in order to give the text position
ce=arrange(mbs33,Plantation,Taxon) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Biomass)-0.5*percent_Biomass)
x11();ggplot(ce, aes(x = Plantation, y = percent_Biomass, fill = Taxon)) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)
# graph by compound
mbs22=melt(mbs2[,c(7:21,77,78)], id.var="Plantation")
colnames(mbs22)[3] <- "weight"
mbs22$Plantation <- factor(mbs22$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(mbs22, aes(x = Plantation, y = weight));c + facet_wrap(~ variable) + geom_bar(stat = "identity")
# a lot of diplopoda in R1Y and not in R10
# coleoptera = oligochaeta ?
# isoptera = hymenoptera
#
## same for functional group
mbsfg2=aggregate(mbs[,c(78:82,75)],by=list(state=mbs$state,age=mbs$age),mean)
mbsfg2$Plantation=paste(mbsfg2$state,mbsfg2$age)
mbsfg2$Plantation=as.factor(mbsfg2$Plantation)
mbsfg2$Plantation <- factor(mbsfg2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"))
#GRAPH
mbsfg3= melt(mbsfg2[,c(3:8)], id.var="Plantation") # melt to get the long format
colnames(mbsfg3)[3] <- "Biomass"
colnames(mbsfg3)[2] <- "Functional_group"
mbsfg3$Plantation <- factor(mbsfg3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
#library(ggplot2)
# with plotted weight
mbsfg3$round=round(mbsfg3$Biomass,3) # round the Biomass for light plotting
#library(plyr)
mbsfg3=ddply(mbsfg3, "Plantation", transform, label_y=cumsum(Biomass)) # add the cumsum in order to give the text position
ce=arrange(mbsfg3,Plantation,Functional_group) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Biomass)-0.5*Biomass)
x11();ggplot(ce, aes(x = Plantation, y = Biomass, fill = Functional_group)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)
# by relative composition
mbsfg33=ddply(mbsfg3, "Plantation", transform,percent_Biomass = Biomass / sum(Biomass) * 100) # add % in the tab
mbsfg33$Plantation <- factor(mbsfg33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mbsfg33$roundp=round(mbsfg33$percent_Biomass,1) # round the Biomass for light plotting
mbsfg33=ddply(mbsfg33, "Plantation", transform, label_y=cumsum(percent_Biomass)) # add the cumsum in order to give the text position
ce=arrange(mbsfg33,Plantation,Functional_group) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Biomass)-0.5*percent_Biomass)
x11();ggplot(ce, aes(x = Plantation, y = percent_Biomass, fill = Functional_group)) + guides(fill=guide_legend(reverse=TRUE))+ geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)
## adjusted p.val  ----
# bonferroni simple sequentialy rejective test
pval=data.frame(test=c("Soil Shannon index","Soil Pielou index","Litter detritivore density","Litter coleoptera density","Litter total density","Litter engineer density","Soil carnivore density","Soil chilopoda density","Soil hymenoptera density","Soil oligochaeta density","Soil engineer density","Soil detritivore density","Soil total density"),type=c("Linear mixed model","Linear mixed model","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function","General linear mixed model, logarithm as link function"),p.val=c(0.0222,0.0042,0.001104,1.168*10^-5,2.258*10^-9,1.769*10^-9,1.904*10^-14,9.105*10^-15,2.2*10^-16,2.2*10^-16,2.2*10^-16,2.2*10^-16,2.2*10^-16),factor=c(1:13))
pval$adjusted.p.val=pval$p.val*pval$factor;pval
## Multivariate analysis ----
# taxon (not significant) ----
## LIT + SOIL
library(ade4)
# CA
coa2<-dudi.coa(mdtot[,c(6:19,75)],scannf=F,nf=5)
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red")

pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(mdtot$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,mdtot$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));par(mfrow = c(1, 2));s.class(bca2$ls,mdtot$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T)
title(xlab="Dim 1 (40%) ",ylab="Dim 2 (33%)")
s.corcircle(bca2$co,box=T,clabel=0.5,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference
x11();plot(pt2,main="Between class inertia")

# SOIL
# CA
coamsoil<-dudi.coa(mds[,c(6:19,76)],scannf=F,nf=4)
# Kaiser-Guttman
barplot(coamsoil$eig, main="Eigenvalues", col="grey")
abline(h=mean(coamsoil$eig), col="red")

pvp=100*coamsoil$eig/sum(coamsoil$eig);pvp;
cumsum(pvp)
#x11();s.class(coamsoil$li,as.factor(mds$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F)#;s.arrow(coamsoil$co,boxes=F,add.plot=T)
# BCA
bcamsoil<-bca(coamsoil,mds$Plantation,scannf=FALSE)
pvp=100*bcamsoil$eig/sum(bcamsoil$eig);pvp;
cumsum(pvp)
#x11();plot(bcamsoil)
x11();par(mfrow = c(1, 2));s.class(bcamsoil$ls,mds$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bcamsoil$co,boxes=F,add.plot=T)
title(xlab="Dim 1 (38%) ",ylab="Dim 2 (35%)")
s.corcircle(bcamsoil$co,clabel=0.5,box=T)
# Permutation test
ptmsoil=rtest(bcamsoil,999)
ptmsoil # non significative difference
x11();plot(ptmsoil,main="Between class inertia")

## LITTER
# CA
coa2<-dudi.coa(mdl[,c(6:19,76)],scannf=F,nf=4)
# Kaiser-Guttman
barplot(coamsoil$eig, main="Eigenvalues", col="grey")
abline(h=mean(coamsoil$eig), col="red")

pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(mdl$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F)#;s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,mdl$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,mdl$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=F,add.plot=T)
title(xlab="Dim 1 (49%) ",ylab="Dim 2 (30%)")
s.corcircle(bca2$co,clabel=0.5,box=T)
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference
x11();plot(pt2,main="Between class inertia")
# functional ----
## LIT + SOIL
library(ade4)
# CA
coa2<-dudi.coa(mdtot[,c(78:82)],scannf=F,nf=2)
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red")

pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(mdtot$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,mdtot$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));par(mfrow = c(1, 2));s.class(bca2$ls,mdtot$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T)
title(xlab="Dim 1 (50%) ",ylab="Dim 2 (28%)",main="p-value: 0.001",font.main=1)
s.corcircle(bca2$co,box=T,clabel=0.5,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # significative difference
x11();plot(pt2,main="Between class inertia")

# SOIL
# CA
coamsoil<-dudi.coa(mds[,c(81:85)],scannf=F,nf=2);
# Kaiser-Guttman
barplot(coamsoil$eig, main="Eigenvalues", col="grey");
abline(h=mean(coamsoil$eig), col="red");
pvp=100*coamsoil$eig/sum(coamsoil$eig);pvp;
cumsum(pvp);
#x11();s.class(coamsoil$li,as.factor(mds$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F)#;s.arrow(coamsoil$co,boxes=F,add.plot=T)
# BCA
bcamsoil<-bca(coamsoil,mds$Plantation,scannf=FALSE);
pvp=100*bcamsoil$eig/sum(bcamsoil$eig);pvp;
cumsum(pvp);
#x11();plot(bcamsoil)
x11();par(mfrow = c(1, 2));s.class(bcamsoil$ls,mds$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);
title(xlab="Dim 1 (49%) ",ylab="Dim 2 (29%)",main="p-value: 0.001",font.main=1);
s.corcircle(bcamsoil$co,clabel=1,box=T,full=F);
# Permutation test
ptmsoil=rtest(bcamsoil,999)
ptmsoil # non significative difference
x11();plot(ptmsoil,main="Between class inertia")
# SOIL II
# CA
a=mds[mds$age=="o",c(81:85,75)]
a$Plantation <- factor(a$Plantation, levels=c("N o","R1 o","R2 o"));
coamsoil<-dudi.coa(a[,c(1:5)],scannf=F,nf=2);
# Kaiser-Guttman
barplot(coamsoil$eig, main="Eigenvalues", col="grey");
abline(h=mean(coamsoil$eig), col="red");
pvp=100*coamsoil$eig/sum(coamsoil$eig);pvp;
cumsum(pvp);
#x11();s.class(coamsoil$li,as.factor(mds$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F)#;s.arrow(coamsoil$co,boxes=F,add.plot=T)
# BCA
bcamsoil<-bca(coamsoil,a$Plantation,scannf=FALSE);
pvp=100*bcamsoil$eig/sum(bcamsoil$eig);pvp;
cumsum(pvp);
#x11();plot(bcamsoil)
x11();par(mfrow = c(1, 2));s.class(bcamsoil$ls,a$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);
title(xlab="Dim 1 (99%) ",ylab="Dim 2 (1%)",main="p-value: 0.009",font.main=1);
s.corcircle(bcamsoil$co,clabel=1,box=T,full=T);
# Permutation test
ptmsoil=rtest(bcamsoil,999)
ptmsoil # non significative difference
x11();plot(ptmsoil,main="Between class inertia")

## LITTER
# CA
coa2<-dudi.coa(mdl[,c(81:85)],scannf=F,nf=2)
# Kaiser-Guttman
barplot(coamsoil$eig, main="Eigenvalues", col="grey")
abline(h=mean(coamsoil$eig), col="red")

pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(mdl$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F)#;s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,mdl$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,mdl$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=F,add.plot=T)
title(xlab="Dim 1 (50%) ",ylab="Dim 2 (42%)",main="p-value: 0.307",font.main=1)
s.corcircle(bca2$co,clabel=0.5,box=T)
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference
x11();plot(pt2,main="Between class inertia")
## indicators ----
library(labdsv)
#tot
iva=indval(mdtot[,6:74],as.numeric(mdtot[,5]))
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(mdtot[,6:74]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
indvalsummary
# soil
mdsind=mds[,c(6:19,76)] # select only the species
mdsind=mdsind[, colSums(abs(mdsind)) != 0] # remove the colums without values
iva=indval(mdsind,as.numeric(mds[,75])) # dont forget to set Plantationment as numeric
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(mdsind>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
indvalsummary
# litter
mdlind=mdl[,c(6:19,76)]
mdlind=mdlind[, colSums(abs(mdlind)) != 0]
iva=indval(mdlind,as.numeric(mdl[,75]))
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(mdlind>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
indvalsummary

### NEMATODES ----
# dataset ----
library(lsmeans);
# density
nemd <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/nematode_densite_0.csv", dec=",");
str(nemd);
nemd$plot=as.factor(nemd$plot);
nemd$age <- factor(nemd$age, levels=c("y","o"));
nemd$Plantation=paste(nemd$state,nemd$age);
nemd$Plantation=as.factor(nemd$Plantation);
nemd$Plantation <- factor(nemd$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));

nemt <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/nematode_taxon_0.csv", dec=",");
str(nemt);
nemt$plot=as.factor(nemt$plot);
nemt$age <- factor(nemt$age, levels=c("y","o"));
nemt$Plantation=paste(nemt$state,nemt$age);
nemt$Plantation=as.factor(nemt$Plantation);
nemt$Plantation <- factor(nemt$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
nemt$shannon=diversity(nemt[,6:59]);
nemt$specnumber=specnumber(nemt[,6:59]);
nemt$pielou = nemt$shannon/log(nemt$specnumber);
nemt$total=apply(nemt[,6:59],1,sum);

nemf <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/nematode_famille_0.csv",dec=",");
str(nemf);
nemf$plot=as.factor(nemf$plot);
nemf$age <- factor(nemf$age, levels=c("y","o"));
nemf$Plantation=paste(nemf$state,nemf$age);
nemf$Plantation=as.factor(nemf$Plantation);
nemf$Plantation <- factor(nemf$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
nemf$shannon=diversity(nemf[,6:46]);
nemf$specnumber=specnumber(nemf[,6:46]);
nemf$pielou = nemf$shannon/log(nemf$specnumber);

nemi <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/nematode_indice_0.csv", dec=",");
str(nemi);
  nemi$plot=as.factor(nemi$plot);
  nemi$age <- factor(nemi$age, levels=c("y","o"));
  nemi$Plantation=paste(nemi$state,nemi$age);
  nemi$Plantation=as.factor(nemi$Plantation);
  nemi$Plantation <- factor(nemi$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
  nemi$BACTERIVORES=nemi$BACTERIVORES_cp1+nemi$BACTERIVORES_AUTRES;
# density ----
#nemd=nemd[nemd$age=="o",]
lm1=lm(number_.of_nematod_per_g~Plantation,data=nemd);
lm0=lm(number_.of_nematod_per_g~1,data=nemd);
anova(lm0,lm1)
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), nemd$Plantation);

lsmeans(lm1,pairwise~Plantation,data=nemd,adjust='bonferroni')
lettre=cld(summary(glht(lm1,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")))
lettre$mcletters
#let=scan(what="")
let=c("a","b","ab","ab","ab")
xx=data.frame(level=levels(nemd$Plantation),letter=let)
x11();ggplot(nemd, aes(x=Plantation, y=number_.of_nematod_per_g))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=4,label = letter)) +
  xlab("Plantation") +
  ylab("Nematode density per gram")
## taxon ----

# richness
lm1=lm(specnumber~Plantation,data=nemt)
lm0=lm(specnumber~1,data=nemt)
anova(lm0,lm1) # difference
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), nemt$Plantation);
lsmeans(lm1,pairwise~Plantation,data=nemt,adjust='bonferroni')
lettre=cld(summary(glht(lm1,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")))
lettre$mcletters
let=c("b","ab","b","a","b")
xx=data.frame(level=levels(nemt$Plantation),letter=let)
x11();ggplot(nemt, aes(x=Plantation, y=specnumber))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=30,label = letter)) +
  xlab("Plantation") +
  ylab("Nematode taxon richness")
# less sp in young
# shannon
lm1=lm(shannon~Plantation,data=nemt)
lm0=lm(shannon~1,data=nemt)
anova(lm0,lm1) # no difference
x11();ggplot(nemt, aes(x=Plantation, y=shannon))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Shannon index (H')")
mean(nemt$shannon)
# pielou
lm1=lm(pielou~Plantation,data=nemt)
lm0=lm(pielou~1,data=nemt)
anova(lm0,lm1) # no difference
x11();ggplot(nemt, aes(x=Plantation, y=pielou))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Pielou index (J')")
mean(nemt$pielou)
## famille ----

# richness
lm1=lm(specnumber~Plantation,data=nemf)
lm0=lm(specnumber~1,data=nemf)
anova(lm0,lm1) # difference
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), nemt$Plantation);
lsmeans(lm1,pairwise~Plantation,data=nemt,adjust='bonferroni')
lettre=cld(summary(glht(lm1,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")))
lettre$mcletters
let=c("a","ab","b","a","b")
xx=data.frame(level=levels(nemd$Plantation),letter=let)
x11();ggplot(nemf, aes(x=Plantation, y=specnumber))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=30,label = letter)) +
  xlab("Plantation") +
  ylab("Nematode family richness")
# shannon
lm1=lm(shannon~Plantation,data=nemf)
lm0=lm(shannon~1,data=nemf)
anova(lm0,lm1) # no difference
x11();ggplot(nemf, aes(x=Plantation, y=shannon))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Shannon index (H')")
# pielou
lm1=lm(pielou~Plantation,data=nemf)
lm0=lm(pielou~1,data=nemf)
anova(lm0,lm1) # no difference
x11();ggplot(nemf, aes(x=Plantation, y=pielou))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Pielou index (J')")
# trophic ----
# bacterivore cp1
lm0=lm(BACTERIVORES_cp1~1,data=nemi);
lm1=lm(BACTERIVORES_cp1~Plantation,data=nemi);
anova(lm0,lm1);
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), nemi$Plantation);
lm0=lm(log(BACTERIVORES_cp1)~1,data=nemi);
lm1=lm(log(BACTERIVORES_cp1)~Plantation,data=nemi);
anova(lm0,lm1);
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), nemi$Plantation);
lsmeans(lm1,pairwise~Plantation,data=nemi,adjust='bonferroni');
lettre=cld(summary(glht(lm1,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")));
lettre$mcletters;
let=c("a","ab","ab","ab","b");
xx=data.frame(level=levels(nemi$Plantation),letter=let);
x11();ggplot(nemi, aes(x=Plantation, y=BACTERIVORES_cp1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=25,label = letter)) +
  xlab("Plantation") +
  ylab("Cp1 bacterivores density per 100g");
# increase in Old plantation: they take advantage on the other grp (cp2) (competitive exclusion ?)
# intensive tillage or fertilisation can increase cp1, but is there this practice in old plantation ?
# + no general increase of all the bacterivore
# bacterivore autres
lm0=lm(BACTERIVORES_AUTRES~1,data=nemi)
lm1=lm(BACTERIVORES_AUTRES~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=BACTERIVORES_AUTRES))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Other bacterivores density per 100g")
# bacterivore total
lm0=lm(BACTERIVORES~1,data=nemi)
lm1=lm(BACTERIVORES~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=BACTERIVORES))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Total bacterivores density per 100g")
# Fungivore
lm0=lm(FONGIVORES~1,data=nemi)
lm1=lm(FONGIVORES~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=FONGIVORES))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Fungivore density per 100g")
# carnivore
lm0=lm(CARNIVORES~1,data=nemi)
lm1=lm(CARNIVORES~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=CARNIVORES))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Carnivore density per 100g")
# omnivore
lm0=lm(OMNIVORES~1,data=nemi)
lm1=lm(OMNIVORES~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=OMNIVORES))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Omnivore density per 100g")
# PHYTOPHAGES
lm0=lm(PHYTOPHAGES~1,data=nemi)
lm1=lm(PHYTOPHAGES~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=PHYTOPHAGES))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Phytophage density per 100g")
# indices ----
# EI opportunist reactivity (food web response to ressource availability)
lm0=lm(EI~1,data=nemi)
lm1=lm(EI~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=EI))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Enrichment Index (EI)")
# ressource availability increase
# SI (more complex food web structure (K stategy))
lm0=lm(SI~1,data=nemi)
lm1=lm(SI~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=SI))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Structural Index (SI)")
# BI
lm0=lm(BI~1,data=nemi)
lm1=lm(BI~Plantation,data=nemi)
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=BI))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Basal Index (BI)")
# MI (maturity index)
lm0=lm(MI~1,data=nemi);
lm1=lm(MI~Plantation,data=nemi);
anova(lm0,lm1);
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), nemi$Plantation);
lsmeans(lm1,pairwise~Plantation,data=nemi,adjust='bonferroni');
lettre=cld(summary(glht(lm1,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")));
lettre$mcletters;
x11();ggplot(nemi, aes(x=Plantation, y=MI,fill=age))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Maturity Index (MI)");
# increasing disturbance (less high cp)
# importante decrease during the y to o phase ?? (it should be the opposite, as the)
# PPI (plant-parasite index)
lm0=lm(PPI~1,data=nemi);
lm1=lm(PPI~Plantation,data=nemi);
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=PPI,fill=age))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Plant-Parasite Index (PPI)");
# increase
# NCR
lm0=lm(NCR~1,data=nemi)
lm1=lm(NCR~Plantation,data=nemi);
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=NCR))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Nematode Channel Ratio (NCR)");
# NCR.
lm0=lm(NCR.~1,data=nemi)
lm1=lm(NCR.~Plantation,data=nemi);
anova(lm0,lm1)
x11();ggplot(nemi, aes(x=Plantation, y=NCR.))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("NCR*");

palette("default")
x11();plot(nemi$EI~nemi$SI,col=as.factor(nemi$Plantation),pch=19,xlim=c(0,100),ylim=c(0,100))
abline(h=50,v=50)
legend("topleft", legend=levels(factor(nemi$Plantation)), text.col=seq_along(levels(factor(nemi$Plantation))))

#p1 <- ggplot(data=nemi, aes(x=SI, y=EI,color=Plantation))+  geom_point()+stat_ellipse(type="euclid");p1
# graph ----
nemi2=aggregate(nemi,by=list(Plantation=nemi$Plantation),mean);
nemi2=nemi2[,-c(2:6,27)];
nemi2$Plantation=as.factor(nemi2$Plantation);
nemi2$Plantation <- factor(nemi2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
nemi3= melt(nemi2[,c(1,4:6,10,22)], id.var="Plantation");
colnames(nemi3)[3] <- "Density";
colnames(nemi3)[2] <- "Trophic_group";
nemi3$Plantation <- factor(nemi3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
#library(ggplot2)
# with plotted weight
nemi3$round=round(nemi3$Density,2);
#library(plyr)
nemi3=ddply(nemi3, "Plantation", transform, label_y=cumsum(Density));
ce=arrange(nemi3,Plantation,Trophic_group);
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Density)-0.5*Density);
x11();ggplot(ce, aes(x = Plantation, y = Density, fill = Trophic_group)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3);
# by relative composition
nemi33=ddply(nemi3, "Plantation", transform,percent_Density = Density / sum(Density) * 100);
nemi33$Plantation <- factor(nemi33$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
nemi33$roundp=round(nemi33$percent_Density,1);
nemi33=ddply(nemi33, "Plantation", transform, label_y=cumsum(percent_Density));
ce=arrange(nemi33,Plantation,Trophic_group);
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Density)-0.5*percent_Density);
x11();ggplot(ce, aes(x = Plantation, y = percent_Density, fill = Trophic_group)) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3);
# graph by compound
nemi22=melt(nemi2[,c(1,4:6,10,22,12)], id.var="Plantation")
colnames(nemi22)[3] <- "Density"
nemi22$Plantation <- factor(nemi22$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
x11();c <- ggplot(nemi22, aes(x = Plantation, y = Density));c + facet_wrap(~ variable) + geom_bar(stat = "identity")
# multivariate ----
library(ade4)
# CA TAXON
# CA
coa2<-dudi.coa(nemt[,c(6:59)],scannf = F, nf = 5)
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red")
pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
x11();par(mfrow = c(1, 2));s.class(coa2$li,as.factor(nemt$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F)
title(xlab="Dim 1 (33%) ",ylab="Dim 2 (20%)")
s.corcircle(coa2$co,box=T,clabel=0.5,full=F)
# BCA
bca2<-bca(coa2,nemt$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,nemt$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T)
title(xlab="Dim 1 (59%) ",ylab="Dim 2 (22%)",main="p-value: 0.079",font.main=1)
s.corcircle(bca2$co,box=T,clabel=0.5,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference

# CA FAMILLE
# CA
coa2<-dudi.coa(nemf[,c(6:46)],scannf = F, nf = 4);
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red")
pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
x11();par(mfrow = c(1, 2));s.class(coa2$li,as.factor(nemt$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F)
title(xlab="Dim 1 (33%) ",ylab="Dim 2 (20%)")
s.corcircle(coa2$co,box=T,clabel=0.5,full=F)
# BCA
bca2<-bca(coa2,nemt$Plantation,scannf=FALSE);
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,nemt$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T);
title(xlab="Dim 1 (61%) ",ylab="Dim 2 (23%)",main="p-value: 0.068",font.main=1);
s.corcircle(bca2$co,box=T,clabel=0.5,full=F);
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference

# CA GRP TROPHIQUE
coa2<-dudi.coa(nemi[,c(6:10,12,13)],scannf = F, nf = 2)
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red")
pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(nemi$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,nemi$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));par(mfrow = c(1, 2));s.class(bca2$ls,nemi$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T)
title(xlab="Dim 1 (63%) ",ylab="Dim 2 (30%)",main="p-value: 0.05",font.main=1)
s.corcircle(bca2$co,box=T,clabel=0.5,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference
x11();plot(pt2,main="Between class inertia")

# PCA INDICE (relevant for indices)
pca2<-dudi.pca(nemi[,c(17:23)],scannf = F, nf = 2)
barplot(pca2$eig, main="Eigenvalues", col="grey")
abline(h=1, col="red")
pvp=100*pca2$eig/sum(pca2$eig);pvp;
cumsum(pvp)
#x11();s.class(pca2$li,as.factor(nemi$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);s.arrow(pca2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(pca2,nemi$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));par(mfrow = c(1, 2));s.class(bca2$ls,nemi$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T)
title(xlab="Dim 1 (55%) ",ylab="Dim 2 (37%)",main="p-value: 0.013",font.main=1)
s.corcircle(bca2$co,box=T,clabel=0.5)
# Permutation test
pt2=rtest(bca2,999)
pt2 # significative difference
x11();plot(pt2,main="Between class inertia")

### MICROORGANISM ----
# required dataset ----
mic <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/microresp_0.csv", dec=",");
mic$plot=as.factor(mic$plot);
mic$inner_replicate=as.factor(mic$inner_replicate);
mic$inner_replicate_2=as.factor(mic$inner_replicate_2);
mic$age <- factor(mic$age, levels=c("y","o"));
mic$Plantation=paste(mic$state,mic$age);
mic$Plantation=as.factor(mic$Plantation);
mic$Plantation <- factor(mic$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mic$bronopol=as.factor(mic$bronopol);
mic[,c(7:14,16:22)]=mic[,c(7:14,16:22)]-mic$H2O;
mic=mic[,-c(20)];
mic[,c(7:21)][mic[,c(7:21)] < 0] <- 0
mic$qCO2=mic$H2O/mic$Glucose;
mic[which(mic$qCO2=="Inf"),]$qCO2=mic[which(mic$qCO2=="Inf"),15]
mic0=mic;
mict=mic[mic$bronopol==0,];
micf=mic[mic$bronopol==1,];
micb=mict;
micb[,c(7:21)]=mict[,c(7:21)]-micf[,c(7:21)];
micb[,c(7:21)][micb[,c(7:21)] < 0] <- 0
mict$specnumber=specnumber(mict[,7:21]);
mict$shannon=diversity(mict[,7:21]);
mict$pielou = mict$shannon/log(mict$specnumber);
mict[which(mict$shannon=="0"),]$pielou=1
micf$specnumber=specnumber(micf[,7:21]);
micf$shannon=diversity(micf[,7:21]);
micf$pielou = micf$shannon/log(micf$specnumber);
micb$specnumber=specnumber(micb[,7:21]);
micb$shannon=diversity(micb[,7:21]);
micb$pielou = micb$shannon/log(micb$specnumber);

mictk=aggregate(mict,by=list(state=mict$state,age=mict$age,plot=mict$plot),mean);
mictk=mictk[,-c(4:9,25)];
mictk$Plantation=paste(mictk$state,mictk$age);
mictk$Plantation=as.factor(mictk$Plantation);
mictk$Plantation <- factor(mictk$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
micbk=aggregate(micb,by=list(state=micb$state,age=micb$age,plot=micb$plot),mean);
micbk=micbk[,-c(4:9,25)];
micbk$Plantation=paste(micbk$state,micbk$age);
micbk$Plantation=as.factor(micbk$Plantation);
micbk$Plantation <- factor(micbk$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
micfk=aggregate(micf,by=list(state=micf$state,age=micf$age,plot=micf$plot),mean);
micfk=micfk[,-c(4:9,25)];
micfk$Plantation=paste(micfk$state,micfk$age);
micfk$Plantation=as.factor(micfk$Plantation);
micfk$Plantation <- factor(micfk$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mict$SIR=apply(mict[,c(7:14,16:21)],1,mean);
micb$SIR=apply(micb[,c(7:14,16:21)],1,mean);
micf$SIR=apply(micf[,c(7:14,16:21)],1,mean);
mictk$SIR=apply(mictk[,c(4:11,13:18)],1,mean);
micfk$SIR=apply(micfk[,c(4:11,13:18)],1,mean);
# prb negative value for bacteria alone
# br ----
#total 
lme1 = lme(fixed= H2O ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);
lme1 = lme(fixed= log(H2O) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); # OK qqplot
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(mict$Plantation),letter=grpplot);
x11();ggplot(mict, aes(x=Plantation, y=H2O))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=10,label = letter)) +
  xlab("Plantation") +
  ylab("Total Basal Respiration (BR)");
# decrease
# bact
lme1 = lme(fixed= H2O ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micb);
summary(lme1);
anova(lme1);
lme1 = lme(fixed= log(H2O+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micb);
anova(lme1);
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); # OK qqplot
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
lm0=lm(H2O~1,data=micbk)
lm1=lm(H2O~Plantation,data=micbk)
anova(lm0,lm1)
#kruskal.test(micbk$H2O~micbk$Plantation)
# fungi
lme1 = lme(fixed= H2O ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(H2O) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);x11();plot(lme1);
a=lm(H2O~1,data=micfk)
b=lm(H2O~Plantation,data=micfk)
anova(a,b); # dif
x11();ggplot(micf, aes(x=Plantation, y=H2O))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Fungal Basal Respiration (BR)")
# BR ratio fungi/(bact+fungi)
mict$BRr=micf$H2O/(mict$H2O)
lme1 = lme(fixed= BRr ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(micb$Plantation),letter=grpplot);
x11();ggplot(micb, aes(x=Plantation, y=BRr))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=1.25,label = letter)) +
  xlab("Plantation") +
  ylab("Total Basal Respiration ratio (BRr)");
# gir ----
# total
lme1 = lme(fixed= Glucose ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(Glucose+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
lm0=lm(mictk$Glucose~1)
lm1=lm(mictk$Glucose~mictk$Plantation)
anova(lm0,lm1) # no diff

p2=ggplot(mict, aes(x=Plantation, y=Glucose))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Total Glucose Induced Respiration (GIR)");
# bact
lme1 = lme(fixed= Glucose ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micb);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micb$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(Glucose+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micb);
summary(lme1);
anova(lme1); # no dif, tranfo to heavy
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
lm0=lm(micbk$Glucose~1)
lm1=lm(micbk$Glucose~micbk$Plantation)
anova(lm0,lm1) # no diff
# fungi
lme1 = lme(fixed= Glucose ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(Glucose+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
a=lm(log(Glucose)~1,data=micfk)
b=lm(log(Glucose)~Plantation,data=micfk)
anova(a,b)
x11();ggplot(micf, aes(x=Plantation, y=Glucose))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Fungi Glucose Induced Respiration (GIR)")

# total
lme1 = lme(fixed= SIR ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(SIR+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
a=lm(SIR~1,data=mict)
b=lm(SIR~Plantation,data=mict)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), mict$Plantation);
a=lm(log(SIR+1)~1,data=mict)
b=lm(log(SIR+1)~Plantation,data=mict)
anova(a,b)
shapiro.test(residuals(a));x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(b), mict$Plantation);

lettre=cld(summary(glht(b,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")))
lettre$mcletters
let=c("ab","a","b","a","ab")
xx=data.frame(level=levels(mict$Plantation),letter=let)
x11();ggplot(mict, aes(x=Plantation, y=SIR))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=40,label = letter)) +
  xlab("Plantation") +
  ylab("Total Substrat Induced Respiration (SIR)");
# fungi
lme1 = lme(fixed= SIR ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(SIR) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
a=lm(log(SIR)~1,data=micf)
b=lm(log(SIR)~Plantation,data=micf)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), mict$Plantation);
kruskal.test(micfk$SIR,micfk$Plantation)
p2=ggplot(micf, aes(x=Plantation, y=Glucose))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Fungi Substrat Induced Respiration (SIR)");
# Urea ----
lme1 = lme(fixed= Urea ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(Urea+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
a=lm(log(Urea)~1,data=mictk)
b=lm(log(Urea)~Plantation,data=mictk)
anova(a,b)
x11();ggplot(mict, aes(x=Plantation, y=Urea))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Total Urea Induced Respiration (GIR)")

lme1 = lme(fixed= Urea ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(Urea+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
a=lm(log(Urea)~1,data=micfk)
b=lm(log(Urea)~Plantation,data=micfk)
anova(a,b)
x11();ggplot(micf, aes(x=Plantation, y=Urea))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Total Urea Induced Respiration (GIR)")

lme1 = lme(fixed= Urea ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micb);
summary(lme1);
anova(lme1); # dif
# qCO2 ----
# total
lme1 = lme(fixed= qCO2 ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict,);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(qCO2) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
a=lm(log(qCO2)~1,data=mictk)
b=lm(log(qCO2)~Plantation,data=mictk)
anova(a,b)
x11();ggplot(mict, aes(x=Plantation, y=qCO2))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Total qCO2")+guides(fill=F)+theme(axis.title=element_text(size=20),axis.text=(element_text(size=15)));
# bact
lme1 = lme(fixed= qCO2 ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micb);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micb$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(qCO2) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micb);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micb$Plantation);#x11();plot(lme1);
a=lm(log(qCO2)~1,data=micbk)
b=lm(log(qCO2)~Plantation,data=micbk)
anova(a,b)
x11();ggplot(micb, aes(x=Plantation, y=qCO2))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Bacterial qCO2")
# fungi
lme1 = lme(fixed= qCO2 ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(qCO2) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=micf);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
a=lm(log(qCO2)~1,data=micfk)
b=lm(log(qCO2)~Plantation,data=micfk)
anova(a,b)
x11();ggplot(micf, aes(x=Plantation, y=qCO2))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=2,label = letter)) +
  xlab("Plantation") +
  ylab("Fungi qCO2")
# div index ----
# shannon
lme1 = lme(fixed= shannon ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(shannon+1) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
a=lm(log(shannon)~1,data=mictk)
b=lm(log(shannon)~Plantation,data=mictk)
anova(a,b)
x11();ggplot(mict, aes(x=Plantation, y=shannon))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Shannon index (H')")
# pielou
lme1 = lme(fixed= pielou ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), micf$Plantation);#x11();plot(lme1);
lme1 = lme(fixed= log(pielou) ~ Plantation,random= ~ 1|plot/inner_replicate/inner_replicate_2,data=mict);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), mict$Plantation);#x11();plot(lme1);
a=lm(log(pielou)~1,data=mictk)
b=lm(log(pielou)~Plantation,data=mictk)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), mictk$Plantationon);
kruskal.test(mictk$pielou,mictk$Plantation)
x11();ggplot(mict, aes(x=Plantation, y=pielou))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Pielou index (J')")
# graph ----
#total
mictk2=aggregate(mictk,by=list(Plantation=mictk$Plantation),mean);
mictk2=mictk2[,-c(2:4,24)];
names(mictk2);
mict3= melt(mictk2[,c(2:16,1)], id.var="Plantation");
colnames(mict3)[3] <- "Biomass";
colnames(mict3)[2] <- "SIRR";
mict3$Plantation <- factor(mict3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
#library(ggplot2)
# with plotted weight
mict3$round=round(mict3$Biomass,2) # round the Biomass for light plotting
#library(plyr)
mict3=ddply(mict3, "Plantation", transform, label_y=cumsum(Biomass)) # add the cumsum in order to give the text position
ce=arrange(mict3,Plantation,SIRR) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Biomass)-0.5*Biomass)
x11();ggplot(ce, aes(x = Plantation, y = Biomass, fill = SIRR,order=SIRR)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))+ylab("Total Substrat Induced Respiration (SIR) (micro g C-CO2/g soil.h)")
# relative
mict4=ddply(mict3, "Plantation", transform,percent_Biomass = Biomass / sum(Biomass) * 100) # add % in the tab
mict4$Plantation <- factor(mict4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mict4$roundp=round(mict4$percent_Biomass,1) # round the GIR for light plotting
mict4=ddply(mict4, "Plantation", transform, label_y=cumsum(percent_Biomass)) # add the cumsum in order to give the text position
ce=arrange(mict4,Plantation,Biomass) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Biomass)-0.5*percent_Biomass)
x11();ggplot(ce, aes(x = Plantation, y = percent_Biomass, fill = SIRR,order=SIRR )) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black")+ylab("Total relative Substrate Induced Respiration (SIR) (micro g C-CO2/g soil.h)")

# fungi
micfk2=aggregate(micfk,by=list(Plantation=micfk$Plantation),mean);
micfk2=micfk2[,-c(2:4,24)];
names(micfk2);
micf3= melt(micfk2[,c(2:16,1)], id.var="Plantation");
colnames(micf3)[3] <- "Biomass";
colnames(micf3)[2] <- "SIRR";
micf3$Plantation <- factor(micf3$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
#library(ggplot2)
# with plotted weight
micf3$round=round(micf3$Biomass,2) # round the Biomass for light plotting
#library(plyr)
micf3=ddply(micf3, "Plantation", transform, label_y=cumsum(Biomass)) # add the cumsum in order to give the text position
ce=arrange(micf3,Plantation,SIRR) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(Biomass)-0.5*Biomass)
x11();ggplot(ce, aes(x = Plantation, y = Biomass, fill = SIRR)) + geom_bar(stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE)) + geom_text(aes(y=label_y, label=round), vjust=1.2, colour="black",size=3)+ylab("Fungi Substrat Induced Respiration (SIR) (micro g C-CO2/g soil.h)")
# relative
micf4=ddply(micf3, "Plantation", transform,percent_Biomass = Biomass / sum(Biomass) * 100) # add % in the tab
micf4$Plantation <- factor(micf4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
micf4$roundp=round(micf4$percent_Biomass,1) # round the GIR for light plotting
micf4=ddply(micf4, "Plantation", transform, label_y=cumsum(percent_Biomass)) # add the cumsum in order to give the text position
ce=arrange(micf4,Plantation,Biomass) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_Biomass)-0.5*percent_Biomass)
x11();ggplot(ce, aes(x = Plantation, y = percent_Biomass, fill = SIRR,order=SIRR)) + guides(fill=guide_legend(reverse=T))+geom_bar(stat = "identity", colour="black")+ylab("Fungal Substrate Induced Respiration (SIR) (micro g C-CO2/g soil.h)")

# tot GIR Vs bronopol
mic=mic0
mic2=aggregate(mic,by=list(Plantation=mic$Plantation,bronopol=mic$bronopol),mean)
mic2=mic2[,-c(3:8,24)]
mic4= melt(mic2[,c(1,8)], id.var="Plantation") # melt to get the long format
colnames(mic4)[3] <- "GIR"
colnames(mic4)[2] <- "glucose"
mic4$Plantation <- factor(mic4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mic4$bronopol=mic2$bronopol
x11();ggplot(mic4, aes(x = Plantation, y = GIR, fill = bronopol)) + geom_bar(position="dodge",stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))
#
# GIR
mic=mic0
micc=mic
micc[1:135,7:21]=micc[1:135,7:21]-micc[136:270,7:21]
mic=micc
mic2=aggregate(mic,by=list(Plantation=mic$Plantation,bronopol=mic$bronopol),mean)
mic2=mic2[,-c(3:8,22)]
mic4= melt(mic2[,c(1,8)], id.var="Plantation") # melt to get the long format
colnames(mic4)[3] <- "GIR"
colnames(mic4)[2] <- "glucose"
mic4$Plantation <- factor(mic4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mic4$type=c("Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Fungi","Fungi","Fungi","Fungi","Fungi")
mic4$type=as.factor(mic4$type)
mic4$type <- factor(mic4$type, levels = c("Bacteria","Fungi"))
x11();ggplot(mic4, aes(x = Plantation, y = GIR, fill = type)) + geom_bar(position="dodge",stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))
# relative
# by relative composition
mic43=ddply(mic4, "Plantation", transform,percent_GIR = GIR / sum(GIR) * 100) # add % in the tab
mic43$Plantation <- factor(mic43$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mic43$roundp=round(mic43$percent_GIR,1) # round the GIR for light plotting
mic43=ddply(mic43, "Plantation", transform, label_y=cumsum(percent_GIR)) # add the cumsum in order to give the text position
ce=arrange(mic43,Plantation,type) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_GIR)-0.5*percent_GIR)
x11();ggplot(ce, aes(x = Plantation, y = percent_GIR, fill = type )) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)

# BR
mic4= melt(mic2[,c(1,11)], id.var="Plantation") # melt to get the long format
colnames(mic4)[3] <- "BR"
colnames(mic4)[2] <- "H2O"
mic4$Plantation <- factor(mic4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mic4$type=c("Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Fungi","Fungi","Fungi","Fungi","Fungi")
mic4$type=as.factor(mic4$type)
mic4$type <- factor(mic4$type, levels = c("Bacteria","Fungi"))
x11();ggplot(mic4, aes(x = Plantation, y = BR, fill = type)) + geom_bar(position="dodge",stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))
# relative
# by relative composition
mic43=ddply(mic4, "Plantation", transform,percent_BR = BR / sum(BR) * 100) # add % in the tab
mic43$Plantation <- factor(mic43$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mic43$roundp=round(mic43$percent_BR,1) # round the BR for light plotting
mic43=ddply(mic43, "Plantation", transform, label_y=cumsum(percent_BR)) # add the cumsum in order to give the text position
ce=arrange(mic43,Plantation,type) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_BR)-0.5*percent_BR)
x11();ggplot(ce, aes(x = Plantation, y = percent_BR, fill = type )) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)
#qCO2
mic4= melt(mic2[,c(1,18)], id.var="Plantation") # melt to get the long format
colnames(mic4)[3] <- "qCO2"
mic4$Plantation <- factor(mic4$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mic4$type=c("Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Fungi","Fungi","Fungi","Fungi","Fungi")
mic4$type=as.factor(mic4$type)
mic4$type <- factor(mic4$type, levels = c("Bacteria","Fungi"))
x11();ggplot(mic4, aes(x = Plantation, y = qCO2, fill = type)) + geom_bar(position="dodge",stat = "identity", colour="black") + guides(fill=guide_legend(reverse=TRUE))
#
mic43=ddply(mic4, "Plantation", transform,percent_qCO2 = qCO2 / sum(qCO2) * 100) # add % in the tab
mic43$Plantation <- factor(mic43$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
mic43$roundp=round(mic43$percent_qCO2,1) # round the qCO2 for light plotting
mic43=ddply(mic43, "Plantation", transform, label_y=cumsum(percent_qCO2)) # add the cumsum in order to give the text position
ce=arrange(mic43,Plantation,type) # to have number in the middle of each bar
ce <- ddply(ce, "Plantation", transform, label_y=cumsum(percent_qCO2)-0.5*percent_qCO2)
x11();ggplot(ce, aes(x = Plantation, y = percent_qCO2, fill = type )) + guides(fill=guide_legend(reverse=TRUE))+geom_bar(stat = "identity", colour="black") + geom_text(aes(y=label_y, label=roundp), vjust=1.2, colour="black",size=3)


# total graph by compound
mictt=melt(mict[,c(7:14,16:21,27,22)], id.var="Plantation");
colnames(mictt)[3] <- "Respiration_total";
mictt$Plantation <- factor(mictt$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
c <- ggplot(mictt, aes(x = Plantation, y = Respiration_total));c + facet_wrap(~ variable) + stat_summary(fun.y="mean", geom="bar")+ylab("Total SIR (micro g C-CO2/g soil.h)");
# fungi graph by compound
micft=melt(micf[,c(7:14,16:21,27,22)], id.var="Plantation")
colnames(micft)[3] <- "Respiration_fungi"
micft$Plantation <- factor(micft$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
c <- ggplot(micft, aes(x = Plantation, y = Respiration_fungi));c + facet_wrap(~ variable) + stat_summary(fun.y="mean", geom="bar")+ylab("Fungi SIR (micro g C-CO2/g soil.h)")
# bacteria graph by compound
micbt=melt(micb[,c(7:14,16:21,24,22)], id.var="Plantation")
colnames(micbt)[3] <- "Respiration_bacteria"
micbt$Plantation <- factor(micbt$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
c <- ggplot(micbt, aes(x = Plantation, y = Respiration_bacteria));c + facet_wrap(~ variable) + stat_summary(fun.y="mean", geom="bar")+ylab("Bacteria SIR (micro g C-CO2/g soil.h)")

# adjusted scale
# total graph by compound
mictt=melt(mict[,c(7:14,16:21,27,22)], id.var="Plantation");
colnames(mictt)[3] <- "Respiration_total";
mictt$Plantation <- factor(mictt$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"));
c <- ggplot(mictt, aes(x = Plantation, y = Respiration_total));c + facet_wrap(~ variable,scales="free_y") + stat_summary(fun.y="mean", geom="bar")+ylab("Total SIR (micro g C-CO2/g soil.h)");
# fungi graph by compound
micft=melt(micf[,c(7:14,16:21,27,22)], id.var="Plantation")
colnames(micft)[3] <- "Respiration_fungi"
micft$Plantation <- factor(micft$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
c <- ggplot(micft, aes(x = Plantation, y = Respiration_fungi));c + facet_wrap(~ variable,scales="free_y") + stat_summary(fun.y="mean", geom="bar")+ylab("Fungi SIR (micro g C-CO2/g soil.h)")
# bacteria graph by compound
micbt=melt(micb[,c(7:14,16:21,24,22)], id.var="Plantation")
colnames(micbt)[3] <- "Respiration_bacteria"
micbt$Plantation <- factor(micbt$Plantation, levels = c("N o","R1 y","R1 o","R2 y","R2 o"))
c <- ggplot(micbt, aes(x = Plantation, y = Respiration_bacteria));c + facet_wrap(~ variable,scales="free_y") + stat_summary(fun.y="mean", geom="bar")+ylab("Bacteria SIR (micro g C-CO2/g soil.h)")

# STAR PLOT
mict6=mic2[1:5,]
lit3=lit3[c(1,4,2,5,3),] # put the row in the right order
palette(rainbow(12, s = 0.6, v = 0.75))
row.names(mict6) <- c("NO","R1O","R2O","R1Y","R2Y") # define row names
x11();stars(mict6[, c(3:18)], len = 0.8, key.loc = c(4.7, 2.3),main = "Litter structure", draw.segments = TRUE)
library("fmsb")
x11();radarchart(mict6[, 3:18],pfcol=c(1:5),plty=1,pdensity=0,axistype=4,seg=4,plwd=2,maxmin=F)
legend("topleft", legend=unique(mict6$Plantation), seg.len=1.5, title="state", pch=1,bty="n" ,lwd=3, y.intersp=0.5, horiz=FALSE, col=c(1:5))
# multi ----
# TOT
# PCA
library(ade4);
coa2<-dudi.pca(mict[,c(7:21)],scannf=F,nf=3);
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red");
coa2$eig;
pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(mdtot$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,mict$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));par(mfrow = c(1, 2));s.class(bca2$ls,mict$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);
title(xlab="Dim 1 (69%) ",ylab="Dim 2 (25%)",main="p-value:0.001",font.main=1);
s.corcircle(bca2$co,box=T,clabel=0.5,full=T);
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference
x11();plot(pt2,main="Between class inertia")
# FUNGI
coa2<-dudi.pca(micf[,c(7:22)],scannf=F,nf=3)
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red")
pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(mdtot$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,mict$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));par(mfrow = c(1, 2));s.class(bca2$ls,mict$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T)
title(xlab="Dim 1 (91%) ",ylab="Dim 2 (8%)",main="p-value:0.001",font.main=1)
s.corcircle(bca2$co,box=T,clabel=0.5,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference
x11();plot(pt2,main="Between class inertia")
# bacteria
coa2<-dudi.pca(micb[,c(7:22)],scannf=F,nf=4)
# Kaiser-Guttman
barplot(coa2$eig, main="Eigenvalues", col="grey")
abline(h=mean(coa2$eig), col="red")
pvp=100*coa2$eig/sum(coa2$eig);pvp;
cumsum(pvp)
#x11();s.class(coa2$li,as.factor(mdtot$Plantation),cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F);s.arrow(coa2$co,boxes=F,add.plot=T)
# BCA
bca2<-bca(coa2,mict$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
#x11();plot(bca2)
x11();par(mfrow = c(1, 2));par(mfrow = c(1, 2));s.class(bca2$ls,mict$Plantation,cell = 1.5, axesell = F, cstar = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);#s.arrow(bca2$co,boxes=T,add.plot=T)
title(xlab="Dim 1 (53%) ",ylab="Dim 2 (42%)",main="p-value:0.001",font.main=1)
s.corcircle(bca2$co,box=T,clabel=0.5,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # non significative difference
x11();plot(pt2,main="Between class inertia")

### QUALITATIVE ----
qua <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/qualitative_0.csv", dec=",")
str(qua)
qua$plot=as.factor(qua$plot)
qua$age <- factor(qua$age, levels=c("y","o"))
qua$Plantation=paste(qua$state,qua$age)
qua$Plantation=as.factor(qua$Plantation)
qua$Plantation <- factor(qua$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"))
qua$slope=as.factor(qua$slope)
qua$light=as.factor(qua$light)
qua$wood=as.factor(qua$wood)
x11();plot(qua[,c(15,3,9:14)])
# more nest
# less dead tree, litter, slope ?
# termite
lm0=lm(termite_nest~1,data=qua)
lm1=lm(termite_nest~Plantation,data=qua)
anova(lm0,lm1)
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), qua$Plantation);
lm0=lm(log(termite_nest+1)~1,data=qua)
lm1=lm(log(termite_nest+1)~Plantation,data=qua)
anova(lm0,lm1)
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), qua$Plantation);
kruskal.test(qua$termite_nest,qua$Plantation) # dif
test=pairwise.wilcox.test(qua$termite_nest,qua$Plantation, p.adj = "bonf");test
x11();ggplot(qua, aes(x=Plantation, y=termite_nest))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Termite nest")+ guides(fill=FALSE)
# dead tree
lm0=lm(dead_tree~1,data=qua)
lm1=lm(dead_tree~Plantation,data=qua)
anova(lm0,lm1)
# litter
lm0=lm(litter_cover_..~1,data=qua)
lm1=lm(litter_cover_..~Plantation,data=qua)
anova(lm0,lm1)
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), qua$Plantation);
lm0=lm(log(litter_cover_..+1)~1,data=qua)
lm1=lm(log(litter_cover_..+1)~Plantation,data=qua)
anova(lm0,lm1)
shapiro.test(residuals(lm1));
bartlett.test(residuals(lm1), qua$Plantation);
kruskal.test(qua$litter_cover_..,qua$Plantation) # dif
test=pairwise.wilcox.test(qua$litter_cover_..,qua$Plantation, p.adj = "bonf");test
x11();ggplot(qua, aes(x=Plantation, y=litter_cover_..))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Litter cover")+ guides(fill=FALSE)


#### CHEMICAL #######################################################################
## required dataset ----
ch <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/chemical_1.csv")
#
ch$Plantation=paste(ch$state,ch$age);
ch$Plantation=as.factor(ch$Plantation);
levels(ch$Plantation);
ch$Plantation <- factor(ch$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
ch$age=as.factor(ch$age)
ch$age <- factor(ch$age, levels=c("y","o")) ;
colnames(ch)[11] <- "SAND";
colnames(ch)[12] <- "SILT";
colnames(ch)[13] <- "CLAY";
#
ch05=ch[1:45,];
ch510=ch[46:90,];
#
ch2=aggregate(ch,by=list(state=ch$state,age=ch$age,plot=ch$plot,depth_cm=ch$depth_cm),mean);
ch2=ch2[,-c(7:14,27)];
ch2$Plantation=paste(ch2$state,ch2$age);
ch2$Plantation=as.factor(ch2$Plantation);
levels(ch2$Plantation);
ch2$Plantation <- factor(ch2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")) ;
ch2$age <- factor(ch2$age, levels=c("y","o")) ;
ch052=ch2[1:15,];
ch5102=ch2[16:30,];
## pH ----
### 05
lme1 = lme(fixed= pH_1.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch05$Plantation),letter=grpplot);
g5=ggplot(ch05, aes(x=Plantation, y=pH_1.1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=6.5,label = letter)) +
  xlab("Plantation") +
  ylab("pH 1:1 0-5 cm")+ guides(fill=FALSE);
# decrease ?
### 10
lme1 = lme(fixed= EC_01.05_dS.m.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation);#x11();plot(lme1); homo
lme1 = lme(fixed= log(EC_01.05_dS.m.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation);#x11();plot(lme1); homo
a=lm(OM_.~1,data=ch05)
b=lm(OM_.~Plantation,data=ch05)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch05$Plantation);
lettre=cld(summary(glht(b,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")))
lettre$mcletters
let=c("c","ab","bc","a","b")
xx=data.frame(level=levels(ch05$Plantation),letter=let)
x11();ggplot(ch05, aes(x=Plantation, y=OM_.))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=2.6,label = letter)) +
  xlab("Plantation") +
  ylab("OM (%) 0-5 cm")+ guides(fill=FALSE);
## OM ----
### 05
lme1 = lme(fixed= OM_. ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
lme1 = lme(fixed= log(OM_.) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
a=lm(OM_.~1,data=ch05);
b=lm(OM_.~Plantation,data=ch05);
anova(a,b);
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch05$Plantation);
lettre=cld(summary(glht(b,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")));
lettre$mcletters;
let=c("c","ab","bc","a","b");
xx=data.frame(level=levels(ch05$Plantation),letter=let);
x11();ggplot(ch05, aes(x=Plantation, y=OM_.))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=2.6,label = letter)) +
  xlab("Plantation") +
  ylab("OM (%) 0-5 cm")+ guides(fill=FALSE);
# 510
lme1 = lme(fixed= OM_. ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch510$Plantation);
lme1 = lme(fixed= log(OM_.) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch510$Plantation);
a=lm(log(OM_.)~1,data=ch510)
b=lm(log(OM_.)~Plantation,data=ch510)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch510$Plantation);
lettre=cld(summary(glht(b,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")))
lettre$mcletters
let=c("b","ab","ab","a","b")
xx=data.frame(level=levels(ch05$Plantation),letter=let)
x11();ggplot(ch510, aes(x=Plantation, y=OM_.))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=2.6,label = letter)) +
  xlab("Plantation") +
  ylab("OM (%) 5-10 cm")+ guides(fill=FALSE);
## Ca ----
### 05
lme1 = lme(fixed= ca.. ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
lme1 = lme(fixed= log(ca..) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch05$Plantation),letter=grpplot);
g1=ggplot(ch05, aes(x=Plantation, y=ca..))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=700,label = letter)) +
  xlab("Plantation") +
  ylab("Ca++ mg.kg-1 0-5 cm")+ guides(fill=FALSE);
# 510
### 10
lme1 = lme(fixed= ca.. ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation); # hetero
lme1 = lme(fixed= log(ca..) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation); # hetero
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(ch510$Plantation),letter=grpplot)
x11();ggplot(ch510, aes(x=Plantation, y=ca..))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=350,label = letter)) +
  xlab("Plantation") +
  ylab("Ca++ mg.kg-1 5-10 cm")
## Mg ----
### 05
lme1 = lme(fixed= mg.._extractNH4OAc_mg_kg.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
summary(lme1);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
lme1 = lme(fixed= log(mg.._extractNH4OAc_mg_kg.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch05$Plantation),letter=grpplot);
g2=ggplot(ch05, aes(x=Plantation, y=mg.._extractNH4OAc_mg_kg.1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=75,label = letter)) +
  xlab("Plantation") +
  ylab("Mg++ mg.kg-1 0-5 cm")+ guides(fill=FALSE);
# 510
lme1 = lme(fixed= mg.._extractNH4OAc_mg_kg.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation); # hetero
lme1 = lme(fixed= log(mg.._extractNH4OAc_mg_kg.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), ch510$Plantation); # homo
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters
xx=data.frame(level=levels(ch510$Plantation),letter=grpplot)
x11();ggplot(ch510, aes(x=Plantation, y=mg.._extractNH4OAc_mg_kg.1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=40,label = letter)) +
  xlab("Plantation") +
  ylab("Mg++ mg.kg-1 5-10 cm")
# EC ----
### 05
lme1 = lme(fixed= EC_01.05_dS.m.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1));
bartlett.test(residuals(lme1), ch510$Plantation); # homo
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch05$Plantation),letter=grpplot);
g6=ggplot(ch05, aes(x=Plantation, y=EC_01.05_dS.m.1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=0.075,label = letter)) +
  xlab("Plantation") +
  ylab("EC 01:05 dS.m-1 0-5 cm")+ guides(fill=FALSE);
### 10
lme1 = lme(fixed= EC_01.05_dS.m.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation);#x11();plot(lme1); homo
lme1 = lme(fixed= log(EC_01.05_dS.m.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation);#x11();plot(lme1); homo
a=lm(log(EC_01.05_dS.m.1)~1,data=ch05)
b=lm(log(EC_01.05_dS.m.1)~Plantation,data=ch05)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch05$Plantation);
## K ----
### 05
lme1 = lme(fixed= K._extractNH4OAc_mg_kg.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), ch05$Plantation); # hetero
lme1 = lme(fixed= log(K._extractNH4OAc_mg_kg.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch05$Plantation); # homo
a=lm(log(K._extractNH4OAc_mg_kg.1)~1,data=ch052);
b=lm(log(K._extractNH4OAc_mg_kg.1)~Plantation,data=ch052);
anova(a,b)
shapiro.test(residuals(b));#x11();qqnorm(residuals(b));qqline(residuals(b));
bartlett.test(residuals(b), ch052$Plantation);

lettre=cld(summary(glht(b,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")));
lettre$mcletters
#let=scan(what="")
let=c("b","b","b","a","a");
xx=data.frame(level=levels(ch05$Plantation),letter=let);
x11();ggplot(ch05, aes(x=Plantation, y=K._extractNH4OAc_mg_kg.1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=150,label = letter)) +
  xlab("Plantation") +
  ylab("K+ 0-5cm (mg/kg)")+guides(fill=F)+theme(axis.title=element_text(size=20),axis.text=(element_text(size=15)));

### 10
lme1 = lme(fixed= K._extractNH4OAc_mg_kg.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), ch510$Plantation); # hetero
lme1 = lme(fixed= log(K._extractNH4OAc_mg_kg.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation); # homo
a=lm(log(K._extractNH4OAc_mg_kg.1)~1,data=ch5102);
b=lm(log(K._extractNH4OAc_mg_kg.1)~Plantation,data=ch5102);
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch5102$Plantation);

lettre=cld(summary(glht(b,linfct=mcp(Plantation="Tukey")), test = adjusted(type = "bonferroni")));
lettre$mcletters
#let=scan(what="")
let=c("b","b","b","a","a");
xx=data.frame(level=levels(ch05$Plantation),letter=let);
x11();ggplot(ch510, aes(x=Plantation, y=K._extractNH4OAc_mg_kg.1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx,size=9, aes(x=c(1:5),y=75,label = letter)) +
  xlab("Plantation") +
  ylab("K+ 5-10cm (mg/kg)")+guides(fill=F)+theme(axis.title=element_text(size=20),axis.text=(element_text(size=15)));
## P ----
#### O5
lme1 = lme(fixed= P_Bray_II_mg_kg.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch05$Plantation); # hetero
lme1 = lme(fixed= log(P_Bray_II_mg_kg.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); normal
bartlett.test(residuals(lme1), ch05$Plantation); # homo
summary(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch05$Plantation),letter=grpplot);
x11();ggplot(ch05, aes(x=Plantation, y=P_Bray_II_mg_kg.1))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=52,label = letter)) +
  xlab("Plantation") +
  ylab("P mg.kg-1 0-5 cm")+guides(fill=F);

### 10
lme1 = lme(fixed= P_Bray_II_mg_kg.1 ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
summary(lme1);
anova(lme1); # dif
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation); #plot(lme1) # hetero
lme1 = lme(fixed= log(P_Bray_II_mg_kg.1) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
shapiro.test(residuals(lme1));#x11();qqnorm(residuals(lme1));qqline(residuals(lme1)); non normal
bartlett.test(residuals(lme1), ch510$Plantation); plot(lme1) # homo
a=lm(log(P_Bray_II_mg_kg.1)~1,data=ch510)
b=lm(log(P_Bray_II_mg_kg.1)~Plantation,data=ch510)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch510$Plantation);
kruskal.test(ch5102$P_Bray_II_mg_kg.1,ch5102$Plantation) # no diff
# SAND ----
#05
lme1 = lme(fixed= SAND ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1)
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
lme1 = lme(fixed= log(SAND) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
a=lm(log(SAND)~1,data=ch05)
b=lm(log(SAND)~Plantation,data=ch05)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch05$Plantation);
kruskal.test(ch052$SAND,ch052$Plantation) # non signif
x11();ggplot(ch05, aes(x=Plantation, y=SAND))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Sand (%)")+ guides(fill=FALSE);
# 510
lme1 = lme(fixed= SAND ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1)
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch510$Plantation);
lme1 = lme(fixed= log(SAND) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1)
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch510$Plantation);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch510$Plantation),letter=grpplot);
x11();ggplot(ch510, aes(x=Plantation, y=SAND))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=90,label = letter)) +
  xlab("Plantation") +
  ylab("Sand (%) (05-10 cm)")+ guides(fill=FALSE);
# CLAY ----
# 05
lme1 = lme(fixed= CLAY ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1)
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch05$Plantation),letter=grpplot);
x11();ggplot(ch05, aes(x=Plantation, y=CLAY))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=20,label = letter)) +
  xlab("Plantation") +
  ylab("Clay (%)")+ guides(fill=FALSE);
# 510
lme1 = lme(fixed= CLAY ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1)
# SILT ----
# 05
lme1 = lme(fixed= SILT ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1)
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
lme1 = lme(fixed= log(SILT) ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch05);
anova(lme1)
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch05$Plantation);
a=lm(log(SILT)~1,data=ch05)
b=lm(log(SILT)~Plantation,data=ch05)
anova(a,b)
shapiro.test(residuals(a));
bartlett.test(residuals(b), ch05$Plantation);
kruskal.test(ch052$SILT,ch052$Plantation)
x11();ggplot(ch05, aes(x=Plantation, y=SILT))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Silt (%)")+ guides(fill=FALSE);
# 5-10
lme1 = lme(fixed= SILT ~ Plantation,random= ~ 1|plot/inner_replicate,data=ch510);
anova(lme1)
shapiro.test(residuals(lme1));
bartlett.test(residuals(lme1), ch510$Plantation);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(ch510$Plantation),letter=grpplot);
x11();ggplot(ch510, aes(x=Plantation, y=SILT))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=25,label = letter)) +
  xlab("Plantation") +
  ylab("Silt (%)")+ guides(fill=FALSE);
# graph texture ----
#install.packages( pkgs = "soiltexture" )
library(soiltexture)
x11();tri=TT.plot( class.sys = "USDA.TT",b.lim=c(0,4))
TT.text(tri.data=ch05[,11:13],geo=tri,labels=ch05$Plantation,cex=0.7,font=2,col=as.numeric(ch05$Plantation))

install.packages("ggtern")
# Load the required libraries
library(ggtern)
library(plyr)
library(grid)
# Load the Data. (Available in ggtern 1.0.3.0 next version)
data(USDA)
# Put tile labels at the midpoint of each tile.
USDA.LAB = ddply(USDA, 'Label', function(df) {
    apply(df[, 1:3], 2, mean)
})
# Tweak
USDA.LAB$Angle = 0
USDA.LAB$Angle[which(USDA.LAB$Label == 'Loamy Sand')] = -35
# Construct the plot.
t8=ggplot(data = USDA, aes(y=Clay, x=Sand, z=Silt)) +
  coord_tern(L="x",T="y",R="z") +
  geom_polygon(aes(fill = Label),
               alpha = 0.75, size = 0.5, color = 'black') +
  geom_text(data = USDA.LAB,
            aes(label = Label, angle = Angle),
            color = 'black',
            size = 3.5) +
  theme_rgbw() +
  theme_showsecondary() +
  theme_showarrows() +
  custom_percent("Percent") +
  theme(legend.justification = c(0, 1),
        legend.position      = c(0, 1),
        axis.tern.padding    = unit(0.15, 'npc')) +
  labs(title = 'USDA Textural Classification Chart',
       fill  = 'Textural Class',
       color = 'Textural Class')+guides(fill=FALSE);  

x11();t8 +geom_point(data=ch510,aes(y=CLAY,x=SAND,z=SILT,colour=Plantation),shape=19,show_guide = TRUE)+scale_color_discrete(name = "Plantation")
## MOISTURE ----
moi<- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/moisture.csv", dec=",");
moi$plot=as.factor(moi$plot);
moi$inner_replicate=as.factor(moi$inner_replicate);
moi$Plantation=paste(moi$state,moi$age);
moi$Plantation=as.factor(moi$Plantation);
moi$Plantation <- factor(moi$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
moi$age <- factor(moi$age, levels=c("y","o"));
lme1 = lme(fixed= X.moisture ~ Plantation,random= ~ 1|plot/inner_replicate,data=moi);
summary(lme1);
anova(lme1);
grp1 = cld(glht(lme1,linfct=mcp(Plantation="Tukey"),adjust='bonferroni'));grp1;
grpplot=grp1$mcletters$Letters;
xx=data.frame(level=levels(moi$Plantation),letter=grpplot);
g4=ggplot(moi, aes(x=Plantation, y=X.moisture))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  geom_text(data=xx, aes(x=c(1:5),y=15,label = letter)) +
  xlab("Plantation") +
  ylab("Moisture content (%)")+ guides(fill=FALSE);
# more moisture in old
# but less in R20

x11();grid.arrange(g1,g2,g3,g4,g5,g6,nrow=2)
## adjusted p.val  ----
# bonferroni simple sequentialy rejective test
pval=data.frame(test=c("Clay 0-5","K 5-10","K 0-5","OM 5-10","pH 0-5","Ca 0-5","Ca 5-10","Mg 0-5","Mg 5-10","EC 0-5","P 0-5","P 5-10","Sand 5-10","Silt 5-10","pH 5-10","OM 0-5"),Type=c("LMM","K","K","LM log","LMM","LMM log","LMM log","LMM log","LMM log","LMM log","LMM log","LM","LMM log","LMM","LM","LM"),p.val=c(0.0378,0.02681,0.02306,0.0206,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,1.545*10^-6,1.545*10^-6),factor=c(1:16));
pval$p.val.adjusted=pval$p.val*pval$factor;pval#----
## Multivariate analysis ----
library(ade4);
## 05
# PCA  

colnames(ch05)[15] <- "pH";
colnames(ch05)[17] <- "EC";
colnames(ch05)[18] <- "OM";
colnames(ch05)[19] <- "P";
colnames(ch05)[20] <- "Ca";
colnames(ch05)[21] <- "Mg";
colnames(ch05)[22] <- "K";

pca2<-dudi.pca(ch05[,c(15,17:22)],scannf=F,nf=2);
barplot(pca2$eig, main="Eigenvalues", col="grey");
abline(h=1, col="red");
# BCA
bca2<-bca(pca2,ch05$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,ch05$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5)
title(xlab="Dim 1 (73%) ",ylab="Dim 2 (17%)",main="p-value: 0.001",font.main=1)
s.corcircle(bca2$co,clabel=0.5,box=T,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # difference
x11();plot(pt2,main="Between class inertia")
# avec texture
pca2<-dudi.pca(ch05[,c(11,12,13,15,17:22)],scannf=F,nf=4);
barplot(pca2$eig, main="Eigenvalues", col="grey")
abline(h=1, col="red")
# BCA
bca2<-bca(pca2,ch05$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,ch05$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5);
title(xlab="Dim 1 (71%) ",ylab="Dim 2 (15%)",main="p-value: 0.001",font.main=1);
s.corcircle(bca2$co,clabel=0.5,box=T,full=F);
# Permutation test
pt2=rtest(bca2,999)
pt2 # difference
x11();plot(pt2,main="Between class inertia")


## 10
# CA
colnames(ch510)[12] <- "pH"
colnames(ch510)[14] <- "EC"
colnames(ch510)[15] <- "OM"
colnames(ch510)[16] <- "P"
colnames(ch510)[17] <- "Ca"
colnames(ch510)[18] <- "Mg"
colnames(ch510)[19] <- "K"

pca2<-dudi.pca(ch510[,c(12,14:19)],scannf=F,nf=2)
barplot(pca2$eig, main="Eigenvalues", col="grey")
abline(h=1, col="red")
# BCA
bca2<-bca(pca2,ch510$Plantation,scannf=FALSE)
pvp=100*bca2$eig/sum(bca2$eig);pvp;
cumsum(pvp)
x11();par(mfrow = c(1, 2));s.class(bca2$ls,ch510$Plantation,cell = 1.5, axesell = F, csta = 1,col=c(1,2,3,4,6),grid=F,clabel=0.5)
title(xlab="Dim 1 (68%) ",ylab="Dim 2 (21%)",main="p-value: 0.001",font.main=1)
s.corcircle(bca2$co,box=T,clabel=0.5,full=F)
# Permutation test
pt2=rtest(bca2,999)
pt2 # difference
x11();plot(pt2,main="Between class inertia")

#### PHYSICAL #######################################################################

## AGGREGATE ----
agg <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/aggregate_1.csv", dec=",");
str(agg);
agg$plot=as.factor(agg$plot);
agg$inner_replicate=as.factor(agg$inner_replicate);
agg$Plantation=paste(agg$state,agg$age);
agg$Plantation=as.factor(agg$Plantation);
agg$Plantation <- factor(agg$Plantation, levels=c("N o","R1 o","R2 o"));
agg$age <- factor(agg$age, levels=c("y","o"));
# degree of aggregation
lme1 = lme(fixed= degree_of_aggregation ~ Plantation,random= ~ 1|plot/inner_replicate,data=agg);
summary(lme1);
anova(lme1);
x11();ggplot(agg, aes(x=Plantation, y=degree_of_aggregation))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Degree of aggregation")
# aggregate stability
agg1=na.exclude(agg)
lme1 = lme(fixed= aggregate_stability ~ Plantation,random= ~ 1|plot/inner_replicate,data=agg1);
summary(lme1);
anova(lme1); # no dif
x11();ggplot(agg, aes(x=Plantation, y=aggregate_stability))  + geom_boxplot(aes(fill=age),notch=F) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") + 
  xlab("Plantation") +
  ylab("Aggregate stability")
# ASI
lme1 = lme(fixed= ASI ~ Plantation,random= ~ 1|plot/inner_replicate,data=agg1);
summary(lme1);
anova(lme1); # no dif
## METEO ----

met <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/meteo.csv")
summary(met)
map=124.36
pdry=(98.4+25.7+6.3+12+76.4)/5
# ??

# probleme d'axe X classe alphabetiquement
met$month <- factor(met$month, ordered = TRUE, levels = as.character(met$month)) # prevent the alphabetic clasification
# ou
#met$month= as.character(met$month); met$month= factor(met$month,levels=unique(met$month))

met$mean.temp=apply(met[,2:3],1,mean)

x11();
par(mar=c(5.1,4.1,4.1,4.1)) # change margin in order: bottom, left, top, right in number of writing line
barplot(met$Total.rainfall,names.arg=met$month,ylab="Average total rainfall (mm)",xlab="Month",ylim=c(0,300),col="blue")
par(new=TRUE)
plot(met$mean.temp, type = c("o"),lty = 'solid',xaxt="n",yaxt="n",xlab="",ylab="", pch=20,col ="red")
points(mean.temp ~ as.numeric(month),col="red", data = met, type = 'l')
axis(4)
mtext("Mean temperature (°C)",side=4,line=3)
legend("topleft",
       c("rainfall","temperature"),
       cex = 0.75, 
       title = "",
       lty=c(1,1),
       lwd=c(2.5,2.5),col=c("blue","red")) # gives the legend lines the correct color and width
       bty = "")#----
## MAP SITES ----
library(ggmap)
palette("Set1")
sites.data <- read.delim("~/Documents/cours/UPMC/Master/S3/stage/data/coord4.csv")
sites.data$label=paste(sites.data$state,sites.data$age);
sites.data$label=as.factor(sites.data$label);
sites.data$label <- factor(sites.data$label, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
map.base<-get_map(location=c(longitude=99.45,latitude=8.92),source = "google",zoom = 13)
colmapping <- unique(sites.data[,10:11])

## scale:
sites.data$colour=as.factor(sites.data$colour)
bb <- attr(map.base,"bb")
distHaversine <- function(long, lat){
  dlong = (long[2] - long[1])*pi/180
  dlat  = (lat[2] - lat[1])*pi/180
  
  # Haversine formula:
  R = 6371;
  a = sin(dlat/2)*sin(dlat/2) + cos(lat[1])*cos(lat[2])*sin(dlong/2)*sin(dlong/2)
  c = 2 * atan2( sqrt(a), sqrt(1-a) )
  d = R * c
  return(d) # in km
}
sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                   lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                   lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                   lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))

sbar$distance = distHaversine(long = c(sbar$lon.start,sbar$lon.end),
                              lat = c(sbar$lat.start,sbar$lat.end))

ptspermm <- 2.83464567  # need this because geom_text uses mm, and themes use pts. Urgh.

map.scale=ggmap(map.base)+ 
  geom_point(aes(x= lon, y=lat,colour=label),shape=19,data=sites.data,show_guide = T)+
  scale_color_discrete(name = "Plantation")+
  theme(legend.text=element_text(size=10))+
  xlab("longitude") +
  ylab("latitude")+
  geom_segment(data = sbar,
               aes(x = lon.start,
                   xend = lon.end,
                   y = lat.start,
                   yend = lat.end)) +
  geom_text(data = sbar,
            aes(x = (lon.start + lon.end)/2,
           y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat),
           label = paste(format(distance, 
                                digits = 4,
                                nsmall = 2),
                         'km')),
           hjust = 0.5,
           vjust = 0,
           size = 8/ptspermm)

# Fix presentation
map.out <- map.scale +  
  theme_bw(base_size = 8)
ggsave(file ="/home/raph/Documents/cours/UPMC/Master/S3/stage/rapport/map4.pdf", 
       plot = map.out,
       width = 4, 
       height = 3)#----
# MAP THAILAND ----
map.base<-get_map(location=c(longitude=99.45,latitude=8.92),source = "google",zoom = 5)
bb <- attr(map.base,"bb")
distHaversine <- function(long, lat){
  dlong = (long[2] - long[1])*pi/180
  dlat  = (lat[2] - lat[1])*pi/180
  
  # Haversine formula:
  R = 6371;
  a = sin(dlat/2)*sin(dlat/2) + cos(lat[1])*cos(lat[2])*sin(dlong/2)*sin(dlong/2)
  c = 2 * atan2( sqrt(a), sqrt(1-a) )
  d = R * c
  return(d) # in km
}
sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                   lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                   lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                   lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))

sbar$distance = distHaversine(long = c(sbar$lon.start,sbar$lon.end),
                              lat = c(sbar$lat.start,sbar$lat.end))

ptspermm <- 2.83464567  # need this because geom_text uses mm, and themes use pts. Urgh.

map.scale=ggmap(map.base)+ 
  xlab("longitude") +
  ylab("latitude")+
  geom_segment(data = sbar,
               aes(x = lon.start,
                   xend = lon.end,
                   y = lat.start,
                   yend = lat.end)) +
  geom_text(data = sbar,
            aes(x = (lon.start + lon.end)/2,
           y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat),
           label = paste(format(distance, 
                                digits = 4,
                                nsmall = 2),
                         'km')),
           hjust = 0.5,
           vjust = 0,
           size = 8/ptspermm)+
  geom_point(aes(x= 99.4, y=8.9,colour="red"),shape=19)+
    scale_color_discrete(name = "District",label="Surat Thani")

# Fix presentation
map.out <- map.scale +  
  theme_bw(base_size = 8)
ggsave(file ="/home/raph/Documents/cours/UPMC/Master/S3/stage/rapport/map5.pdf", 
       plot = map.out,
       width = 4, 
       height = 3)#----

# altitude ----
sites.data <- read.delim("~/Documents/cours/UPMC/Master/S3/stage/data/coord4.csv")
rownames(sites.data) = sites.data[,1]

googEl <- function(locs)  {
  require(RJSONIO)
  locstring <- paste(do.call(paste, list(locs[, 2], locs[, 1], sep=',')),
                     collapse='|')
  u <- sprintf('http://maps.googleapis.com/maps/api/elevation/json?locations=%s&sensor=false',
               locstring)
  res <- fromJSON(u)
  out <- t(sapply(res[[1]], function(x) {
    c(x[['location']]['lat'], x[['location']]['lng'], 
      x['elevation'], x['resolution']) 
  }))    
  rownames(out) <- rownames(locs)
  return(out)
}

alt=googEl(sites.data[,c(9,8)]) # longitude have to be in the first column
#write.table(alt, "/home/raph/Documents/cours/UPMC/Master/S3/stage/data/alt.txt",row.names = T) 
alt <- read.table("~/Documents/cours/UPMC/Master/S3/stage/data/alt.txt", header=T, quote="\"")

plot3d(alt$lat,alt$long,alt$elevation)


library(rgl)
y <- 2 * volcano # Exaggerate the relief
x <- 10 * (1:nrow(y)) # 10 meter spacing (S to N)
z <- 10 * (1:ncol(y)) # 10 meter spacing (E to W)
ylim <- range(y)
ylen <- ylim[2] - ylim[1] + 1
colorlut <- terrain.colors(ylen) # height color lookup table
col <- colorlut[ y-ylim[1]+1 ] # assign colors to heights
rgl.open()
rgl.surface(x, z, y, color=col, back="lines") 

# MCOA ----
## required dataset ----
ch <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/chemical_1.csv")
#
ch$Plantation=paste(ch$state,ch$age);
ch$Plantation=as.factor(ch$Plantation);
levels(ch$Plantation);
ch$Plantation <- factor(ch$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
ch$age=as.factor(ch$age)
ch$age <- factor(ch$age, levels=c("y","o")) ;
colnames(ch)[11] <- "SAND";
colnames(ch)[12] <- "SILT";
colnames(ch)[13] <- "CLAY";
#
ch05=ch[1:45,];
ch510=ch[46:90,];
#
ch2=aggregate(ch,by=list(state=ch$state,age=ch$age,plot=ch$plot,depth_cm=ch$depth_cm),mean);
ch2=ch2[,-c(7:14,27)];
ch2$Plantation=paste(ch2$state,ch2$age);
ch2$Plantation=as.factor(ch2$Plantation);
levels(ch2$Plantation);
ch2$Plantation <- factor(ch2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")) ;
ch2$age <- factor(ch2$age, levels=c("y","o")) ;
ch052=ch2[1:15,];
ch5102=ch2[16:30,];
colnames(ch05)[15] <- "pH";
colnames(ch05)[17] <- "EC";
colnames(ch05)[18] <- "OM";
colnames(ch05)[19] <- "P";
colnames(ch05)[20] <- "Ca";
colnames(ch05)[21] <- "Mg";
colnames(ch05)[22] <- "K";

md <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_density_0.csv");
md[md$strata=="Litter",c(6:84)]=md[md$strata=="Litter",c(6:84)]*16;
md[md$strata=="Soil",c(6:84)]=md[md$strata=="Soil",c(6:84)]*160;
md$oligochaeta=md$oligochaeta+md$oligochaeta_diapause; # sum oligochaeta
names(md[c(6,12,14,22,24,80:84)]);
md=md[,-c(6,12,14,22,24,80:84)]; # remove mesofauna & egg & oligochaeta_diapause
str(md);
md$plot=as.factor(md$plot);
md$inner_replicate=as.factor(md$inner_replicate);
md$age <- factor(md$age, levels=c("y","o"));
md$Plantation=paste(md$state,md$age);
md$Plantation=as.factor(md$Plantation);
md$Plantation <- factor(md$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
md$sp=apply(md[,20:74],1,sum); # sum sp
md$total=apply(md[,6:74],1,sum);
md$specnumber=specnumber(md[,6:74]);
md$shannon=diversity(md[,6:74]);
md[md$specnumber==0,]$shannon=NA; # non sense
md$pielou = md$shannon/log(md$specnumber); # NaN problem:
which(md$shannon==0);md$pielou[c(1,5,9,31,41,55,61,77,83)]=0;
# functionnal group
# engineer: hymenoptera, isoptera, oligochaeta
# carnivore: araneae, chilopoda
# omnivore: blattodea
# herbivore: gastropoda, hemiptera, lepidoptera, orthoptera 12131719
# detritivore: coleoptera, coleoptera_larva, diplopoda, isopoda 9101115
# + sp
fg=data.frame(md[,c(1:5)],engineer=md[,14]+md[,16]+md[,18],carnivore=md[,6]+md[,8],omnivore=md[,7],herbivore=md[,12]+md[,13]+md[,17]+md[,19],detritivore=md[,9]+md[,10]+md[,11]+md[,15]);
md=data.frame(md,fg[,c(6:10)]);
mdl=subset(md,strata=="Litter");
mds=subset(md,strata=="Soil");
mdtot=data.frame(mds[,c(1:4,75)],mds[,c(6:74,76:78,81:85)]+mdl[,c(6:74,76:78,81:85)]);
mdtot$shannon=diversity(mdtot[,6:74]);
mdtot$pielou= mdtot$shannon/log(mdtot$specnumber);
# sum inner replicate
mdtot2=aggregate(mdtot[,c(6:19,75:84)],by=list(state=mdtot$state,age=mdtot$age,plot=mdtot$plot),median); # create a tab with the mean of inner-replicate
mdtot2$Plantation=paste(mdtot2$state,mdtot2$age); # create a variable with gathered state and age
mdtot2$Plantation=as.factor(mdtot2$Plantation); # define it as factor
levels(mdtot2$Plantation); # give the current level of $Plantation
mdtot2$Plantation <- factor(mdtot2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mdl2=aggregate(mdl[,c(6:19,74,76:85)],by=list(state=mdl$state,age=mdl$age,plot=mdl$plot),median); # create a tab with the mean of inner-replicate
mdl2$Plantation=paste(mdl2$state,mdl2$age); # create a variable with gathered state and age
mdl2$Plantation=as.factor(mdl2$Plantation); # define it as factor
levels(mdl2$Plantation); # give the current level of $Plantation
mdl2$Plantation <- factor(mdl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mds2=aggregate(mds[,c(6:19,74,76:85)],by=list(state=mds$state,age=mds$age,plot=mds$plot),median); # create a tab with the mean of inner-replicate
mds2$Plantation=paste(mds2$state,mds2$age); # create a variable with gathered state and age
mds2$Plantation=as.factor(mds2$Plantation); # define it as factor
levels(mds2$Plantation); # give the current level of $Plantation
mds2$Plantation <- factor(mds2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mdol <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_density_0.csv");
mdol$plot=as.factor(mdol$plot);
mdol$inner_replicate=as.factor(mdol$inner_replicate);
mdol$age <- factor(mdol$age, levels=c("y","o"));
mdol$Plantation=paste(mdol$state,mdol$age);
mdol$Plantation=as.factor(mdol$Plantation);
mdol$Plantation <- factor(mdol$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mdoll=subset(mdol,strata=="Litter");
mdols=subset(mdol,strata=="Soil");
mdtol=data.frame(mdols[,c(1:4,85)],mdols[,c(21,22)]+mdoll[,c(21,22)]);
#
mdtol2=aggregate(mdtol[,6:7],by=list(state=mdtot$state,age=mdtot$age,plot=mdtot$plot),median); # create a tab with the mean of inner-replicate
mdtol2$Plantation=paste(mdtol2$state,mdtol2$age); # create a variable with gathered state and age
mdtol2$Plantation=as.factor(mdtol2$Plantation); # define it as factor
levels(mdtol2$Plantation); # give the current level of $Plantation
mdtol2$Plantation <- factor(mdtol2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mds2=aggregate(mds,by=list(state=mds$state,age=mds$age,plot=mds$plot),mean);
mds2=mds2[,-c(4:8)];
mds2$Plantation=paste(mds2$state,mds2$age); # create a variable with gathered state and age
mds2$Plantation=as.factor(mds2$Plantation); # define it as factor
levels(mds2$Plantation); # give the current level of $Plantation
mds2$Plantation <- factor(mds2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
#
mb <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_biomass_0.csv");
mb[mb$strata=="Litter",c(6:84)]=mb[mb$strata=="Litter",c(6:84)]*16;
mb[mb$strata=="Soil",c(6:84)]=mb[mb$strata=="Soil",c(6:84)]*160;

str(mb);
mb$plot=as.factor(mb$plot);
mb$inner_replicate=as.factor(mb$inner_replicate);

mb$oligochaeta=mb$oligochaeta+mb$oligochaeta_diapause; # sum oligochaeta
names(mb[c(6,12,14,22,24,80:84)]);
mb=mb[,-c(6,12,14,22,24,80:84)]; # remove mesofauna & egg & oligochaeta_diapause
mb$age <- factor(mb$age, levels=c("y","o"));
mb$Plantation=paste(mb$state,mb$age);
mb$Plantation=as.factor(mb$Plantation);
mb$Plantation <- factor(mb$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mb$sp=apply(mb[,20:74],1,sum); # sum sp
mb$biomass=apply(mb[,6:74],1,sum);
#
## functionnal group
# engineer: hymenoptera, isoptera, oligochaeta
# carnivore: araneae, chilopoda
# omnivore: blattodea
# herbivore: gastropoda, hemiptera, lepidoptera, orthoptera 12131719
# detritivore: coleoptera, coleoptera_larva, diplopoda, isopoda 9101115
fg=data.frame(mb[,c(1:5)],engineer=mb[,14]+mb[,16]+mb[,18],carnivore=mb[,6]+mb[,8],omnivore=mb[,7],herbivore=mb[,12]+mb[,13]+mb[,17]+mb[,19],detritivore=mb[,9]+mb[,10]+mb[,11]+mb[,15],sp=mb[,76]);
mb=data.frame(mb,fg[,c(6:10)]);
mbl=subset(mb,strata=="Litter");
mbs=subset(mb,strata=="Soil");
mbtot=data.frame(mbs[,c(1:4,75)],mbs[,c(6:74,76:82)]+mbl[,c(6:74,76:82)]);
#
mbl2=aggregate(mbl,by=list(state=mbl$state,age=mbl$age,plot=mbl$plot),mean);
mbl2=mbl2[,-c(4:8)];
mbl2$Plantation=paste(mbl2$state,mbl2$age); # create a variable with gathered state and age
mbl2$Plantation=as.factor(mbl2$Plantation); # define it as factor
levels(mbl2$Plantation); # give the current level of $Plantation
mbl2$Plantation <- factor(mbl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mdl2=aggregate(mdl,by=list(state=mdl$state,age=mdl$age,plot=mdl$plot),mean);
mdl2=mdl2[,-c(4:8)];
mdl2$Plantation=paste(mdl2$state,mdl2$age); # create a variable with gathered state and age
mdl2$Plantation=as.factor(mdl2$Plantation); # define it as factor
levels(mdl2$Plantation); # give the current level of $Plantation
mdl2$Plantation <- factor(mdl2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
#
mbs2=aggregate(mbs,by=list(state=mbs$state,age=mbs$age,plot=mbs$plot),mean); # create a tab with the mean of inner-replicate
mbs2=mbs2[,-c(4:8)]; # delete columns with is no more needed
mbs2$Plantation=paste(mbs2$state,mbs2$age); # create a variable with gathered state and age
mbs2$Plantation=as.factor(mbs2$Plantation); # define it as factor
levels(mbs2$Plantation); # give the current level of $Plantation
mbs2$Plantation <- factor(mbs2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#
mbol <- read.csv("~/Documents/cours/UPMC/Master/S3/stage/data/fauna_biomass_0.csv");
mbol$plot=as.factor(mbol$plot);
mbol$inner_replicate=as.factor(mbol$inner_replicate);
mbol$age <- factor(mbol$age, levels=c("y","o"));
mbol$Plantation=paste(mbol$state,mbol$age);
mbol$Plantation=as.factor(mbol$Plantation);
mbol$Plantation <- factor(mbol$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o"));
mboll=subset(mbol,strata=="Litter");
mbols=subset(mbol,strata=="Soil");
mbtol=data.frame(mbols[,c(1:4,85)],mbols[,c(21,22)]+mboll[,c(21,22)]);
#
mbtol2=aggregate(mbtol,by=list(state=mbtol$state,age=mbtol$age,plot=mbtol$plot),mean); # create a tab with the mean of inner-replicate
mbtol2=mbtol2[,-c(4:8)]; # delete columns with is no more needed
mbtol2$Plantation=paste(mbtol2$state,mbtol2$age); # create a variable with gathered state and age
mbtol2$Plantation=as.factor(mbtol2$Plantation); # define it as factor
levels(mbtol2$Plantation); # give the current level of $Plantation
mbtol2$Plantation <- factor(mbtol2$Plantation, levels=c("N o","R1 y","R1 o","R2 y","R2 o")); # define the right order of each Plantation
#

# ----
require(ade4)
# e.g.
#dudi.pca(Table1)->acp1
#dudi.pca(Table2)->acp2
#dudi.pca(Table3)->acp3
#list<-ktab.list.dudi(list(acp1, acp2, acp3), tabnames=c("table1", "table2", "table3"))
#mcoa(list)->MCOA
#plot(MCOA)

macro<-dudi.pca(mdtot[,c(78:82)],scannf=F,nf=2)
macro$eig/sum(macro$eig)
chem<-dudi.pca(ch05[,c(15,17:22)],scannf=F,nf=2);
chem$eig/sum(chem$eig)
all.equal(macro$lw,chem$lw)
coi=coinertia(macro,chem,scan=F,nf=2)
coi
coi$eig[1]/sum(coi$eig)
summary(coi)
# RV coinertia ratio
randtest(coi,nrepet=999)
x11();plot(coi)
list<-ktab.list.dudi(list(macro, chem), tabnames=c("macrofauna", "chemical"))
mcoa(list)->MCOA
x11();plot(MCOA)

#### covington ----
x0=c(0:30)
y0=-5.25*x0^(1.24)*exp(-0.0649*x0^1.063)+86.75
y0[30]
x1=c(31:60)
y1=-5.25*x1^(1.24)*exp(-0.0649*x1^1.063)+53.40873
x2=c(61:90)
y1[30]
y2=-5.25*x2^(1.24)*exp(-0.0649*x2^1.063)+47.96021
x=c(x0,x1,x2)
y=c(y0,y1,y2)
a=data.frame(x,y)

ggplot(a,aes(x=x,y=y))+geom_line()+xlim(0,90)+ylim(0,100)

x=c(0:86)
y0=-5.25*x[c(0:29)]^(1.24)*exp(-0.0800*x[c(0:29)]^1.063)+100
y0[29]
y1=-5.25*x[c(0:29)]^(1.24)*exp(-0.0800*x[c(0:29)]^1.063)+79.36643
y1[29]
y2=-5.25*x[c(0:29)]^(1.24)*exp(-0.0800*x[c(0:29)]^1.063)+58.73286
y2[y2<0]=0
y=c(y0,y1,y2)
a=data.frame(x,y)
p2=ggplot(a,aes(x=x,y=y))+geom_line()+xlim(0,100)+ylim(0,100)+geom_vline(aes(xintercept=29), linetype="dashed")+geom_vline(aes(xintercept=58), linetype="dashed")

x=c(0:86)
y=-5.25*x^(1.24)*exp(-0.0800*x^1.063)+100
a=data.frame(x,y)
p1=ggplot(a,aes(x=x,y=y))+geom_line()+xlim(0,100)+ylim(0,100)

grid.arrange(p1,p2,nrow=1);

no=mean(ch05[ch05$Plantation=="N o",]$OM_.)
r1o=mean(ch05[ch05$Plantation=="R1 o",]$OM_.)
r2o=mean(ch05[ch05$Plantation=="R2 o",]$OM_.)
