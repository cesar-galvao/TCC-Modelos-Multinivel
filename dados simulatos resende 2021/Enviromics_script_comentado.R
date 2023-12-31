#devtools::install_github("variani/lme4qtl")
library(lme4); library(lme4qtl); 
library(sommer)
library(randomForest);#library(nnet)
library(ggplot2)
library(reshape2); library(scales)
library(nadiv) #para lidar com pedigrees
library(psych) #para funcao de traco de matrizes
library(dplyr);library(tidyr)
#rm(list=ls())

#FUNÇÃO PARA PADRONIZAR VARIÁVEIS(μ=0,σ²=1):
padr <- function(x){(x-mean(x))/sd(x)}

##Lendo dados ambientípicos:
envdat<-read.table("~/r_projects/curso_ambientomica/envdat.txt",header=TRUE)
dim(envdat) #10000 pixels * 100 covars
envdat[1:5,1:5]

ggplot(envdat) + #dando uma olhada em uma cov ambiental
  geom_tile(aes(Longitude, Latitude, fill = EA003)) +
  scale_fill_gradientn(colours = rev(jet.colors(7))) +
  geom_point(data = envdat[!is.na(envdat$Trial),],aes(Longitude,Latitude), size = 2) 

##Lendo dados fenotípicos:
pheno<-read.table("~/r_projects/curso_ambientomica/pheno.txt",header=TRUE); dim(pheno)
pheno$Trial <- as.factor(pheno$Trial)
pheno$Genotype <- as.factor(pheno$Genotype)
head(pheno);hist(pheno$Yield,breaks=100)#BLUE

genos <- levels(pheno$Genotype)#salvando codigos de genotipos

##Lendo dados pedigree:
ped<-read.table("~/r_projects/curso_ambientomica/pedigree.txt",header=TRUE);
head(ped); dim(ped)

##########-----------------------##########

##Filtrando dados ambientípicos para os fenotípicos
env <- envdat[envdat$Trial%in%levels(pheno$Trial),] #filtra so os 50 casos com trial
rownames(env)<-env$Trial; 
phe <- data.frame(pheno[,-2],env[pheno$Trial,]); rownames(phe) <- NULL #join com dados de fenotipo
phe[1:5,1:10];dim(phe)

par(mfrow=c(2,1))#Verificando representatividade
hist(env$EA003,xlim=0:1);hist(envdat$EA003,xlim=0:1)
dev.off()

# Geração de múltiplos marcadores ambientômicos: ----

#target population of environments
# sem yield, indica se a area foi fenotipada ou nao
tpe <- rbind(data.frame("type"="phen",env),
             data.frame("type"="area",envdat))
tpe$type<-as.factor(tpe$type)
tpe[,5:104]<-apply(tpe[,5:104],2,padr)
#confere distribuição da variavel EA003 em area fenotipada e nao fenotipada
ggplot(tpe, aes(x=type, y=EA003))+
  geom_boxplot()#confere EA

#desejavel que o boxplot de uma variavel ambiental
#nos pixels que tem fenótipo esteja dentro do boxplot da area

fitness<-c(); ENVMkrs <- tpe[,1:4]
for(i in 1:150){ 
  #sorteando um número de amostras:
  nsamp <- sample(2:10,1);nsamp #sorteando de 2 a 10 genotipos
  #subset contendo apenas os genotipos sorteados no passo anterior:
  # a partir da tabela phe
  boots<-droplevels(phe[phe$Genotype%in%sort(sample(genos,nsamp)),])
  
  ##Ajuste Linear:
  poptr<-sort(sample(1:nrow(boots),nrow(boots)/2)) #amostra com metade do tamanho do boot p/ treino
  popvl<-(1:nrow(boots))[-poptr] #resto da amostra para validacao
  
  # em boots, usa amostra de treino e vars ambientais para treino
  treino<-lm(Yield ~ ., 
             data=droplevels(boots[poptr,c(2,6:105)])) 
  #rgg é a correlacao entre valores preditos na amostra na amostra de validacao e seus valores reais
  rgg<-cor(predict(treino,boots[popvl,c(2,6:105)]),boots[popvl,]$Yield)
  fitness<-c(fitness,rgg)
  
  #ENVFIT é um fit para todos os casos do bootstrap, sem dividir entre amostra de treino e val.
  ENVfit<-lm(Yield ~ ., 
               data=droplevels(boots[,c(2,6:105)]))
  #R é a correlação entre os valores ajustados do ENVFIT e os valores reais do bootstrap
  R<-cor(predict(ENVfit),boots$Yield)
  
  ##Ajuste não-linear:
  #ENVfit<-randomForest(Yield ~ ., size=5,
  #             data=droplevels(boots[,c(2,6:105)]))
  #R<-cor(predict(ENVfit),boots$Yield);R

  message("process= ",i,", nsamp= ",nsamp,", R= ",round(R,3),", rgg= ",round(rgg,3))  
  
  #Usando o modelo gerado do bootstrap, gerar previsões para todo o grid  
  ENVMkr <- predict(ENVfit,tpe[,-(1:4)]); #ENVMkr[ENVMkr<0]<-0; 
  #ENVMkr<-rescale(ENVMkr); #hist(ENVMkr)
  
  ENVMkrs <- data.frame(ENVMkrs,ENVMkr)
};
ENVnames<-paste("ENV",sprintf("%03d", 1:150),sep="")#se alterar a quantidade de marcadores, alterar o 150 aqui
colnames(ENVMkrs)[5:154] <- ENVnames#se alterar a quantidade de marcadores, alterar o 154 aqui
names(fitness)<-ENVnames

#rgg é considerado uma medida de goodness of fit
hist(fitness,breaks=50); abline(v=0,lty=2,lwd=3,col="red")
ENVMkrs[45:55,1:10]; dim(ENVMkrs) #conferida rapida
ENVMkrs[,ENVnames]<-apply(ENVMkrs[,ENVnames],2,padr)


# Selecao dos marcadores ambientomicos com rgg superior a 0.5
ENVsel<-names(fitness[fitness>0.5]);fitness[fitness>0.5]
#ENVMkrs<-ENVMkrs[,c("type", "Trial","Longitude","Latitude",names(fitness[fitness>0.4]))]



ggplot(ENVMkrs, aes(x=type, y=ENV005))+
  geom_boxplot()#confere ENV



##Criando Kernels:

# Resgata a matriz resultado de todos os fit de bootstrap
# 

W <- as.matrix(ENVMkrs[ENVMkrs$type=="phen",ENVsel])
hist(cor(W)[upper.tri(cor(W))],breaks=100)
LDenv<- mean(cor(W)[upper.tri(cor(W))]); LDenv

Ecor <- as.matrix(cor(t(W)))# matriz correlacao entre ambientes/trial
E <- (W%*%t(W)) / (tr(W%*%t(W))/nrow(W))
plot(Ecor[upper.tri(Ecor)],E[upper.tri(E)])

heatmap(E,cexRow = 0.5, cexCol = 0.5)

# matriz A é a que dá a estrutura de covariancia dos ensaios
A <- as.matrix(makeA(prepPed(ped)))
heatmap(A,cexRow = 0.5, cexCol = 0.5)

#not run: Muito pesado!
#modG<-mmer(Yield ~ 1,random=~Genotype, data=phe)
#modE<-mmer(Yield ~ 1,random=~vs(Trial, Gu=E), data=phe)

#desbalanceando os dados (20% = 1000/5000)
set.seed(2023);fdesb<-sort(sample(1:5000,1000))
reduced_phe<-droplevels(phe[fdesb,])#Enxugando os dados

#Conferindo balanceamento dos Genotipos:
apply(table(reduced_phe$Genotype,reduced_phe$Trial),1,sum); message("ngen= ",length(apply(table(reduced_phe$Genotype,reduced_phe$Trial),1,sum))," | OK!")
#Conferindo balanceamento dos Trials:
apply(table(reduced_phe$Genotype,reduced_phe$Trial),2,sum); message("ntrl= ",length(apply(table(reduced_phe$Genotype,reduced_phe$Trial),2,sum))," | OK!")

#system.time(modG0<-mmer(Yield ~ Trial,random=~ vsr(Genotype,Gu=A), data=reduced_phe))
#system.time(modG1<-mmer(Yield ~ Trial,random=~ vsr(Genotype,Gu=A), data=phe, init=modG0$sigma_scaled))
#summary(modG0); #summary(modG1)
#EGV <- unlist(modG0$U)+modG0$Beta[1:1,]$Estimate; #EGV = estimated genetic values
#media<-aggregate(Yield~Genotype,FUN=mean,phe)$Yield
#plot(media,EGV);abline(0,1,lty=2,col="gray")

##Não usamos as covariáveis ambientais, e sim os marcadores:
dat <- phe[,-(6:106)]
#usa apenas os marcadores selecionados no bootstrap
dat <- data.frame(dat,apply(ENVMkrs[phe$Trial,ENVsel],2,padr))
dat <- droplevels(dat[fdesb,])#Desbalanceando os dados

############
# MODELOS ---- 
############

vldat <- phe[,-(6:106)]
vldat<-droplevels(vldat[!(1:5000%in%fdesb),])#Desbalanceando os dados
dim(vldat)

vlmod0<-lmer(Yield ~ 1+ (1|Genotype),vldat)
vlblup0<-ranef(vlmod0)$Genotype+fixef(vlmod0)

#####MODELO BASICO
#MOD0 <- lmer(Yield ~ 1 + (1|Genotype)+(1|Genotype:Trial),phe)  
# MOD0 <- lmer(Yield ~ 1 + (1|Genotype),dat)
# summary(MOD0)
# blup0<-ranef(MOD0)$Genotype+fixef(MOD0) #EGV = estimated genetic values
# plot(unlist(vlblup0),unlist(blup0));abline(0,1,lty=2,col="darkgrey")
# cor(unlist(vlblup0),unlist(blup0))
# 
# CP_mod0<-c()
# for(i in unique(dat$Trial)){
#   temp<-vldat[vldat$Trial%in%i,]
#   gcorresp<-intersect(temp$Genotype,rownames(blup0))
#   CP_mod0<-c(CP_mod0,cor(temp[temp$Genotype%in%gcorresp,]$Yield,blup0[gcorresp,]))
# }

#####MODELO KERNEL
#EBLUP - Kernel ambientomico(compatível com GBLUP)
### ISSO ENTREGA OS =VALORES PREDITOS= POR GENÓTIPO EM CADA AMBIENTE
# MOD1<-mmer(Yield ~ 1,
#              random=~vsr(Genotype,Trial, Gu=E),
#              data=dat)
# load("MOD1.RData")
# summary(MOD1)
# 
# blup1<-melt(unlist(MOD1$U));head(blup1)
# blup1$factors<-rownames(blup1)
# blup1<-as.data.frame(blup1%>%separate_wider_delim(factors,":Trial.Yield.",names=c("Genotype","Trial")))
# blup1<-blup1[order(blup1$Trial),]
# 
# CP_mod1<-c()
# for(i in unique(dat$Trial)){
#   temp1<-vldat[vldat$Trial%in%i,]
#   temp2<-blup1[blup1$Trial%in%i,]
#   gcorresp<-intersect(temp1$Genotype,temp2$Genotype)
#   CP_mod1<-c(CP_mod1,cor(temp1[temp1$Genotype%in%gcorresp,]$Yield,temp2[temp2$Genotype%in%gcorresp,]$value))
# }

#"RRBLUP" Ambientomico
### ISSO ENTREGA OS =EFEITOS= POR GENÓTIPO EM CADA AMBIENTE

ENVsel
rrmod<-formula(paste("Yield ~ 1 + (",paste(ENVsel,collapse = "+"),"|Genotype)"))
rrmod

MOD2 <- lmer(rrmod,dat)
ranef(MOD2)

blup2<-data.frame()
for(g in genos){
  blup2<-rbind(blup2,data.frame("Trial"=rownames(ENVMkrs[1:50,ENVsel]),"Genotype"=g,ENVMkrs[1:50,ENVsel]))
}
  
blup2$value<-predict(MOD2, blup2)

CP_mod2<-c()
for(i in unique(dat$Trial)){
  temp1<-vldat[vldat$Trial%in%i,]
  temp2<-blup2[blup2$Trial%in%i,]
  gcorresp<-intersect(temp1$Genotype,temp2$Genotype)
  CP_mod2<-c(CP_mod2,cor(temp1[temp1$Genotype%in%gcorresp,]$Yield,temp2[temp2$Genotype%in%gcorresp,]$value))
}

data.frame("Modelo"=c("Basico","Kernel","RandReg"),
           "CP"=round(c(mean(CP_mod0),mean(CP_mod1),mean(CP_mod2)),3))

#not run:
#paste(paste("ENV",sprintf("%03d", 1:100),sep=""),collapse="+")

MOD2A<-relmatLmer(Yield~1+(ENV001+ENV002+ENV003|Genotype),dat,relmat=list(Genotype=A))
ranef(MOD2A)

## Cross-validation 4-fold ##
#not run...

reduced_phe <- data.frame("Fold"=NA,reduced_phe)
reduced_phe[reduced_phe$Longitude%in%1:50&reduced_phe$Latitude%in%1:50,]$Fold <- "C"
reduced_phe[reduced_phe$Longitude%in%1:50&reduced_phe$Latitude%in%51:100,]$Fold <- "A"
reduced_phe[reduced_phe$Longitude%in%51:100&reduced_phe$Latitude%in%1:50,]$Fold <- "D"
reduced_phe[reduced_phe$Longitude%in%51:100&reduced_phe$Latitude%in%51:100,]$Fold <- "B"

ggplot(reduced_phe,aes(Longitude, Latitude,label=Fold)) + 
  geom_text()+theme_bw()+
  geom_hline(yintercept=50,lty=2,color="blue") + 
  geom_vline(xintercept=50,lty=2,color="blue")
