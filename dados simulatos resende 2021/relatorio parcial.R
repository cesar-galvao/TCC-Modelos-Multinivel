# Pacotes ----
if(!require("pacman")){install.packages("pacman")}
pacman::p_load(tidyverse, scales, sommer, car, corrplot, lme4, magrittr)

# Pedigree ----

ped <- read.csv("pedigree.txt", sep = "\t")

# quantos genotipos
ped$id |> length() #sao 100

# maes
ped %>%
  group_by(dam) %>%
  count()%>%
  ungroup()%>%
  mutate(conta = cumsum(n),
         sex = "F")

#pais
ped %>%
  group_by(sire) %>%
  count() %>%
  ungroup()%>%
  mutate(conta = cumsum(n),
         sex = "M")

# casais
# ped %>%
#   group_by(sire, dam) %>%
#   count() %>%
#   arrange(sire)

### -------------- COMENTARIOS ----------------- ####
# Os primeiros 20 genótipos são pais ou maes.
# Quantidade de filiações é desbalanceada quanto a progenitores.
# São 20 progenitores e 80 F1
### ------------------------------------------- ####


# Fenotipo ----
feno <- read.table("pheno.txt")

# trials
trials <- feno %>%
  group_by(Trial) %>%
  summarise(
    n = n(),
    mean_yield = mean(Yield),#media de cada grupo
    sigma_i = var(Yield)
  ) %>%
  ungroup() %>%
  mutate(mu = mean(mean_yield), #media geral
         LI = mean_yield - qnorm(.975)*sqrt(sigma_i),#tem n = 100, da para usar normal ja
         LS = mean_yield + qnorm(.975)*sqrt(sigma_i)) %>%
  pivot_longer(cols = c(LS, LI), names_to = "IC")

# #distribuicao das medias e media geral
# ggplot(trials, aes(Trial, mean_yield))+
#   geom_point()+
#   geom_abline(intercept = trials$mu, slope = 0, color = "red")

#IC de cada trial, ordenado, comparados à média geral
ggplot(trials, aes(reorder(Trial, mean_yield), mean_yield))+
  geom_point()+
  geom_line(data = trials, aes(Trial,value, group = Trial))+
  geom_abline(intercept = trials$mu, slope = 0, color = "red", linetype = "dashed")+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x = "Célula experimental (Trial)", y = "Produção (Yield)")


# ggplot(feno, aes(Yield, group = Trial))+
#   geom_boxplot()+
#   coord_flip()

# nao sei se vale a pena usar boxplot


### -------------- COMENTARIOS ----------------- ####
# Experimentos são balanceados, todos os 100 genotipos ocorrem em cada trial.
# Média dos trials parecem estar distribuidas em torno da média, mas com variancias diferentes
# TESTAR VARIANIAS DIFERENTES
# Com um IC para cada Trial, nenhum deles aparenta ser significativamente diferente da média
# 50 trials, 100 replicacoes em cada, 5.000 experimentos
### ------------------------------------------- ####

# Famílias ----

mu_y <- mean(feno$Yield)

familia <- feno %>%
  inner_join(ped, by = c("Genotype" = "id")) %>% 
  filter(!is.na(dam)) %>%
  group_by(dam,sire) %>%
  mutate(
    mean_yield = mean(Yield),
    n = n(),
    sigma_i = var(Yield),
    LI = mean_yield - qnorm(.975)*sqrt(sigma_i),
    LS = mean_yield + qnorm(.975)*sqrt(sigma_i) 
  ) %>%
  pivot_longer(cols = c(LS, LI), names_to = "IC")%>%
  ungroup() %>%
  mutate(prog = paste(dam, sire, sep = " x "))

ggplot(familia, aes(reorder(prog, mean_yield), mean_yield))+
  geom_point()+
  geom_line(data = familia, aes(prog,value, group = prog))+
  geom_abline(intercept = mu_y, slope = 0, color = "red", linetype = "dashed")+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x = "Progenitores", y = "Produção (Yield)")

# testando diferença entre as variâncias ----
leveneTest(Yield ~ prog, data = familia)
#TukeyHSD(aov(Yield ~ prog, data = familia))
pairwise.t.test(familia$Yield, familia$prog, p.adjust.method = "bonferroni")

### -------------- COMENTARIOS ----------------- ####
# Alguns casais parecem produzir filhos com maior producao, mas nenhum realmente diferente
## da média
# Teste para diferença de variâncias não parece bom. Talvez por causa de uma amostra grande?
### ------------------------------------------- ####
  

# Ambiente ----
amb <- read.table("envdat.txt") %>% 
  mutate(experimento = ifelse(is.na(Trial), 0, 1),
         text = ifelse(!is.na(Trial),str_sub(Trial, start = 3), ""))

# mapa de onde teve trial
amb %>% 
  ggplot(aes(x = Longitude, y = Latitude, 
             fill = factor(experimento,
                           levels = c(0,1),
                           labels = c("Sem experimento", "Com experimento")))) +
  geom_tile(color = "lightgray") +
  scale_fill_manual(values = c("white", "black")) +
  geom_text(label = amb$text, vjust = -.5)+
  theme_minimal()+
  theme(axis.line = element_blank(),
        panel.grid = element_blank())+
  labs(fill = "Experimento")

# ilustrar heat map de uma das variaveis
ggplot(amb) +
  geom_tile(aes(Longitude, Latitude, fill = EA005)) +
  scale_fill_gradientn(colours = rev(jet.colors(7))) +
  geom_point(data = amb[!is.na(amb$Trial),],aes(Longitude,Latitude), size = 2,
             shape = 3)+
  theme_minimal()+
  theme(axis.line = element_blank(),
        panel.grid = element_blank())

# clusterizar e ver média de yield por região 
matriz <- amb %>%
  dplyr::select(starts_with("EA")) %>%
  as.matrix() %>%
  scale()


cluster_table <- tibble(cent = 1:20)%>%
  mutate(withinss = map_vec(cent, function(x){kmeans(matriz, centers = x)$tot.withinss}))

cluster_table %>%
  ggplot(aes(cent, withinss))+
  geom_point()+
  theme_bw()


# gera clusters
amb <- amb %>%
  mutate(cluster = factor((kmeans(matriz, centers = 5))$cluster))
rm(matriz)

# plot dos clusters
ggplot(amb) +
  geom_tile(aes(Longitude, Latitude, fill = cluster)) +
  scale_fill_brewer()+
  labs(fill = "Cluster")+
  theme_minimal()+
  theme(axis.line = element_blank(),
        panel.grid = element_blank())

# dataset de yield e clusters para avaliar IC  
desc_clusters <- amb %>%
  dplyr::select(Trial, cluster) %>%
  inner_join(feno) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_yield = mean(Yield),
    sigma_i = var(Yield)
  ) %>%
  ungroup() %>%
  mutate(mu = mean(mean_yield), #tem n = 100, da para usar normal ja
         LI = mean_yield - qnorm(.975)*sqrt(sigma_i),
         LS = mean_yield + qnorm(.975)*sqrt(sigma_i)) %>%
  pivot_longer(cols = c(LS, LI), names_to = "IC")


#IC de cada trial, ordenado, comparados à média geral
ggplot(desc_clusters, aes(reorder(cluster, mean_yield), mean_yield))+
  geom_point()+
  geom_line(data = desc_clusters, aes(cluster,value, group = cluster))+
  geom_abline(intercept = desc_clusters$mu, slope = 0, color = "red", linetype = "dashed")+ 
  theme_bw()


### -------------- COMENTARIOS ----------------- ####
# Comportamento dos clusters varia muito de acordo com número (...)
### ------------------------------------------- ####


# Modelo linear não hierarquizado ----
data_fit_lm <- amb %>%
  dplyr::select(Trial, starts_with("EA")) %>% #seleciona trial e vars ambientais
  filter(!is.na(Trial))%>% #filtra apenas para casos em que Trial nao é NA
  inner_join(feno) %>% #junta com o banco que 
  dplyr::select(Genotype, Yield, starts_with("EA")) %>%
  mutate(across(starts_with("EA"), ~scale(.)))

# um modelo linear por genótipo
lm_fits <- data_fit_lm %>%
  group_by(Genotype) %>%
  nest(data = c(Yield, EA001:EA100)) %>%
  mutate(lm_result = map(data, ~lm(Yield ~ (.), data = .)))

rm(data_fit_lm, lm_fits)

### -------------- COMENTARIOS ----------------- ####
# Não funciona porque acusa sigularidades
### ------------------------------------------- ####

data_fit_lm2 <- amb %>%
  dplyr::select(Trial, starts_with("EA")) %>% #seleciona trial e vars ambientais
  filter(!is.na(Trial))%>% #filtra apenas para casos em que Trial nao é NA
  inner_join(familia) %>% #junta com o banco que 
  dplyr::select(prog, Yield, starts_with("EA")) %>%
  mutate(across(starts_with("EA"), ~scale(.)))

lm_fits2 <- data_fit_lm2 %>%
  group_by(prog) %>%
  nest(data = c(Yield, EA001:EA100)) %>%
  mutate(lm_result = map(data, ~lm(Yield ~ (.), data = .)))

rm(data_fit_lm2, lm_fits2)

### -------------- COMENTARIOS ----------------- ####
# Não funciona porque acusa sigularidades
# uma saída pode ser justamente o bootstrap do Tassinari
### ------------------------------------------- ####

matriz <- amb %>%
  dplyr::select(starts_with("EA")) %>%
  as.matrix() %>%
  scale()

colnames(matriz) = NULL

cor_matrix <- cor(matriz[,])

corrplot(cor_matrix, method = "circle", type = "lower", order = "hclust", tl.col = "transparent")

# Modelo linear hierarquizado sem marcador ambientomico ----

# sem a criação de marcadores ambientomicos

data_fit_lmer <- amb %>%
  dplyr::select(Trial, starts_with("EA")) %>% #seleciona trial e vars ambientais
  filter(!is.na(Trial))%>% #filtra apenas para casos em que Trial nao é NA
  inner_join(feno) %>% #junta com o banco que 
  dplyr::select(Genotype, Yield, starts_with("EA")) %>%
  mutate(across(starts_with("EA"), ~ as.vector(scale(.))))

formula_lmer <- paste0("EA", str_pad(1:100, width = 3, side = "left", pad = "0")) %>%
  str_c(.,collapse = " + ") %>%
  paste0("Yield ~ (1 + ", ., "|Genotype)") %>%
  formula()

hfit <- lmer(formula_lmer, data = data_fit_lmer)

rm()

### -------------- COMENTARIOS ----------------- ####
# sem estrutura de covariancia

# Sem transformar as EA em ENVMrkr nao funciona. Tem mais
# efeito aleatorio que observacao.
### ------------------------------------------- ####

# Modelo linear hierarquizado com marcador ambientomico ----

# montagem das bases usadas
# target population of environments com covars ambientais, tipo e trial
tpe <- read.table("envdat.txt") %>%
  mutate(type = factor(
              ifelse(!is.na(Trial),"phen", "area")), #cria a variavel type
         across(starts_with("EA"), ~ as.vector(scale(.))) # padroniza as EA
    ) %>%
  dplyr::select(type, everything())

# dados das células que contém experimento
phe <- tpe %>%
  filter(!is.na(Trial)) %>%
  dplyr::select(-type) %>%
  inner_join(feno) %>% 
  dplyr::select(Genotype, Yield, Trial, Longitude, Latitude, everything())
  


#essa tabela receberá os marcadores ambientômicos gerados no loop/bootstrap

ENVMkrs <- tpe %>% 
  dplyr::select(type, Trial, Longitude, Latitude) %>%
  mutate(
    #cria os vetores vazios para reservar memoria e acelerar o bootstrap
    as.data.frame(replicate(150, numeric(10000)), col.names = ENVnames)
    )

#renomeia as demais variaveis
names(ENVMkrs)[-c(1:4)]<-paste("ENV",
                               str_pad(1:150, width = 3, side = "left", pad = "0"),
                               sep="")

#lista de valores únicos de genótipo
genos <- phe$Genotype %>% unique()

#vetor vazio para aceitar os valores de fitness
fitness<-numeric(150)
names(fitness)<-ENVnames

### Construcao dos ENVMRKR - FOR LOOP ----

### -------------- COMENTARIOS ----------------- ####
# DA PRA PARALELIZAR
### ------------------------------------------- ####

for(i in 1:150){ 
  #sorteando um número de amostras:
  nsamp <- sample(2:10,1);nsamp #sorteando de 2 a 10 genotipos
  
  #subset contendo apenas os genotipos sorteados no passo anterior:
  boots <- phe %>%
    filter(Genotype %in% sample(genos,nsamp)) %>%
    droplevels()
  
  # montagem das amostras de treino e de validacao
  poptr<-sample(1:nrow(boots),nrow(boots)/2) %>%
                sort() #amostra de individuos com metade do tamanho do boot p/ treino
  popvl<-(1:nrow(boots))[-poptr] #resto da amostra para validacao
  
  # em dados de boots TREINO, usa amostra de treino e vars ambientais para treino
  treino<-lm(Yield ~ ., data = droplevels(boots[poptr,c(2, #Yield
                                                        6:105)])) #EA
  
  #rgg é a correlacao entre valores preditos na amostra na amostra de validacao e seus valores reais
  rgg<-cor(
    predict(treino,boots[popvl,c(2,6:105)]), #predict do modelo na particao de TREINO
           boots[popvl,]$Yield)
  
  fitness[i]<-rgg
  
  #ENVFIT é um fit para todos os casos do bootstrap, sem dividir entre amostra de treino e val.
  ENVfit<-lm(Yield ~ ., data=droplevels(boots[,c(2,6:105)]))
  
  #R é a correlação entre os valores ajustados do ENVFIT e os valores reais do bootstrap
  R<-cor(predict(ENVfit),boots$Yield)
  
  
  message("process= ",i,", nsamp= ",nsamp,", R= ",round(R,3),", rgg= ",round(rgg,3))  
  
  #Usando o modelo gerado do bootstrap, gerar previsões para todo o grid
  #do TPE. Esse novo predict são os valores do marcador para todos os 10.000 pontos
  ENVMkr <- predict(ENVfit,tpe[,-(1:4)])
  
  # Atribui os valores do marcador à variavel correspondente
  ENVMkrs[,i+4] <- ENVMkr
}



#rgg é considerado uma medida de goodness of fit
#é a correlacao entre os fitted values do modelo (modelo de treino sobre 
#os dados de validacao) com os dados reais da validacao
hist(fitness,breaks=50); 
abline(v=0,lty=2,lwd=3,col="red")


#padroniza os marcadores ambientais
ENVMkrs[,ENVnames]<-apply(ENVMkrs[,ENVnames],2,padr)


# Selecao dos marcadores ambientomicos com rgg superior a 0.5
ENVsel<-names(fitness[fitness>0.5]);fitness[fitness>0.5]

### -------------- COMENTARIOS ----------------- ####
# Usando o bootstrap proposto pelos autores
# Procurar fundamentacao teorica
### ------------------------------------------- ####

## Modelo sem estrutura de covariancia ----

# NOTAS ----
# Experimentar composição de ENVMRKR com análise fatorial