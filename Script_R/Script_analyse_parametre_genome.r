#Script permettant de realiser differents graphiques d'analyser de different parametres de binning mais sur les génomes#

setwd('C:\\Users\\Robin\\Desktop\\pasteur\\3c')

#J'appelle mes differentes librairies"
library('ggplot2')
library("ggrepel")
library("dplyr")
library("ggpubr")

#Graphique global pour voir nos genomes complet non contamines #

comple_vs_conta <- read.csv(file ="binning/metator_20_80/miComplete2.csv", sep ="\t", header=T) #Je prend le tableau de sortie de miComplete

comple_vs_conta$Completeness <- df$Completeness*100 #Je transforme la completion en %
comple_vs_conta$Redundancy <- (df$Redundancy*100)-100 #Je transforme la contamination en %

#Je plot la contamination en fonction de la completion
ggplot(comple_vs_conta, aes(x = Redundancy, y= Completeness)) +
  geom_point()+
  xlab("% de contamination")+
  ylab("Completion du genome")+
  ggtitle("Completion et Contamination de nos 25 genomes avec miComplete")+
  theme_minimal()+  
  geom_vline(xintercept=10, col="red") + #Seuil vertical
  geom_hline(yintercept=90, col="red") #Seuil horizontal

#Je vais faire des boxplot pour voir l'impact des iteration sur la taille, la completion et la conta#

df_1 <- read.csv(file ="binning/metator_20_80/miComplete2_1_80.csv", sep ="\t", header=T) #Rapport miComplet pour -O 80 et -i 1
df_5 <- read.csv(file ="binning/metator_20_80/miComplete2_5_80.csv", sep ="\t", header=T) #Rapport miComplet pour -O 80 et -i 5
df_20 <- read.csv(file ="binning/metator_20_80/miComplete2_20_80.csv", sep ="\t", header=T) #Rapport miComplet pour -O 80 et -i 20
df_50 <- read.csv(file ="binning/metator_20_80/miComplete2_50_80.csv", sep ="\t", header=T) #Rapport miComplet pour -O 80 et -i 50

#Jutilise dyplyr pour coller vers le base les differente dataframe pour grouper nos datas dans l'ordre
summary <- bind_rows(df_1,df_5)
summary <- bind_rows(summary,df_20)
summary <- bind_rows(summary,df_50)

#je cree un vecteur numerique pour les iteration avec le nombre de genome dans chaque fichier
iter<- c((rep(c(1), times=32)),(rep(c(5), times=34)),(rep(c(20), times=33)),(rep(c(50), times=35)))
summary$iteration <- as.factor(iter) #je creer colonne interation et stock le vecteur

summary$Completeness <- summary$Completeness*100 #Je transforme la completion en %
summary$Redundancy <- (summary$Redundancy*100)-100 #Je transforme la contamination en %


######################################
#je fais les trois graphique (taille, completion et conta) en fonction des differentes iterations en boxplot
lenght <-ggplot(summary, aes(x= iteration, y= Length, color=iteration, group=iteration)) +
  geom_boxplot()+
  theme_minimal()

completeness <-ggplot(summary, aes(x= iteration, y= Completeness, color=iteration, group=iteration)) +
  geom_boxplot()+
  theme_minimal()

contamination <-ggplot(summary, aes(x= iteration, y= Redundancy, color=iteration, group=iteration)) +
  geom_boxplot()+
  theme_minimal()

######################################
##ggarange permet d'arranger differents plot en une figure et permet de faire de tres belle figure article-ready
figure <-ggarrange(lenght,completeness,contamination, ncol=3, nrow=1, common.legend = TRUE, legend = "right", align="hv")

#annotate_figure permet de donner un beau titre à notre figure globale
annotate_figure(figure, top = text_grob("La taille, la complÃ©tion et la contamination en fonction des itÃ©ration 1 ,5 , 20 et 50 ", 
                                        color = "red", face = "bold", size = 10))
