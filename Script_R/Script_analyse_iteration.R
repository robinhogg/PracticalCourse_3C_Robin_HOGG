#Script permettant de realiser differents graphiques d'analyser de different parametres de binning#

setwd('C:\\Users\\Robin\\Desktop\\pasteur\\3c')

#J'appelle mes differentes librairies"
library('ggplot2')
library("ggrepel")
library('ggpubr')

#Graphique des tests de differentes iterations et overlapping#

inte_over <- read.csv(file = "temp.csv", sep = "\t", header= F) #Je prend le tableau et nomme les colonnes
colnames(inte_over)=c("bins", "iteration", "overlapping")

#Je plot le nb de bins >100Kb en fonction des différents iterations et en fonction (couleur) des overlapping
ggplot(inte_over, aes(x= iteration, y=bins, group=overlapping)) + 
  geom_line(position=position_dodge(width =1), aes(color=overlapping)) + 
  geom_point(position=position_dodge(width =1), aes(color=overlapping)) +
  theme_classic()+
  xlab("Nombre d'iteration")+
  ylab("nb de bins >100Kb")+
  labs(fill = "Overlaping")+
  ggtitle("Graphique recapitulatif de nos tests avec differentes iteration et overlaping en montrant les bins > 100Kb")+
  ylim(0,42)


#Graphique global des differents intervals #

####################################################################
d10_x_100 <- read.csv(file = "temp_1.csv", sep = "\t", header= F) #Table 10kb>x<100Kb, somme des tailles bins, it {1-50}
colnames(d10_x_100)=c("sum", "it", "over")


d10_x_100$it <- factor(d10_x_100$it,      #Je met en facteur et defini des levels pour les iterations                             
                 levels = c("i1", "i2","i3","i4","i5","i10","i20","i30","i40","i50"))

#je stock le plot du contenue en sequence en fonction des differentes interation dans une variable
d10_x_100_plot <- ggplot(d10_x_100, aes(x=it, y=sum, color=it)) + 
  geom_bar(stat = "identity")+
  xlab("10kb_x_100kb")+
  ylab("Contenu en sequence (bp)")+
  labs(fill ="Nombre d'iterations")+
  ylim(0,180000000)+
  theme(axis.text.x = element_text(angle = 90))

####################################################################
d100_x_500 <- read.csv(file = "temp_2.csv", sep = "\t", header= F) #Table 100kb>x<500Kb, somme des tailles bins, it {1-50}
colnames(d100_x_500)=c("sum", "it", "over")

d100_x_500$it <- factor(d100_x_500$it,        #Je met en facteur et defini des levels pour les iterations                           
                 levels = c("i1", "i2","i3","i4","i5","i10","i20","i30","i40","i50"))

#je stock le plot du contenue en sequence en fonction des differentes interation dans une variable
d100_x_500_plot <- ggplot(d100_x_500, aes(x=it, y=sum, color=it)) + 
  geom_bar(stat = "identity")+
  xlab("100kb_x_500kb")+
  ylim(0,180000000)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90))

####################################################################
d500_x_1000 <- read.csv(file = "temp_3.csv", sep = "\t", header= F) #Table 500kb>x<1000Kb, somme des tailles bins, it {1-50}
colnames(d500_x_1000)=c("sum", "it", "over")

d500_x_1000$it <- factor(d500_x_1000$it,      #Je met en facteur et defini des levels pour les iterations                          
                 levels = c("i1", "i2","i3","i4","i5","i10","i20","i30","i40","i50"))

#je stock le plot du contenue en sequence en fonction des differentes interation dans une variable
d500_x_1000_plot <- ggplot(d500_x_1000, aes(x=it, y=sum, color=it)) + 
  geom_bar(stat = "identity")+
  xlab("500kb_x_1000kb")+
  ylim(0,180000000)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90))

####################################################################
dx_1000 <- read.csv(file = "temp_4.csv", sep = "\t", header= F) #Table x>1000Kb, somme des tailles bins, it {1-50}
colnames(dx_1000)=c("sum", "it", "over")

dx_1000$it <- factor(dx_1000$it,                     #Je met en facteur et defini des levels pour les iterations                
                 levels = c("i1", "i2","i3","i4","i5","i10","i20","i30","i40","i50"))

#je stock le plot du contenue en sequence en fonction des differentes interation dans une variable
dx_1000_plot <-ggplot(dx_1000, aes(x=it, y=sum, color=it)) + 
  geom_bar(stat = "identity")+
  xlab("x>1000kb")+
  ylim(0,180000000)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90))

####################################################################
total <- read.csv(file = "temp2.csv", sep = "\t", header= F)#Table total, somme des tailles bins, it {1-50}
colnames(total)=c("sum", "it", "over")


total$it <- factor(total$it,                        #Je met en facteur et defini des levels pour les iterations                        
                 levels = c("i1", "i2","i3","i4","i5","i10","i20","i30","i40","i50"))

#je stock le plot du contenue en sequence en fonction des differentes interation dans une variable
total_plot <-ggplot(total, aes(x=it, y=sum, color=it)) + 
  geom_bar(stat = "identity")+
  xlab("total")+
  ylim(0,180000000)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90))
####################################################################
####################################################################

##ggarange permet d'arranger differents plot en une figure et permet de faire de tres belle figure article-ready
figure <-ggarrange(d10_x_100_plot,d100_x_500_plot,d500_x_1000_plot,dx_1000_plot,total_plot, ncol=5, nrow=1, common.legend = TRUE, legend="right", align="hv")

#annotate_figure permet de donner un beau titre à notre figure globale
annotate_figure(figure, top = text_grob("Analyse de l'Ã©volution des groupes de contigs en fonction du nombre d'itÃ©rations de l'algorithme de Louvain ", 
                                      color = "red", face = "bold", size = 10))
