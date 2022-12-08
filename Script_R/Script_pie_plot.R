#Script permettant de realiser differents graphiques d'analyser nos genomes sur le plan taxonomique
setwd('C:\\Users\\Robin\\Desktop\\pasteur\\3c')

#J'appelle mes differentes librairies"
library('ggplot2')
library("ggrepel")
library("dplyr")
library("ggpubr")
library("tidyr")


summary <- read.table(file = "bin_summary.txt", sep ="\t", header= TRUE) #je lis le résultat summary des bins

#Je filtre pour prendre que les genomes avec une completion>50% et une contamination<20%)
summary_filtred <- subset(summary, completness>=50 & contamination<=20)

#je sépare la colonne avec toute la taxonmy en fonction de chaque element
summary_filtred <- separate(data=summary_filtred, col ="taxonomy", into =c("k", "p","c", "o","f", "g","s"), sep =";")

######################################################
#Je met dans une nouvelle df les phylum et leur count
phylum <- summary_filtred  %>%
  group_by(p) %>%
  summarise(n = n())%>%
  arrange(desc(n))

#je fais les levels.. (a ameliorer clairement)
phylum$p <- factor(phylum$p,                        #Je met en facteur et defini des levels pour les iterations                        
                   levels = c("p__Firmicutes", "p__Proteobacteria","p__Actinobacteria","p__Bacteroidetes",
                              "p__Spirochaetes","p__Crenarchaeota","p__Cyanobacteria",
                              "p__Euryarchaeota","p__Tenericutes","p__Verrucomicrobia"))

#je fais le pie plot taxo distribution
phylum_d <- ggplot(phylum, aes(x="", y=n, fill=p))+
  geom_col(color = "black") +
  geom_label(aes(label = n),
             color = "white",
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y") + theme_void()+
  ggtitle(" MAGs taxonimic distribution (Phylum)")+
  labs(fill = "Phylum")+
  theme(plot.title = element_text(hjust = 0.5))


#Je met dans une nouvelle df les orders et leur count
order <- summary_filtred  %>%
  group_by(o) %>%
  summarise(n = n())%>%
  arrange(desc(n))

#je fais les levels.. (a ameliorer clairement)
order$o <- factor(order$o,                        #Je met en facteur et defini des levels pour les iterations                        
                   levels = c(
                      "o__Clostridiales",         
                      "o__Actinomycetales",       
                      "o__Lactobacillales",      
                      "o__Rhizobiales",       
                      "o__Bacillales",     
                      "o__Bacillales_3",    
                      "o__Bacteroidales" ,   
                      "o__Burkholderiales",
                      "o__Enterobacteriales",
                      "o__Rhodobacterales",  
                      "o__Spirochaetales",  
                      "o__Campylobacterales",  
                      "o__Chroococcales",
                      "o__Erysipelotrichales",   
                      "o__Methanobacteriales",   
                      "o__Pseudomonadales" ,    
                      "o__Selenomonadales",
                      "o__Sulfolobales",    
                      "o__Verrucomicrobiales",  
                      "o__Vibrionales"
                                  ))

#je fais le pie plot taxo distribution
order_d <-ggplot(order, aes(x="", y=n, fill=o))+
  geom_col(color = "black") +
  geom_label(aes(label = n),
             color = "white",
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y") + theme_void()+
  ggtitle(" MAGs taxonimic distribution (Order)")+
  labs(fill = "Order")+
  theme(plot.title = element_text(hjust = 0.5))

######################################################

#je creer une colonne avec une normaltion du HiC_cov 
summary_filtred$abondance = summary_filtred$HiC_Coverage/sum(summary_filtred$HiC_Coverage)

#on créer une DF avec que l'abondance est le phylum
phylum_abon <- summary_filtred[,-c(1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18)]

#je groupe par phylum et sum les abondance
phylum_abon2 <- phylum_abon %>%
    group_by(p) %>%
    summarise(abondance2 = sum(abondance))%>%
    arrange(desc(abondance2)) %>% mutate(across(where(is.numeric), ~ round(., 3)))

#je met en %
phylum_abon2$abondance2 <- phylum_abon2$abondance2*100

#je fais les levels.. (a ameliorer clairement)
phylum_abon2$p <- factor(phylum_abon2$p,                        #Je met en facteur et defini des levels pour les iterations                        
                   levels = c("p__Proteobacteria",     
                              "p__Firmicutes",  
                              "p__Spirochaetes",
                              "p__Actinobacteria",
                              "p__Crenarchaeota",
                              "p__Bacteroidetes",
                              "p__Euryarchaeota",
                              "p__Verrucomicrobia",
                              "p__Cyanobacteria",
                              "p__Tenericutes"))

#je fais le pie plot taxo abondance
phylum_a <-ggplot(phylum_abon2, aes(x="", y=abondance2, fill=p))+
  geom_col(color = "black") +
  geom_label(aes(label = abondance2),
             color = "white",
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y") + theme_void()+
  ggtitle(" MAGs taxonimic abondance (Phylum)")+
  labs(fill = "Phylum")+
  theme(plot.title = element_text(hjust = 0.5))


 
#on créer une DF avec que l'abondance est le order
order_abon <- summary_filtred[,-c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18)]

#je groupe par order et sum les abondance
order_abon2 <- order_abon %>%
  group_by(o) %>%
  summarise(abondance2 = sum(abondance))%>%
  arrange(desc(abondance2)) %>% mutate(across(where(is.numeric), ~ round(., 3)))

#je met en %
order_abon2$abondance2 <- order_abon2$abondance2*100

#je fais les levels.. (a ameliorer clairement)
order_abon2$o <- factor(order_abon2$o,                        #Je met en facteur et defini des levels pour les iterations                        
                         levels = c("o__Burkholderiales",
                                    "o__Spirochaetales",
                                    "o__Actinomycetales",
                                    "o__Sulfolobales",
                                    "o__Rhizobiales",
                                    "o__Bacillales_3",
                                    "o__Vibrionales",
                                    "o__Rhodobacterales",
                                    "o__Campylobacterales",
                                    "o__Enterobacteriales",
                                    "o__Selenomonadales",
                                    "o__Lactobacillales",
                                    "o__Bacteroidales",
                                    "o__Clostridiales",
                                    "o__Bacillales",
                                    "o__Pseudomonadales",
                                    "o__Methanobacteriales",
                                    "o__Verrucomicrobiales",
                                    "o__Chroococcales",
                                    "o__Erysipelotrichales"))


#je fais le pie plot taxo abondance
order_a <-ggplot(order_abon2, aes(x="", y=abondance2, fill=o))+
  geom_col(color = "black") +
  geom_label(aes(label = abondance2),
             color = "white",
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y") + theme_void()+
  ggtitle(" MAGs taxonimic abondance (order)")+
  labs(fill = "order")+
  theme(plot.title = element_text(hjust = 0.5))

####################################################
####################################################

##ggarange permet d'arranger differents plot en une figure et permet de faire de tres belle figure article-ready
ggarrange(phylum_d,phylum_a,order_d, order_a, ncol=2, nrow=2, align="hv")%>%
  ggexport(filename = "figure_pie_plot.png", width = 1800, height = 900)
                                   


