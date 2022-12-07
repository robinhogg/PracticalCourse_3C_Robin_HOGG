# Session 2: Annotation

L'annotation de notre métagénome répondra à deux objectifs principaux : 
 * Caractériser les différents organismes présente dans notre communauté articifiel (annotation taxonomique); 
 * Caractériser les différentes fonctions (annotation fonctionnelle), présentes dans l'échantillon.

## A. Recherche des phases ouvertes de lecture

Les séquences protéiques évoluent moins vite que les séquences nucléotidique. On peut utiliser des séquences nucléotidiques pour comparer des espèces proches. Mais, pour des metagénome, on regarder plutot les séquence protéique. Ici, on va regarder les phases ouverte de lecture (ORF) avec **Prodigal**. Il se base sur du machine learning. Il va scanner les séquences et nous proposer les gènes qu'ils prédits.

```
prodigal -p meta -a annotations/prodigal/assembly_prot.fa -o annotations/prodigal/assembly.gene -d annotations/prodigal/assembly_gene.fa -i  assemblage/assembly_all.fa > log_files/prodigal.log  2>&1
```
où
> * -p = pour définir que on travaille sur du métagénome;
> * -a : fichier d'output où seront écrites les séquences de protéines putatives;
> * -i : fichier d'input de l'assemblage global;
> * -o : Output pour ORFs détectées;
> * -d : Output ORFs putatives mais les séquences nucléotidiques.

Le programme Prodigal va nous rendre un fichier avec des séquences de ce format :
```
>contig1_1 # 2 # 232 # -1 # ID=1_1 ; partial=10 ; start_type=ATG ; rbs_motif=None ; rbs_spacer=None ; gc_cont=0.424
SEQUENCE
```
L'en-tête correspond à : >numéro_du_contig # début # fin # sens # ID_du_gène ; partial ; type_codon_start ; motif_site_liaison_ribosome ; distance_espacement_binding_site_ribosome ; contenu_en_GC

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 17 : Combien de gènes putatifs détectez-vous ? combien de protéines ?**
>


**Question 18 : Combien de gènes complets ? incomplets ?**



**Question 19: Quelle est la longueur moyenne des gènes détectés ?**



**Question 20 : Quel est le plus long gène détecté ?**



**Question 21 : Quel est la longueur totale des gènes détectés ?**



**Question 22 : Quelle est la densité en séquences codantes de votre assemblage ? cette valeur vous semble-t-elle cohérente ?**
-----------------------------------------------------------------------------------------------------------------------------------------------------------
L'étape d'après est de confronter ces séquences de gènes putatif contre des bases de données afin de voir des homologies entre ces gènes putatifs et des gènes déjà décrits. On peut aussi utiliser des algorithmes de type HMM afin de voir des motifs caractéristiques de différent élements génétique.

## B. Caractérisation des ORFs putatives

### a. Diamonds
