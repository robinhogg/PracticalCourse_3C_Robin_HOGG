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
Pour les gènes putatifs :
```
grep -c "^>" annotations/prodigal/assembly_gene.fa
```
Pour les protéines putatifs :
```
grep -c "^>" annotations/prodigal/assembly_prot.fa
```
>On trouve 177 970 gènes et protéines.

**Question 18 : Combien de gènes complets ? incomplets ?**
On sait que les gènes complets sont annotés dans l'en-tête par : partial=00 :
```
grep -c "partial=00" annotations/prodigal/assembly_gene.fa 
```
>On trouve 153 432 gènes complets et donc par logique : 177 970 - 153 432 = 24 538 gènes imcomplets.

**Question 19: Quelle est la longueur moyenne des gènes détectés ?**
On utilise awk :
```
cat annotations/prodigal/assembly_gene.fa | grep "^>" | awk '{print $5-$3+1}'|awk '{sum+=$1} END {print sum/NR}' | awk '{print "en moyennne les gènes font : "$1"bp"}'
```
>On a une longueur moyenne des gènes de 862.814 bp.

**Question 20 : Quel est le plus long gène détecté ?**
On utilise awk :
```
cat annotations/prodigal/assembly_gene.fa | grep "^>" | awk '{print $5-$3+1}'| sort -k 1,1 -g -r |head -1| awk '{print "le plus long gène fait : "$1"bp"}'
```
>Le gène le plus long fait : 19 529 bp.

**Question 21 : Quel est la longueur totale des gènes détectés ?**
On utilise awk :
```
cat annotations/prodigal/assembly_gene.fa | grep "^>" | awk '{print $5-$3}'|awk '{sum+=$1} END {print sum}' | awk '{print "la longueur total des gènes est de : "$1"bp"}'

```
>Cumulés, les gènes font : 153 377 084 bp

**Question 22 : Quelle est la densité en séquences codantes de votre assemblage ? cette valeur vous semble-t-elle cohérente ?**
On utilise awk :
```
cat annotations/prodigal/assembly_gene.fa | grep "^>" | awk '{print $5-$3}'|awk '{sum+=$1} END {print sum}' | awk '{print $1/171321468*100}'|  awk '{print "la densité codante est de : "$1"%"}'
```
>On trouve une densité codante de 89.52 %
-----------------------------------------------------------------------------------------------------------------------------------------------------------
L'étape d'après est de confronter ces séquences de gènes putatif contre des bases de données afin de voir des homologies entre ces gènes putatifs et des gènes déjà décrits. On peut aussi utiliser des algorithmes de type HMM afin de voir des motifs caractéristiques de différent élements génétique.

## B. Caractérisation des ORFs putatives

### a. Diamonds
