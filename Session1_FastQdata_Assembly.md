# Session 1: Contrôle qualité et traitement des séquences brutes issues du séquençage

Avant de traiter les données, ont les récupèrent :
  Récupérer le jeu de données 7 : 
 ```
  scp rhogg@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/2022-2023/ANALYSE_DES_GENOMES_2022_2023/TP_Meta3C/fastq/lib9_* fastq/
   ```
  Récupérer les séquences des adaptateurs illumina : 
   ```
  scp -r rhogg@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/2022-2023/ANALYSE_DES_GENOMES_2022_2023/TP_Meta3C/database/ ./
   ```
  Récuperer des software : 
   ```
  scp -r rhogg@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/2022-2023/ANALYSE_DES_GENOMES_2022_2023/TP_Meta3C/software/ ./
   ```

Avant de procéder à l'analyse ou à l'exploitation d'un ensemble de données de séquençage, il est important de réaliser des contrôles de qualité des séquences brutes et d'appliquer des filtres si nécessaires. Cette opération permettra de s'assurer qu'il n'y a pas de problèmes au niveaux des reads (qualité des reads, %GC, répartition des bases homogènes, contamination, etc....)

## A. On visualise les fichiers fastq

On visualise les données avec 
```
zcat (nom du fichier.gz) | head
```
-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 1 : Combien de lignes un read occupe-t-il ?**

>Dans un ficher .fastq, il y a 4 lignes par read.

**Question 2 : A quoi correspond chaque ligne ?**
> * La première ligne (@) donne la description de la séquence/reads (ici : @nom_du_sequenceur:flowcell:coordonée:sens_du_read(2 = reverse):index)
> * La deuxième ligne est la séquence.
> * La troisième ligne est séparatrice (+).
> * La dernière ligne est une ligne de qualité en format ASCII (X (qualité) + 33).	
					
**Question 3 : Combien de reads forward et reverse avez-vous dans vos jeux de données ?**

Façon d'y répondre laborieuse :
```
zgrep -c "^@" nom_du_fichier
```
Façon d'y répondre avec awk :
```
for type in 3C SG; do for sens in for rev; do echo "lib in progress:""$type""$sens"; zcat fastq/lib9_"$type"_"$sens".fastq.gz | wc -l | awk '{print $1/4}' | awk '{print $1, "total reads"}' ; zcat fastq/lib9_"$type"_"$sens".fastq.gz | head -2 | tail -1 | wc | awk '{print $3-1,"read lenght"}' ; echo ""; done ; done 
```


|Nom du fichier (.fastq.gz) |Nombre de reads|Longueur des reads (bp)|
|--------|--------|--------|
|    SG_for    |    10 000 000   |    150    |
|    SG_rev    |    10 000 000     |    150    |
|    3C_for   |    5 000 000    |    35   |
|    3C_rev    |    5 000 000    |    35    |

**Question 4 : Quelle est la longueur des reads (SG et 3C) ?**

>Les reads 3C sont beaucoup plus petits. Pour le 3C on veut du réseau d'interaction = pas besoin de grand reads qui coûte chère. Pour l'assemblage, il faut des grands reads (SG)

**Question 5 : Quels "Tags" sont associés à vos librairies ?**

>Tags = index de l'adaptateur. SG = GAATTCGT+AGGCTATA , 3C = AGCGATAG+AGGCTATA

**Question 6 : Quelles différences observez-vous entre les Reads SG et les Reads 3C ?**
			
>Ils possèdent une taille différente, un tag différent et un nombre de reads total différent.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
## B. On visualise la qualité de nos fichiers .fastq
On va réaliser le contrôle de qualité à l'aide du programme **FastQC**

On lance **FastQC** avec cette commande pour chaque fichier : 
```
/Formation_AdG/FastQC/fastqc -t 2 --nogroupe -o fastq/rapport_qualite/ fastq/lib9_SG_for.fastq.gz > log_files/fastqc_raw_SG_for.log 2>&1
```
où 
```
-t (nombre de CPU offert pour la tâche 
--nogroupe (graph pour chaque base)
```
**Résultats observés :**	
>La qualité globale des reads 3C sont excellent. Comparer à un small RNA seq, on n'a pas de séquence répété alarmante. On voit par contre, dans les >séquences surrepressenté les adaptateur TruSeq adaptater index 7 (mais vraiment léger).
-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 7 : En analysant et comparant les rapports de qualité, quelles différences observez vous entre vos différentes banques ? Quelle est l’enzyme que vous avez utilisée pour faire votre banque 3C ?**
>Pour le contenue en base de nos reads, pour SG, on voit que c'est homogène tout le long. Pour 3C, en début de reads, on voit un enrichissement de la >séquence : GATC. Cela proviens des enzymes de restriction que l'on as utilisées. On a utilisé une enzyme équivalente à DpnII (GATC). Donc, on enrichie >notre librairie en bout de séquence par ce motif. Cela permet de voir aussi que la banque 3C à bien marché : c'est le cas ici.
-----------------------------------------------------------------------------------------------------------------------------------------------------------

## C. Cutadapt : détection et retrait des séquences d’adaptateurs

On va nettoyer/grignoter nos reads. On peut par exemple enlever des séquences "enzyme de restriction" en début de reads, grignoter les reads si la qualité est mauvaise dans les bouts, supprimer les adaptateurs etc.
Pour ce faire on va utilisé **Cutadapt** (accepte missmatch dans les adaptateurs jusqu'à un certain niveau).

Les commandes que j'ai utilisé pour cut sont les suivantes :
```
cutadapt -q 20 -m 45 -j 2 -a file:database/adaptateur.fasta -A file:database/adaptateur.fasta -o fastq/lib9_filtre_SG_for.fastq.gz -p fastq/lib9_filtre_SG_rev.fastq.gz fastq/lib9_SG_for.fastq.gz fastq/lib9_SG_rev.fastq.gz > log_files/cutadapt_SG.log 2>&1
```
```
cutadapt -q 20 -m 33 -j 2 -a file:database/adaptateur.fasta -A file:database/adaptateur.fasta -o fastq/lib9_filtre_3C_for.fastq.gz -p fastq/lib9_filtre_3C_rev.fastq.gz fastq/lib9_3C_for.fastq.gz fastq/lib9_3C_rev.fastq.gz > log_files/cutadapt_3C.log 2>&1
```
Où
```
-q représente la qualité minimal accepté
-m est la longueur minimal du read 
-j nombre de CPU pour la tâche
-a séquence adaptateur que on veut cut
```
On relance ensuite **FastQC** sur les fichiers filtrés.
Pour chaque fichier :
```
/Formation_AdG/FastQC/fastqc -t 2 --nogroup -o fastq/rapport_qualite/ fastq/lib9_filtre_SG_for.fastq.gz > log_files/fastqc_filter_SG_for.log 2>&1
```
**Résultats observés :**	
>On a plus d'adaptateur sur-représenter dans mes données sauf dans le 3C rev où on retrouve quand même une séquence intriguante (GGCTACAGCATCTGCGTTTGTCTATGGTGATGAAT) dans les séquences sur-représentées du rapport fastQC n'ayant aucun hit dans leur banque et dans le blastn.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 8 : Combien de reads avez-vous gardé après cette étape de filtration ?**
|Nom du fichier (.fastq.gz) |Nombre de reads au départ|Nombre de read gardés| Nombre de reads perdu |
|--------|--------|--------|--------|
|    SG_for    |    10 000 000   |    9 971 770    |    28 230 (0.28 %)    |
|    SG_rev    |    10 000 000     |    9 971 7700    |    28 230 (0.28 %)    |
|    3C_for   |    5 000 000    |    4 664 225   |    335 775 (6.7%)    |
|    3C_rev    |    5 000 000    |    4 664 2255    |    335 775 (6.7%)    |
-----------------------------------------------------------------------------------------------------------------------------------------------------------
## D. Assemblage : Megahit

On va lancer un programme d'assemblage au nom de **MegaHit* qui va donc assembler nos reads sous forme de contigs. On fera ce megahit sur SG (pas 3C car, petit reads et pas beaucoup de reads et les for et rev peuvent être éloignés de plusieurs Mb)
```
/Formation_AdG/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -t 2 -1 fastq/lib9_filtre_SG_for.fastq.gz -2 fastq/lib9_filtre_SG_rev.fastq.gz -o assemblage/lib9/ > log_files/megahit_lib9_log 2>&1
```
`Par manque de temps, on va sauter cette étape et récupérer directement un assemblage complet réaliser par l'encadrant :`
```
scp rhogg@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/2022-2023/ANALYSE_DES_GENOMES_2022_2023/TP_Meta3C/assemblage/assembly_all.fa assemblage/ 
```
-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 9: Combien de contigs avez vous dans votre assemblage ?**
Pour répondre à cette question on peut faire simplement :
```
grep -c "^>" assembly_all.fa
```
>On trouve alors que on à 20566 contigs.

**Question 10: Quelle est la taille de votre plus grand contig ?**
J'ai fais :
```
cat assembly_all.fa | sed 's/_\(.*\)_\(.*\)_\(.*\)/ \1 \2 \3 /' | sort -k 4,4 -g -r | head -1 | awk '{print "le plus grand contig fait : "$4"bp"}'
```
Une version épurée :
```
cat assembly_all.fa | grep "^>" | sed 's/_/ /g' | sort -k 4,4 -g -r | head -1 | awk '{print "le plus grand contig fait : "$4"bp"}'
```
>La taille du plus grand contig est de 1 041 876 bp.

**Question 11: Quelle est la taille moyenne de vos contigs ?**
Pour ce faire, pas le choix, on passe par awk et son NR (nombre de ligne d'un fichier) :
```
cat assembly_all.fa | grep "^>" | awk -F "_" '{sum+=$4} END {print sum/NR}' | awk '{print "en moyennne les contigs font : "$1"bp"}'
```
>La taille moyenne de nos contigs est de 8 473.06 bp.
-----------------------------------------------------------------------------------------------------------------------------------------------------------
On va ensuite chercher à obtenir des informations statistique sur notre assemblage. Pour ce faire, on va utiliser le programme **Quast** tel que :
```
/Formation_AdG/quast-5.1.0rc1/quast.py assemblage/assembly_all.fa -o assemblage/rapport_assemblage/all_lib > log_files/quast_all.log 2>&1
```
-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 14 : Quelles sont les données fournies par Quast ?**
>Les données sont assez diverses. On a par exemple la taille totale de notre assemblage : 171 321 468 bp qui nous permet d'estimer à environs 50, le nombre d'espèces attendus (3-4 Mb par génome). On a aussi des information sur la qualité de l'assemblage (avec les N et L).  On a aussi la couverture en %GC et encore bon nombre d'informations sur l'assemblage.

**Question 15 : Donnez une définition du N50 ?**
> * N50 =information sur la fragmentation de assemblage : taille du plus petit conting permettant de couvrir 50% de l'assemblage. (tout les contigs au dessus (avec >lui) permette de faire 50% de l'assemblage)
> * L50 = rang de ce contigs dans tout les contigs.

**Question 16 : Quelle est la valeur théorique du N100 ?**
>C'est la taille du plus petits contigs de notre assemblage.
-----------------------------------------------------------------------------------------------------------------------------------------------------------

Texte normal suivi d’un [color=#26B260]texte coloré en vert[/color] dans un paragraphe.
FastQC en fasta : sed -n '1~4s/^@/>/p;2~4p' fastq_dir/reads.LegPneuPar3X.fastq | fold -w 80 > fasta_dir/reads.LegPneuPar3X.fasta
