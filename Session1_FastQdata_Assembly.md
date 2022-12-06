# Session 1: Contrôle qualité et traitement des séquences brutes issues du séquençage

Avant de traiter les données ont les récupèrent :
  Récuperer le jeu de données 7 : 
 ```
  scp rhogg@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/2022-2023/ANALYSE_DES_GENOMES_2022_2023/TP_Meta3C/fastq/lib9_* fastq/
   ```
  Récuperer séq. adaptateur illumina : 
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

Façon d'y repondre laborieuse :
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

>Les reads 3C sont beaucoup plus petits. Pour le 3C on veut du réseau d'intéraction = pas besoin de grand reads qui coute chère. Pour l'assemblage, il faut des grands reads (SG)

**Question 5 : Quels "Tags" sont associés à vos librairies ?

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

FastQC en fasta : sed -n '1~4s/^@/>/p;2~4p' fastq_dir/reads.LegPneuPar3X.fastq | fold -w 80 > fasta_dir/reads.LegPneuPar3X.fasta
