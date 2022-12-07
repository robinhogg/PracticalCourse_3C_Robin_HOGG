# Session 3: Binning des contigs

## A. MetaTOR - Metagenomic Tridimensional Organisation-based Reassembly

L'idée est de se baser sur les banques 3C et l'assemblage pour faire l'alignement de 3C sur l'assemblage pour construire le réseau d'intéraction. On va normaliser ce réseau et faire tourner l'algo de Louvain pour trouver des sous-réseau et évaluer la qualité de nos génomes bactériens avec des marqueurs taxonomiques présent en une copie dans les génomes connues (43 gènes présente en une copie dans 90% des génomes bactéries connus). On garde les résultats correct et on refais tourner l'algo de Louvain sur les communautés contaminées etc.

### a. Génération du réseau d’interactions inter-contig

On va d'abors créer ce réseau d'intération :
```
metator network -t 2 -1 fastq/lib9_filtre_3C_for.fastq.gz -2 fastq/lib9_filtre_3C_rev.fastq.gz -a assemblage/assembly_all.fa -o binning/metator/
```
Où
```
-t = le nombre de threads;
-1 = fastq file 3C for;
-2 = fastq file 3C rev;
-a = l'assemblage;
-o output.
```
En 3C, on aligne le for et rev de manière indépendante et on merge ensuite le résultat. La qualité du mapping est à 30. Le programme utilise l'aligneur Bowtie2.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 27 : Combien de nœuds (ou contigs) contient votre réseau global ?**




**Question 28 : Combien de paires de reads ont été alignées ?**




**Question 29 : Combien de paires de reads ont été alignées sur deux contigs différents ?**




**Question 30 : déduisez en le 3D ratio (nb de reads liant 2 contigs différent par rapport au nombre total de reads alignés)**




-----------------------------------------------------------------------------------------------------------------------------------------------------------
