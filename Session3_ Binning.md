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
En 3C, on aligne le for et rev de manière indépendante et on merge ensuite le résultat. La qualité du mapping est à 30. Le programme utilise l'aligneur Bowtie2. Il faut réevaluer les données en normalisant avec les bruits de fond. La normalisation dite empirique est celle qui apporte les meilleurs résultats actuellement. On ne garde pas aussi les intéractions d'un contig sur lui-même.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 27 : Combien de nœuds (ou contigs) contient votre réseau global ?**
On passe par awk :
```
cat binning/metator/network.txt | awk '{print $1"\n"$2}' | sort -u | wc -l | awk '{print "il y a "$1" noeud"}'
```
>On a : 9 369 noeuds (contigs) (contre 20 566 de base).

**Question 28 : Combien de paires de reads ont été alignées ?**
On peut lire directement dans le fichier log créer par metator :
> On a : 3 948 736 pairs aligned.

**Question 29 : Combien de paires de reads ont été alignées sur deux contigs différents ?**
On peut lire directement dans le fichier log créer par metator :
>On a 2274943 contacts inter-contigs in the library.

**Question 30 : déduisez en le 3D ratio (nb de reads liant 2 contigs différent par rapport au nombre total de reads alignés)**
Le 3D ratio est une mesure de qualité. On peut lire directement dans le fichier log créer par metator :
>Le 3D ratio est : 0.5761192948832234
On a donc environ un read sur deux qui map sur deux contigs différents.
-----------------------------------------------------------------------------------------------------------------------------------------------------------

### b. Partitionnement du réseau d'interaction

Maintenant que on a le réseau, on va faire tourner l'algo de Louvain :
```
export LOUVAIN_PATH=software/gen-louvain/
```
Puis :
```
metator partition -i 1 -O 100 -t 2 -n binning/metator/network.txt -c binning/metator/contig_data_network.txt -a assemblage/assembly_all.fa -o binning/metator/
```
Où :
> * -i = itération de louvain
> * -t = threads
> * -n = le network
> * -c = le fichier contig
> * -a = fichier assemblage
> * -o = output
> * -O = seuil de overlap (%)

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 31 : Combien de bins détectez-vous ?**
On peut regarder la sortie d'erreur standard du programme :
>On a 692 core bins.


**Question 32 : Combien de contigs ne sont associés à aucun autre (ou combien de communautés ne comprennent qu'un seul contig) ?**
Manière longue :
```
grep -c "^.*       .*      .*      .*      .*      .*      .*      -" binning/metator/contig_data_partition.txt
```
Manière awk :
```
cat binning/metator/contig_data_partition.txt | awk '{print $2,$8}' | grep "-" | wc -l | awk '{print "Il y a "$1" contigs orphelins"}'
```
>On a 11 197 contigs qui ne sont associés à aucun autre.


**Question 33 : Combien de bin contiennent plus de 10 Kb, 100 Kb, 500 Kb et 1 Mb de séquences ?**
On utilise awk en replançant le chiffre après $10>= par 10 000, 100 000, 500 000 et 1 000 000 :
```
cat binning/metator/contig_data_partition.txt | awk '$10>=10000 {print $8,$10}'|sort -u | wc -l | awk '{print "Il y a "$1" contigs >=10kb"}'
```
>On a 59 bins >= 10 Kb, 42 >= 100 Kb, 34 >= 500 Kb et 34 >= 1 Mb. On peut faire l'approximation que on a bien entre 30 et 40 génomes dans notre assemblage initial.
-----------------------------------------------------------------------------------------------------------------------------------------------------------
On refais tourner l'algo (-F supprimer et remplace) :
```
metator partition -i 1 -O 100 -F -t 4 -n binning/metator/network.txt -c binning/metator/contig_data_network.txt -a assemblage/assembly_all.fa -o binning/metator/
```
-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 34 : Détectez-vous le même nombre de communautés que précédemment ? Ces communautés sont-elles de la même taille ? Qu'en déduisez-vous ?**
>Non, cette algo n'est pas déterministe (=stochastique). On va avoir des résultats différents entre chaque lancement.
-----------------------------------------------------------------------------------------------------------------------------------------------------------
On va faire tourner l'algo de Louvain avec plusieurs -i et -O différent pour voir comment le binning évolue.
```
for iter in 1 2 3 4 5 10 20 30 40 50; do for o in 50 60 70 80 90 100; do metator partition -i "$iter" -O "$o" -t 4 -n binning/metator/network.txt -c binning/metator/contig_data_network.txt -a assemblage/assembly_all.fa -o binning/metator_"$iter"_"$o"; done; done
```
Ensuite on va lancer une boucle pour compiler des données dans un fichier texte :
```
for iter in 1 2 3 4 5 10 20 30 40 50; do for o in 50 60 70 80 90 100; do cat binning/metator_"$iter"_"$o"/contig_data_partition.txt | awk '$10>=100000 {print $8,$10}'|sort -u | wc -l | awk '{print $1}' >> temp1.txt; echo "$iter" >> temp1.txt;  echo "o$o" >> temp1.txt ; done; done
```
On paste ensuite pour avoir une beau tableau :
```
cat temp1.txt |paste - - - > temp.csv
```
Ensuite on passr sur R : AJOUTER FICHIER SCRIPT R ICI

On a :
![prodigal](/pictures/Graph1.png)

```
for iter in 1 2 3 4 5 10 20 30 40 50; do for o in 80; do cat binning/metator_"$iter"_"$o"/contig_data_partition.txt | awk '$10>=100000 {print $8,$10}'|sort -u | awk '{sum+= $2} END {print sum}' | awk '{print $1}' >> temp2.txt; echo "$iter" >> temp2.txt;  echo "o$o" >> temp2.txt ; done; done
```
