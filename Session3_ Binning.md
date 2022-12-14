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
On paste ensuite pour avoir un beau tableau :
```
cat temp1.txt |paste - - - > temp.csv
```
Ensuite on passe sur R : [Script](Script_R/Script_analyse_iteration.R)

On a :
![prodigal](/pictures/Graph1.png)

>On voit que à 20 de itérations, on a quelque chose de stable et avec un -O de 80/90.

Ensuite, on veut afficher pour différente range de Kb, la somme des tailles par itérations :

```
for iter in 1 2 3 4 5 10 20 30 40 50; do for o in 80; do cat binning/metator_"$iter"_"$o"/contig_data_partition.txt | awk '$10>=100000 && $10<=500000 {print $8,$10}'|sort -u | awk '{sum+= $2} END {print sum}' | awk '{print $1}' >> temp_2.txt; echo "i$iter" >> temp_2.txt;  echo "o$o" >> temp_2.txt ; done; done
```
On fais ça pour 4 range et le total et avec un -O de 80.
Ensuite on passe sur R : [Script](Script_R/Script_analyse_iteration.R)

On as :
![prodigal](/pictures/Graph2.png)

## B. Validation des bins

Les logiciels se basent sur les 43 gènes présents en une seule copie dans un génome bactérien pour évaluer la qualité des bins.

On va installer miComplete :
```
sudo pip install micomplete
```
On prend -i 20 et -O 80 pour la suite de l'analyse. Metator copie les génomes qu'il détecte dans le dossier overlapping_bin. C'est les fichiers de se dossier que on va donner à miComplete :
```
var=$(ls -l binning/metator_20_80/overlapping_bin/ | sed '1d' | awk '{print $9}' |grep ".fa" | awk -F "." '{print $1}') ; for i in $var; do cp binning/metator_20_80/overlapping_bin/"$i".fa binning/metator_20_80/overlapping_bin/"$i".fna; done
```
On lui créer une liste aussi :
```
find binning/metator_20_80/overlapping_bin/ -maxdepth 1 -type f -name "*.fna" | miCompletelist.sh > binning/metator_20_80/overlapping_bin/listbins.tab
```
Puis, on lance l'analyse :
```
miComplete binning/metator_20_80/overlapping_bin/listbins.tab --threads 8 --hmms Bact105 -o binning/metator_20_80/miComplete.txt 
```
Dans le fichier de miComplete.txt  le completeness c'est 100 + % de contamination (1.17 = 17% de contamination) et Markers = completetion/100 donc 0.72 de markers = 72% completion.
On peut alors représenter la qualité de nos génomes ici (un utilisant R : [Script](Script_R/Script_analyse_parametre_genome.r)) :
![prodigal](/pictures/Graph3.png)
(Les lignes rouges représentes les seuils d'un "bon génome" (>90% de complétion et <10% de contamination)

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 35 : Combien de génome(s) reconstruit(s) et complet(s) avez-vous ? Quelle proportion en terme de séquence cela représente t il ?**
>On voit alors des génomes très contaminer et d'autres incomplets. On peut retenir alors 16 génomes sur les 25 (5 très contaminés et 3 trop incomplets)

-----------------------------------------------------------------------------------------------------------------------------------------------------------
On fais ça avec itérations = 1/5/20/50. FICHIER : [Script](Script_R/Script_analyse_parametre_genome.r)
![prodigal](/pictures/Graph4.png)

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 36 : faites une analyse comparative de votre binning en fonction des différentes itérations**
>On voit que, pour notre communauté artificielle, l'effet de l'itération est très faible sur les différents paramètres testés (taille, couverture, contamination).

-----------------------------------------------------------------------------------------------------------------------------------------------------------

## C. Analyse des bins obtenus

On télécharge le dossier :
```
scp -r rhogg@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/2022-2023/ANALYSE_DES_GENOMES_2022_2023/TP_Meta3C/metator_final/ binning/ .
```
Dans ce dossier, on trouve pas mal de fichier donc un bin_summary.txt qui va nous servir pour réaliser plusieurs analyses. D'abors, on peut l'utiliser pour réaliser des picharts de notre communauté artificielle [Script](Script_R/Script_pie_plot.R)
![prodigal](/pictures/Graph5.png)

(L'abondance est en %). On voit que donc une bactérie présente plusieurs fois peut avoir une abondance faible et inversement.

## D. Matrice d'intéraction

On va créer cette matrice d'intéraction à l'aide de métator (et sa fonction contact map) et **hicstuff** tel que :
```
metator contactmap -t 8 -a assemblage/assembly_all.fa -c metator_final/contig_data_final.txt -e DpnII,HinfI -n "MetaTOR_22_2" -p metator_final/alignment_0.pairs -f -o matrices/MetaTOR_22_2/ -O "final_bin"
```
```
hicstuff view -n -b 10kb -f matrices/MetaTOR_22_2/MetaTOR_22_2.frags.tsv -o matrices/MetaTOR_22_2/mat_10kb_norm.pdf matrices/MetaTOR_22_2/MetaTOR_22_2.mat.tsv
```
Voila le genre de rendu que l'on a : ![prodigal](/pictures/Graph6.png)

Le petit carré tout seul en bas a droit représente le phage Spp1. Les lignes blanches sont les séquences répétées.

On regarde ensuite les contigs les plus couverts pour faire leur map :
```
cat metator_final/contig_data_final.txt | sed '1d' | awk '$3>= 30000 {print $2,($5*35)/$3,$17}' | sort -k 2,2 -g -r | head
```
On va maintenant partir à la recherche du plus grand phage circulaire detecté par :
```
cat database/VirSorter_results_table.txt | grep "circular"
```
On trouve donc un énorme phage de plus de 200 Kb qui après des grep infecte cette organisme (Metator_49_0) :
Hote du phage : k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas_aeruginosa

On fais aussi la matrice d'intéraction de ce contig de 200 kb du phage où l'on voit bien que c'est circulaire (1px : 10Kb) :
[Matrice](/pictures/Graph8.png)

On peut partir de la séquence du contig de se phage pour le blaster sur le NCBI avec le taxis bactial virus qui nous donnera que c'est le phage PhikZ ([Report_Blast](blast_report/T61EHBT201N-Alignment.txt)) ! On regarde ensuite la converture et %GC de Metator_49_0 [Script](Script_bash/bin_analysis.sh) :

[Graphique](/pictures/Graph7.png)

On voit bien la tache ne haut à droite avec un log(couverture) et %GC très différent de son hôte en base a droite.

## D. Conclusion générale

AU début, on avais mis : 6 phages, une dizaine de plasmide, 44 microorganisme et des crottes de souris avec 12 bactéries connues. On voit seulement 32 génomes nous ! On arrive donc pas a caractériser toutes les espèces que on a mis à l'origine ! Aussi, avec le 3C et nos algos, on detecte les bactéries à plusieurs chromosomes que les algo de binning ont du mal à voir.
