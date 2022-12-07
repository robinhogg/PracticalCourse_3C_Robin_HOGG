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
>On trouve une densité codante de 89.52 %. Le fait de trouver une bons nombres de gènes, complet, avec une taille moyenne et une densité codante cohérente avec les procaryotes, nous indiique que notre assemblage est de plutôt bonne qualité.
-----------------------------------------------------------------------------------------------------------------------------------------------------------
L'étape d'après est de confronter ces séquences de gènes putatif contre des bases de données afin de voir des homologies entre ces gènes putatifs et des gènes déjà décrits. On peut aussi utiliser des algorithmes de type HMM afin de voir des motifs caractéristiques de différent élements génétique.

## B. Caractérisation des ORFs putatives

### a. Diamond (blast-like)

C'est un aligneur comme blast du NCBI avec les protéines et l'ADN traduites, mais plus rapide avec les gros jeux de données. L'ensemble des paramètres du programme **Diamond** est conséquent.

On va d'abors créer une database (DB) avec les gènes de resistances à l'antibiotique :
```
diamond makedb -d database/Res_aa --in database/ARmeta-polypeptides.fa > log_files/blastdb_diamond.log 2>&1
```

On lance ensuite le programme :
```
diamond blastp -p 2 --db database/Res_aa.dmnd -o annotations/blast_output/prot_vs_AMR.txt --outfmt 6 -q annotations/prodigal/assembly_prot.fa
```
Le outpout fmt6 est le même que celui du blast NCBI.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 23 : Combien de gènes/protéines candidats obtenez vous ?**
```
cat annotations/blast_output/prot_vs_AMR.txt | awk '{print $1}' | sort -u | wc -l | awk '{print "Le nombre de gènes uniques qui s alignent est de : "$1}'
```
>On trouve 7 586 gènes qui s'alignent à l'aide de diamond. 

**Question 24: combien de vos gènes répondent aux critères :" on estime que l'on a un vrai homologue lorsque l'on a une identité d'au moins 80% sur 80% de la longueur du gène." ?**
On refais un diamond :
```
diamond blastp -p 8 --db database/Res_aa.dmnd -o annotations/blast_output/prot_vs_AMR_2.txt --outfmt 6 qseqid sseqid pident slen length -q annotations/prodigal/assembly_prot.fa
```
On passe ensuite par awk :
```
cat annotations/blast_output/prot_vs_AMR_2.txt | awk '{print $1,$3,($5/$4)*100}' | awk '$2>=80 && $3>=80 {print $1}' | sort -u | wc -l | awk '{print "le nombre de gènes candidats est de : "$1}'
```
>On trouve 210 gènes candidats avec ces seuils assez stricts.
-----------------------------------------------------------------------------------------------------------------------------------------------------------
### b. HMM (Modèle de Markov caché)

Le programme le plus utilisé pour les HMM est **HMMER**. C'est une boîte à outil qui peut faire pleins d'actions. Pour rappel, les modèles HMM se base sur la reconaissance de forme.
```
hmmsearch --tblout annotations/hmm_output/prot_vs_resfam.txt --cpu 2 database/Resfams.hmm  annotations/prodigal/assembly_prot.fa  >  annotations/hmm_output/prot_vs_resfam.out
```
-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 25 : Combien de candidats obtenez vous avec cette méthode ?**

On utiliser la foncton :
```
cat annotations/hmm_output/prot_vs_resfam.txt | grep "NODE" | awk '{print $1}' | sort -u | wc -l | awk '{print "le nombre de gènes candidats est de : "$1}'
```
>On trouve 14 826 candidats trouvés avec HMMER sans modifier les paramètres.
On relance le HMMER avec un nouveau paramètre (-E) permettant de sélectionner la E-Value. (Les auteurs de la base de données ResFams préconisent un seuil de détection de 1e-50 au niveau de la protéine entière.)
```
hmmsearch --tblout annotations/hmm_output/prot_vs_resfam.txt -E 1e-50 --cpu 2 database/Resfams.hmm  annotations/prodigal/assembly_prot.fa  >  annotations/hmm_output/prot_vs_resfam.out
```
>On relance la fonction pour chercher le nombre de candidats et on trouve 2 392 candidats avec cette -E de  1e-50.
-----------------------------------------------------------------------------------------------------------------------------------------------------------

A disposition, on a des fichiers résultats de programmes spécifiques qui va chercher directement des gènes candidats avec des caractéristique précises. Par exemple : **Resfinder** chercher des gènes de résistance aux antibiotiques.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 26 : Combien de gènes candidats le programme Resfinder détecte t il ? faites une comparaison avec les méthodes précédentes.**
On utilise awk :
```
cat database/Resfinder_results_table.txt  | grep "NODE" | awk '{print $7}' | sort -u | wc -l | awk '{print "le nombre de gènes candidats résistance aux antibiotiques : "$1}'
```
>On trouve 50 gènes candidats aux gènes de résistance aux antibiotiques
On créer 3 fichiers pour le diagramme de veine au format NODE_X_NODE_X :
```
cat database/Resfinder_results_table.txt  | grep "NODE" | awk '{print $6}' | sort -u > annotations/Resfinder_vein.txt
```
```
cat annotations/hmm_output/prot_vs_resfam.txt | grep "NODE" | awk '{print $1}' | sort -u | sed 's/_/ /g' | awk '{print $1"_"$2"_"$3"_"$4}' | sort -u > annotations/HMM_resfinder_vein.txt
```
```
cat annotations/blast_output/prot_vs_AMR_2.txt | awk '{print $1,$3,($5/$4)*100}' | awk '$2>=80 && $3>=80 {print $1}' | sort -u | sed 's/_/ /g' | awk '{print $1"_"$2"_"$3"_"$4}' | sort -u > annotations/blast_vein.txt
```
On fait le diagramme de vein sur : http://bioinformatics.psb.ugent.be/webtools/Venn/

![prodigal](/pictures/venn_result29840.png){:height="36px" width="36px"}.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
