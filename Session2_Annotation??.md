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


-----------------------------------------------------------------------------------------------------------------------------------------------------------
**Question 1 : Combien de lignes un read occupe-t-il ?**

>Dans un ficher .fastq, il y a 4 lignes par read.
