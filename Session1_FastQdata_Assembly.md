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

>Ceci est une **zone en retrait**.
>La zone continue ici

>Ceci est une autre **zone de retrait**.
Cette zone continue également dans la ligne suivante.
Cependan, cette ligne n’est plus en retrait
