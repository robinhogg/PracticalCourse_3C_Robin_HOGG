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


```
cutadapt -q 20 -m 33 -j 2 -a file:database/adaptateur.fasta -A file:database/adaptateur.fasta -o fastq/lib9_filtre_3C_for.fastq.gz -p fastq/lib9_filtre_3C_rev.fastq.gz fastq/lib9_3C_for.fastq.gz fastq/lib9_3C_rev.fastq.gz > log_files/cutadapt_3C.log 2>&1
```
