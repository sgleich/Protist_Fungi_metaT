**PARAGON MetaT Bioinformatic Pipeline**  
**By: Samantha Gleich & Syrena Whitner**   
**Last modified: 10/4/23**

![](static/protist.png)
![](static/fungi.tiff)

## Trim Sequences - Trimmomatic
Trim sequences using trimmomatic version 0.40
```
java -jar /paragon/trimmomatic_2023/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 16 -phred33 -trimlog [sample.log] [sample_R1.fastq] [sample_R2.fastq] [sample_R1_trimmed.fastq] [sample_R1_orphan.fastq] [sample_R2_trimmed.fastq] [sample_R2_orphan.fastq] ILLUMINACLIP:/path/to/directory/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50
```
## Separate rRNA from the read pool - SortMeRNA
We will use the [PR2 database] (https://pr2-database.org) version 5.0.0 and SortMeRNA version 4.3.6
```
sortmerna --ref pr2_version_5.0.0_SSU_mothur.fasta --reads sample_R1_trimmed.fastq --reads sample_R2_trimmed.fastq
```
