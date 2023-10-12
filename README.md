# PARAGON MetaT Bioinformatic Pipeline
## By: Samantha Gleich & Syrena Whitner  
## Last modified: 10/12/23

![](static/protist.png)
![](static/fungi.tiff)

## Trim Sequences - Trimmomatic
Trim sequences using trimmomatic version 0.40.
```
java -jar /paragon/trimmomatic_2023/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 16 -phred33 -trimlog [sample.log] [sample_R1.fastq] [sample_R2.fastq] [sample_R1_trimmed.fastq] [sample_R1_orphan.fastq] [sample_R2_trimmed.fastq] [sample_R2_orphan.fastq] ILLUMINACLIP:/path/to/directory/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50
```
## Separate rRNA from the read pool - SortMeRNA
Use the [PR2 database](https://pr2-database.org) version 5.0.0 and SortMeRNA version 4.3.6 to separate rRNA reads from non-rRNA reads (i.e., other).
```
sortmerna --ref pr2_version_5.0.0_SSU_mothur.fasta --reads sample_R1_trimmed.fastq --reads sample_R2_trimmed.fastq --sam --fastx --aligned aligned --other other --paired_in --out2
```
## miTag (rRNA) analysis - QIIME2 VSEARCH
Use qiime2 (version 2020.6) vsearch and PR2 version 5.0.0 to assign taxonomic annotations to rRNA reads (i.e., reads in the "aligned" files from SortMeRNA).
```
code
```
## Metatrascriptome assembly - rnaSPAdes
Concatenated all non-rRNA reads into single R1 and R2 files.
```
code
```
Use the SPAdes assembly program (version X) to assemble the metatranscriptome.
```
code
```
