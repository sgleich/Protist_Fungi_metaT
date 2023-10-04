**PARAGON MetaT Bioinformatic Pipeline**  
**By: Samantha Gleich & Syrena Whitner**   
**Last modified: 10/4/23**

## Trim Sequences - Trimmomatic
Trim sequences using trimmomatic version 0.40
```
java -jar /paragon/trimmomatic_2023/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE [sample-laneNum-readDir.fastq] ILLUMINACLIP:/path/to/directory/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50
```
