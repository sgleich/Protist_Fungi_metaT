# PARAGON MetaT Bioinformatic Pipeline
## By: Samantha Gleich & Syrena Whitner  
## Last modified: 1/5/24

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
## miTag (rRNA) analysis - Make QIIME2 classifier using the PR2 18S rRNA database
Use qiime2 (version 2022.2) to make a PR2 (version 5.0.0) 18S rRNA classifier.  
```
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads pr2_version_5.0.0_seq.qza --i-reference-taxonomy pr2_version_5.0.0_tax.qza --o-classifier pr2_version_5.0.0_classifier.qza
```
## miTag (rRNA) analysis - QIIME2 VSEARCH
Use qiime2 (version 2022.2) vsearch and PR2 version 5.0.0 to assign taxonomic annotations to rRNA reads (i.e., reads in the "aligned" files from SortMeRNA).  
  
First upload the rRNA reads into qiime2 format. The manifest.txt file is a comma-separated manifest file in the qiime2-specified format (sample-id,absolute-filepath,direction). This manifest.txt file lists all of the R1 and R2 rRNA read files obtained using SortMeRNA.
```
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --output-path demux.qza --input-format PairedEndFastqManifestPhred33
```
Then, join the paired-end reads.
```
qiime vsearch join-pairs --i-demultiplexed-seqs demux.qza --o-joined-sequences demux-joined.qza --p-minmergelen 200
```
Dereplicate the sequences (i.e., collapse identical sequences).
```
qiime vsearch dereplicate-sequences --i-sequences demux-joined.qza --o-dereplicated-table table.qza --o-dereplicated-sequences rep-seqs.qza
```
Cluster sequences at 97% similarity - create 97% OTUs.
```
qiime vsearch cluster-features-de-novo --i-table table.qza --i-sequences rep-seqs.qza --p-perc-identity 0.97 --o-clustered-table table-dn-97.qza --o-clustered-sequences rep-seqs-dn-97.qza
```
Assign taxonomy to the 97% OTUs using the previously made PR2 classifier. 
```
qiime feature-classifier classify-sklearn --i-classifier pr2_version_5.0.0_classifier.qza --i-reads rep-seqs-dn-97.qza --o-classification tax_sklearn.qza
```
Export OTU table and taxonomy table for rRNA analysis/visualization. 
```
qiime tools export --input-path tax_sklearn.qza --output-path output

qiime tools export --input-path table-dn-97.qza --output-path output

cd ./output
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
```
## Metatrascriptome assembly - rnaSPAdes
Concatenated all non-rRNA reads into single R1 and R2 files.
```
cat sample1_other_R1.fq sample2_other_R1.fq sample3_other_R1.fq > all_R1.fq  
cat sample1_other_R2.fq sample2_other_R2.fq sample3_other_R2.fq > all_R2.fq
```
Re-pair reads using bbmap repair.sh
```
repair.sh in=all_R1.fq in2=all_R2.fq out=all_new_R1.fq out2=all_new_R2.fq
```
Use the SPAdes assembly program (v 3.15.5) to assemble the metatranscriptome.
```
rnaspades.py -t 16 -m 1100 -o rnaspades_out -1 all_new_R1.fq -2 ./all_new_R2.fq 
```
## Assembly statistics - rnaQUAST
Calculate assembly statistics using rnaQUAST v. 2.2.2
```
rnaQUAST.py -c hard_filtered_transcripts.fasta -o rnaspades_quast_out
```
## Map reads to assembled contigs - salmon
Use salmon v. 1.10.1 to map transcripts to rnaSPAdes assembled contigs. First make a salmon index.
```
salmon index -t hard_filtered_transcripts.fasta -i salmon_index
```
Then map transcripts from each sample to this index. 
```
salmon quant -i salmon_index -l A -1 sample1_other_R1.fq -2 sample1_other_R2.fq -o sample1_salmon
```
## Identify putative protein-coding regions - Transdecoder
Use transdecoder v. 5.7.1 to identify putative protein coding regions.
```
TransDecoder.LongOrfs -t hard_filtered_transcripts.fasta -O transdecoder_rnaspades
```
## Functional annotations - eggnog-mapper
Use eggnog-mapper v. 2.0.1 to assign functional annotations to the putative protein-coding regions. 
```
emapper.py -i longest_orfs.pep --output eggnog_rnaspades -m diamond
```
## Taxonomic annotations - EUKulele
```
code coming.
```
