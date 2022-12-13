# Metagenomics Analysis
Metagenomics analysis of 280 faecal samples from swine treated with different antibiotic treatments and sampled at four different timepoints.

## 1. Data Pre-processing

### 1.1 Raw Read Quality Analysis
Raw reads are renamed with ids and relevant metadata for initial raw read quality analysis with FastQC and MultiQC software:
```
conda activate fastqc
fastqc {file}.fastq.gz -o {output_dir}/ -t 8
cd {output_dir}/
multiqc {file}.zip
```

### 1.2 Host Decontamination and Trimming 
After initial analysis of read quality, sequence counts, duplication levels, and adapter content, reads are host decontaminated using Bowtie2 and trimmed using Trimmomatic. 
First, a Bowtie2 index must be created including the potential contaminant genomes. In this case, *Sus scrofa* (GCA_000003025.6_Sscrofa11.1) and *PhiX* bacteriophage genoems are downloaded and concatenated for index generation:
```
$PATH = /scratch/jguitar/raw_data/all_files/
conda activate kneaddata
cat GCA_000003025.6_Sscrofa11.1.fasta phix.fasta > sscrofa11.1_phix.fasta # in a new directory named /contaminant_genomes/
bowtie2-build sscrofa11.1_phix.fasta > sscrofa11.1_phix.index
```
Kneaddata software is used for host decontamination and trimming of low quality and adapter sequences with the following options:
```
conda activate kneaddata
kneaddata --remove-intermediate-output -t 16 
          --input {file_name}_R1.fastq.gz --input {file_name}_R2.fastq.gz --output /$PATH/kneaddata_output 
          --reference-db /$PATH/contaminant_genomes/Sscrofa11.1_phiX.index 
          --bowtie2-options "--very-sensitive --dovetail" 
          --trimmomatic /mnt/beegfs/scratch/jguitar/.conda/envs/kneaddata/share/trimmomatic-0.39-2/ 
          --trimmomatic-options "ILLUMINACLIP:adaptors.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50"
kneaddata_read_count_table --input /$PATH/kneaddata_output/ --output kneaddata_read_counts.txt
```
*kneaddata_read_counts.txt* is a summary file that contains the total reads included in each category of kneaddata (paired, unmatched, index_paired, index_unmatched) which allows for analysis of contamination levels.

### 1.3 Trimmed Read Quality Analysis 
After host decontamination and trimming, quality of trimm reads is again analysed with FastQC and MultiQC software: 
```
conda activate fastqc
cd /$PATH/kneaddata_output/
fastqc *kneaddata_paired* -o {output_dir}/ -t 8
cd {output_dir}/
multiqc {file}.zip
```




