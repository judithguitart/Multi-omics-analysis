# Metagenomics Analysis
Metagenomics analysis of 280 faecal samples from swine treated with different antibiotic treatments and sampled at four different timepoints.

## 1. Data Pre-processing

### 1.1 Raw Read Quality Analysis
Raw reads are renamed with ids and relevant metadata for initial raw read quality analysis with FastQC v0.11.9 and MultiQC v1.13 software:
```
conda activate fastqc
fastqc *.fastq.gz -o raw_fastqc_output/ -t 8
cd raw_fastqc_output/
multiqc *.zip
```

### 1.2 Host Decontamination and Trimming with kneaddata v0.10.0
After initial analysis of read quality, sequence counts, duplication levels, and adapter content, reads are host decontaminated using Bowtie2 v2.4.4 and trimmed using Trimmomatic v0.39.2. 
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
name = $(echo $file | sed -e 's/_1.fastq.gz//g')
kneaddata --remove-intermediate-output -t 16 
          --input {$name}_R1.fastq.gz --input {$name}_R2.fastq.gz --output /$PATH/kneaddata_output 
          --reference-db /$PATH/contaminant_genomes/Sscrofa11.1_phiX.index 
          --bowtie2-options "--very-sensitive --dovetail" 
          --trimmomatic /mnt/beegfs/scratch/jguitar/.conda/envs/kneaddata/share/trimmomatic-0.39-2/ 
          --trimmomatic-options "ILLUMINACLIP:adaptors.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50"
kneaddata_read_count_table --input /$PATH/kneaddata_output/ --output kneaddata_read_counts.txt
```
*kneaddata_read_counts.txt* is a summary file that contains the total reads included in each category of kneaddata (paired, unmatched, index_paired, index_unmatched) which allows for analysis of contamination levels. *Kneaddata_paried* files are saved in a separate folder for further analyses.

### 1.3 Trimmed Read Quality Analysis 
After host decontamination and trimming, quality of trimm reads is again analysed with FastQC and MultiQC software: 
```
conda activate fastqc
cd /$PATH/kneaddata_output/
fastqc *kneaddata_paired* -o trim_reads_fastqc/ -t 8
cd trim_reads_fastqc/
multiqc *.zip
```


## 2. Taxonomy Analyses 

### 2.1 Taxonomic Classification of Reads
Kraken2 v2.1.2 classifier is used for read taxonomic classification using the Maxikraken2 database from the trimmed fastq reads. This process is run twice to obtain output in MPA style and kraken report style:
```
conda activate kraken2
cd /$PATH/kneaddata_output/kneaddata_paired/
name=$(echo $file | sed -e 's/_R1_kneaddata_paired_1.fastq//g')
# mpa-style report run
kraken2 --db /prod/apps/kraken/maxikraken2_1903_140GB/ --confidence 0.01 --threads 16 --use-names 
        --report kraken_reports_mpa/${name}_report.txt --report-zero-counts --use-mpa-style --output kraken_outputs_mpa/${name}_output.txt 
        --paired ${name}_R1_kneaddata_paired_1.fastq ${name}_R1_kneaddata_paired_2.fastq --unclassified-out unclassified_mpa/$name#.fq
# kraken2-style report run
kraken2 --db /prod/apps/kraken/maxikraken2_1903_140GB/ --confidence 0.01 --threads 16 --use-names 
        --report kraken_reports/${name}_report.txt --report-zero-counts --output kraken_outputs/${name}_output.txt 
        --paired ${name}_R1_kneaddata_paired_1.fastq ${name}_R1_kneaddata_paired_2.fastq --unclassified-out unclassified/$name#.fq
```
From run .err file, percentages of sequences classified and unclassified are extracted to analyse Kraken2 results.

### 2.2 Abundance and Diversity Data Processing with Microeco 
Taxonomy reports generated with kraken2 in a mpa-style format are used to generate abundance and taxonomy tables required for Microeco input.
Here, there is a summary of the main functions of Microeco, including rarefaction (complete output uploaded in *microeco_rarefaction_ed.docx*):

```
# In R:
library(microeco)
library(magrittr)
library(mecodev)
library(tidyverse)

# Upload abundance, taxonomy, and metadata tables
microeco_tax_bacteria %<>% tidy_taxonomy()
dataset_bacteria <- microtable$new(otu_table = microeco_abu_bacteria, sample_table = metadata, tax_table = microeco_tax_bacteria)
dataset_bacteria_new
dataset_bacteria_new$tidy_dataset()
dataset_bacteria_new$sample_sums() %>%  range
rarefaction <- trans_rarefy$new(dataset_bacteria, alphadiv = "Shannon", depth = c(0, 10, 50, 500, 2000, 4000, 6000, 8000))
rarefaction$plot_rarefy(color_values = rep("grey", 279), show_point = TRUE, add_fitting = FALSE, show_legend = FALSE)

# Calculation of taxa abundance and alpha and beta divesity
dataset_bacteria$cal_abund()
dataset_bacteria$save_abund(dirpath = "taxa_abund_bacteria")
dataset_bacteria$cal_alphadiv(PD = FALSE)
dataset_bacteria$save_alphadiv(dirpath = "alpha_diversity_bacteria")
dataset_bacteria$cal_betadiv()
dataset_bacteria$save_betadiv(dirpath = "beta_diversity_bacteria")

# Comparison between groups and visits and plotting of alpha and beta diversity indexes
ta_bacteria <- trans_alpha$new(dataset = dataset_bacteria, group = "Group/Visit")
ta_bacteria$cal_diff(method = "anova")
ta_bacteria$res_alpha_diff[1:6, ]
ta_bacteria$plot_alpha(add_letter = T, measure = "Shannon", use_boxplot = TRUE)
tb_bacteria <- trans_beta$new(dataset = dataset_bacteria, group = "Group/Visit", measure = "bray")
tb_bacteria$cal_ordination(ordination = "PCoA")
tb_bacteria$plot_ordination(plot_color = "Group/Visit", plot_shape = "Group/Visit", plot_group_ellipse = TRUE)

# Statistical analysis
tb_bacteria$cal_manova(cal_manova_all = TRUE)
tb_bacteria$res_manova$aov.tab
tb_bacteria$cal_manova(cal_manova_paired = TRUE)
tb_bacteria$res_manova
tb_bacteria$cal_manova(cal_manova_set = "Group + Visit")
tb_bacteria$res_manova$aov.tab
tb_bacteria$cal_betadisper()
tb_bacteria$res_betadisper

# More plots
g1 <- tb_bacteria$plot_group_distance(distance_pair_stat = TRUE)
subset_bacteria_phylum <- clone(dataset_bacteria)
subset_bacteria_phylum$tax_table %<>% subset(Phylum != "p__Unassigned")
tphylum_bacteria_subset <- trans_abund$new(dataset = subset_bacteria_phylum, taxrank = "Phylum", ntaxa = 12)
tphylum_bacteria_subset$plot_bar(others_color = "grey70", facet = "Visit", facet2 = "Group",  xtext_keep = FALSE, legend_text_italic = FALSE) # exclude facet2 if necessary
tspecies_bacteria_subset$plot_box(group = "Group/Viisit")
tspecies_mean_bacteria_subset <- trans_abund$new(dataset = subset_bacteria_species, taxrank = "Species", ntaxa = 12, groupmean = "Group/Visit")
tspecies_mean_bacteria_subset$plot_bar(xtext_keep = FALSE, legend_text_italic = FALSE) + theme_classic() + theme(axis.title.y = element_text(size = 18))
tspecies_bacteria_subset$plot_heatmap(facet = "Visit/Group", xtext_keep = FALSE, withmargin = FALSE)
```

### 2.3 Statistical Analysis of Alpha Diversity Comparisons 
From *alpha_diversity.csv* file obtained with Microeco, statistical analysis are performed to compare visits within groups and groups within visits. Here are the main tests performed, after subsetting data (complete output uploaded *group_visit_comparisons_ed.docx*)

```
# In R:
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
diversity  <- data.frame(metadata, index)

# Wilcoxon non-parametric t-test for all data
pairwise.wilcox.test(x = diversity$Shannon, g = diversity$Group, p.adjust.method = "BH")
pairwise.wilcox.test(x = diversity$Shannon, g = diversity$Visit, p.adjust.method = "BH")

# One-way non-parametric ANOVA per group or visit
kruskal.test(G1$Shannon ~ G1$Visit) # example with G1
# Wilcoxon non-parametric t-test per group or visit
pairwise.wilcox.test(x = G1$Shannon, g= G1$Visit, p.adjust.method = "BH")
ggplot(G1, aes(x = Visit, y = Shannon)) + geom_boxplot() + ylim(3.30,4.30) 

# Further analysis: NMDS plot from ANOSIM test and PCA plots from beta-diversity
```

### 2.4 Multivariate Association Analysis with Maaslin2
From abundance, taxonomy, and metadata tables from Microeco: 

```
# In R:
library(Maaslin2)
abundance <- t(abundance_table_bateria)
abundance <- abundance[rownames(metadata2),]

# Function to perform pre-filtering to keep only the feautures representing the 0.01% of the total abundance
low.count.removal = function(data, # OTU count data frame of size n (sample) x p (OTU)
                        percent=0.0001 # cutoff chosen) {
    keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
    data.filter = data[,keep.otu]
    return(list(data.filter = data.filter, keep.otu = keep.otu)) 
    
# call the function then apply on the offset data 
result.filter <- low.count.removal(abundance, percent=0.0001)
input_data <- result.filter$data.filter

# Example for G1 within visits
input_metadataG1 <- subset(input_metadata, Group=='1')
input_dataG1 <- input_data[rownames(input_metadataG1),]
fit_data <- Maaslin2(input_dataG1, input_metadataG1, 'output_genus_groupG1',
    fixed_effects = c('Visit'),
    reference=c("Visit, V1"),
    min_abundance=0,
    #random_effects = c('Visit'),
    analysis_method="LM",
    max_significance = 0.1,
    min_prevalence=0.15,
    cores = 3,
    normalization = 'CLR',
    transform = 'NONE',
    standardize = FALSE)
```
From *significant_results* files obtained from each group and each visit, matrixes are created in Excel. These matrixes are compared for significant differences observed between visits compared against all groups and between groups contrasting all visits. 

**FALTA PHINCH???**


# 3. Co-assembly and Pangenome Analyses (not finished)

## 3.1 Co-assembly of reads with Megahit
Assembly of trimmed reads from the same group and visit together with Megahit v1.0.2. Example for G1_V1, ids in '-1' and '-2' must be separated individually per commas:
```
conda activate assembly
megahit -1 *_G1_V1_1_kneaddata_paired_1.fastq.gz -2 *_G1_V1_1_kneaddata_paired_2.fastq.gz
        --min-contig-len 1000 -t 40 --k-min 127 -o coassembly/coassembly_k_G1_V1 --presets meta-large 
```
Main output is in *final.contigs.fa* which must be renamed for further analyses.

## 3.2 Pangenome analysis (with anvi'o?)
**FALTA AQUESTA PART**


# 4. Assembly and Binning for MAG Generation 

## 4.1 Assembly with SPAdes assembler
The metaspades.py script from SPAdes version 3.15.5 is used to assembly fastq files of paired reads with the following options: 
```
conda activate assembly
metaspades.py --only-assembler -1 ${name}_R1_kneaddata_paired_1.fastq -2 ${name}_R1_kneaddata_paired_2.fastq -o assembly/${name}_assembly
for i in `dir *_assembly/scaffolds.fasta`; do name=$(echo $i | sed "s/assembly\///"); cp $i /$PATH/scaffolds/$name; done
```
After running all assemblies, the for loop modifies each *scaffolds.fasta* file with the name id and copies them to a new directory containing all scaffolds. In this directory all fasta files are checked for quality with metaQUAST and CheckM using the following scripts:
```
conda activate recycler
metaquast.py --threads 20 --label=name_1,$name_2,$name_3.... $name_1_scaffolds.fasta $name_2_scaffolds.fasta $name_3_scaffolds.fasta 
```
```
conda activate checkm
checkm lineage_wf -t 36 -x fasta /$PATH/kneaddata_paired/scaffolds/{input} /$PATH/kneaddata_paired/scaffolds/checkm_results/{output}
```

## 4.2 Binning and Binning Refinement with MetaWrap




# 5. Functional Annotation and Bin Taxonomy





# 6. Resistance, Virulence, and Plasmid Identification



