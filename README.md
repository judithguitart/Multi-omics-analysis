# Multi-omics analysis of swine gut microbiome
Analysis of shotgun metagenomics data of 280 faecal samples and shotgun metatranscriptomics data of 147 faecal samples. Seven different groups of piglets treated with common post-weaning diarrhoea treatments at the time of the study were sampled longitudinally one day before weaning (ST1), three days after treatment (ST2), and two (ST3) and four (ST4) weeks after treatment.
Data availability in BioProject accession number PRJNA101070.
Related publications: https://doi.org/10.1186/s42523-024-00307-6


## 1. Pre-processing of metagenomics data

### 1.1 Raw Read Quality Analysis
Raw reads are renamed with ids and relevant metadata for initial raw read quality analysis with FastQC v0.11.9 and MultiQC v1.13 software:
```bash
conda activate fastqc
fastqc *.fastq.gz -o raw_fastqc_output/ -t 8
cd raw_fastqc_output/
multiqc *.zip
```

### 1.2 Host Decontamination and Trimming with kneaddata v0.10.0
After initial analysis of read quality, sequence counts, duplication levels, and adapter content, reads are host decontaminated using Bowtie2 v2.4.4 and trimmed using Trimmomatic v0.39.2. 
First, a Bowtie2 index must be created including the potential contaminant genomes. In this case, *Sus scrofa* (GCA_000003025.6_Sscrofa11.1) and *PhiX* bacteriophage genomes are downloaded and concatenated for index generation:
```bash
$PATH = ~/raw_data/all_files/
conda activate kneaddata
# Move output file into a new "contaminant_genomes/" folder
cat GCA_000003025.6_Sscrofa11.1.fasta phix.fasta contaminant_genomes/sscrofa11.1_phix.fasta 
bowtie2-build sscrofa11.1_phix.fasta > sscrofa11.1_phix.index
```
Kneaddata software is used for host decontamination and trimming of low quality and adapter sequences with the following options:
```bash
conda activate kneaddata
name=$(echo $file | sed -e 's/_1.fastq.gz//g')
kneaddata --remove-intermediate-output -t 16 
          --input1 ${name}_R1.fastq.gz --input2 ${name}_R2.fastq.gz --output /$PATH/kneaddata_output 
          --reference-db /$PATH/contaminant_genomes/Sscrofa11.1_phiX.index 
          --bowtie2-options "--very-sensitive --dovetail" 
          --trimmomatic ~/.conda/envs/kneaddata/share/trimmomatic-0.39-2/ 
          --trimmomatic-options "ILLUMINACLIP:adaptors.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50"
kneaddata_read_count_table --input /$PATH/kneaddata_output/ --output kneaddata_read_counts.txt
```
*kneaddata_read_counts.txt* is a summary file that contains the total reads included in each category of kneaddata (paired, unmatched, index_paired, index_unmatched) which allows for analysis of contamination levels. *Kneaddata_paried* files are saved in a separate folder for further analyses.

### 1.3 Trimmed Read Quality Analysis 
After host decontamination and trimming, quality of trimm reads is again analysed with FastQC and MultiQC software: 
```bash
conda activate fastqc
cd /$PATH/kneaddata_output/
fastqc *kneaddata_paired* -o trim_reads_fastqc/ -t 8
cd trim_reads_fastqc/
multiqc *.zip
```


## 2. Taxonomic classification and microbial diversity analyses

### 2.1 Taxonomic Classification of Reads
Kraken2 v2.1.2 classifier is used for read taxonomic classification using the Maxikraken2 database from the trimmed fastq reads. This process is run twice to obtain output in MPA style, by adding the `--use-mpa-style` flag, and kraken report style:
```bash
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
cd kraken_reports_mpa
python combine_mpa.py -i {all_mpa_reports_separated_by_spaces} -o all_mpa_reports.txt  #from KrakenTools
```
From run .err file, percentages of sequences classified and unclassified are extracted to analyse Kraken2 results.

### 2.2 Abundance and Diversity Data Processing with Microeco 
Taxonomy reports generated with kraken2 in a mpa-style format are used to generate abundance and taxonomy tables required for Microeco input.
Here, there is a summary of the main functions of Microeco, including rarefaction (complete output uploaded in *microeco_rarefaction_ed.docx*):
```R
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
rarefaction <- trans_rarefy$new(dataset_bacteria, alphadiv = "Shannon", depth = c(0, 10, 50, 500, 2000, 4000, 6000, 8000, 665830))
rarefaction$plot_rarefy(color_values = rep("grey", 279), show_point = TRUE, add_fitting = FALSE, show_legend = FALSE)
# Calculation of taxa abundance and alpha and beta divesity
dataset_bacteria$cal_abund()
dataset_bacteria$save_abund(dirpath = "taxa_abund_bacteria")
dataset_bacteria$cal_alphadiv(PD = FALSE)
dataset_bacteria$save_alphadiv(dirpath = "alpha_diversity_bacteria")
dataset_bacteria$cal_betadiv()
dataset_bacteria$save_betadiv(dirpath = "beta_diversity_bacteria")
# Comparison between groups and sampling times and plotting of alpha and beta diversity indexes
ta_bacteria <- trans_alpha$new(dataset = dataset_bacteria, group = "Group/ST")
ta_bacteria$cal_diff(method = "anova")
ta_bacteria$res_alpha_diff[1:6, ]
ta_bacteria$plot_alpha(add_letter = T, measure = "Shannon", use_boxplot = TRUE)
tb_bacteria <- trans_beta$new(dataset = dataset_bacteria, group = "Group/ST", measure = "bray")
tb_bacteria$cal_ordination(ordination = "PCoA")
tb_bacteria$plot_ordination(plot_color = "Group/ST", plot_shape = "Group/ST", plot_group_ellipse = TRUE)
# Statistical analysis
tb_bacteria$cal_manova(cal_manova_all = TRUE)
tb_bacteria$res_manova$aov.tab
tb_bacteria$cal_manova(cal_manova_paired = TRUE)
tb_bacteria$res_manova
tb_bacteria$cal_manova(cal_manova_set = "Group + ST")
tb_bacteria$res_manova$aov.tab
tb_bacteria$cal_betadisper()
tb_bacteria$res_betadisper
# More plots
g1 <- tb_bacteria$plot_group_distance(distance_pair_stat = TRUE)
subset_bacteria_phylum <- clone(dataset_bacteria)
subset_bacteria_phylum$tax_table %<>% subset(Phylum != "p__Unassigned")
tphylum_bacteria_subset <- trans_abund$new(dataset = subset_bacteria_phylum, taxrank = "Phylum", ntaxa = 12)
tphylum_bacteria_subset$plot_bar(others_color = "grey70", facet = "ST", facet2 = "Group",  xtext_keep = FALSE, legend_text_italic = FALSE) # exclude facet2 if necessary
tspecies_bacteria_subset$plot_box(group = "Group/ST")
tspecies_mean_bacteria_subset <- trans_abund$new(dataset = subset_bacteria_species, taxrank = "Species", ntaxa = 12, groupmean = "Group/ST")
tspecies_mean_bacteria_subset$plot_bar(xtext_keep = FALSE, legend_text_italic = FALSE) + theme_classic() + theme(axis.title.y = element_text(size = 18))
tspecies_bacteria_subset$plot_heatmap(facet = "ST/Group", xtext_keep = FALSE, withmargin = FALSE)
```

### 2.3 Statistical Analysis of Alpha Diversity Comparisons 
From *alpha_diversity.csv* file obtained with Microeco, statistical analysis are performed to compare sampling times within groups and groups within sampling times. Here are the main tests performed, after subsetting data (complete output uploaded *group_ST_comparisons_ed.docx*):
```R
# In R:
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
diversity  <- data.frame(metadata, index)
# Wilcoxon non-parametric t-test for all data
pairwise.wilcox.test(x = diversity$Shannon, g = diversity$Group, p.adjust.method = "BH")
pairwise.wilcox.test(x = diversity$Shannon, g = diversity$ST, p.adjust.method = "BH")
# One-way non-parametric ANOVA per group or ST
kruskal.test(G1$Shannon ~ G1$ST) # example with G1
# Wilcoxon non-parametric t-test per group or ST
pairwise.wilcox.test(x = G1$Shannon, g= G1$ST, p.adjust.method = "BH")
ggplot(G1, aes(x = ST, y = Shannon)) + geom_boxplot() + ylim(3.30,4.30) 
```

### 2.4 Multivariate Association Analysis with Maaslin2
From abundance, taxonomy, and metadata tables from Microeco: 
```R
# In R:
library(Maaslin2)
abundance <- t(abundance_table_bateria)
abundance <- abundance[rownames(metadata2),]
# Function to perform pre-filtering to keep only the feautures representing the 0.01% of the total abundance
# @param data - OTU count data frame of size n (sample) x p (OTU)
# @param percent - cutoff chosen 
low.count.removal = function(data, percent=0.0001) {
    keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
    data.filter = data[,keep.otu]
    return(list(data.filter = data.filter, keep.otu = keep.otu))
}
# Call the function then apply on the offset data 
result.filter <- low.count.removal(abundance, percent=0.0001)
input_data <- result.filter$data.filter
# Example for G1 within sampling times
input_metadataG1 <- subset(input_metadata, Group=='1')
input_dataG1 <- input_data[rownames(input_metadataG1),]
fit_data <- Maaslin2(input_dataG1,
	input_metadataG1,
	'output_genus_groupG1',
    	fixed_effects = c('ST'),
	reference = c("ST, V1"),
    	min_abundance = 0,
    	analysis_method = "LM",
    	max_significance = 0.1,
    	min_prevalence = 0.15,
    	cores = 3,
    	normalization = 'CLR',
    	transform = 'NONE',
    	standardize = FALSE)
```
From *significant_results* files obtained from each group and each sampling time, matrixes are created in Excel. These matrixes are compared for significant differences observed between sampling times compared against all groups and between groups contrasting all sampling times. 


## 3. Resistome characterization
From trimmed metagenomic reads, Resfinder is used for resistome characterization. Resfinder v4.2.5 software is chosen because it offers the possibility to obtain quantification of reads mapping against the Resfinder database with the output json file. 
```bash
conda activate resfinder
for file in *; do
	if [[ "$file" == *_fixed_1.fastq ]] ; then
	          name=$(echo $file | sed -e 's/_fixed_1.fastq//g')
	          echo $name
	          python3 -m resfinder -ifq ${name}_fixed_1.fastq ${name}_fixed_2.fastq -o ~/final_dataset/repeat_resfinder/${name}_repeat_resfinder -j ~/final_dataset/repeat_resfinder/${name}_resfinder.json -acq
	fi
done
```
To obtain sequencing depth of each ARG identified from the each sample's json file, the following JavaScript code was developed to obtain a matrix of sequencing depths per ARG and sample:
```JS
const fs = require('fs');
const globalDictionaryResult = {};
const globalAllGenes = [];
let globalCount = 0;

function init() {
    fs.readdir('./json_files/', (err, listOfFiles) => {
        for (file of listOfFiles) {
            if (file.endsWith('.json')) {
                const fileWithoutSufix = file.replace('.json', '');

                readFileAndWriteOutput(fileWithoutSufix, () => {
                    globalCount++;

		    if (globalCount === listOfFiles.length) {
                        const finalMatrix = [];
                        globalAllGenes.sort((genA, genB) => genA.localeCompare(genB))
                        finalMatrix.push(['Name', ...globalAllGenes].join('\t'));
                        Object.entries(globalDictionaryResult).forEach(([key, value]) => {
                            const matrixLine = [key];
                            for (gene of globalAllGenes) {
                                if (value[gene]) {
                                    matrixLine.push(value[gene])
                                } else {
                                    matrixLine.push('-')
                                }
                            }
                            finalMatrix.push(matrixLine.join('\t'));
                        })
                        writeInFile('superMatrix_metaG.tsv', finalMatrix.join('\n'))
                    }
                });
            }
        }
    })
}

function writeInFile(outputFileName, text) {
    fs.open(outputFileName, "a", (err, fd)=>{
        if (err){
            console.error('Error writing in' + outputFileName + '\n' + err.message);
            return;
        }
        fs.write(fd, text, (err) => {
            if (err) {
                console.error('Error writing in' + outputFileName + '\n' + err.message)
            }
        })
    })
}

function readFileAndWriteOutput(inputFilename, cb) {
    fs.readFile('./json_files/' + inputFilename + '.json', 'utf8', (error, data) => {
        if (error) {
            console.error(error);
            return;
        }

        const fileID = inputFilename.replace('_resfinder', '');
        const parsedJSON = JSON.parse(data).seq_regions;

        globalDictionaryResult[fileID] = {};
    
        Object.values(parsedJSON).forEach((value) => {
            if (globalDictionaryResult[fileID][value.name]) {
                globalDictionaryResult[fileID][value.name] += value.depth;
            } else {
                globalDictionaryResult[fileID][value.name] = value.depth;
            }
            if (globalAllGenes.indexOf(value.name) === -1) {
                globalAllGenes.push(value.name);
            }
        })
        cb();
    })
}
init();
```
The *'superMatrix_metaG.tsv'* file is then analysed in Excel for normalization by the total number of reads per sample and calculated in parts per million (PPM) and summed per antibiotic class for further statistical analyses. 


## 4. Binning stratey

### 4.1 Single sample assemblies (SA)
The metaspades.py script from SPAdes version 3.15.5 is used to assembly fastq files of trimmed paired reads with the following options: 
```bash
conda activate assembly
metaspades.py --only-assembler -1 ${name}_R1_kneaddata_paired_1.fastq -2 ${name}_R1_kneaddata_paired_2.fastq -o assembly/${name}_assembly
for i in `dir *_assembly/scaffolds.fasta`; do
	name=$(echo $i | sed "s/assembly\///");
	cp $i /$PATH/scaffolds/$name;
done
```
After running all assemblies, the for loop modifies each *scaffolds.fasta* file with the name id and copies them to a new directory containing all scaffolds. In this directory all fasta files are checked for quality with metaQUAST using the following script:
```bash
conda activate recycler
metaquast.py --threads 20 --label=name_1,$name_2,$name_3.... $name_1_scaffolds.fasta $name_2_scaffolds.fasta $name_3_scaffolds.fasta 
```

### 4.2 Co-assembly per treatment group and sampling time
Assembly of trimmed reads from the same group and sampling time together with Megahit v1.0.2. Example for G1_V1, ids in '-1' and '-2' must be separated individually per commas:
```bash
conda activate assembly
megahit -1 *_G1_V1_1_kneaddata_paired_1.fastq.gz -2 *_G1_V1_1_kneaddata_paired_2.fastq.gz
        --min-contig-len 1000 -t 40 --k-min 127 -o coassembly/coassembly_k_G1_V1 --presets meta-large 
```
Main output is in *final.contigs.fa* which must be renamed for further analyses.

### 4.3 Binning and Binning refinement
For this section, the MetaWRAP software is implemented in its version v1.3. The first implemented module is the MetaWRAP Binning module with CONCOCT, MaxBin, and metaBAT2 together. Before running the binning, it is necessary to fix the trimmed fastq files to make them readable for the MetaWRAP software. 
```bash
conda activate metawrap-env
repair.sh in1=${name}_1_kneaddata_paired_1.fastq in2=${name}_1_kneaddata_paired_2.fastq 
          out1=fixed_fastq/${name}_fixed_1.fastq out2=fixed_fastq/${name}_fixed_2.fastq 
          outsingle=fixed_fastq/${name}_single.fastq
metawrap binning -o binning_results/${name}_binning -t 96 -a ${name}_scaffolds.fasta 
                 --metabat2 --maxbin2 --concoct               
                 fixed_fastq/${name}_fixed_1.fastq fixed_fastq/${name}_fixed_2.fastq
```
Tip: for this module, create or move binning_results/ and fixed_fastq/ directories before running inside the directory where the *_scaffolds.fasta* files are saved. 
The same script is performed for co-assemblies, but only including MaxBin2 and MetaBAT2 tools due to computational requirements. Then, the refinement module of MetaWRAP is applied to all bins obtained from single assemblies first, and refined bins from single assemblies are then refined with the two sets of bins from co-assemblies frome each tool (MaxBin2 and MetaBAT2) with a minimum completeness of 70% and a maximum contamination of 10%:
```bash
conda activate metawrap-env
metawrap bin_refinement -o binning_results/${name}_binning/${name}_refinement -t 96 -A binning_results/${name}_binning/concoct_bins -B binning_results/${name}_binning/maxbin2_bins -C binning_results/${name}_binning/metabat2_bins -c 70 -x 10
# Single assemblies refinement results are then moved to a new folder per treatment group and sampling time named *'GX_STX_binning/GX_STX_bins/'* to be refined with the two sets of bins from co-assemblies:
metawrap bin_refinement -o binning_results/*GX_STX*_binning/*GX_STX*_refinement -t 96 -A binning_results/GX_STX_binning/GX_STX_bins/ -B binning_results/GX_STX_binning/coassembly_maxbin2_bins -C binning_results/GX_STX_binning/coassembly_metabat2_bins -c 70 x 10
```
From *'metawrap_70_10.stats'* files, completeness and contamination means are analysed, as well as taxonomic assginments. Bins quality is also checked with the CheckM software v1.2.0:
```bash
conda activate checkm 
checkm lineage_wf -t 36 -x fasta /$PATH/kneaddata_paired/scaffolds/{input} /$PATH/kneaddata_paired/scaffolds/checkm_results/{output}
```


## 5. Metagenome-Assembled Genomes (MAG) as a reference database

### 5.1 De-replication and quality analysis
From the total 17698 bins generated in the previous step, the dRep v3.4.5 software is used to remove repeated genomes and CheckM and CoverM v0.7.0 are used for quality analysis.
```bash
conda activate drep
dRep dereplicate drep_output/ -p 16 -comp 70 -con 10 -g {all_bins}.fa
conda deactivate
conda activate checkm
checkm lineage_wf -t 16 -x fa $PATH/dereplicated_bins/ $PATH/dereplicated_bins_checkm
conda deactivate
conda activate activate coverM
coverm genome --threads 4 -x fa --genome-fasta-directory $PATH/dereplicated_bins/ --methods relative_abundance mean length variance count --min-covered-fraction 0
conda deactivate
```

### 5.3 Phylogeny 
From final MAG catalogue, a phylogenetic tree is generated with PhyloPhlan v3.1.68 with the phylophlan database and following options: 
```bash
conda activate phylophlan
phylophlan -i dereplicated_bins -d phylophlan -t a --databases_folder phylophlan_bins/databases/ --diversity high -f supermatrix_aa.cfg --nproc 8 --verbose --force_nucleotides --genome_extension ".fa"
```

### 5.4 Taxonomy and functional annotation
For taxonomy and functional annotation of final MAGs, GTDB-Tk v2.0.1 and DRAM v1.5.0 are used, respectively:
```bash
conda activate GTDBTK-1.5.0
gtdbtk classify_wf --genome_dir dereplicated_bins/ --out_dir MAGs_gtdbk.dir -x *.fa --cpus X
```
```bash
conda activate DRAM_15
prokka --outdir prokka_bins_annotation *.fa
cd prokka_bins_annotation/
cat *.faa > annotation_file.faa
DRAM.py annotate --use_uniref -i annotation_file.faa --threads 20 -o dram_annotation
```
From the output annotation table, gene ids per fasta file of each bin is filtered together with each KO_id and KEGG_id. This 'MAGs_KO_filtered' table is used to transfrom annotations into multiple Genome Inferred Functional Traits (GIFTs) using the distillR package from R:
```R
library(tidyverse)
library(devtools)
library(distillR)
library(Rtsne)
library(patchwork)
library(viridis)
library(Polychrome)
library(broom)
library(dplyr)
library(tibble)
library(Rtsne)
library(ggplot2)
library(viridis)

# Load the Annotation table
SMGAtoGIFTs <- MAGs_KO_filtered %>%
  select(GENE,fasta,ko_id)
GIFTs <- distill(SMGAtoGIFTs,GIFT_db,genomecol=2,annotcol=3)
# Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs,GIFT_db)
# Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements,GIFT_db)
# Aggregate function-level GIFTs into overall Biosynthesis, Degradation and Structural GIFTs
GIFTs_domains <- to.domains(GIFTs_functions,GIFT_db)
avg_mci <- rowMeans(GIFTs_functions) %>%
  enframe(name = "Genome",value = "MCI")

Metadata <- metadata %>%
  select(Genome,phylum,order,tax) %>%
  distinct(Genome,.keep_all = T)
# Join the tables
avg_mci_Meta <- avg_mci %>%
  left_join(Metadata,by="Genome")

####Get the heatmap
Elements <- GIFT_db %>%
  select(Code_element,Domain, Function)
GIFTs_elements %>%
  pheatmap::pheatmap(cluster_cols = F)
#Add taxonomy to this table
GIFTs_W_Tax <- as_tibble(GIFTs_elements,rownames = NA) %>% 
  rownames_to_column("Genome")  %>%
  left_join(Metadata,by="Genome")
GIFTSannotCol <- Elements %>%
  select(Code_element,Function) %>%
  distinct() %>%
  column_to_rownames("Code_element")
GIFTSTaxAnnot <- GIFTs_W_Tax %>%
  select(Genome,phylum) %>%
  column_to_rownames("Genome")
ColorAnnot <- list(phylum=(Metadata %>%
                            select(phylum) %>%
                            distinct() %>%
                            mutate(Color=c(paletteer_dynamic("cartography::pastel.pal", 20),"#F7B6D2","#C7C7C7","#DBDB8D","#9EDAE5")) %>%
                            deframe()),
                   Function=(Elements %>%
                               select(Function) %>%
                               distinct() %>%
                               mutate(Color=rev(c(paletteer_dynamic("cartography::pastel.pal", 20),"#A174A2"))) %>%
                               deframe()))
Color <- rev(paletteer_c("grDevices::Purple-Blue", 10))
GIFTsElementsPH_p <- GIFTs_elements %>%
  pheatmap::pheatmap(cluster_cols = F,
                     color = Color,
                     annotation_col = GIFTSannotCol, annotation_colors=ColorAnnot,
                     annotation_row = GIFTSTaxAnnot,
                     show_colnames = F,
                     show_rownames = F)
ggsave(GIFTsElementsPH_p, file="GIFTsElementsPH_p2.pdf", width = 12, height = 12)
```
ARGs from final MAGs are also analysed with Resfinder v4.3.0:
```bash
conda activate resfinder
for file in *; do
	if [[ "$file" == *.fa ]] ; then
    		name=$(echo $file | sed -e 's/.fa//g')
    		echo $name 
		python3 -m resfinder -ifa ${name}.fa -o ~/metatranscriptomics/bins_annotation/dereplicated_bins/resfinder_bins/${name}_resfinder -j ~/metatranscriptomics/bins_annotation/dereplicated_bins/resfinder_bins/${name}_resfinder.json -acq
	fi
done
```


## 6. Pre-processing of metatranscriptomics data and mapping
Raw metatranscriptomic reads are renamed with ids and relevant metadata and filtered with the Kneaddata v0.12.0 software. Again, *Sus scrofa* (GCA_000003025.6_Sscrofa11.1) and *PhiX* bacteriophage genomes are filtered as contaminant genomes:
```bash
conda activate kneaddata
for file in *; do
	if [[ "$file" == *_1P.fastq.gz ]] ; then
		name=$(echo $file | sed -e 's/_1P.fastq.gz//g')
		echo $name
		kneaddata --remove-intermediate-output -t 16 --input ${name}_1P.fastq.gz --input ${name}_2P.fastq.gz --output $PATH/kneaddata_output --reference-db $PATH/contaminant_genomes/ --bowtie2-options "--very-sensitive --dovetail" --trimmomatic ~/.conda/envs/kneaddata/share/trimmomatic-0.39-2/ --trimmomatic-options "ILLUMINACLIP:adaptors.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50"
	fi
done
kneaddata_read_count_table --input $PATH/kneaddata_output/ --output kneaddata_read_count_table.tsv
column -t -s $'\t' kneaddata_read_count_table.tsv > kneaddata_read_counts.txt
```
The *kneaddata_read_count_table* utility allows for quality analysis of trimmed, decontaminated, and final paired reads.
Then, ribosomal RNA is removed with SortMeRNA v4.3.6 with Rfam and Silva databases with the following options:
```bash
conda activate sortmerna
for file in *; do
	if [[ "$file" == *_1P.fastq.gz ]] ; then
		name=$(echo $file | sed -e 's/_1P.fastq.gz//g')
		echo $name
		sortmerna --workdir kneaddata_output/ --out2 --paired_out --fastx --other /out/$input.noSalmon.noRNA --threads 12 --ref rRNA_databases/silva-bac-16s-id90.fasta --ref rRNA_databases/silva-arc-16s-id95.fasta --ref rRNA_databases/silva-bac-23s-id98.fasta --ref rRNA_databases/silva-arc-23s-id98.fasta --ref rRNA_databases/silva-euk-18s-id95.fasta --ref rRNA_databases/silva-euk-28s-id98.fasta --ref rRNA_databases/rfam-5.8s-database-id98.fasta --ref rRNA_databases/rfam-5s-database-id98.fasta --reads ${name}_1.fastq.gz --reads ${name}_2.fastq.gz
	fi
done
```


## 7. Transcript abundances and differential expression analyses
Then, trimmed and decontaminated reads are pseudoaligned for transcript abundance quantification with the Kallisto v0.50.1 software, using the annotation file from the MAGs catalogue:
```bash
conda activate kallisto
kallisto index --index=kallisto_index annotation_file.ffn -t 32
for file in *; do
	if [[ "$file" == *_1P.fastq.gz ]] ; then
		name=$(echo $file | sed -e 's/_1P.fastq.gz//g')
		echo $name
		kallisto quant -i kallisto_index.indx --plaintext ${name}_1P.fastq.gz ${name}_2P.fastq.gz -t 32 -o ${name}.kallisto_index.indx.dir
	fi
done
# Join abundance kallisto tables from all samples and extract tpms (sequence of 5,10,15,20...) or counts (sequence of 4,9,14,19...). Example for tpms: 
paste *.kallisto_index.indx.dir/abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv
```
For differential expression analyses, the table containing counts is used as DESeq2 uses raw counts and performs normalization:
```R
library("readr")
library("DESeq2")
library("ggplot2")
library("ggrepel")

cts_4out <- read_csv("cts_4out.csv")
coldata_4out <- read_csv("coldata_4out.csv")
head(coldata_4out)
gene_ids <- read_delim("gene_ids3.tab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cts_4out <- data.frame(cts_4out, row.names = 1, check.names = FALSE)
cts_4out <- round(cts_4out)
coldata_4out <- data.frame(coldata_4out, row.names = 1)
row.names(coldata_4out)
colnames(cts_4out)
all(rownames(coldata_4out) == colnames(cts_4out))
# Analysis with 'Group_Visit' variable:
dds_4out <- DESeqDataSetFromMatrix(countData = cts_4out, colData = coldata_4out, design = ~ Group_Visit + ID_mother)
keep_4out <- rowSums(counts(dds_4out)) >= 10
dds_4out <- dds_4out[keep_4out,]
dds_4out <- DESeq(dds_4out)
save(dds_4out, file = "DESeqData_4out.rda")
vsd_4out <- vst(dds_4out)
save(vsd_4out, file = "vsd_4out.rda")
pcaData_4out <- plotPCA(vsd_4out, intgroup=c("Group","Visit"), returnData = TRUE)
percentVar_4out <- round(100*attr(pcaData_4out, "percentVar"))
ggplot(pcaData_4out, aes(PC1,PC2, color=Group, shape=Visit)) + 
  geom_point(size = 3, stroke = 1.2) + scale_shape_manual(values=c(2,4,0)) +
  xlab(paste0("PC1: ",percentVar_4out[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_4out[2],"% variance")) +
  coord_fixed()
ggsave("PCA_counts_4out.png")
write.csv(pcaData_4out, "pcaData_4out.csv")

#Example of differential expression analyses between G2 and G1 at ST1:
resST1_G2vsG1 <- results(dds_4out, contrast=c("Group_Visit","G2_ST1","G1_ST1"))
resST1_G2vsG1_ord <- resST1_G2vsG1[order(resST1_G2vsG1$pvalue),]
resST1_G2vsG1_sig <- subset(resST1_G2vsG1_ord, padj < 0.05)
summary(resST1_G2vsG1_sig)
head(resST1_G2vsG1_sig)
resST1_G2vsG1_df <- data.frame(genes=rownames(resST1_G2vsG1_sig), resST1_G2vsG1_sig, row.names = NULL)
modify_ids <- match(resST1_G2vsG1_df$genes, gene_ids$TAXID)
resST1_G2vsG1_final <- transform(resST1_G2vsG1_df, genes = ifelse(!is.na(modify_ids), gene_ids$GENEID[modify_ids], genes))
head(resST1_G2vsG1_final)
#Exporting and visualizing the results:
write.csv(as.data.frame(resST1_G2vsG1_final), "ST1/resST1_G2vsG1_final.csv")
volcano_plot_ST1_G2vsG1 <- ggplot(resST1_G2vsG1_final, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.01)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "purple")) +
  labs(title = "Volcano Plot ST1 G2vsG1",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Significantly\nDifferentially Expressed") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")
volcano_plot_ST1_G2vsG1 + geom_text_repel(data = subset(resST1_G2vsG1_final, padj < 0.05), aes(label = genes), size = 2)
ggsave("volcano_plot_ST1_G2vsG1.jpg", path = "ST1")
```
From all results, the total number of DEGs between conditions can be identified, as well as specific comparisons. 

Then, to find the differentially expressed functional traits (GIFTs) between groups or sampling times, the *distillR* package can be used again from the results obtained from these specific comparisons. For this step, prokka_ids and KO_ids can be extracted from the DRAM output of the MAGs annotation. Then, KO_ids are added to the counts table used as input for DESeq2 and specific DEGs are extracted using the prokka_id from each specific comparison. Then, from the resulting filtered table of counts per DEG per sample, the following script extracts specific KO_ids present in each sample: 
```JS
const fs = require('fs');
const GLOBAL = {};

function init() {
    fs.readFile('./KO_comparison_X.csv', 'utf8', (error, data) => {
        const allRows = data.split('\n');
        const samplesID = allRows[0].split(';').slice(2, allRows.length);
        for (let index = 0; index < samplesID.length; index++) {
            const sampleID = samplesID[index];
            GLOBAL[sampleID] = [];
            for (let rowIndex = 1; rowIndex < allRows.length; rowIndex++) {
                const rowValues = allRows[rowIndex].split(';');
                const koId = rowValues[1];
                const value = rowValues[2 + index];

                if (value > 0) {
                    GLOBAL[sampleID].push(koId);
                }
            }
        }

        let finalData = '';
        Object.entries(GLOBAL).forEach(([key, values]) => {
            console.log(key)
            values.forEach((value) => {
                finalData += `${key.trim()};${value}\n`;
            });
        });

        fs.writeFile('./annotations_comparison_X.csv', finalData, () => {});
    });
}

init();
```
The output from this script can be uploaded in R to follow the same script as in *5.4.* and obtain the functional traits of the DEGs identified between specific conditions. 


## 8. Expressed Resistome characterization
For the expressed resistome, two different strategies can be followed. On one hand, the same analysis performed in section *3*, but using the trimmed metatranscriptomic reads as input. On the other hand, resistance-related terms can be filtered from the DESeq2 output results and further analysed. 
```bash
for file in *.csv; do
    name=$(echo $file | cut -d '_' -f1-2)
    grep -E "baseMean|resistan|antibiotic|antimicrobial|colistin|cillin|efflux|pump|drug|multi|lactam|cyclin|macrolid|midazol|penem|cepha|quinolon|mycin|polypeptid|sulfonamid|phenicol" ${name}_final.csv > ${name}_filtered.csv
done
```


## Other analyses not included in the thesis

### Antibiotic Resistance and Virulence genes identification with ABRicate
From assembled scaffolds obtained using metaSPAdes, ABRicate v1.0.1 is used to screen for antimicrobial resistance and virulence genes. The chosen databases are CARD, NCBI, and VFDB. Before launching ABRicate, *.fasta* extension must be changed to *.fa* to make files readable by ABRicate software.
```bash
rename .fasta.fa .fa *.fasta
conda activate abricate
abricate ${name}_scaffolds.fa --db card > /$PATH/kneaddata_paired/scaffolds_oct/abricate_scaffolds/abricate_card/${name}_card.tab
cd abricate_card/
abricate --summary *.tab > card_report.tab
abricate ${name}_scaffolds.fa --db ncbi > /$PATH/kneaddata_paired/scaffolds_oct/abricate_scaffolds/abricate_ncbi/${name}_ncbi.tab
cd abricate_ncbi/
abricate --summary *.tab > ncbi_report.tab
abricate ${name}_scaffolds.fa --db vfdb > /$PATH/kneaddata_paired/scaffolds_oct/abricate_scaffolds/abricate_vfdb/${name}_vfdb.tab
cd abricate_vfdb/
abricate --summary *.tab > vfdb_report.tab
```
These reports are downloaded and analysed in Excel to generate a presence/absence matrix with the percentage of coverage. 


### Plasmid assembly and identification with metaplasmidSPAdes
The chosen software for plasmid assembly from *.fastq* files is SPAdes version 3.15.5. It appeared to be more efficient than SCAPP, and web-based software like PLSDB did not support large metagenomics files. This is the script implemented with the following identification of plasmid scaffolds using BLAST software:
```bash
conda activate assembly
metaplasmidspades.py -t 48 --only-assembler 
                     -1 ${name}_R1_kneaddata_paired_1.fastq -2 ${name}_R1_kneaddata_paired_2.fastq 
                     -o /$PATH/kneaddata_paired/plasmid_assembly/${name}_plasmid
for i in `dir *_plasmid/scaffolds.fasta`; do name=$(echo $i | sed "s/\//_/"); cp $i /$PATH/plasmid_scaffolds/$name; done
conda deactivate
conda activate blast/2.13
blastn -task megablast -db nt -query ${name}_plasmid_scaffolds.fasta -out ${name}_plasmidspades.blast -evalue 1e-06 
       -outfmt "6 qseqid sseqid slen qstart qend length evalue score pident nident mismatch qcovs stitle sseq" 
       -max_target_seqs 5 -perc_identity 95 -num_threads $SLURM_CPUS_PER_TASK
```
