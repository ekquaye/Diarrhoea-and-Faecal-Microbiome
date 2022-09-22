################################################################################
################## Script for analysis in QIIME2 (v2021.4) #####################
################################ DON'T RUN IN R ################################

# changing working directory
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro ~ % cd Desktop/Seqs


# listing files in wd
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % ls
Undetermined_S0_L001_I1_001.fastq	metadata.tsv
Undetermined_S0_L001_R1_001.fastq	sample-metadata-thesis.tsv
Undetermined_S0_L001_R2_001.fastq	sample-metadata_ORIGINAL.xlsx


# creating .gz compressed files of sequences while maintaining the original files
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % gzip -c Undetermined_S0_L001_I1_001.fastq > barcodes.fastq.gz

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % gzip -c Undetermined_S0_L001_R1_001.fastq > forward.fastq.gz

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % gzip -c Undetermined_S0_L001_R2_001.fastq > reverse.fastq.gz


# creating folders
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % mkdir bacteria

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % mkdir emp-paired-end-sequences


# importing sequences into QIIME2
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime tools import \
> --type EMPPairedEndSequences \
> --input-path emp-paired-end-sequences \
> --output-path emp-paired-end-sequences.qza
Imported emp-paired-end-sequences as EMPPairedEndDirFmt to emp-paired-end-sequences.qza


# demultiplexing
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime demux emp-paired \
> --m-barcodes-file metadata.tsv \
> --m-barcodes-column Barcode_Sequence_16S \
> --i-seqs emp-paired-end-sequences.qza \
> --p-no-golay-error-correction \
> --o-per-sample-sequences bacteria/demux-bac.qza \
> --o-error-correction-details bacteria/demux-bac-details
Saved SampleData[PairedEndSequencesWithQuality] to: bacteria/demux-bac.qza
Saved ErrorCorrectionDetails to: bacteria/demux-bac-details.qza


# summarising the demultiplexing step
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime demux summarize \
> --i-data bacteria/demux-bac.qza \
> --o-visualization bacteria/demux-bac.qzv
Saved Visualization to: bacteria/demux-bac.qzv


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime tools view bacteria/demux-bac.qzv


# dada2 denoising step
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime dada2 denoise-paired \
> --i-demultiplexed-seqs bacteria/demux-bac.qza \
> --p-trunc-len-f 240 \
> --p-trunc-len-r 200 \
> --o-table bacteria/table-bac.qza \
> --o-representative-sequences bacteria/rep-seqs-bac.qza \
> --o-denoising-stats bacteria/denoising-stats-bac.qza
Saved FeatureTable[Frequency] to: bacteria/table-bac.qza
Saved FeatureData[Sequence] to: bacteria/rep-seqs-bac.qza
Saved SampleData[DADA2Stats] to: bacteria/denoising-stats-bac.qza


# feature table
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime feature-table summarize \
--i-table bacteria/table-bac.qza \
--o-visualization bacteria/table-bac.qzv \
--m-sample-metadata-file metadata.tsv
Saved Visualization to: bacteria/table-bac.qzv


# 16S rRNA sequences
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime feature-table tabulate-seqs \
--i-data bacteria/rep-seqs-bac.qza \
--o-visualization bacteria/rep-seqs-bac.qzv

Saved Visualization to: bacteria/rep-seqs-bac.qzv


# denoising stats
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime metadata tabulate \
--m-input-file bacteria/denoising-stats-bac.qza \
--o-visualization bacteria/denoising-stats-bac.qzv
Saved Visualization to: bacteria/denoising-stats-bac.qzv


# downloading Greengenes SEPP reference phylogenetic tree

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % conda install wget


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % wget \
-O "sepp-refs-gg-13-8.qza" \
"https://data.qiime2.org/2021.4/common/sepp-refs-gg-13-8.qza"
--2021-06-15 17:08:18--  https://data.qiime2.org/2021.4/common/sepp-refs-gg-13-8.qza
Resolving data.qiime2.org... 54.200.1.12
Connecting to data.qiime2.org|54.200.1.12|:443... connected.
HTTP request sent, awaiting response... 302 FOUND
Location: https://s3-us-west-2.amazonaws.com/qiime2-data/2021.4/common/sepp-refs-gg-13-8.qza [following]
--2021-06-15 17:08:20--  https://s3-us-west-2.amazonaws.com/qiime2-data/2021.4/common/sepp-refs-gg-13-8.qza
Resolving s3-us-west-2.amazonaws.com... 52.218.176.248
Connecting to s3-us-west-2.amazonaws.com|52.218.176.248|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 50161069 (48M) [binary/octet-stream]
Saving to: 'sepp-refs-gg-13-8.qza'

sepp-refs-gg-13-8.qza   100%[============================>]  47.84M   834KB/s    in 37s     

2021-06-15 17:08:58 (1.30 MB/s) - 'sepp-refs-gg-13-8.qza' saved [50161069/50161069]


# building phylogenetic tree from the Greengages SEPP tree backbone

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime fragment-insertion sepp \
> --i-representative-sequences bacteria/rep-seqs-bac.qza \
> --i-reference-database sepp-refs-gg-13-8.qza \
> --o-tree bacteria/insertion-tree-bac.qza \
> --o-placements bacteria/insertion-tree-placements.qza \
> --p-threads 3
Saved Phylogeny[Rooted] to: bacteria/insertion-tree-bac.qza
Saved Placements to: bacteria/insertion-tree-placements.qza


# filtering features that did not fit into reference tree

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime fragment-insertion filter-features \
--i-table bacteria/table-bac.qza \
--i-tree bacteria/insertion-tree-bac.qza \ 
--o-filtered-table filtered-table-bac.qza \
--o-removed-table removed-table-bac.qza
Saved FeatureTable[Frequency] to: filtered-table-bac.qza
Saved FeatureTable[Frequency] to: removed-table-bac.qza


######################################################################################################
################## disregard this section ###########################
#### I rebuilt the phylogenetic tree for a second time because I could not import the previous one into R as a phyloseq object… 
# Had to be sure this was not the problem. 
# Issue had to do with the ape/qiime2R package import step ####
# Files were succesfully imported on Windows PC 


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro bacteria % qiime fragment-insertion sepp \
--i-representative-sequences rep-seqs-bac.qza \
--i-reference-database sepp-refs-gg-13-8.qza \
--o-tree insertion-tree-two.qza \
--o-placements insertion-placements-two.qza \
--p-threads 3


Saved Phylogeny[Rooted] to: insertion-tree-two.qza
Saved Placements to: insertion-placements-two.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro bacteria % qiime fragment-insertion filter-features \
> --i-table table-bac.qza \
> --i-tree insertion-tree-two.qza \
> --o-filtered-table filtered-table-two.qza \
> --o-removed-table removed-table-two


Saved FeatureTable[Frequency] to: filtered-table-two.qza
Saved FeatureTable[Frequency] to: removed-table-two.qza

####################################################################################################


######## continued

# training feature classifier de novo, for classification of my sequences…Greengages 13.8 99 version

Last login: Sun Jun 20 12:46:03 on console


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro ~ % cd Desktop/Seqs 


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % mkdir training-feature-classifiers


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % cd training-feature-classifiers


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro training-feature-classifiers % ls
gg_13_8_otus		gg_13_8_otus.tar.gz


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro training-feature-classifiers % ls
99_otu_taxonomy.txt	gg_13_8_otus
99_otus.fasta		gg_13_8_otus.tar.gz


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro training-feature-classifiers % qiime tools import \
> --type 'FeatureData[Sequence]' \
> --input-path 99_otus.fasta \
> --output-path 99_otus.qza 

Imported 99_otus.fasta as DNASequencesDirectoryFormat to 99_otus.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro training-feature-classifiers % qiime tools import \
> --type 'FeatureData[Taxonomy]' \
> --input-format HeaderlessTSVTaxonomyFormat \
> --input-path 99_otu_taxonomy.txt \
> --output-path ref-taxonomy.qza



Imported 99_otu_taxonomy.txt as HeaderlessTSVTaxonomyFormat to ref-taxonomy.qza

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro training-feature-classifiers % qiime feature-classifier extract-reads \
> --i-sequences 99_otus.qza \
> --p-f-primer GTGCCAGCMGCCGCGGTAA \
> --p-r-primer GGACTACHVGGGTWTCTAAT \
> --p-min-length 200 \
> --p-max-length 400 \
> --o-reads ref-seqs.qza

Saved FeatureData[Sequence] to: ref-seqs.qza



(qiime2-2021.4) welcome@Welcomes-MacBook-Pro training-feature-classifiers % qiime feature-classifier fit-classifier-naive-bayes \
> --i-reference-reads ref-seqs.qza \
> --i-reference-taxonomy ref-taxonomy.qza \
> --o-classifier classifier.qza

Saved TaxonomicClassifier to: classifier.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro training-feature-classifiers % cd ..



(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime feature-classifier classify-sklearn \
> --i-reads bacteria/rep-seqs-bac.qza \
> --i-classifier training-feature-classifiers/classifier.qza \
> --o-classification bacteria/taxonomy-bac.qza

Saved FeatureData[Taxonomy] to: bacteria/taxonomy-bac.qza



(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime metadata tabulate \
> --m-input-file bacteria/taxonomy-bac.qza \
> --m-input-file bacteria/rep-seqs-bac.qza \
> --o-visualization bacteria/taxonomy-bac.qzv

Saved Visualization to: bacteria/taxonomy-bac.qzv



(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % qiime tools view bacteria/taxonomy-bac.qzv
Press the 'q' key, Control-C, or Control-D to quit. This view may no longer be accessible or work correctly after quitting.




####################################################################################################
# Verifying observed misclassification in Greengages for Alloiococcus with Silva 138


Last login: Thu May 12 07:52:23 on console
(base) welcome@Welcomes-MacBook-Pro ~ % conda activate qiime2-2021.4

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro ~ % pwd
/Users/welcome

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro ~ % ls 
Applications				Movies
CytoscapeConfiguration			Music
Desktop					Pictures
Documents				Public
Downloads				Zotero
Heritage.R				miniconda3
Library					qiime2-2021.4-py38-osx-conda.yml
MiCoNE

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro ~ % cd Desktop/Seqs

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % ls
16S:ITS import-demux script.rtf		fungi
Undetermined_S0_L001_I1_001.fastq	gg-13-8-99-515-806-nb-classifier.qza
Undetermined_S0_L001_R1_001.fastq	metadata.tsv
Undetermined_S0_L001_R2_001.fastq	removed-table-bac.qza
bacteria				sample-metadata-thesis.tsv
emp-paired-end-sequences		sample-metadata_ORIGINAL.xlsx
emp-paired-end-sequences.qza		sepp-refs-gg-13-8.qza
empress-tree-trial.qzv			silva-138-99-515-806-nb-classifier.qza
filtered-table-bac.qza			training-feature-classifiers

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % mkdir silva-alternate-tax-classification

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Seqs % cd silva-alternate-tax-classification

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro silva-alternate-tax-classification % ls
rep-seqs-bac.qza			silva-138-99-515-806-nb-classifier.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro silva-alternate-tax-classification % qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-515-806-nb-classifier.qza \
--i-reads rep-seqs-bac.qza \
--o-classification taxonomy-bac-silva.qza
Saved FeatureData[Taxonomy] to: taxonomy-bac-silva.qza

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro silva-alternate-tax-classification % qiime metadata tabulate \
> --m-input-file taxonomy-bac-silva.qza \
> --m-input-file rep-seqs-bac.qza \
> --o-visualization taxonomy-bac-silva.qzv
Saved Visualization to: taxonomy-bac-silva.qzv

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro silva-alternate-tax-classification % qiime tools view taxonomy-bac-silva.qzv 
Press the 'q' key, Control-C, or Control-D to quit. This view may no longer be accessible or work correctly after quitting.
#############################################################################################################




################################################################################

#############################################################
# in QIIME 2 for network analysis

# combined data (all samples)
# importing genera feature table 
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime tools import \
--input-path otu_biom.biom \
--type 'FeatureTable[Frequency]' \ 
--input-format BIOMV100Format \
--output-path feature-table-genera.qza


# importing taxonomy (up to genera) file
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_analyses % qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path tax-genera.txt \
--output-path taxonomy-genera.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime metadata tabulate \
--m-input-file taxonomy-genera.qza \
--o-visualization taxonomy-genera.qzv
Saved Visualization to: taxonomy-genera.qzv

# SCNIC (q2-scnic package) workflow

# using the sparcc filter
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime SCNIC sparcc-filter \
--i-table feature-table-genera.qza \
--o-table-filtered feature-table-genera-filtered.qza
Saved FeatureTable[Frequency] to: feature-table-genera-filtered.qza


# calculating correlations
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime SCNIC calculate-correlations \
--i-table feature-table-genera-filtered.qza \
--p-method sparcc \
--o-correlation-table correlations.qza

# plug in error
Plugin error from SCNIC:
  
  [Errno 2] No such file or directory: 'fastspar'

Debug info has been saved to /var/folders/65/rbbb7x4j3m779h2gn8r3jl_r0000gn/T/qiime2-q2cli-err-5zn7cgzi.log

# solving using the solution on the qiime 2 forum
# https://forum.qiime2.org/t/scnic-calculate-correlations-fastspar-not-found/18379
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % conda install -c bioconda -c conda-forge --override-channels fastspar


# calculating correlations (after solving the plugin error)
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime SCNIC calculate-correlations \
--i-table feature-table-genera-filtered.qza \
--p-method sparcc \
--o-correlation-table correlations.qza
Saved PairwiseFeatureData to: correlations.qza


# making a correlation network
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime SCNIC build-correlation-network-r \
--i-correlation-table correlations.qza \
--p-min-val 0.35 \
--o-correlation-network network.qza
Saved Network to: network.qza


## detecting and summarising modules of features (genera)
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime SCNIC make-modules-on-correlations \
--i-correlation-table correlations.qza \
--i-feature-table feature-table-genera-filtered.qza \
--p-min-r 0.35 \
--o-collapsed-table feature-table-genera-collapsed.qza \
--o-correlation-network modules.qza \
--o-module-membership membership.qza
Saved FeatureTable[Frequency] to: feature-table-genera-collapsed.qza
Saved Network to: modules.qza
Saved ModuleMembership to: membership.qza


# saving and visualising the membership file in .qzv format
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime metadata tabulate \
--m-input-file membership.qza \
--o-visualization membership.qzv
Saved Visualization to: membership.qzv


# exporting network file for visualisation outside qiime 2
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % qiime tools export \
--input-path modules.qza \                   
--output-path module-exported-file
Exported modules.qza as GraphModelingLanguageDirectoryFormat to directory module-exported-file


# NETWORK for heathy controls

mkdir Network_AGE_Healthy
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro R_to_QIIME2 % cd Network_AGE_Healthy


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime tools import \
--input-path otu_healthy_biom.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path feature-table-healthy.qza
Imported otu_healthy_biom.biom as BIOMV100Format to feature-table-healthy.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path tax-healthy-genera.txt \
--output-path taxonomy-genera-healthy.qza
Imported tax-healthy-genera.txt as HeaderlessTSVTaxonomyFormat to taxonomy-genera-healthy.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime metadata tabulate \
--m-input-file taxonomy-genera-healthy.qza \
--o-visualization taxonomy-genera-healthy.qzv                       
Saved Visualization to: taxonomy-genera-healthy.qzv


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC sparcc-filter \
--i-table feature-table-healthy.qza \
--o-table-filtered feature-table-healthy-filtered.qza
Saved FeatureTable[Frequency] to: feature-table-healthy-filtered.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC calculate-correlations \
--i-table feature-table-healthy-filtered.qza \
--p-method sparcc \
--o-correlation-table correlations-healthy.qza
Saved PairwiseFeatureData to: correlations-healthy.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC build-correlation-network-r \
--i-correlation-table correlations-healthy.qza \
--p-min-val 0.35 \
--o-correlation-network network-healthy.qza
Saved Network to: network-healthy.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC make-modules-on-correlations \
--i-correlation-table correlations-healthy.qza \
--i-feature-table feature-table-healthy-filtered.qza \
--p-min-r 0.35 \
--o-collapsed-table feature-table-healthy-collapsed.qza \
--o-correlation-network modules-healthy.qza \
--o-module-membership membership-healthy.qza
Saved FeatureTable[Frequency] to: feature-table-healthy-collapsed.qza
Saved Network to: modules-healthy.qza
Saved ModuleMembership to: membership-healthy.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime metadata tabulate \
--m-input-file membership-healthy.qza \
--o-visualization membership-healthy.qzv
Saved Visualization to: membership-healthy.qzv


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime tools export \
--input-path modules-healthy.qza \ 
--output-path modules-healthy-exported-file
Exported modules-healthy.qza as GraphModelingLanguageDirectoryFormat to directory modules-healthy-exported-file



# NETWORK for AGE cases

(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime tools import \
--input-path otu_AGE_biom.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path feature-table-AGE.qza
Imported otu_AGE_biom.biom as BIOMV100Format to feature-table-AGE.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path tax-AGE-genera.txt \
--output-path taxonomy-genera-AGE.qza
Imported tax-AGE-genera.txt as HeaderlessTSVTaxonomyFormat to taxonomy-genera-AGE.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime metadata tabulate \
--m-input-file taxonomy-genera-AGE.qza \
--o-visualization taxonomy-genera-AGE.qzv
Saved Visualization to: taxonomy-genera-AGE.qzv


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC sparcc-filter \
--i-table feature-table-AGE.qza \
--o-table-filtered feature-table-AGE-filtered.qza
Saved FeatureTable[Frequency] to: feature-table-AGE-filtered.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC calculate-correlations \
--i-table feature-table-AGE-filtered.qza \
--p-method sparcc \
--o-correlation-table correlations-AGE.qza
Saved PairwiseFeatureData to: correlations-AGE.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC build-correlation-network-r \
--i-correlation-table correlations-AGE.qza \
--p-min-val 0.35 \
--o-correlation-network network-AGE.qza
Saved Network to: network-AGE.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime SCNIC make-modules-on-correlations \
--i-correlation-table correlations-AGE.qza \
--i-feature-table feature-table-AGE-filtered.qza \
--p-min-r 0.35 \
--o-collapsed-table feature-table-AGE-collapsed.qza \
--o-correlation-network modules-AGE.qza \
--o-module-membership membership-AGE.qza
Saved FeatureTable[Frequency] to: feature-table-AGE-collapsed.qza
Saved Network to: modules-AGE.qza
Saved ModuleMembership to: membership-AGE.qza


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime metadata tabulate \
--m-input-file membership-AGE.qza \
--o-visualization membership-AGE.qzv
Saved Visualization to: membership-AGE.qzv


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro Network_AGE_Healthy % qiime tools export \
--input-path modules-AGE.qza \
--output-path modules-AGE-exported-file
Exported modules-AGE.qza as GraphModelingLanguageDirectoryFormat to directory modules-AGE-exported-file






################################################################################
##### Extracting sequences for submission in a public repository
(qiime2-2021.4) welcome@Welcomes-MacBook-Pro ~ % cd Desktop/Seqs/bacteria


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro bacteria % ls


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro bacteria % mkdir sequences        


(qiime2-2021.4) welcome@Welcomes-MacBook-Pro bacteria % qiime tools extract \
--input-path demux-bac.qza \
--output-path sequences
Extracted demux-bac.qza to directory sequences/a5f13fe3-3a7c-407b-beea-1f8a097d8329
