library(missForest)
library(glmnet)
library(caret)
library(biomaRt)
library(igraph)
library(parallel)
library(doParallel)
source('COSMO_function.R')


#################### Integrated COSMO ####################
rna_raw_file <- 'Example_data/rna.tsv'
pro_raw_file <- 'Example_data/pro.tsv'
cli_input_file <- 'Example_data/cli.tsv'

run_COSMOR(rna_raw_file, pro_raw_file, cli_input_file, impute_missing = F, cor_cutoff=0.5)







#################### Step-by-step Manual  ####################
rna_raw_file <- 'Example_data/rna.tsv'
pro_raw_file <- 'Example_data/pro.tsv'
cli_input_file <- 'Example_data/cli.tsv'

########## Preprocessing ##########
## Read Raw Data Input
raw_rnaseq <- read.table(rna_raw_file, sep='\t', header=T)     # featureN x sampleN
raw_proteome <- read.table(pro_raw_file, sep='\t', header=T)   # featureN x sampleN
clinical <- read.table(cli_input_file, sep='\t', header=T)     # sampleN x attributes


## Annotate chromosomes information
chr_annotate <- prpc_annotate(raw_rnaseq, raw_proteome)

## Preprocess RNAseq & Proteomics data
rnaseq <- prpc_omics(raw_rnaseq, chr_annotate, impute_missing = T)
proteome <- prpc_omics(raw_proteome, chr_annotate, impute_missing=T)




########## Pairwise Alignment ##########
align_outputs <- align_omics(rnaseq, proteome)

corsample    <- align_outputs$sampleCor_1
matcher      <- align_outputs$matching_2
nonmatch     <- which(matcher$mismatch_status == 1)




########## Attribute Prediction ##########
predict_outputs <- predict_gender(clinical, chr_annotate, rnaseq, proteome, matcher)

traincli <- predict_outputs$prediction_2
cli_suspect <- which(traincli$misannotate == 1)




#################### Label Correction ####################
correct_outputs <- determine_error(corsample, matcher, traincli)

swapped_sample <- correct_outputs$swapped_sample
shifted_sample <- correct_outputs$shifted_sample


## generate label corrected-table
final_tab <- generateCorrectedTable(traincli, swapped_sample, shifted_sample)


## visualize mislabeled samples
errorDir <- 'error_sample/'
dir.create(errorDir, recursive = TRUE, showWarnings = FALSE)
visualizeMismatch(corsample, nonmatch, errorDir)

write.table(swapped_sample, paste0(errorDir, 'swapped_samples.tsv'), sep='\t', row.names=F)
write.table(shifted_sample, paste0(errorDir, 'shifted_samples.tsv'), sep='\t', row.names=F)

visualizeSwapping(corsample, swapped_sample, errorDir)
visualizeShifting(corsample, shifted_sample, traincli, errorDir)






