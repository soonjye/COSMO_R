library(missForest)
library(glmnet)
library(caret)
library(biomaRt)
library(igraph)
library(parallel)
library(doParallel)

source('COSMO_function.R')


#### input files
rna_raw_file <- 'Example_data/rna.tsv'
pro_raw_file <- 'Example_data/pro.tsv'
cli_input_file <- 'Example_data/cli.tsv'


#### output files
interDir <- 'intermediate_result/'
dir.create(interDir, recursive = TRUE, showWarnings = FALSE)

chr_annot_file <- '1_chr_annotations.tsv'
prpc_rna_file  <- '1_prpc_rna.tsv'
prpc_pro_file  <- '1_prpc_pro.tsv'

geneCor_file         <- '2_Gene_correlation.tsv'
geneCor_histoImg     <- '2_Gene_correlation_histogram.png'
sampleCorMatrix_file <- '2_Sample_correlation_matrix.tsv'
sampleCor_heatmapImg <- '2_Sample_correlation_heatmap.png'
matching_result_file <- '2_Pairwise_matching.tsv'

attribute_prediction_file <- '3_Attribute_prediction.tsv'

final_table_file <- '4_Final_table.tsv'


#################### Preprocessing ####################
## Read Raw Data Input
raw_rnaseq <- read.table(rna_raw_file, sep='\t', header=T)    # featureN x sampleN
raw_proteome <- read.table(pro_raw_file, sep='\t', header=T)  # featureN x sampleN
clinical <- read.table(cli_input_file, sep='\t', header=T)    # sampleN x attributes


## Annotate chromosomes information
chr_annotate <- prpc_annotate(raw_rnaseq, raw_proteome)

## Preprocess RNAseq & Proteomics data
rnaseq <- prpc_omics(raw_rnaseq, chr_annotate, impute_missing = F)
proteome <- prpc_omics(raw_proteome, chr_annotate, impute_missing = F)


## Write into files
write.table(chr_annotate, paste0(interDir, chr_annot_file), sep='\t', col.names=T, row.names=F)
write.table(rnaseq, paste0(interDir, prpc_rna_file), sep='\t', col.names=T, row.names=T)
write.table(proteome, paste0(interDir, prpc_pro_file), sep='\t', col.names=T, row.names=T)




#################### Pairwise Alignment ####################
align_outputs <- align_omics(rnaseq, proteome, cutoff = 0.5)

gene_cor  <- align_outputs$geneCor_2
corsample <- align_outputs$sampleCor_2
matcher   <- align_outputs$matching_2
nonmatch  <- which(matcher$mismatch_status == 1)
cat('Found', length(nonmatch), 'omics mismatched samples! \n\n')


#### write into files
write.table(gene_cor, paste0(interDir, geneCor_file), sep='\t')

png(paste0(interDir, geneCor_histoImg), width=400, height=400)
align_outputs$geneCorPlot_2
dev.off()

write.table(corsample, paste0(interDir, sampleCorMatrix_file), sep='\t', row.names=T, col.names=T)

png(paste0(interDir, sampleCor_heatmapImg), width = 600, height = 600)
heatmap(corsample, Rowv=NA, Colv = NA, scale='none')
dev.off()

write.table(matcher, paste0(interDir, matching_result_file), sep='\t', row.names=F)




#################### Attribute Prediction ####################
rna_sex <- partitionSexMatrix(rnaseq, chr_annotate)
pro_sex <- partitionSexMatrix(proteome, chr_annotate)

predict_outputs <- predict_gender(clinical, rna_sex, pro_sex, matcher)

traincli <- predict_outputs$prediction_2
cli_suspect <- which(traincli$misannotate == 1)
cat(length(cli_suspect), 'is suspected to have gender mislabeled! \n\n')



#### write into file
write.table(traincli, paste0(interDir, attribute_prediction_file), sep='\t')




#################### Label Correction ####################
## Perform Label Correction
correct_outputs <- determine_error(corsample, matcher, traincli)

swapped_sample <- correct_outputs$swapped_sample
shifted_sample <- correct_outputs$shifted_sample


#### generate label-corrected table
final_tab <- generateCorrectedTable(traincli, swapped_sample, shifted_sample)
write.table(final_tab, paste0(interDir, final_table_file), sep='\t', row.names=F)


## visualize mislabeled samples
errorDir <- 'error_sample'
dir.create(errorDir, recursive = TRUE, showWarnings = FALSE)
#visualizeMismatch(corsample, nonmatch)

write.table(swapped_sample, paste0(errorDir, 'swapped_samples.tsv'), sep='\t', row.names=F)
write.table(shifted_sample, paste0(errorDir, 'shifted_samples.tsv'), sep='\t', row.names=F)

visualizeSwapping(corsample, swapped_sample, errorDir)
visualizeShifting(corsample, shifted_sample, traincli, errorDir)





#################### Integrated COSMO ####################



