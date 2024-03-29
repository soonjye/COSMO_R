#################### Functions ####################
#### Preprocessing: annotate genes with chromosomes position
prpc_annotate <- function(omics1, omics2) {
  geneSyms <- union(rownames(omics1), rownames(omics2))
  cat('Total number of unique genes =', length(geneSyms), '... \n')
  cat('Annotating with chromosomes information... \n')
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  chr_annotate <- getBM(attributes=c('hgnc_symbol', 'chromosome_name'), filters='hgnc_symbol', values=geneSyms, mart=mart)
  
  return(chr_annotate)
}

#### Load: gene chromosome annotation and extract sex chromosomes genes
getSexGenes <- function(chr_annotate, id_type="hgnc_symbol"){
  sexgenes <- which(chr_annotate$chromosome_name %in% c('X', 'Y'))
  sexgenes <- unique(chr_annotate[,id_type][sexgenes])
  return(sexgenes)
}

#### Preprocessing: clean omics data
prpc_omics <- function(omics, chr_annotate, impute_missing = FALSE) {
  cat('Preprocessing data matrix, dimension =', dim(omics), '\n')
  
  ## filter variables with no variance (all observation has the same or missing value)
  #rowvar <- apply(omics, 1, function(x) var(x, na.rm=T))
  #cat(' ', sum(rowvar == 0 | is.na(rowvar)), 'genes with no variance >> Removed... ')
  #omics <- omics[rowvar != 0 & !(is.na(rowvar)), ]
  #cat('Remaining =', nrow(omics), '\n')
  
  ## partition sex genes and autosomal genes
  sexgenes <- getSexGenes(chr_annotate)
  omics_atsm <- omics[setdiff(rownames(omics), sexgenes), ]
  omics_sex  <- omics[intersect(rownames(omics), sexgenes), ]
  
  
  #### For sex genes, replace missing values with 0
  cat('  Sex genes = ', nrow(omics_sex), '\n')
  
  rowvar <- apply(omics_sex, 1, function(x) var(x, na.rm=T))
  cat(sprintf('      %d genes with no variance >> Removed... \n', sum(rowvar == 0 | is.na(rowvar))))
  omics_sex <- omics_sex[rowvar != 0 & !(is.na(rowvar)), ]
  
  rowmiss <- apply(omics_sex, 1, function(x) sum(is.na(x)))
  cat(sprintf('      %d genes with missing value >> Replaced with 0 \n', sum(rowmiss > 0)))
  omics_sex[is.na(omics_sex)] <- 0
  cat('    Remaining sex gene =', nrow(omics_sex), '\n\n')
  
  
  #### For autosomal genes, filter variables with missing values in > 50% of sample size
  cat('  Autosomal gene =', nrow(omics_atsm), '\n')
  
  rowvar <- apply(omics_atsm, 1, function(x) var(x, na.rm=T))
  cat(sprintf('      %d genes with no variance >> Removed... \n', sum(rowvar == 0 | is.na(rowvar))))
  omics_atsm <- omics_atsm[rowvar != 0 & !(is.na(rowvar)), ]
  
  rowmiss <- apply(omics_atsm, 1, function(x) sum(is.na(x)))
  cat(sprintf('      %d genes with missing value in > 50%% (%.1f) of samples >> Removed... \n', sum(rowmiss > ncol(omics_atsm)/2), ncol(omics_atsm)/2))
  omics_atsm <- omics_atsm[rowmiss <= ncol(omics_atsm)/2, ]
  
  rowmiss <- apply(omics_atsm, 1, function(x) sum(is.na(x)))
  cat(sprintf('      %d genes with missing value in <= 50%% (%.1f) of samples >> ', sum(rowmiss <= ncol(omics_atsm)/2 & rowmiss > 0), ncol(omics_atsm)/2))
  
  #### For remaining autosomal genes, if impute_missing == TRUE, then impute missing values
  if (impute_missing == TRUE) {
    cat('Impute missing values... \n')
    omics_atsm <- t(omics_atsm)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    auto_imp <- missForest(omics_atsm, parallelize="variables")
    stopCluster(cl)
    omics_atsm <- auto_imp$ximp
    cat('        Imputation Out-Of-Bag error =', auto_imp$OOBerror, '\n')
    omics_atsm <- t(omics_atsm)
  } else {
    cat('Removed... \n')
    omics_atsm <- omics_atsm[rowmiss == 0, ]
  }
  cat('    Remaining autosomal genes =', nrow(omics_atsm), '\n')
  
  
  preprocessed <- rbind(omics_atsm, omics_sex)
  cat('Dimension of preprocessed data matrix =', dim(preprocessed), '\n')
  return(preprocessed)
}



#### Alignment: clean data input
cleanData <- function(omics) {
  cat('Dimension of data matrix = ', dim(omics), '\n')
  rowmiss <- apply(omics, 1, function(x) sum(is.na(x)))
  cat('Number of genes with missing value =', sum(rowmiss != 0), '>> Removed... \n')
  omics <- omics[rowmiss == 0, ]
  cat('Dimension = ', dim(omics), '\n')
  
  omics <- na.omit(omics)
  omics <- scale(t(omics))
  cat('Scaling and Centering... Completed \n')
  cat('Dimension = ', dim(omics), '\n\n')
  return(omics)
}

#### Alignment: Extract correlated genes
computeCorGene <- function(rnamatrix, promatrix){
  rnamatrix <- as.matrix(rnamatrix)
  promatrix <- as.matrix(promatrix)
  
  intergene <- intersect(colnames(rnamatrix), colnames(promatrix))
  cat(sprintf('Compute Pearson Correlation of %d genes (present in both dataset)... ', length(intergene)))
  gene_cor  <- rep(0, length(intergene))
  names(gene_cor) <- intergene
  for (gene in intergene) {
    gene_cor[gene] <- cor(rnamatrix[, gene], promatrix[, gene], method='pearson')
  }
  
  gene_cor <- gene_cor[order(gene_cor, decreasing = T)]
  gene_cor <- as.data.frame(gene_cor)
  colnames(gene_cor) <- 'Pearson_cor'
  cat('Completed \n')
  return(gene_cor)
}


#### Alignment: Invert sample correlation matrix into probability matrix
computeRankdist <- function(corsample){
  rnadistpro <- corsample
  prodistrna <- corsample
  for (i in 1:nrow(corsample)){
    rnadistpro[i,] <- exp(scale(rnadistpro[i,])) / sum(exp(scale(rnadistpro[i,])))
    prodistrna[,i] <- exp(scale(prodistrna[,i])) / sum(exp(scale(prodistrna[,i])))
  }
  rankdist <- sqrt(rnadistpro * prodistrna)
  
  return(rankdist)
}

#### Alignment: Perform stable matching based on probability matrix
stableMarriage <- function(cormatrix, mode = 'rnapro'){
  probmatrix <- computeRankdist(cormatrix)
  sampleN <- nrow(probmatrix)
  
  rankrow <- t(apply(probmatrix, 1, function(x) rank(-x, ties.method='first')))
  if (mode == 'clinic') {
    rankcol <- t(rankrow)
  } else {
    rankcol <- apply(probmatrix, 2, function(x) rank(-x, ties.method='first'))
  }
  
  matches <- data.frame(omics1=1:nrow(probmatrix), omics1_label=NA, omics2=NA, omics2_label=NA)
  matches$self_cor   <- diag(cormatrix)
  matches$self_score <- diag(rankcol + rankrow)
  matches$match_cor   <- NA
  matches$match_score <- NA
  singlemales <- 1:nrow(probmatrix)
  femalerings <- rep(NA, ncol(probmatrix))
  
  while (length(singlemales) != 0){
    for (ppsmale in singlemales){
      propose <- 1
      single <- TRUE
      while (single == TRUE) {
        ppsfmle <- which(rankrow[ppsmale,] == propose)
        engaged <- femalerings[ppsfmle]
        if (is.na(engaged) || rankcol[engaged, ppsfmle] > rankcol[ppsmale, ppsfmle]) {
          matches$omics2[ppsmale] <- ppsfmle
          femalerings[ppsfmle]    <- ppsmale
          matches$match_score[ppsmale] <- rankrow[ppsmale, ppsfmle] + rankcol[ppsmale, ppsfmle]
          matches$match_cor[ppsmale]   <- cormatrix[ppsmale, ppsfmle]
          singlemales             <- setdiff(singlemales, ppsmale)
          if (!(is.na(engaged))) singlemales <- c(singlemales, engaged) 
          single <- FALSE
        } else {
          propose <- propose + 1
        }
      }
    }
  }
  
  matches$omics1_label <- rownames(probmatrix)[matches$omics1]
  matches$omics2_label <- colnames(probmatrix)[matches$omics2]
  
  nonmatch <- which(matches$omics1 != matches$omics2)
  nonmatch <- c(nonmatch, which(matches$omics1 == matches$omics2 & matches$match_score > max(2, sampleN*0.1)))
  matches$mismatch_status <- 0
  matches$mismatch_status[nonmatch] <- 1
  return(matches)
}

#### Alignment: Perform Stable Marriage only on mislabeled instances
attentionStableMarriage <- function(cormatrix, threshold = 0.80){
  pervertical <- cormatrix
  perhorizontal <- cormatrix
  for (i in 1:nrow(cormatrix)){
    j <- i
    pervertical[i, ] <- (cormatrix[i, ] - min(cormatrix[i, ], cormatrix[, i]) ) / ( max(cormatrix[i, ], cormatrix[, i]) - min(cormatrix[i, ], cormatrix[, i]) )
    perhorizontal[, j] <- (cormatrix[, j] - min(cormatrix[, j], cormatrix[j, ]) ) / ( max(cormatrix[, j], cormatrix[j, ]) - min(cormatrix[, j], cormatrix[j, ]) )
  }
  
  #f1area <- 2 * pervertical * perhorizontal / (pervertical + perhorizontal)
  f1area <- sqrt(pervertical * perhorizontal)
  sampleN <- nrow(f1area)
  mismatches <- which(diag(f1area) < threshold)
  
  rankrow <- t(apply(f1area, 1, function(x) rank(-x, ties.method='first')))
  rankcol <- apply(f1area, 2, function(x) rank(-x, ties.method='first'))
  match_score <- rankrow+rankcol
  
  matches <- data.frame(omics1=1:nrow(f1area), omics1_label=NA, omics2=1:ncol(f1area), omics2_label=NA)
  matches$self_cor <- diag(cormatrix)
  matches$self_prob <- diag(f1area)
  matches$self_score <- diag(match_score)
  
  matches$match_cor <- NA
  matches$match_prob <- NA
  matches$omics2[mismatches] <- NA
  
  singlemales <- mismatches
  femalerings <- matches$omics2
  #femalerings[mismatches] <- NA
  
  while (length(singlemales) != 0){
    for (ppsmale in singlemales){
      propose <- 1
      single <- TRUE
      while (single == TRUE) {
        ppsfmle <- which(rankrow[ppsmale,] == propose)
        engaged <- femalerings[ppsfmle]
        if (ppsfmle %in% mismatches && (is.na(engaged) || f1area[engaged, ppsfmle] < f1area[ppsmale, ppsfmle])) {
          matches$omics2[ppsmale] <- ppsfmle
          femalerings[ppsfmle]    <- ppsmale
          singlemales             <- setdiff(singlemales, ppsmale)
          if (!(is.na(engaged))) singlemales <- c(singlemales, engaged)
          single <- FALSE
        } else {
          propose <- propose + 1
        }
      }
    }
  }
  
  for (i in 1:sampleN) {
    matches$match_prob[i] <- f1area[matches$omics1[i], matches$omics2[i]]
    matches$match_cor[i]  <- cormatrix[matches$omics1[i], matches$omics2[i]]
    matches$match_score[i]<- match_score[matches$omics1[i], matches$omics2[i]]
  }
  matches$omics1_label <- rownames(f1area)[matches$omics1]
  matches$omics2_label <- colnames(f1area)[matches$omics2]
  
  matches$mismatch_status <- 0
  matches$mismatch_status[mismatches] <- 1
  return(matches)
}

#### Alignment: Comprehensive function to generate every table and plot
align_omics <- function(rnaseq, proteome, cutoff = 0.5) {
  align_output <- list()            # initialize output with a list of intermediate results
  rnaseq <- cleanData(rnaseq)       # Remove feature with missing values
  proteome <- cleanData(proteome)   # Remove feature with missing values
  
  ## Extract overlapped features and samples
  sample_label <- intersect(rownames(rnaseq), rownames(proteome))
  sampleN      <- length(sample_label)
  cat('Number of samples present in both omics data =', sampleN, '\n')
  
  feature_label <- intersect(colnames(rnaseq), colnames(proteome))
  cat('Number of variables present in both omics data =', length(feature_label), '\n\n')
  
  rnaseq <- rnaseq[sample_label, feature_label]
  proteome <- proteome[sample_label, feature_label]
  
  
  ## First iteration of matching
  cat('First iteration: Compute Gene Correlation ... \n')
  gene_cor <- computeCorGene(rnaseq, proteome)
  align_output$geneCor_1 <- gene_cor
  align_output$geneCorPlot_1 <- ggplot(gene_cor, aes(x=Pearson_cor)) +
    geom_histogram(color='black', fill='white') +
    geom_vline(xintercept = cutoff, linetype="dashed") +
    labs(x='Pearson Correlation', title='Distribution of Gene-wise Pearson Correlation',
         subtitle=sprintf('Mean = %.3f; SD = %.3f', mean(gene_cor$Pearson_cor, na.rm=T), sd(gene_cor$Pearson_cor, na.rm=T)))
  
  corgene <- rownames(gene_cor)[gene_cor$Pearson_cor > cutoff & !(is.na(gene_cor$Pearson_cor))]
  cat(sprintf('Selected %d genes (correlation > %.2f) \n', length(corgene), cutoff))
  
  corsample <- cor(t(rnaseq[, corgene]), t(proteome[, corgene]))
  align_output$sampleCor_1 <- corsample
  cat('First iteration: Perform data matching ... \n')
  
  matcher <- stableMarriage(corsample)
  align_output$matching_1 <- matcher
  
  nonmatch <- which(matcher$mismatch_status == 1)
  cat('First iteration: Found', length(nonmatch), 'pair mismatch\n\n')
  
  
  ## Second iteration of matching
  if (length(nonmatch) > 0) {
    cat('Second iteration: Compute Gene Correlation without mismatch sample... \n')
    gene_cor <- computeCorGene(rnaseq[-nonmatch, ], proteome[-nonmatch, ])
    align_output$geneCor_2 <- gene_cor
    align_output$geneCorPlot_2 <- ggplot(gene_cor, aes(x=Pearson_cor)) +
      geom_histogram(color='black', fill='white') +
      geom_vline(xintercept = cutoff, linetype="dashed") +
      labs(x='Pearson Correlation', title='Distribution of Gene-wise Pearson Correlation',
           subtitle=sprintf('Mean = %.3f; SD = %.3f', mean(gene_cor$Pearson_cor, na.rm=T), sd(gene_cor$Pearson_cor, na.rm=T)))
    
    corgene <- rownames(gene_cor)[gene_cor$Pearson_cor > cutoff & !(is.na(gene_cor$Pearson_cor))]
    cat(sprintf('Selected %d genes (correlation > %.2f) \n', length(corgene), cutoff))
    
    corsample <- cor(t(rnaseq[, corgene]), t(proteome[, corgene]))
    align_output$sampleCor_2 <- corsample
    
    cat('Second iteration: Perform data matching ... \n')
    matcher <- attentionStableMarriage(corsample)
    align_output$matching_2 <- matcher
    
    nonmatch <- which(matcher$mismatch_status == 1)
    cat('Second iteration: Found', length(nonmatch), 'pair mismatch\n\n')
  } else {
    align_output$geneCor_2 <- align_output$geneCor_1
    align_output$geneCorPlot_2 <- align_output$geneCorPlot_1
    align_output$sampleCor_2 <- align_output$sampleCor_1
    align_output$matching_2 <- align_output$matching_1
  }
  
  return(align_output)
}




#### Prediction: Partition sex chromosomes genes
partitionSexMatrix <- function(omics, chr_annotate) {
  sexgenes <- getSexGenes(chr_annotate)
  omics_sex  <- omics[intersect(rownames(omics), sexgenes), ]
  return(t(omics_sex))
}

#### Prediction: Get weight of training instance, distributed equally by class
getClassWeight <- function(labels){
  weight <- rep( 1, length(labels) )
  labelclass <- unique(labels)
  weight[labels == labelclass[1]] <- length(labels)/ 2 / sum(labels == labelclass[1])
  weight[labels == labelclass[2]] <- length(labels)/ 2 / sum(labels == labelclass[2])
  return(weight)
}

#### Prediction: Training Regularized Logistic Regression model
trainGLM <- function(msiLabel, inputmtx, alpha){
  weight <- getClassWeight(msiLabel)
  if (sum(weight) - length(msiLabel) > .Machine$double.eps^0.9){
    cat(sprintf('Error: sum(classweight) = %.20f does not equal to length(msiLabel) = %.20f! \n', sum(weight), length(msiLabel)))
  }   # should be nrow(trainset)
  
  # perform cross validation of elasticnet to determine optimum lambda
  cv.glm <- cv.glmnet(as.matrix(inputmtx), msiLabel, family="binomial", weights=weight, alpha=alpha)
  (best_lambda <- cv.glm$lambda.1se)
  fit <- glmnet(as.matrix(inputmtx), msiLabel, family="binomial", weights=weight, alpha=alpha, lambda=best_lambda)
  
  return(fit)
}

#### Prediction: Perform k-fold Cross Validation to predict clinical attributes
trainGLMcv <- function(msiLabel, inputmtx, alpha = 0.3, k = 5){
  inputmtx <- as.matrix(inputmtx)
  flds <- createFolds(msiLabel, k = k, list = TRUE, returnTrain = FALSE)
  predoutput <- data.frame(gender=msiLabel, gender_prob=0)
  
  for (f in 1:k) {
    testidx <- flds[[f]]
    numvar <- 0
    iter   <- 0
    while (numvar < 4 && iter < 50){
      fit1 <- trainGLM(msiLabel[-testidx], inputmtx[-testidx, ], 0.3)
      numvar1 <- sum(coef(fit1) > 0)
      if (numvar1 >= numvar) {
        fit <- fit1
        numvar <- numvar1
      }
      iter <- iter + 1
    }
    
    trainpred <- predict(fit, inputmtx[-testidx, ], type='class')[, 1]
    accuracy  <- sum(trainpred == msiLabel[-testidx]) / length(trainpred)
    cat('  Fold', f, '- Iter', iter, '- Gender predictor model - Selected', numvar, 'variable; ')
    cat(sprintf('Accuracy = %.4f \n', accuracy))
    
    predoutput$gender[testidx] <- predict(fit, inputmtx[testidx, ], type='class')[, 1]
    predoutput$gender_prob[testidx] <- predict(fit, inputmtx[testidx, ], type='response')[, 1]
  }
  return(predoutput)
}

#### Prediction: Comprehensive function for prediction gender
predict_gender <- function(clinical, rna_sex, pro_sex, matcher, omictype1 = 'omics1', omictype2 = 'omics2') {
  predict_outputs <- list()                               # initialize a list for intermediate results

  sample_label <- intersect(rownames(rna_sex), rownames(pro_sex))
  sampleN      <- length(sample_label)
  cat('Number of samples present in both omics data =', sampleN, '\n\n')
  rna_sex <- rna_sex[sample_label, ]
  pro_sex <- pro_sex[sample_label,]
  
  nonmatch <- which(matcher$mismatch_status == 1)
  
  
  ## Perform prediction
  cli_exists <- sample_label %in% clinical$sample
  if (sum(!(cli_exists)) > 0) {
    cat('Error: Clinical information not found for some samples! \n\n')
    return()
  }
  
  rownames(clinical) <- clinical$sample
  traincli <- clinical[sample_label, c('sample', 'gender')]
  clinical$gender <- as.factor(clinical$gender)
  traincli$gender_prob <- apply(traincli, 1, function(x) if (x['gender'] == levels(clinical$gender)[1]) 0 else 1)
  rownames(traincli) <- 1:sampleN
  
  ## First round of prediction (after removing RNA/PRO mismatch samples)
  matchingidx <- setdiff(1:sampleN, nonmatch)
  
  cat('First Iteration: Cross Validation for', omictype1, 'samples ... \n')
  predoutput <- trainGLMcv(traincli$gender[matchingidx], rna_sex[matchingidx, ], 0.3)
  traincli$omics1_gender <- traincli$gender
  traincli$omics1_gender_prob <- traincli$gender_prob
  traincli$omics1_gender[matchingidx] <- predoutput$gender
  traincli$omics1_gender_prob[matchingidx] <- predoutput$gender_prob
  
  cat('First Iteration: Cross Validation for', omictype2, 'samples ... \n')
  predoutput <- trainGLMcv(traincli$gender[matchingidx], pro_sex[matchingidx, ], 0.3)
  traincli$omics2_gender <- traincli$gender
  traincli$omics2_gender_prob <- traincli$gender_prob
  traincli$omics2_gender[matchingidx] <- predoutput$gender
  traincli$omics2_gender_prob[matchingidx] <- predoutput$gender_prob
  
  traincli$avg_omic_gender <- (traincli$omics1_gender_prob + traincli$omics2_gender_prob) / 2
  cli_suspect <- which(abs(traincli$gender_prob - traincli$avg_omic_gender) > 0.6)
  cat('First iteration: Found', length(cli_suspect), 'potential clinical swapping.. \n\n')
  predict_outputs$prediction_1 <- traincli
  
  ## Second round of prediction (after removing RNA/PRO mismatch samples & clinical swapping suspect)
  matchingidx <- setdiff(matchingidx, cli_suspect)
  
  cat('Second Iteration: Training Model for', omictype1, 'gender prediction ... \n')
  numvar <- 0
  iter <- 0
  while (numvar < 4 && iter < 50){
    fit1 <- trainGLM(traincli$gender[matchingidx], rna_sex[matchingidx, ], 0.3)
    numvar1 <- sum(coef(fit1) > 0)
    if (numvar1 > numvar) {
      fit <- fit1
      numvar <- numvar1
    }
    iter <- iter + 1
  }
  cat('  Iter', iter, '- Gender predictor model has', numvar, 'variable --> ')
  traincli$omics1_gender <- predict(fit, as.matrix(rna_sex), type='class')[, 1]
  traincli$omics1_gender_prob <- predict(fit, as.matrix(rna_sex), type='response')[, 1]
  
  trainpred <- predict(fit, rna_sex[matchingidx, ], type='class')[, 1]
  accuracy  <- sum(trainpred == traincli$gender[matchingidx]) / length(trainpred)
  predict_outputs$omics1_model_accuracy <- accuracy
  cat(sprintf('Training Accuracy = %.4f \n', accuracy))
  
  cat('Second Iteration: Training Model for', omictype2, 'gender prediction ... \n')
  numvar <- 0
  iter   <- 0
  while (numvar < 4 && iter < 50){
    fit1 <- trainGLM(traincli$gender[matchingidx], pro_sex[matchingidx, ], 0.3)
    numvar1 <- sum(coef(fit1) > 0)
    if (numvar1 > numvar) {
      fit <- fit1
      numvar <- numvar1
    }
    iter <- iter + 1
  }
  cat('  Iter', iter, '- Gender predictor model has', numvar, 'variable --> ')
  traincli$omics2_gender <- predict(fit, as.matrix(pro_sex), type='class')[, 1]
  traincli$omics2_gender_prob <- predict(fit, as.matrix(pro_sex), type='response')[, 1]
  
  trainpred <- predict(fit, pro_sex[matchingidx,], type='class')[, 1]
  accuracy  <- sum(trainpred == traincli$gender[matchingidx]) / length(trainpred)
  predict_outputs$omics2_model_accuracy <- accuracy
  cat(sprintf('Training Accuracy = %.4f \n', accuracy))
  
  
  traincli$avg_omic_gender <- (traincli$omics1_gender_prob + traincli$omics2_gender_prob)/2
  traincli$difference      <- abs(traincli$gender_prob - traincli$avg_omic_gender)
  cli_suspect <- which(traincli$difference > 0.6)
  cli_suspect <- setdiff(cli_suspect, nonmatch)
  traincli$misannotate <- 0
  traincli$misannotate[cli_suspect] <- 1
  cat('Second iteration: Found', length(cli_suspect), 'clinical swapping.. \n\n')
  predict_outputs$prediction_2 <- traincli
  
  return(predict_outputs)
}




#### Correction: Comprehensive function for label correction
determine_error <- function(corsample, matcher, traincli, omictype1 = 'omics1', omictype2 = 'omics2') {
  sample_label <- rownames(corsample)
  sampleN   <- length(sample_label)
  rankdist  <- computeRankdist(corsample)
  nonmatch  <- which(matcher$mismatch_status == 1)
  
  correct_outputs <- list()
  
  
  #### Determine swapping
  cat('Determine swapping cases ... \n')
  swapped_id <- c()
  for (r in nonmatch){
    if (r %in% swapped_id)    next
    p <- matcher$omics2[r]
    
    if (matcher$omics2[p] != r || matcher$omics2[r] == r)       next
    if (matcher$match_score[r] + matcher$match_score[p] > log(sampleN)*2)   next
    
    swapped_id <- c(swapped_id, r, p)
  }
  correct_outputs$swapped_id <- swapped_id
  
  swapped_sample <- data.frame(id=swapped_id, traincli[swapped_id, c(1:7)], swap_type=rep(0, length(swapped_id)))
  swapped_sample$omics1_error <- abs(swapped_sample$gender_prob - swapped_sample$omics1_gender_prob)
  swapped_sample$omics2_error <- abs(swapped_sample$gender_prob - swapped_sample$omics2_gender_prob)
  
  i <- 1
  while (i < length(swapped_id)) {
    j <- i + 1
    if (swapped_sample$gender[i] == swapped_sample$gender[j]) {
      cat(sprintf('Same sex: %d <--> %d (both %s)... unable to infer which get swapped \n', swapped_sample$id[i], swapped_sample$id[j], swapped_sample$gender[i]))
    } else {
      omics1_error <- swapped_sample$omics1_error[i] + swapped_sample$omics1_error[j]
      omics2_error <- swapped_sample$omics2_error[i] + swapped_sample$omics2_error[j]
      if (omics1_error > omics2_error) {
        cat(sprintf('Swap in %s: %d <--> %d (Error rate = %s swapping (%.3f) vs %s swapping (%.3f) \n', omictype1, swapped_sample$id[i], swapped_sample$id[j], omictype1, omics1_error, omictype2, omics2_error))
        swapped_sample$swap_type[i:j] <- 1
      } else if (omics2_error > omics1_error){
        cat(sprintf('Swap in %s: %d <--> %d (Error rate = %s swapping (%.3f) vs %s swapping (%.3f) \n', omictype2, swapped_sample$id[i], swapped_sample$id[j], omictype1, omics1_error, omictype2, omics2_error))
        swapped_sample$swap_type[i:j] <- 2
      } else {
        cat('Error: Inconsistent sex swapping but same error rate!')
      }
    }
    i <- i + 2
  }
  correct_outputs$swapped_sample <- swapped_sample
  
  
  #### Determine duplication and shifting
  cat('\n')
  cat('Determine Shifting cases ... \n')
  
  dup_shift <- setdiff(nonmatch, swapped_id)
  chains <- list()
  if (length(dup_shift) == 0) {
    cat('No shifting & duplication suspect! Skip... \n')
  } else {
    cat(length(dup_shift), 'samples suspected for duplication and shifting =', dup_shift, '\n')
    shiftdist <- matcher[dup_shift,]
    shiftdist <- shiftdist[order(-shiftdist$match_score), ]
    
    lose_starts <- shiftdist$omics1[shiftdist$match_score > max(2, log(sampleN))]
    lose_ends   <- shiftdist$omics2[shiftdist$match_score > max(2, log(sampleN))]
    
    ### chain identification
    shifttype <- c()
    i <- 1
    for (start in lose_starts){
      chain <- c(start)
      cnext <- start
      while (!(cnext %in% lose_ends)) {
        cnext <- shiftdist$omics1[shiftdist$omics2 == cnext]
        chain <- c(chain, cnext)
      }
      
      clast <- which.max(rankdist[, chain[length(chain)]])
      if (clast != chain[length(chain)])     chain <- c(chain, clast)
      cfirst <- which.max(rankdist[chain[1],])
      if (cfirst != chain[1])    chain <- c(cfirst, chain)
      
      chains[[i]] <- chain
      i <- i + 1
    }
    
    if (sum(!(dup_shift %in% unlist(chains))) == 0){
      cat('All suspected samples are found in chain! \n\n')
    } else {
      cat('Warning! Not all suspected samples are found in chain. Circular shifting suspected! \n\n')
    }
  }
  correct_outputs$shifted_chains <- chains
  
  shifted_sample <- data.frame()
  c <- 0
  for (chain in chains) {
    c <- c + 1
    cat('Chain', c, ': ', paste(chain, collapse= ' '), '\n')
    shifted <- traincli[chain, 1:7]
    shifted <- cbind(data.frame(id=chain, row.names=NULL), shifted)
    shifted$chain_num <- c
    shifted$omics1_error <- abs(shifted$gender_prob - shifted$omics1_gender_prob)
    shifted$omics2_error <- abs(shifted$gender_prob - shifted$omics2_gender_prob)
    shifted$shift_type <- 0       # 0 = unknown; 1 = omics1 swap; 2 = omics2 swap
    
    lenchain <- length(chain)
    omics1_error <- sum(shifted$omics1_error)
    omics2_error <- sum(shifted$omics2_error)
    distfront <- corsample[chain[2], chain[1]]
    distback  <- corsample[chain[lenchain], chain[lenchain-1]]
    if (length(unique(shifted$gender)) == 1) {
      cat('  Same attribute, Inspecting correlation: \n')
      if (distback > distfront) {
        cat(sprintf('    Correlation: %s_%d <-> %s_%d = %.4f \t %s_%d <-> %s_%d = %.4f \n', omictype1, chain[2], omictype2, chain[1], distfront, omictype1, chain[length(chain)], omictype2, chain[length(chain)-1], distback))
        cat(' ', omictype2, 'shift: ', chain[1], paste(chain[3:lenchain], collapse=' '), chain[length(chain)], '\n')
        shifted$shift_type <- 2
      } else {
        cat(sprintf('    Correlation: %s_%d <-> %s_%d = %.4f \t %s_%d <-> %s_%d = %.4f \n', omictype1, chain[2], omictype2, chain[1], distfront, omictype1, chain[length(chain)], omictype2, chain[length(chain)-1], distback))
        cat(' ', omictype1, 'shift: ', chain[1], paste(chain[1:(length(chain)-2)], collapse = ' '), chain[length(chain)], '\n')
        shifted$shift_type <- 1
      }
    } else {
      cat('  Inspecting Attribute: Label =', paste(shifted$gender, sep=' '), '\t', omictype1, '=', paste(shifted$omics1_gender, sep=' '), '\t', omictype2, '=', paste(shifted$omics2_gender, sep=' '), '\n')
      cat(sprintf('    Error rate of %s shifting (%.4f) vs %s shifting (%.4f) \n', omictype1, omics1_error, omictype2, omics2_error))
      if (omics2_error > omics1_error) {          # supposedly distfront < distback
        cat(sprintf('    Correlation: %s_%d <-> %s_%d = %.4f \t %s_%d <-> %s_%d = %.4f \n', omictype1, chain[2], omictype2, chain[1], distfront, omictype1, chain[length(chain)], omictype2, chain[length(chain)-1], distback))
        cat(' ', omictype2, 'shift: ', chain[1], paste(chain[3:lenchain], collapse=' '), chain[length(chain)], '\n')
        shifted$shift_type <- 2
      } else if (omics1_error > omics2_error) {   # supposedly distfront > distback
        cat(sprintf('    Correlation: %s_%d <-> %s_%d = %.4f \t %s_%d <-> %s_%d = %.4f \n', omictype1, chain[2], omictype2, chain[1], distfront, omictype1, chain[length(chain)], omictype2, chain[length(chain)-1], distback))
        cat(' ', omictype1, 'shift: ', chain[1], paste(chain[1:(length(chain)-2)], collapse = ' '), chain[length(chain)], '\n')
        shifted$shift_type <- 1
      }
    }
    
    if (nrow(shifted_sample) == 0) {
      shifted_sample <- shifted
    } else {
      shifted_sample <- rbind(shifted_sample, shifted)
    }
    cat('\n')
  }
  correct_outputs$shifted_sample <- shifted_sample
  return(correct_outputs)
}

#### Correction: Generate Plots
generateIndPlot <- function(probmatrix, pt, omictype1 = 'omics1', omictype2 = 'omics2') {
  sample_label <- rownames(probmatrix)
  plot(probmatrix[pt,], probmatrix[, pt], xlab=sprintf('Correlation of %s (%s) to other %s', omictype1, sample_label[pt], omictype2),
       ylab=sprintf('Correlation of %s (%s) to other %s', omictype2, sample_label[pt], omictype1))
  points(probmatrix[pt, pt], probmatrix[pt, pt], col=4, pch=19)
  text(probmatrix[pt, pt], probmatrix[pt, pt], labels=sample_label[pt], pos=1, col=4)
}

visualizeMismatch <- function(corsample, nonmatch, outputDir, omictype1 = 'omics1', omictype2 = 'omics2') {
  for (nm in nonmatch) {
    png(paste0(outputDir, '/mismatch_', nm, '.png'), width = 400, height = 400)
    generateIndPlot(corsample, nm, omictype1, omictype2)
    dev.off()
  }  
}


generateSwapPlot <- function(probmatrix, pt, swap, omictype1 = 'omics1', omictype2 = 'omics2') {
  sample_label <- rownames(probmatrix)
  plot(probmatrix[pt,], probmatrix[, pt], xlab=sprintf('Correlation of %s (%s) to other %s', omictype1, sample_label[pt], omictype2),
       ylab=sprintf('Correlation of %s (%s) to other %s', omictype2, sample_label[pt], omictype1))
  points(probmatrix[pt, pt], probmatrix[pt, pt], col=4, pch=19)
  text(probmatrix[pt, pt], probmatrix[pt, pt], labels=sample_label[pt], pos=1, col=4)
  points(probmatrix[pt, swap], probmatrix[swap, pt], col=3, pch=19)
  text(probmatrix[pt, swap], probmatrix[swap, pt], labels=sample_label[swap], pos=2, col=3)
}

visualizeSwapping <- function(corsample, swapped_sample, outputDir, omictype1 = 'omics1', omictype2 = 'omics2') {
  idx <- 1
  while (idx < nrow(swapped_sample)) {
    swap_pair <- swapped_sample$id[idx:(idx+1)]
    
    png(paste0(outputDir, '/swap_', paste(swap_pair, collapse='_'), '.png'), width = 400, height = 400)
    generateSwapPlot(corsample, swap_pair[1], swap_pair[2], omictype1, omictype2)
    dev.off()
    png(paste0(outputDir, '/swap_', paste(swap_pair[2:1], collapse='_'), '.png'), width = 400, height = 400)
    generateSwapPlot(corsample, swap_pair[2], swap_pair[1], omictype1, omictype2)
    dev.off()
    
    idx <- idx + 2
  }
}


generateShiftPlot <- function(probmatrix, pt, shift1, shift2, omictype1 = 'omics1', omictype2 = 'omics2') {
  sample_label <- rownames(probmatrix)
  plot(probmatrix[pt,], probmatrix[, pt], xlab=sprintf('Correlation of %s (%s) to other %s', omictype1, sample_label[pt], omictype2),
       ylab=sprintf('Correlation of %s (%s) to other %s', omictype2, sample_label[pt], omictype1))
  points(probmatrix[pt, pt], probmatrix[pt, pt], col=4, pch=19)
  text(probmatrix[pt, pt], probmatrix[pt, pt], labels=sample_label[pt], pos=1, col=4)
  points(probmatrix[pt, shift1], probmatrix[shift1, pt], col=3, pch=19)
  text(probmatrix[pt, shift1], probmatrix[shift1, pt], labels=sample_label[shift1], pos=2, col=3)
  points(probmatrix[pt, shift2], probmatrix[shift2, pt], col=3, pch=19)
  text(probmatrix[pt, shift2], probmatrix[shift2, pt], labels=sample_label[shift2], pos=2, col=3)
}

visualizeShifting <- function(corsample, shifted_sample, traincli, outputDir, omictype1 = 'omics1', omictype2 = 'omics2') {
  for (c in unique(shifted_sample$chain_num)) {
    chain <- shifted_sample$id[shifted_sample$chain_num == c]
    
    clen <- length(chain)
    
    png(paste0(outputDir, '/shift_c', c, '_scatter.png'), width = 400*clen, height = 400, pointsize = 20)
    par(mfrow=c(1,clen))
    for (i in 1:clen){
      if (i == 1) {
        generateSwapPlot(corsample, chain[i], chain[i+1], omictype1, omictype2)
      } else if (i == clen) {
        generateSwapPlot(corsample, chain[i], chain[i-1], omictype1, omictype2)
      } else {
        generateShiftPlot(corsample, chain[i], chain[i-1], chain[i+1], omictype1, omictype2)
      }
    }
    dev.off()
    
    subcor <- corsample[chain, chain]
    edge_yes <- 1:clen
    edge_yes[clen] <- clen^2
    if (shifted_sample$shift_type[shifted_sample$chain_num == c][1] == 1){
      edge_yes[2:(clen-1)] <- (2:(clen-1)-1) * clen + (2:(clen-1)-1)
    } else {
      edge_yes[2:(clen-1)] <- 2:(clen-1) * clen + 2:(clen-1) 
    }
    
    png(paste0(outputDir, '/shift_c', c, '_edge.png'), width = 400, height = 400)
    g <- graph.incidence(subcor, weighted = T)
    V(g)$gender <- as.vector(as.matrix(traincli[chain, c('omics1_gender', 'omics2_gender')]))
    V(g)$color <- c("tomato", "cyan")[as.factor(V(g)$gender)]
    V(g)$shape <- c("square", "circle")[V(g)$type+1]
    V(g)$label.color <- c('blue', 'black')[V(g)$type+1]
    V(g)$label.cex=1
    E(g)$color <- "gray70"
    E(g)$color[edge_yes] <- 'gray30'
    plot(g, edge.width=E(g)$weight*20, layout=layout.bipartite)
    dev.off()
  }
}


generateCorrectedTable <- function(traincli, swapped_sample, shifted_sample, omictype1 = 'omics1', omictype2 = 'omics2') {
  sample_label <- traincli$sample
  sampleN      <- length(sample_label)
  final_tab    <- data.frame(id=1:sampleN, sample=sample_label, Clinical=1:sampleN, omics1=1:sampleN, omics2=1:sampleN)
  colnames(final_tab)[4:5] <- c(omictype1, omictype2)
  
  ## flag gender mislabeling
  cli_suspect <- which(traincli$misannotate == 1)
  final_tab$Clinical[cli_suspect] <- -1
  
  ## flag label for swapping
  idx <- 1
  while (idx < nrow(swapped_sample)) {
    swap_pair <- swapped_sample$id[idx:(idx+1)]
    if (swapped_sample$swap_type[idx] == 0) {
      # if unable to infer, label -1 for both RNA & PRO
      final_tab[swap_pair, omictype1] <- -1
      final_tab[swap_pair, omictype2] <- -1
    } else if (swapped_sample$swap_type[idx] == 1) {
      # if RNAseq get mismatch, switch label of RNAseq
      final_tab[swap_pair, omictype1] <- swap_pair[2:1]
    } else if (swapped_sample$swap_type[idx] == 2) {
      # if RNAseq get mismatch, switch label of RNAseq
      final_tab[swap_pair, omictype2] <- swap_pair[2:1]
    }
    idx <- idx + 2
  }
  
  ## flag label for shifting
  if (nrow(shifted_sample)) {
    for (cnum in 1:max(shifted_sample$chain_num)) {
      chain <- shifted_sample$id[shifted_sample$chain_num == cnum]
      clen  <- length(chain)
      
      if (shifted_sample$shift_type[shifted_sample$chain_num == cnum][1] == 1) {
        # if RNAseq get mismatch, switch label of RNA
        final_tab[chain[2:(clen-1)], omictype1] <- chain[1:(clen-2)]
      } else if (shifted_sample$shift_type[shifted_sample$chain_num == cnum][1] == 2) {
        # if Proteomics get mismatch, switch label of PRO
        final_tab[chain[2:(clen-1)], omictype2] <- chain[3:(clen)]
      }
    }
  }
  
  return(final_tab)
}


#################### Integrated COSMO ####################
run_COSMOR <- function(rna_raw_file, pro_raw_file, cli_input_file, cor_cutoff = 0.5, impute_missing = F, omictype1 = 'omics1', omictype2 = 'omics2') { 
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
  rnaseq <- prpc_omics(raw_rnaseq, chr_annotate, impute_missing = impute_missing)
  proteome <- prpc_omics(raw_proteome, chr_annotate, impute_missing = impute_missing)
  
  #### Write into files
  write.table(chr_annotate, paste0(interDir, chr_annot_file), sep='\t', col.names=T, row.names=F)
  write.table(rnaseq, paste0(interDir, prpc_rna_file), sep='\t', col.names=T, row.names=T)
  write.table(proteome, paste0(interDir, prpc_pro_file), sep='\t', col.names=T, row.names=T)
  
  
  
  
  #################### Pairwise Alignment ####################
  align_outputs <- align_omics(rnaseq, proteome, cutoff = cor_cutoff)
  
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
  
  predict_outputs <- predict_gender(clinical, rna_sex, pro_sex, matcher, omictype1, omictype2)
  
  traincli <- predict_outputs$prediction_2
  cli_suspect <- which(traincli$misannotate == 1)
  cat(length(cli_suspect), 'is suspected to have gender mislabeled! \n\n')
  
  #### write into file
  write.table(traincli, paste0(interDir, attribute_prediction_file), sep='\t')
  
  
  
  
  #################### Label Correction ####################
  ## Determine types of error
  correct_outputs <- determine_error(corsample, matcher, traincli, omictype1, omictype2)
  
  swapped_sample <- correct_outputs$swapped_sample
  shifted_sample <- correct_outputs$shifted_sample
  
  
  ## generate label-corrected table
  final_tab <- generateCorrectedTable(traincli, swapped_sample, shifted_sample, omictype1, omictype2)
  write.table(final_tab, paste0(interDir, final_table_file), sep='\t', row.names=F)
  
  
  ## visualize mislabeled samples
  errorDir <- 'error_sample/'
  dir.create(errorDir, recursive = TRUE, showWarnings = FALSE)
  #visualizeMismatch(corsample, nonmatch, errorDir)
  
  write.table(swapped_sample, paste0(errorDir, 'swapped_samples.tsv'), sep='\t', row.names=F)
  write.table(shifted_sample, paste0(errorDir, 'shifted_samples.tsv'), sep='\t', row.names=F)
  
  visualizeSwapping(corsample, swapped_sample, errorDir, omictype1, omictype2)
  visualizeShifting(corsample, shifted_sample, traincli, errorDir, omictype1, omictype2)
  
  return(final_tab)
}

