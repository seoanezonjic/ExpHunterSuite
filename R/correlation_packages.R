#' @importFrom DESeq2 DESeqDataSetFromMatrix rlog
#' @importFrom WGCNA allowWGCNAThreads
#' @importFrom diffcoexp diffcoexp
#' @importFrom utils write.table
#' @importFrom SummarizedExperiment assay
analysis_diffcoexp <- function(data, path, target) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, 
                                          colData = target["treat"], 
                                          design = ~ treat)
    rld <- DESeq2::rlog(dds, blind=FALSE)
    mat_counts <- SummarizedExperiment::assay(rld)
    treat <- mat_counts[,target$treat=="Treat"]
    control <- mat_counts[,target$treat=="Ctrl"]

    WGCNA::allowWGCNAThreads()
    res <- diffcoexp::diffcoexp(exprs.1 = control,
                                exprs.2 = treat, 
                                r.method = "spearman", 
                                q.dcgth=1)

    utils::write.table(res$DCLs, file=file.path(path, "DCLs.txt"), sep="\t", 
                       quote=FALSE, row.names=FALSE)
    utils::write.table(res$DCGs, file=file.path(path, "DCGs.txt"), sep="\t", 
                       quote=FALSE, row.names=FALSE)
    return(res)
}


build_design_for_WGCNA <- function(target, 
                                   string_factors=NULL, 
                                   numeric_factors=NULL){
# FOR WGCNA: Check that the appropriate factor columns can be found in the 
# target file and makes a data frame with the specified factor
    if(!is.null(string_factors)){ 

        target_string_factors <- target[,colnames(target) %in% string_factors,
                                                    drop = FALSE]
        target_string_factors <- data.frame(sapply(target_string_factors, 
                                                   as.factor), 
                                            stringsAsFactors=TRUE)
        return(target_string_factors)
    }

    if(!is.null(numeric_factors)){
        if(sum(nchar(numeric_factors) == 0)) return("") 

        target_numeric_factors <- target[,colnames(target) %in% 
                                              numeric_factors,
                                              drop = FALSE]
        return(target_numeric_factors)
    }
    
    return(NULL)
}



perform_WGCNA_combinations <- function(WGCNA_all=FALSE, 
                    WGCNA_input, 
                    index_treatmn_cols, 
                    index_control_cols, 
                    path, 
                    target_numeric_factors, 
                    target_string_factors, 
                    WGCNA_memory, 
                    WGCNA_deepsplit, 
                    WGCNA_detectcutHeight, 
                    WGCNA_mergecutHeight, 
                    WGCNA_min_genes_cluster, 
                    WGCNA_blockwiseNetworkType, 
                    WGCNA_blockwiseTOMType, 
                    WGCNA_minCoreKME, 
                    WGCNA_minCoreKMESize, 
                    WGCNA_minKMEtoStay,
                    corType){
    results <- list()
    if(WGCNA_all == TRUE) { #TODO => Este bloque de código es repetitivo. 
        WGCNA_input_treatment <- WGCNA_input[, index_treatmn_cols]
        WGCNA_input_control <- WGCNA_input[, index_control_cols]
        WGCNA_treatment_path <- file.path(path, "Treatment_only_data")
        dir.create(WGCNA_treatment_path)
        message('Performing WGCNA correlation analysis for treated samples')
        results[['WGCNA_treatment']] <- analysis_WGCNA(
                       data=WGCNA_input_treatment,
                       path=WGCNA_treatment_path,
                       target_numeric_factors=target_numeric_factors,
                       target_string_factors=target_string_factors,
                       WGCNA_memory=WGCNA_memory,
                       WGCNA_deepsplit=WGCNA_deepsplit,
                       WGCNA_detectcutHeight=WGCNA_detectcutHeight,
                       WGCNA_mergecutHeight=WGCNA_mergecutHeight,
                       WGCNA_min_genes_cluster=WGCNA_min_genes_cluster,
                       cor_only=TRUE, 
                       blockwiseNetworkType = WGCNA_blockwiseNetworkType, 
                       blockwiseTOMType = WGCNA_blockwiseTOMType,
                       WGCNA_minCoreKME = WGCNA_minCoreKME,
                       WGCNA_minCoreKMESize = WGCNA_minCoreKMESize,
                       WGCNA_minKMEtoStay = WGCNA_minKMEtoStay,
                       corType = corType
        )

        WGCNA_control_path <- file.path(path, "Control_only_data")
        dir.create(WGCNA_control_path)

        message('Performing WGCNA correlation analysis for control samples\n')
        results[['WGCNA_control']] <- analysis_WGCNA(data=WGCNA_input_control,
                       path=WGCNA_control_path,
                       target_numeric_factors=target_numeric_factors,
                       target_string_factors=target_string_factors,
                       WGCNA_memory=WGCNA_memory,
                       WGCNA_deepsplit=WGCNA_deepsplit,
                       WGCNA_detectcutHeight=WGCNA_detectcutHeight,
                       WGCNA_mergecutHeight=WGCNA_mergecutHeight,
                       WGCNA_min_genes_cluster=WGCNA_min_genes_cluster,
                       cor_only=TRUE, 
                       blockwiseNetworkType = WGCNA_blockwiseNetworkType, 
                       blockwiseTOMType = WGCNA_blockwiseTOMType,
                       WGCNA_minCoreKME = WGCNA_minCoreKME,
                       WGCNA_minCoreKMESize = WGCNA_minCoreKMESize,
                       WGCNA_minKMEtoStay = WGCNA_minKMEtoStay,
                       corType = corType)
    }
        

    message('Performing WGCNA correlation analysis for all samples\n')
    results[['WGCNA_all']] <- analysis_WGCNA(data=WGCNA_input,
                    path=path,
                    target_numeric_factors=target_numeric_factors,
                    target_string_factors=target_string_factors,
                    WGCNA_memory=WGCNA_memory,
                    WGCNA_deepsplit=WGCNA_deepsplit,
                    WGCNA_detectcutHeight=WGCNA_detectcutHeight,
                    WGCNA_mergecutHeight=WGCNA_mergecutHeight,
                    WGCNA_min_genes_cluster=WGCNA_min_genes_cluster,
                    cor_only=FALSE, 
                    blockwiseNetworkType = WGCNA_blockwiseNetworkType, 
                    blockwiseTOMType = WGCNA_blockwiseTOMType,
                    WGCNA_minCoreKME = WGCNA_minCoreKME,
                    WGCNA_minCoreKMESize = WGCNA_minCoreKMESize,
                    WGCNA_minKMEtoStay = WGCNA_minKMEtoStay,
                    corType = corType)        
    return(results)
}





#' @importFrom utils assignInNamespace write.table
#' @importFrom WGCNA pickSoftThreshold blockwiseModules labels2colors 
#' plotDendroAndColors moduleEigengenes corPvalueStudent 
#' binarizeCategoricalVariable orderMEs
#' @importFrom grDevices pdf dev.control recordPlot dev.off
#' @importFrom graphics par text abline
#' @importFrom stats cor
#' @importFrom ggplot2 ggplot scale_x_continuous scale_y_continuous 
#' geom_count aes_string
analysis_WGCNA <- function(data,
                          path, 
                          target_numeric_factors, 
                          target_string_factors, 
                          WGCNA_memory, WGCNA_deepsplit, 
                          WGCNA_detectcutHeight, 
                          WGCNA_mergecutHeight, 
                          WGCNA_min_genes_cluster, 
                          cor_only, 
                          blockwiseNetworkType, 
                          blockwiseTOMType, 
                          WGCNA_minCoreKME, 
                          WGCNA_minCoreKMESize, 
                          WGCNA_minKMEtoStay,
                          corType) {

    if(is.null(WGCNA_minCoreKMESize)){
        WGCNA_minCoreKMESize <- WGCNA_min_genes_cluster/3
    } 

    data <- t(data)
    nSamples <- nrow(data)

    ####################################################################
    ## THRESHOLDING - EFFECTS OF BETA ON TOPOLOGY AND AUTO SELECTION
    ####################################################################
    utils::assignInNamespace(x="..minNSamples", 
                            value=3, #Overwrite harcoded limit from 4 to 3
                            ns="WGCNA")
    powers <- c(seq(10), seq(from = 12, to=30, by=2))

    bicor <- WGCNA::bicor

    corFnc <- cor
    if(corType == "bicor")
        corFnc <- bicor

    sft <- WGCNA::pickSoftThreshold(data, powerVector = powers, verbose = 5, corFnc = corFnc)

    # Calculate Power automatically
    sft_mfs_r2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
    min_pow_ind <- which(sft_mfs_r2 > 0.9)[1] # First time the value passes 0.9

    in_advance <- 2
    if(is.na(min_pow_ind)) {
        # Plan B to calculate power
        for(i in seq(length(sft_mfs_r2)-in_advance)) {
            vals <- sft_mfs_r2[i:(i+in_advance)]
            diffs <- abs(diff(vals))
            # Condition to meet
            if(sum(diffs < 0.01) == in_advance & sft_mfs_r2[i] > 0.8) {
              min_pow_ind <- i
                 break
            }
        }
    }

    if(is.na(min_pow_ind)) {
        warning("Could not obtain a valid power (beta) value for WGCNA",
            " so the default of 30 will be used - proceed with caution")
        # assumes 30 will be the largest testable power
        min_pow_ind <- length(sft_mfs_r2) 
    }
    pow <- sft$fitIndices[min_pow_ind, "Power"]
    # cor <- WGCNA::cor # TO CORRECT FUNCTION OVERRIDE DUE TO OTHER PACKAGES

    grDevices::pdf(file.path(path, "thresholding.pdf"))
        grDevices::dev.control(displaylist="enable")
        graphics::par(mfrow = c(1,2))
        cex1 <- 0.9
    # Scale-free topology fit index as function of the soft-thresholding power
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",
             ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence"));
        graphics::text(sft$fitIndices[,1], 
                      -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                      labels=powers,cex=cex1,col="red");
        # this line corresponds to using an R^2 cut-off of h
        graphics::abline(h=0.90, col="red")
        graphics::abline(h=0.80, col="red", lty="dashed")
        graphics::abline(v=pow, col="black", lty="dotted")

        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
        graphics::text(sft$fitIndices[,1], sft$fitIndices[,5], 
                       labels=powers, cex=cex1,col="red")
        graphics::abline(v=pow, col="black", lty="dotted")

        power_threshold_effects <- grDevices::recordPlot()
    grDevices::dev.off()

    ####################################################################
    ## CLUSTER SAMPLES TO GENERATE MODULES
    ####################################################################

    # Have to improve this as generating fork.
    tom_file_base <- "generatedTOM"
    if(length(list.files(path=path, 
                         pattern=paste0(tom_file_base, ".*\\.RData"))) > 0) {
        loadTOM_TF <- TRUE
        saveTOM_TF <- FALSE
    } else {
        loadTOM_TF <- FALSE
        saveTOM_TF <- TRUE
    }

    cor <- WGCNA::cor # This issue has been reported to the WGCNA author
    net <- WGCNA::blockwiseModules(data, power = pow,
        # Increase to memory limit in order to obtain more realistic results
                            maxBlockSize = WGCNA_memory, 
                            minModuleSize = WGCNA_min_genes_cluster,
                            deepSplit = WGCNA_deepsplit, 
                            detectCutHeight = WGCNA_detectcutHeight,
                            reassignThreshold = 0, 
                            mergeCutHeight = WGCNA_mergecutHeight,
                            numericLabels = TRUE, 
                            pamRespectsDendro = FALSE,
                            loadTOM = loadTOM_TF,    
                            saveTOM = saveTOM_TF,
                            saveTOMFileBase = file.path(path, tom_file_base), 
                            verbose = 5, 
                            networkType = blockwiseNetworkType,
                            TOMType = blockwiseTOMType,
                            minCoreKME = WGCNA_minCoreKME,
                            minCoreKMESize = WGCNA_minCoreKMESize, 
                            minKMEtoStay = WGCNA_minKMEtoStay,
                            corType = corType)
    cor<-stats::cor

    moduleColors <- WGCNA::labels2colors(net$colors)

    # Plot the dendrogram and the module colors underneath
    grDevices::pdf(file.path(path, 'wcgnaModules.pdf'))
        WGCNA::plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
    grDevices::dev.off()

    MEs <- WGCNA::moduleEigengenes(data, net$colors)$eigengenes

    gene_module_cor <- as.data.frame(WGCNA::cor(data, MEs, use = "p"));
    gene_module_cor_p <- as.data.frame(WGCNA::corPvalueStudent(
                                        as.matrix(gene_module_cor), nSamples));
    colnames(gene_module_cor_p) <- colnames(gene_module_cor) <- gsub("ME", 
                                      "Cluster_", colnames(gene_module_cor_p) )

    utils::write.table(gene_module_cor, 
                       file=file.path(path, "gene_MM.txt"), 
                       sep="\t", quote=FALSE)
    utils::write.table(gene_module_cor_p, 
                       file=file.path(path, "gene_MM_p_val.txt"), 
                       sep="\t", quote=FALSE)
    utils::write.table(MEs, 
                      file=file.path(path, "eigen_values_per_samples.txt"), 
                      sep="\t", quote=FALSE)

    if(cor_only == TRUE) {
        return("cor_only")
    }

    ####################################################################
    ## PROCESS TRAIT DATA
    ####################################################################

    trait <- data.frame(row.names=rownames(data))

    if(is.data.frame(target_numeric_factors)) {
        trait <- data.frame(trait, target_numeric_factors)
    }
    if(is.data.frame(target_string_factors)) {
    # Code to convert string factors to numeric (1 vs. 0 for each category).
    # NOTE YOU NEED A MINIMUM COUNT OF 2 FOR EACH CATEGORY IN THE FACTORS
        binarized_string_factors <- lapply(names(target_string_factors), 
            function(factor_name) {
                string_factor <- target_string_factors[[factor_name]]
                binarized_table <- WGCNA::binarizeCategoricalVariable(
                                  string_factor,
                                minCount = 1,
                                includePairwise = FALSE,
                                includeLevelVsAll = TRUE)
                colnames(binarized_table) <- paste(factor_name, 
                                               levels(string_factor),sep="_")
                return(binarized_table)
        })
        trait <- data.frame(trait, binarized_string_factors)
    }

    ###################################################################
    ## CLUSTER DATA
    ###################################################################

    ####################################################################
    ## CLUSTER MODULES WITH TRAITS AND PRODUCE HEATMAP
    ####################################################################

    MEs_colors <- MEs <- WGCNA::moduleEigengenes(data, net$colors)$eigengenes
    ME_numeric <- as.numeric(gsub("ME", "", colnames(MEs)))
    colnames(MEs_colors) <- paste0("ME", WGCNA::labels2colors(ME_numeric))

    moduleTraitCor <- WGCNA::cor(MEs_colors, trait, use = "p")
    moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
 
    textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)
    colnames(moduleTraitCor) <- names(trait)
    rownames(moduleTraitCor) <- names(MEs)
    moduleTraitCor_df <- as.data.frame(as.table(moduleTraitCor))
    colnames(textMatrix) <- names(trait)
    rownames(textMatrix) <- names(MEs)
    textMatrix_df <- as.data.frame(as.table(textMatrix))
    colnames(moduleTraitCor_df) <- c("Module","Trait","Correlation")
    moduleTraitCor_df$Text_correlation <- textMatrix_df[["Freq"]]
    ####################################################################
    # Produce trait AND module correlations and p-values
    ####################################################################
    ME_numeric_traits <- WGCNA::orderMEs(cbind(MEs, trait))

    trait_and_module_cor <- WGCNA::cor(ME_numeric_traits, use = "p")
    trait_and_module_cor_p <- WGCNA::corPvalueStudent(trait_and_module_cor, 
                                                      nSamples)
    trait_module_cor_val <- as.data.frame(as.table(trait_and_module_cor))
    trait_module_p_val <- as.data.frame(as.table(trait_and_module_cor_p))

    trait_module_all_vals <- data.frame(trait_module_cor_val, 
                                        trait_module_p_val[,3])
    names(trait_module_all_vals) <- c("A", "B", "correlation", "p-value")
    utils::write.table(trait_module_all_vals, 
                       file=file.path(path, "trait_and_module_all_vals.txt"), 
                       sep="\t", quote=FALSE, row.names=FALSE)

    ####################################################################
    ## Produce Tables
    ####################################################################
    # Report tables: 
    # Genes per module 
    gene_module_cor <- as.data.frame(WGCNA::cor(data, MEs, use = "p"))
    gene_module_cor_p <- as.data.frame(WGCNA::corPvalueStudent(
                                        as.matrix(gene_module_cor), nSamples))
    colnames(gene_module_cor_p) <- colnames(gene_module_cor) <- gsub("ME", 
                                       "Cluster_", colnames(gene_module_cor_p))
    # Genes per trait
    gene_trait_cor <- as.data.frame(WGCNA::cor(data, trait, use = "p"));
    gene_trait_cor_p <- as.data.frame(WGCNA::corPvalueStudent(
                                          as.matrix(gene_trait_cor), nSamples))
    # Module per trait 
    # (also produced above for the plot - should give smae results.)
    module_trait_cor <- WGCNA::cor(MEs, trait, use = "p")
    module_trait_cor_p <- WGCNA::corPvalueStudent(module_trait_cor, nSamples)
    row.names(module_trait_cor) <- row.names(module_trait_cor_p) <- gsub("ME", 
                                     "Cluster_", row.names(module_trait_cor_p))

    utils::write.table(trait, 
       file=file.path(path, "sample_trait.txt"), 
       sep="\t", quote=FALSE)
    utils::write.table(gene_trait_cor, 
           file=file.path(path, "gene_trait.txt"), sep="\t", quote=FALSE)
    utils::write.table(gene_trait_cor_p, 
           file=file.path(path, "gene_trait_p_val.txt"), sep="\t", quote=FALSE)
    utils::write.table(module_trait_cor, 
           file=file.path(path, "module_trait.txt"), sep="\t", quote=FALSE)
    utils::write.table(module_trait_cor_p, 
           file=file.path(path, "module_trait_p_val.txt"), 
           sep="\t", quote=FALSE)

    cluster_ID<- net$colors
    Cluster_MM <- unlist(lapply(names(cluster_ID), function(x){
        gene_module_cor[x, paste0("Cluster_", cluster_ID[x])]
    }))
    Cluster_MM_pval <- unlist(lapply(names(cluster_ID), function(x){
        # Get the MM value
        gene_module_cor_p[x, paste0("Cluster_", cluster_ID[x])] 
    }))

    # Plot cluster ID vs. ID of cluster with lowest MM p-value for each gene
    MM_Cluster_ID <- apply(gene_module_cor_p, 
                           1, 
                           function(x) which(x == min(x))) - 1 
                           # Cluster ID with minimum MM p value
    cluster_vs_MM <- ggplot2::ggplot(as.data.frame(cbind(cluster_ID, 
                                                         MM_Cluster_ID)), 
                              ggplot2::aes_string(x = "cluster_ID" ,
                                                  y = "MM_Cluster_ID")) +
                              ggplot2::geom_count() + 
                              ggplot2::scale_x_continuous("Cluster ID", 
                                labels = function(x) paste0("Cluster_", x)) + 
                   ggplot2::scale_y_continuous("Cluster with highest MM value", 
                                labels = function(x) paste0("Cluster_", x))

    # NEW - get the p value for the cluster id, not the lowest p-value

    return(list(
          gene_cluster_info = data.frame(ENSEMBL_ID = names(cluster_ID), 
                                         Cluster_ID = as.numeric(cluster_ID),
                                         Cluster_MM = Cluster_MM,
                                         Cluster_MM_pval = Cluster_MM_pval),
          package_objects=list(gene_module_cor = gene_module_cor, 
                               gene_module_cor_p = gene_module_cor_p,
                               gene_trait_cor = gene_trait_cor, 
                               gene_trait_cor_p = gene_trait_cor_p, 
                               module_trait_cor = module_trait_cor, 
                               module_trait_cor_p = module_trait_cor_p),
          plot_objects=list(sorted_colours = WGCNA::labels2colors(ME_numeric),
                            trait_vs_module = moduleTraitCor_df,
                            trait_and_module = ME_numeric_traits,
                            power_threshold_effects = power_threshold_effects,
                            cluster_vs_MM = cluster_vs_MM)
            )
    )
}

#' corM2igraph
#'
#' `corM2igraph` Creates an igraph out of a correlation matrix.
#'
#' @param corM correlation matrix.
#' @param cor_abs_thr Correlation threshold (absolute value)
#' @returns An Igraph.'
#' @examples
#'  \dontrun{
#'      corM2igraph(matrix, 0.75)
#'  }
#' @export
corM2igraph <- function(corM, cor_abs_thr = 0.75){
  cor_df <- data.frame(row=rownames(corM)[row(corM)[upper.tri(corM)]], 
           col=colnames(corM)[col(corM)[upper.tri(corM)]], 
           corr=corM[upper.tri(corM)])
  cor_df$corr_type <- ifelse(c(cor_df$corr) > 0, "corr","anticorr")
  cor_df$corr <- abs(cor_df$corr)
  cor_df_fil <- cor_df[cor_df$corr >= cor_abs_thr,]
  nodes <- data.frame(nodes = unique(c(cor_df_fil$row, cor_df_fil$col)))
  return(list(nodes = nodes, edges = cor_df_fil))
}
    
#' @importFrom stats cor
analysis_diff_correlation <- function(all_genes_stats, data, control, treat, PCIT_filter=FALSE){
    # DE and PIF calculation, see "Microarray data processing, normalization and differential expression" in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    # DE => Diff exp, PIF =>Phenotypic Impact Factor
    nonZeroVal <- 10^-5
    if(sum(data==0) > 0) data[data == 0] <- replicate(sum(data == 0), jitter(nonZeroVal))
    if(sum(control==0) > 0) control[control == 0] <- replicate(sum(control == 0), jitter(nonZeroVal))
    if(sum(treat==0) > 0) treat[treat == 0] <- replicate(sum(treat == 0), jitter(nonZeroVal))
    metrics <- all_genes_stats['logFC_DESeq2']
    names(metrics) <- 'de'
    data <- log(data, 2)
    control <- log(control, 2)
    treat <- log(treat, 2)

    average <- apply(data, 1, mean) # Ai in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    treat_average <- apply(treat, 1, mean) # Ai in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    control_average <- apply(control, 1, mean) # Ai in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    metrics$average <- average
    metrics$treat_average <- treat_average
    metrics$control_average <- control_average
    metrics$pif <- average * metrics$de

    if(!PCIT_filter) {
        control <- stats::cor(t(control)) # Compute correlations
        treat <- stats::cor(t(treat)) # Compute correlations
    }else{
        # JRP Removed CeTF as it depended on GenomeTools which has been deprecated
        # control <- CeTF::PCIT(control, tolType = "mean")$adj_sig
        # treat <- CeTF::PCIT(treat, tolType = "mean")$adj_sig
    }

    #diferential_wiring: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    dw <- control - treat
    prevalent_gene_stats <- all_genes_stats[all_genes_stats$genes_tag == 'PREVALENT_DEG', ]
    if(nrow(prevalent_gene_stats) == 0){
        message('Using POSSIBLE_DEG instead PREVALENT_DEG')
        prevalent_gene_stats <- all_genes_stats[all_genes_stats$genes_tag == 'POSSIBLE_DEG', ]
    }
    metrics$rif1 <- scale(get_rif(1, metrics, prevalent_gene_stats, dw, control, treat)) # Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    metrics$rif2 <- scale(get_rif(2, metrics, prevalent_gene_stats, dw, control, treat)) # Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    #Differential hubbing (differential conectivity): https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
    metrics$ctrl_cn <- apply(control,1,function(x) length(which(abs(x) > 0.8)))
    metrics$treat_cn <- apply(treat,1,function(x) length(which(abs(x) > 0.8)))
    metrics$diff_cn <- metrics$ctrl_cn - metrics$treat_cn
    # As Horvath (WGCNA): DiffK (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1998880/pdf/335_2007_Article_9043.pdf)
    k_control <- metrics$ctrl_cn/max(metrics$ctrl_cn)
    k_treat <- metrics$treat_cn/max(metrics$treat_cn)
    metrics$diffK <- k_control - k_treat
    return(metrics)
}

get_rif <- function(rif_version, metrics, deg_stats, dw, control_r, treat_r){
    genes <- row.names(metrics)
    rif <- NULL
    if(rif_version == 1){
        rif <- unlist(lapply(genes, compute_rif1, data=metrics, deg_stats=deg_stats, dw=dw))
    }else if(rif_version == 2){
        rif <- unlist(lapply(genes, compute_rif2, data=metrics, deg_stats=deg_stats, control_r=control_r, treat_r=treat_r))
    }
    return(rif)
}

# Regulatory Impact Factor (RIF):  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000382
compute_rif1 <- function(g, data, deg_stats, dw){ #eq 4
    #adapted: It will be computed for all genes and not only TF. Also, all genes DW will be take into account instead of only significant DE genes
    # When implemented significant RNAseq DEgenes change to DE only.
    # Also, implement eq 5 and perform the mean between two eqs. 
    # A posterior paper, saids that the two measures are very different.
    deg_pif <- data[row.names(deg_stats),"pif" ]
    rif <- sum(deg_pif * dw[row.names(deg_stats), g]^2)/length(deg_pif)
    return(rif)
}

compute_rif2 <- function(g, data, deg_stats, control_r, treat_r){ #eq 5
    #adapted: It will be computed for all genes and not only TF. Also, all genes DW will be take into account instead of only significant DE genes
    # When implemented significant RNAseq DEgenes change to DE only.
    # Also, implement eq 5 and perform the mean between two eqs. 
    # A posterior paper, saids that the two measures are very different.
    rif_treat <- (data[row.names(deg_stats),"treat_average" ] * treat_r[row.names(deg_stats), g])^2
    rif_control <- (data[row.names(deg_stats),"control_average" ] * control_r[row.names(deg_stats), g])^2
    rif <- sum(rif_control - rif_treat)/length(rif_treat)
    return(rif)
}
