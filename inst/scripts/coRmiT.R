#!/usr/bin/env Rscript

################### INITILAIZE
options(scipen = 0.001, 
        digits = 3)

if( Sys.getenv('DEGHUNTER_MODE') == 'DEVELOPMENT' ){

    full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  
                   error=function(e) # works when using R CMD
                  normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', 
                    commandArgs())], '='))[2]))
    main_path_script <- dirname(full.fpath)
    root_path <- file.path(main_path_script, '..', '..')
    custom_libraries <- c("plotting_functions.R", "write_report.R", 
        "general_functions.R",  "statistics_functions.R", 
        "miRNA_RNA_functions.R", "functional_analysis_library.R", 
        "main_cormit.R", "mkinfer_modified.R")
    for (lib in custom_libraries){
        source(file.path(root_path, 'R', lib))
      }
    template_folder <- file.path(root_path, 'inst', 'templates')
    organism_table_path <- file.path(root_path,"inst","external_data", 
        "organism_table.txt")

} else {
    require('ExpHunterSuite')
    root_path <- find.package('ExpHunterSuite')
    template_folder <- file.path(root_path, 'templates')
    organism_table_path <- file.path(root_path, "inst","external_data",
        "organism_table.txt")
}
 
 #### NOTAS PARA EL SCRIPT

################### OPTIONS
option_list <- list(
    optparse::make_option(c("-r", "--RNAseq_folder"), type="character", 
        default="",
        help="ExpHunterSuite RNAseq execution folder"),
    optparse::make_option(c("-m", "--miRNAseq_folder"), type="character", 
        default=NULL,
        help="ExpHunterSuite miRNAseq execution folder"),
    optparse::make_option(c("--selected_targets"), type="character", 
        default=NULL,
        help="One column file indicating external selected targets"),
    optparse::make_option(c("-o", "--output_files"), type="character",
     default=".", 
        help = "Output folder"),
    optparse::make_option(c("-S", "--strategies"), type="character", 
        default="EE,Eh,Ed,hd,hE",
        help = paste0("EE = Anticorrelation between RNAseq modules Eigengene",
            " and miRNAseq modules Eigengene, Eh = Anticorrelation between ",
            "RNAseq modules Eigengene and miRNAseq hub genes, Ed = ",
            "Anticorrelation between RNAseq modules Eigengene and miRNAseq ",
            "DEG expression profiles, hd = Anticorrelation between RNAseq hub",
            " genes and miRNAseq DEG expression profiles, hE = Anticorrelation",
            " between RNAseq hub genes and miRNAseq modules Eigengene. ",
            "Default : %default")),
    optparse::make_option(c("-d", "--databases"), type="character", 
        default=paste0("targetscan,mirdb,diana_microt,elmmo,microcosm,",
            "miranda,pictar,pita"),
        help = paste0("Set prediction databases included on multiMiR",
            " to use as gold standard.")),
    optparse::make_option("--tag_filter target_tag,miRNA_tag", type="character", 
        default=paste0("putative,putative"),
        help = paste0("Set filter type for RNAseq and miRNAseq input data by comma separated string. Available options are: 'prevalent' to use only PREVALENT_DEG. 'all_possible' to use all PREVALENT_DEG and POSIBLE_DEG. 'putative' to use PREVALENT_DEG, POSIBLE_DEG and miRNA or RNAs that correlates with their expression profile.")),
    optparse::make_option(c("-C","--corr_coef"), type="character", 
        default="pearson",
        help = paste0("Correlation method: chose between pearson, spearman and kendall")),
    optparse::make_option(c("--organism"), type ="character", default="hsa",
        help= paste0("Reference organism to use. Available 'hsa' for human",
            " and 'mmu' for mouse.")),
    optparse::make_option(c("-M", "--multimir_db"), type ="character", 
        default=NULL,
        help= "Indicate .RData file with parsed multiMiR data."),
    optparse::make_option(c("-p", "--p_val_cutoff"), type="double", 
        default=0.05,
        help="Correlation P value threshold . Default=%default"),
    optparse::make_option("--corr_type TYPE", type="character", 
        default=paste0("lower"),
        help = paste0("Set if correlations are [lower] or [higher] than the --p_val_cutoff. Default=%default")),
    optparse::make_option(c("-c", "--corr_cutoff"), type="double", 
        default=-0.7,
        help="Correlation threshold . Default=%default"),
    optparse::make_option(c("-s", "--sample_proportion"), type="double", 
        default=0.01,
        help="Score distribution sample proportion. Default=%default"),
    optparse::make_option(c("-P", "--permutations"), type="double", 
        default=50,
        help="Permutations of random tests. Default=%default"),
    optparse::make_option(c("-u", "--filter_db_theshold"), type="double", 
        help=paste0("Minimun predicted databases that must ",
            "support a pair. Default=%default"),
        default=0),
    optparse::make_option(c("-f", "--database_to_filter"), type="character", 
        default=paste0("targetscan,mirdb,diana_microt,elmmo,microcosm,",
            "miranda,pictar,pita"),
        help=paste0("Predicted databases that must ",
            "support a pair on --filter_db_threshold. Default=%default")),
    optparse::make_option(c("--module_membership_cutoff"), type="double", 
        default=0.7,
        help=paste0("This script reject genes with lower module membership to",
            " their modules. Default=%default")),
    optparse::make_option(c("-R", "--report"), ,type = "character", 
        default="miRNA_RNA_comparison.html",
        help="Name of the html file. Default : %default"),
    optparse::make_option(c("-t", "--translation_file"), type = "character", 
        default = NULL,
        help = paste0('Two columns (\t) file with miRNA names translation: ',
            'Same IDs as input on first column, and miRBase ID on second',
            ' column (MIMAT000000)')),
    optparse::make_option(c("-T", "--translate_ensembl"), type="logical", 
        default=FALSE, action = "store_true",
        help = paste0('Translate mRNA Ensembl ID to ENTREZ ID and GENESYMBOL',
            ' and include the translation in output.')),
    optparse::make_option(c("--compare_pred_scores"), type="logical", 
        default=FALSE, action = "store_true",
        help = 'Compute predition scores comparison. Default=%default'),
    optparse::make_option(c("--mapping_output"), type="character", 
        default="Target_log2FC",
        help = "Select the output column to show in functional report: Predicted_DB_count, Validated_DB_count, Correlation, Target_log2FC, miRNA_log2FC. Default=%default")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


check_and_create_dir(opt$output_files)

miRNA_cor_results <-  coRmiT( 
    RNAseq_folder=opt$RNAseq_folder,
    miRNAseq_folder=opt$miRNAseq_folder,
    output_files=opt$output_files,
    strat_names=unlist(strsplit(opt$strategies, ",")),
    organism=opt$organism,
    multimir_db=opt$multimir_db,
    sample_proportion = opt$sample_proportion,
    p_val_cutoff=opt$p_val_cutoff,
    corr_cutoff=opt$corr_cutoff,
    MM_cutoff=opt$module_membership_cutoff,
    permutations = opt$permutations,
    report=opt$report,
    databases = opt$databases,
    translate_ensembl = opt$translate_ensembl,  
    mc_cores = opt$mc_cores,
    tag_filter = opt$tag_filter,
    filter_db_theshold = opt$filter_db_theshold,
    database_to_filter = opt$database_to_filter,
    translation_file = opt$translation_file,
    organism_table_path = organism_table_path, 
    corr_type = opt$corr_type,
    corr_coef = opt$corr_coef,
    selected_targets_file = opt$selected_targets,
    template_folder = template_folder,
    compare_pred_scores = opt$compare_pred_scores    )
miRNA_cor_results$all_pairs$RNAseq_mod <- miRNA_cor_results$RNAseq$DH_results[match(miRNA_cor_results$all_pairs$RNAseq, miRNA_cor_results$RNAseq$DH_results$gene_name),"Cluster_ID"]

miRNA_cor_results$weighted_c_table <- NULL
miRNA_cor_results <- c(miRNA_cor_results, 
    list(report_name =opt$report,
         template_folder =template_folder,
         output_files =normalizePath(opt$output_files),
         p_val_cutoff =opt$p_val_cutoff,
         corr_cutoff =opt$corr_cutoff,
         sample_proportion =opt$sample_proportion,
         mapping_output = opt$mapping_output))


do.call("write_miRNA_cor_report", miRNA_cor_results)