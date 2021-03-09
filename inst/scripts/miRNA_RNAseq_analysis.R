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
        "main_miRNA_RNAseq_analysis.R")
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
        default=NULL,
        help="ExpHunterSuite RNAseq execution folder"),
    optparse::make_option(c("-m", "--miRNAseq_folder"), type="character", 
        default=NULL,
        help="ExpHunterSuite miRNAseq execution folder"),
    optparse::make_option(c("-o", "--output_files"), type="character",
     default=".", 
        help = "Output folder"),
    optparse::make_option(c("-a", "--approaches"), type="character", 
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
    optparse::make_option(c("--organism"), type ="character", default="hsa",
        help= paste0("Reference organism to use. Available 'hsa' for human",
            " and 'mmu' for mouse.")),
    optparse::make_option(c("-M", "--multimir_db"), type ="character", 
        default=NULL,
        help= "Indicate .RData file with parsed multiMiR data."),
    optparse::make_option(c("--only_known"), type="logical", default=FALSE, 
        action = "store_true",
        help= paste0("Perform analysis only using known miRNAs. Faster but ",
            "you lose information of novel miRNAs. Only available when a ",
            "translation file is given")),
    optparse::make_option(c("-p", "--p_val_cutoff"), type="double", 
        default=0.05,
        help="Correlation P value threshold . Default=%default"),
    optparse::make_option(c("-P", "--permutations"), type="double", default=10,
        help="Permutations of random tests. Default=%default"),
    optparse::make_option(c("-c", "--corr_cutoff"), type="double", default=-0.7,
        help="Correlation threshold . Default=%default"),
    optparse::make_option(c("-u", "--unmaintained_filter"), type="double", 
        default=1,
        help=paste0("Minimun predicted and unmaintained databases that must ",
            "support a pair. Default=%default")),
    optparse::make_option(c("-C", "--mc_cores"), type="double", default=1,
        help="Dist. Default=%default"),
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
            ' and include the translation in output.'))
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))



dir.create(opt$output_files, recursive = TRUE)


miRNA_cor_results <- miRNA_RNAseq_analysis(
    RNAseq_folder=opt$RNAseq_folder,
    miRNAseq_folder=opt$miRNAseq_folder,
    output_files=opt$output_files,
    approaches=opt$approaches,
    organism=opt$organism,
    multimir_db=opt$multimir_db,
    p_val_cutoff=opt$p_val_cutoff,
    corr_cutoff=opt$corr_cutoff,
    module_membership_cutoff=opt$module_membership_cutoff,
    permutations = opt$permutations,
    report=opt$report,
    databases = opt$databases,
    translate_ensembl = opt$translate_ensembl,
    mc_cores = opt$mc_cores,
    filter_unmaintained = opt$unmaintained_filter,
    translation_file = opt$translation_file,
    #ABAJO file.path(root_path, "R", "organism_table.txt"),
    organism_table_path = organism_table_path, 
    template_folder = template_folder #file.path(root_path, "inst", "templates")
    )
write_miRNA_cor_report(miRNA_cor_results = miRNA_cor_results, 
                        template_folder = template_folder, 
                        output_files = opt$output_files, 
                        report_name = "miRNA_RNA_comparison.html"
                        )
