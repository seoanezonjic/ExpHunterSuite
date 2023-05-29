#!/usr/bin/env Rscript
library(dplyr)

filter_and_rank <- function(cv_cont_table){
    # cv_cont_table <- cv_cont_table[cv_cont_table$Pvalue <= 0.05,]
    if (nrow(cv_cont_table) == 0 ) return(NULL)

    cv_cont_table <- cv_cont_table[,c("miRNA","strategy", "corr_cutoff", "crossval_status", "Odds_ratio", "TP")]
    cv_cont_table <- reshape(cv_cont_table, idvar = c("miRNA","strategy", "corr_cutoff"), timevar = "crossval_status", direction = "wide")
    if (sum(c("Odds_ratio.test", "Odds_ratio.train") %in% colnames(cv_cont_table)) != 2) return(NULL)
    # cv_cont_table <- cv_cont_table  %>% arrange(desc(Odds_ratio.test), desc(TP.test), corr_cutoff)
    # cv_cont_table$OR_test_rank <- seq(1, nrow(cv_cont_table))    
    # str(cv_cont_table)
    cv_cont_table$OR_test_rank <- rank(-cv_cont_table$Odds_ratio.test, ties.method ="min")
    cv_cont_table$OR_train_rank <- rank(-cv_cont_table$Odds_ratio.train,ties.method ="min")
    return(cv_cont_table)

}
get_test_rank <- function(cv_cont_table){
    top_rank_train_val <- min(cv_cont_table$OR_train_rank)
    top_rank_train <- cv_cont_table[cv_cont_table$OR_train_rank ==top_rank_train_val,]
    lower_rank_test <- min(top_rank_train$OR_test_rank)
   # print(top_rank_train)
    result <- data.frame(test_rank = top_rank_train[top_rank_train$OR_test_rank == lower_rank_test, c("miRNA","strategy","corr_cutoff","OR_test_rank")])
    return(result)
}


################### OPTIONS
option_list <- list(
    optparse::make_option(c("-i", "--input_cormit"), type="character", 
        default="",
        help="coRmit folder. Note that this script needs an execution of comit coRmiT with '--save_temp' option activated."),
    optparse::make_option(c("-s", "--test_size"), type="double", 
        help="Proportion of input taken as test. Default=%default",
        default=0.25),
   optparse::make_option(c("-o", "--output_files"), type="character",
     default=".", 
        help = "Output folder")

        )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

################### INITIALIZE
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



cormit_path <- file.path(opt$input_cormit, "temp.RData")
load(cormit_path)

p_val_cutoff <- 0.05
all_pairs <- miRNA_cor_results$all_pairs
miRNA_strat_result <- miRNA_cor_results$miRNA_cont_tables
miRNA_strat_result <-  miRNA_strat_result[miRNA_strat_result$Pvalue <= 0.05 & miRNA_strat_result$db_group == "multimir",]
best_miRNA_strategies <- as.data.frame(select_best_strategy(miRNA_strat_result))
strat_combinations <- as.data.frame(unique(miRNA_strat_result[, c("corr_cutoff", "strategy")]))
miRNA_cont_tables <- data.frame()



for (combination in seq(1,nrow(strat_combinations))) {
    corr_cutoff <- strat_combinations[combination, "corr_cutoff"]
    strategy <- strat_combinations[combination, "strategy"]

    contingency_tables <- prepare_for_stats(all_pairs, strategy, corr_cutoff,
    db_groups = "multimir", p_val_cutoff = p_val_cutoff , crossval = TRUE, test_sample = opt$test_size)
    miRNA_cont_tables <- rbind(miRNA_cont_tables, contingency_tables[["strat_miRNA_ct"]])
}


message("Stats prepared")
miRNA_cont_tables <- v_get_stats(miRNA_cont_tables, 
                              selected_stats = c("v.fisher.test","odds_ratio"))

miRNA_cont_tables$iterations <- paste(miRNA_cont_tables$miRNA, miRNA_cont_tables$crossval_it, sep = "_")

splitted_cont_tables <- split(miRNA_cont_tables, miRNA_cont_tables$iterations)


cross_rank <- lapply(splitted_cont_tables, filter_and_rank)
full_ranking <-as.data.frame(data.table::rbindlist(cross_rank, 
                                     use.names = TRUE, 
                                     idcol = "iters"))
cross_rank <- lapply(cross_rank, get_test_rank)
message("Train and test datasets ordered")
cross_rank <-as.data.frame(data.table::rbindlist(cross_rank, 
                                     use.names = TRUE, 
                                     idcol = "iters"))
colnames(cross_rank) <- c("iteration", "miRNA", "strategy","corr_cutoff","test_rank")
best_miRNA_strategies$acc_ratio <- 0
best_miRNA_strategies$total_strats <- 0
# save(cross_rank,best_miRNA_strategies,full_ranking, file = "test.Rdata")
# load("test.Rdata")
train_rank_dist <- data.frame()
for (miRNA in best_miRNA_strategies$miRNA) {
    best_miRNA_strat <- unlist(best_miRNA_strategies[best_miRNA_strategies$miRNA == miRNA,])

    best_rank_in_train <- full_ranking[full_ranking$miRNA == miRNA &
                                        full_ranking$strategy == best_miRNA_strat[7] &
                                        full_ranking$corr_cutoff == best_miRNA_strat[8],]
    train_rank_dist <- rbind(train_rank_dist, best_rank_in_train)
}

str(full_ranking)
message("Computing report")
dir.create(opt$output_files)
outf <- file.path(normalizePath(opt$output_files),"coRmiT_cv.html")
    rmarkdown::render(file.path(template_folder, 'cormit_cv.Rmd'),
                      output_file = outf, intermediates_dir = opt$output_files)