#############################################
### QC GRAPH FUNCTIONS FROM DROP 
#############################################

# This uses the data.table package. Not only functions, also
# data table manipulation syntax.
# This function here uses counts(ods), which is merely the count matrix (gene-counts)
# and the merged bam stats, which is total counts per sample. Need to see how
# Hunter works first, but I can prepare whatever in the meantime.
# By summing counts(ods) by columns we get a similar result, but the numbers
# are much smaller (counts undergo filtering before becoming unfitted ods).
# The order is maintained. If viable, just the counts matrix should be enough.
# Maybe turn it into an outrider dataset.

plotCountRank <- function(CountsFile, ReadsFile)
{
    if(!missing(ReadsFile)){
        message("Reads file detected. Using total reads for sample ranking")
        ReadsFile = ReadsFile
    } else {
        message("Reads file not detected. Using total counts for sample ranking")
        ReadsFile = NA
    }
    coverage_dt <- get_counts(CountsFile = CountsFile, ReadsFile = ReadsFile)
    count_plot <- countRank(coverage_dt = coverage_dt, ReadsFile = ReadsFile)
    writePlot(plot = count_plot, name = "Counts Rank Plot")
}

get_counts <- function(CountsFile, ReadsFile = ReadsFile)
{
    if (missing(CountsFile)){
        stop('Please provide counts matrix')
    }
    # Load files. We do not need gene names from counts matrix
	cnts_mtx <- data.table::fread(CountsFile)

    # Rewriting in base R until I get help importing special data.table functions (such as :=)
    #bam_coverage$sampleID <- as.character(bam_coverage)

    if (!is.na(ReadsFile)){
        total_counts <- data.table::fread(ReadsFile)
        colnames(total_counts) <- c("sample", "attr", "val")
        total_counts <- total_counts[total_counts$attr=="initial_total_sequences",-2]
        # Total reads might have been counted without taking into account ExpHunterSuite blacklist,
        # which would lead to errors. This next line removes blacklisted samples from total reads table.
        total_counts <- total_counts[total_counts$sample %in% colnames(cnts_mtx[,-1])]
        # Remove redundant sampleID column
        total_counts <- as.numeric(total_counts$val)

    } else {
        total_counts = colSums(cnts_mtx[,-1])
    }

    coverage_dt <- data.table::data.table(sampleID = sort(colnames(cnts_mtx[,-1])),
                                total_counts = total_counts)

    coverage_dt <- coverage_dt[order(coverage_dt$total_counts)]

    coverage_dt$count_rank <- c(1:nrow(coverage_dt))

    return(coverage_dt)
}

countRank <- function(coverage_dt, ReadsFile = ReadsFile)
{
    if (!is.na(ReadsFile)) {
        xtitle="Sample Rank (total reads)"
    } else {
        xtitle="Sample Rank (total counts)"
    }
    plot <- ggplot2::ggplot(coverage_dt, aes(x = count_rank, y = total_counts)) +
        ggplot2::geom_point(size = 3) +
        cowplot::theme_cowplot() +
        cowplot::background_grid() +
        ggplot2::labs(title = "Read Counts", x=xtitle, y = "Reads Counted") +
        ggplot2::ylim(c(0,NA)) +
        ggplot2::scale_color_brewer(palette="Dark2") +
        ggplot2::ggtitle("Total counts sorted by rank")
}

writePlot <- function(plot, name, dir='/mnt/home/users/bio_267_uma/aestebanm/test/hunterTest/')
{
    png(filename=paste0(dir,name,".png"))
    plot(plot)
    dev.off()
}


  # #' # Filtering
  # #' **local**: A pre-filtered summary of counts using only the local (from BAM) counts. Omitted if no external counts  
  # #' **all**: A pre-filtered summary of counts using only the merged local (from BAM) and external counts  
  # #' **passed FPKM**: Passes the user defined FPKM cutoff in at least 5% of genes  
  # #' **min 1 read**: minimum of 1 read expressed in 5% of genes  
  # #' **min 10 reads**: minimum of 10 reads expressed in 5% of genes  

  # quant <- .95

  # if(has_external){
  #     filter_mtx <- list(
  #       local = cnts_mtx_local,
  #       all = cnts_mtx,
  #       `passed FPKM` = cnts_mtx[rowData(ods)$passedFilter,],
  #       `min 1 read` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 1, ],
  #       `min 10 reads` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 10, ]
  #     )
  #     filter_dt <- lapply(names(filter_mtx), function(filter_name) {
  #       mtx <- filter_mtx[[filter_name]]
  #       data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
  #     }) %>% rbindlist
  #     filter_dt[, filter := factor(filter, levels = c('local', 'all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
  # } else {
  #     filter_mtx <- list(
  #       all = cnts_mtx,
  #       `passed FPKM` = cnts_mtx[rowData(ods)$passedFilter,],
  #       `min 1 read` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 1, ],
  #       `min 10 reads` = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 10, ]
  #     )
  #     filter_dt <- lapply(names(filter_mtx), function(filter_name) {
  #       mtx <- filter_mtx[[filter_name]]
  #       data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
  #     }) %>% rbindlist
  #     filter_dt[, filter := factor(filter, levels = c('all', 'passed FPKM', 'min 1 read', 'min 10 reads'))]
  # }

  # binwidth <- .2
  # p_hist <- ggplot(filter_dt, aes(x = median_counts, fill = filter)) +
  #   geom_histogram(binwidth = binwidth) +
  #   scale_x_log10() +
  #   facet_wrap(.~filter) +
  #   labs(x = "Mean counts per gene", y = "Frequency", title = 'Mean Count Distribution') +
  #   guides(col = guide_legend(title = NULL)) +
  #   scale_fill_brewer(palette = "Paired") +
  #   theme_cowplot() +
  #   theme(legend.position = "none")

  # p_dens <- ggplot(filter_dt, aes(x = median_counts, col = filter)) +
  #   geom_density(aes(y=binwidth * ..count..), size = 1.2) +
  #   scale_x_log10() +
  #   labs(x = "Mean counts per gene", y = "Frequency") +
  #   guides(col = guide_legend(title = NULL)) +
  #   scale_color_brewer(palette = "Paired") +
  #   theme_cowplot() +
  #   theme(legend.position = "top",
  #         legend.justification="center",
  #         legend.background = element_rect(color = NA))

  # #+ meanCounts, fig.height=6, fig.width=12
  # plot_grid(p_hist, p_dens)

  # #' ### Expressed Genes
  # exp_genes_cols <- c(Rank = "expressedGenesRank",`Expressed\ngenes` = "expressedGenes", 
  #                     `Union of\nexpressed genes` = "unionExpressedGenes", 
  #                     `Intersection of\nexpressed genes` = "intersectionExpressedGenes", 
  #                     `Genes passed\nfiltering` = "passedFilterGenes", `Is External` = "isExternal")

  # expressed_genes <- as.data.table(colData(ods)[,exp_genes_cols])
  # colnames(expressed_genes) <- names(exp_genes_cols)

  # #+ expressedGenes, fig.height=6, fig.width=8
  # plotExpressedGenes(ods) + 
  #   theme_cowplot() +
  #   background_grid(major = "y") +
  #   geom_point(data =melt(expressed_genes,id.vars = c("Rank","Is External")),
  #              aes(x = Rank, y = value, col = variable, shape = `Is External`),show.legend = has_external)

  # if(has_external){
  #     DT::datatable(expressed_genes[order(Rank)],rownames = F)
  # } else{
  #     DT::datatable(expressed_genes[order(Rank),-"Is External"],rownames = F)
  # }

  # write.table(expressed_genes, paste0(opt$dataset, "_expressed_genes.tsv"), sep="\t", row.names= F)

