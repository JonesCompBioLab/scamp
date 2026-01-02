# suppressMessages(library(ArchR, quietly=T))

# ------------ UTILITIES ------------- #
average_windows_containing_region <- function(cna_data, query_gr, copy_number_array) {
    gr <- rowRanges(cna_data)
    windows <- subsetByOverlaps(gr, query_gr, maxgap=1e5, minoverlap=0L, type='any')

    if (length(windows) == 0) {
        return(NULL)
    }

    cns_subset <- copy_number_array[,windows$name,drop=F]
    return(apply(cns_subset, 1, mean))
}


# ------------- ARGUMENTS -------------- #
args <- commandArgs(trailingOnly=TRUE)
CNA_DIR <- args[[1]]
OUTPUT_DIR <- args[[2]]
REFERENCE <- args[[3]]
ROOT_DIR <- args[[4]]

# ------------- PIPELINE -------------- #

if (REFERENCE == 'hg19') {
    gene_ranges <- readRDS(paste0(
        ROOT_DIR, "/reference/geneAnnoHg19.rds")
    )$genes
} else if (REFERENCE == 'hg38') {
   gene_ranges <- readRDS(paste0(
        ROOT_DIR, "/reference/geneAnnoHg38.rds")
    )$genes
} else if (REFERENCE == 'mm10') {
    gene_ranges <- readRDS(paste0(
        ROOT_DIR, "/reference/geneAnnoMm10.rds")
    )$genes
} else {
    stop('Choose one of the following referenge genomes: `hg38`, `hg19`, `mm10`')
}

gene_ranges <- gene_ranges[!is.na(gene_ranges$symbol),]

files <- list.files(CNA_DIR)
for (file in files) {
    message(paste0("Processing ", file, '...'))

    message("Loading in data...")
    cna_windows <- readRDS(paste0(CNA_DIR, "/", file))

    # creating copy-number array
    cns = assays(cna_windows)$CNs
    rownames(cns) <- rowData(cna_windows)$name
    cns = as.data.frame(t(cns))

    message("Computing copy numbers for gene bodies...")
    cn_list = list()
    x <- 0
    for (gene in gene_ranges$symbol) {
        if(x %% 1000 == 0){
                message(sprintf("%s of %s", x, length(gene_ranges$symbol)))
        }
        x <- x + 1

        copy_number <- average_windows_containing_region(cna_windows, subset(gene_ranges, symbol == gene), cns)
        if (is.null(copy_number)) {
            next
        }
        
        cn_list[[gene]] <- copy_number
    }

    cn_dataframe <- data.frame(cn_list)

    file.stub <- strsplit(file, '.rds', fixed=T)[[1]]
    agg.name <- paste0(OUTPUT_DIR, "/", file.stub, ".all.tsv")

    write.table(cn_dataframe, agg.name, sep='\t')
    message(paste0("Output written to ", agg.name))

}

message("Done.")




