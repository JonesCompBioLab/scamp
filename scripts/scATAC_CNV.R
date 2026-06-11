#!/usr/bin/env Rscript

suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(R.utils))
set.seed(1)

"%ni%" <- Negate("%in%")

args <- commandArgs(trailingOnly=TRUE)
FRAGS_DIR <- args[[1]]
WINDOW_SIZE <- args[[2]]
STEP_SIZE <- args[[3]]
NEIGHBORS <- args[[4]]
WHITELIST <- args[[5]]
OUTPUT_DIR <- args[[6]]
REFERENCE <- args[[7]]
ROOT_DIR <- args[[8]]

# load in utilities
source(paste0(
    ROOT_DIR, "/scripts/scATAC_CNV_utilities.R"
))

frag_files <- readFragFiles(FRAGS_DIR)

if (grepl("atac", names(frag_files)[[1]], fixed=T)) {
    names(frag_files) <- names(frag_files) %>% gsub('_atac.*','',.) %>% gsub('.*/','',.)
} else {
    names(frag_files) <- names(frag_files) %>% gsub('.fragments.*','',.) %>% gsub('.*/','',.)
}

whitelist <- read.table(WHITELIST, sep='\t', header=F, comment.char='')[,1]

if (REFERENCE == 'hg19') {
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    blacklist <- rtracklayer::import.bed(
        paste0(ROOT_DIR, '/reference/hg19-blacklist.v2.bed.gz')
    )
    genome <- BSgenome.Hsapiens.UCSC.hg19
} else if (REFERENCE == 'hg38') {
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    blacklist <- rtracklayer::import.bed(
        paste0(ROOT_DIR, '/reference/hg38.blacklist.bed.gz')
    )
    genome <- BSgenome.Hsapiens.UCSC.hg38
} else if (REFERENCE == 'mm10') {
    suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
    blacklist <- rtracklayer::import.bed(
        paste0(ROOT_DIR, '/reference/mm10-blacklist.v2.bed.gz')
    )
    genome <- BSgenome.Mmusculus.UCSC.mm10
} else {
    stop('Choose one of the following referenge genomes: `hg38`, `hg19`, `mm10`')
}

windows <- makeWindows(genome = genome,
                           blacklist = blacklist %>% subset(width(.) > 1000),
                           windowSize = as.numeric(WINDOW_SIZE),
                           slidingSize = as.numeric(STEP_SIZE))

fragFiles2cnaObj(frag_files = frag_files, 
                     windows = windows,
                     whitelist = whitelist,
                     dir = OUTPUT_DIR,
                     neighbors = as.numeric(NEIGHBORS))

