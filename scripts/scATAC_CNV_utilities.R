# A set of utilities for processing scATACseq adata and computing copy-number
# profiles.

"%ni%" <- Negate("%in%")


#' A function for reading fragment files from a directory.
#' Args:
#       data_dir: Contains .tsv.gz and .tsv.gz.tbi files in subdirectories named
#'          by sample.
readFragFiles <- function(data_dir) {
    data_dir <- normalizePath(data_dir)
    frag_files <- list.files(data_dir, pattern = "*fragments.tsv.gz$", recursive = TRUE, full.names = TRUE) 
    names(frag_files) <- frag_files %>% gsub(paste0(data_dir,'/'), "", .) %>% gsub('/atac.*','',.)
    return(frag_files)
}

getArrowFiles <- function(archr_proj_dir) { #archr_proj_dir should be an ArchR project.

    archr_proj_dir <- normalizePath(archr_proj_dir)
    arrowFiles <- list.files(paste0(archr_proj_dir, "/ArrowFiles"), recursive = TRUE, full.names = TRUE)
    names(arrowFiles) <- arrowFiles %>% gsub(paste0(archr_proj_dir, "/ArrowFiles/"), "", .) %>% gsub(".arrow", "", .)
    return(arrowFiles)

}

#' A function for reading fragment files.
#' Args:
#'      file: A ATAC fragment file from a 10X experiment.
#'      convertArrowFile: Convert ArchR ArrowFile into a fragment file
#'      minFrags: Minimum number of fragments for a given cell.
readFragments <- function(file, convertArrowFile = FALSE, minFrags = 100) {
    message("Reading in fragment files...")
        if (convertArrowFile) {
                fragments <- getFragmentsFromArrow(file)
                mcols(fragments)$N <- 1
        } else {
                fragments <- data.frame(fread(file, header=FALSE))
                if (ncol(fragments) == 4) {
                        fragments[,5] <- 1
                }
                fragments <- GRanges(
                        seqnames = fragments[,1],
                        IRanges(fragments[,2]+1, fragments[,3]),
                        RG = fragments[,4],
                        N = fragments[,5]
                )
        }
    message("Filtering Lowly Represented Cells...")
    tabRG <- table(fragments$RG)
    keep <- names(tabRG)[which(tabRG >= minFrags)]
    fragments <- fragments[fragments$RG %in% keep,]
    fragments <- sort(sortSeqlevels(fragments))
}

#' A helper to chunkify data.
#' Args:
#'      x: Data to chunkify.
#'      n: Number of chunks.
divideIntoChunks <- function(x,n) {
    if (n == 1) { return( list(x) ) } 
    else { split(x, cut(seq_along(x), n, labels = FALSE)) }
}


addCNs <- function(cnaObj, bgd_CN = 2){
    assays(cnaObj)$CNs <- bgd_CN * ( 2 ^ assays(cnaObj)$log2FC )
    return(cnaObj)
}

fragFiles2cnaObj <- function(frag_files, 
                             windows,
                             whitelist, 
                             dir,
                             neighbors = 100,
                             force = TRUE,
                             remove = c("chrM","chrX","chrY"),
                             convertArrowFile = FALSE,
                             addRG = FALSE,
                             max_granges_length = 4e+7){
    
    for (sample_name in names(frag_files)) {

        message(paste0("CNA object being generated for ", sample_name))
        fragments <- readFragments(frag_files[sample_name], convertArrowFile=convertArrowFile)
        if (addRG) {
            fragments$RG <- paste0(sample_name, "#", fragments$RG)
        }
	    if (!is.null(whitelist)) {
            fragments <- fragments[mcols(fragments)[,"RG"] %in% whitelist]
	    }
                
        split_into <- ceiling( length(fragments) / max_granges_length )
        chunks <- divideIntoChunks(unique(fragments$RG), split_into)
        cnaObj <- lapply(chunks, function(chunk){
            fragments_chunk <- fragments %>% subset(RG %in% chunk)
            cnaObj <- scCNA(windows, 
                            fragments_chunk, 
                            neighbors = neighbors, 
                            LFC = 1.5, 
                            FDR = 0.1, 
                            force = force, 
                            remove = remove)
            return(cnaObj)
        }) %>% do.call(cbind, .)
        cnaObj <- addCNs(cnaObj)
        saveRDS(cnaObj, file = paste0(dir, '/', 'cnaObj_', sample_name, '.rds'))

    }    
    
    cnaObj_files <- list.files(dir, pattern = c("cnaObj_.*.rds$"), full.names = TRUE) %>%
        subset( grepl(paste(names(frag_files), collapse = '|'), .) )
    return(cnaObj_files)
}


labelCNAWindows <- function(cnaObj, 
                            min_fraction_cells_CNA){ #minimum fraction of cells with CNA in any given window
    min_n_cells_CNA <- min_fraction_cells_CNA * ncol(cnaObj)
    n_cells_CNA <- assays(cnaObj)$CNA %>% apply(1, function(row){ length(which(row == 1)) })
    rowData(cnaObj)$CNA <- n_cells_CNA >= min_n_cells_CNA    
    return(cnaObj)
}

countInsertions <- function(query, fragments, by = "RG"){
  #Count By Fragments Insertions
  seqlevelsStyle(fragments) <- 'UCSC' #convert seqnames to UCSC, change to seqlevelsStyle(query)
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), 
            RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), 
            RG = mcols(fragments)[,by])
  )
#  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

makeWindows <- function(genome, blacklist, windowSize = 10e6, slidingSize = 2e6){
    chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
    chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
    windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% 
        .[which(width(.)==windowSize),]
    mcols(windows)$wSeq <- as.character(seqnames(windows))
    mcols(windows)$wStart <- start(windows)
    mcols(windows)$wEnd <- end(windows)
    message("Subtracting Blacklist...")
    windowsBL <- lapply(seq_along(windows), function(x){
            if(x %% 100 == 0){
                message(sprintf("%s of %s", x, length(windows)))
            }
            gr <- GenomicRanges::setdiff(windows[x,], blacklist)
            mcols(gr) <- mcols(windows[x,])
            return(gr)
        })
    names(windowsBL) <- paste0("w",seq_along(windowsBL))
    windowsBL <- unlist(GRangesList(windowsBL), use.names = TRUE)
    mcols(windowsBL)$name <- names(windowsBL)
    message("Adding Nucleotide Information...")
    windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
    windowNuc <- lapply(seq_along(windowSplit), function(x){
        message(sprintf("%s of %s", x, length(windowSplit)))
        idx <- which( as.character(seqnames(chromSizes)) == names(windowSplit)[x] )
        chrSeq <- Biostrings::getSeq(genome,chromSizes[idx])
        grx <- windowSplit[[x]]
        aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
        mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
        mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
        return(grx)
      }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
    windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
    windowNuc
}

scCNA <- function(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, 
                  remove = c("chrM","chrX","chrY")){
	
	#Keep only regions in filtered chromosomes
	windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
	fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
	windows <- windows[seqnames(windows) %ni% remove]
	fragments <- fragments[seqnames(fragments) %ni% remove]

	#Count Insertions in windows
	message("Getting Counts...")
	counts <- countInsertions(windows, fragments, by = "RG")[[1]]
	message("Summarizing...")
	windowSummary <- GRangesList()
	countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
	for(x in seq_along(unique(mcols(windows)$name))){
		if(x %% 100 == 0){
			message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
		}
		idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
		wx <- windows[idx,]
		wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
		mcols(wo)$name <- mcols(wx)$name[1]
		mcols(wo)$effectiveLength <- sum(width(wx))
		mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
		mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
		mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
		mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
		countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
		windowSummary[[x]] <- wo
	}
	windowSummary <- unlist(windowSummary)
	
	#Keep only regions with less than 0.1% N
	keep <- which(windowSummary$N < 0.001) 
	windowSummary <- windowSummary[keep,]
	countSummary <- countSummary[keep,]
	
	#Now determine the nearest neighbors by GC content
	message("Computing Background...")
	bdgMean <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	bdgSd <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	log2FC <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	z <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
	pval <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))

	for(x in seq_len(nrow(countSummary))){
		if(x %% 100 == 0){
			message(sprintf("%s of %s", x, nrow(countSummary)))
		}
		#Get Nearest Indices
		idxNN <- head(order(abs(windowSummary$GC[x] - windowSummary$GC)), neighbors + 1)
		idxNN <- idxNN[idxNN %ni% x]
		#Background
		if(any(colMeans(countSummary[idxNN, ])==0)){
			if(force){
				message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
			}else{
				stop("Background Mean = 0!")
			}
		}
		bdgMean[x, ] <- colMeans(countSummary[idxNN, ])
		bdgSd[x, ] <- matrixStats::colSds(countSummary[idxNN, ])
		log2FC[x, ] <- log2((countSummary[x, ]+1e-5) / (bdgMean[x, ]+1e-5))
		z[x, ] <- (countSummary[x,] - bdgMean[x, ]) / bdgSd[x, ]
		pval[x, ] <- 2*pnorm(-abs(z[x, ]))
	}
	padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
	CNA <- matrix(0, nrow=nrow(countSummary), ncol=ncol(countSummary))
	CNA[which(log2FC >= LFC & padj <= FDR)] <- 1

	se <- SummarizedExperiment(
		assays = SimpleList(
				CNA = CNA,
				counts = countSummary,
				log2FC = log2FC,
				padj = padj,
				pval = pval,
				z = z,
				bdgMean = bdgMean,
				bdgSd = bdgSd
			),
		rowRanges = windowSummary
	)
	colnames(se) <- colnames(counts)

	return(se)
}
