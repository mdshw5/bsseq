read.gemstone <- function(files, sampleNames, rmZeroCov = FALSE, verbose = TRUE){
    ## Argument checking
    if (anyDuplicated(files)){
        stop("duplicate entries in 'files'")
    }
    if (length(sampleNames) != length(files) | anyDuplicated(sampleNames)){
        stop("argument 'sampleNames' has to have the same length as argument 'files', without duplicate entries")
    }
    ## Process each file
    idxes <- seq_along(files)
    names(idxes) <- sampleNames
    allOut <- lapply(idxes, function(ii){
        if (verbose) {
            cat(sprintf("[read.gemstone] Reading file '%s' ... ", files[ii]))
        }
        ptime1 <- proc.time()
        raw <- read.gemstoneFileRaw(thisfile = files[ii])
        M <- matrix(elementMetadata(raw)[, "M2"], ncol = 1)
        Cov <- matrix(elementMetadata(raw)[, "Cov"], ncol = 1)
        elementMetadata(raw) <- NULL
        out <- BSseq(gr = raw, M = M, Cov = Cov,
                     sampleNames = sampleNames[ii], rmZeroCov = rmZeroCov)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if (verbose) {
            cat(sprintf("done in %.1f secs\n", stime))
        }
        out
    })
    if (verbose) {
        cat(sprintf("[read.gemstone] Joining samples ... "))
    }
    ptime1 <- proc.time()
    allOut <- combineList(allOut)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }
    allOut
}

read.gemstoneFileRaw <- function(thisfile, verbose = TRUE){
    ## Set up the 'what' argument for scan()
    columnHeaders <- c("rname", "mt", "pos", "M2", "M1", "M", "C", "N")
    what0 <- replicate(length(columnHeaders), character(0))
    names(what0) <- columnHeaders
    int <- which(columnHeaders %in% c("pos", "mt", "M2", "M1", "M", "C", "N"))
    what0[int] <- replicate(length(int), integer(0))
    ## Read in the file
    if (grepl("\\.gz$", thisfile))
        con <- gzfile(thisfile)
    else
        con <- file(thisfile, open = "r")
    out <- scan(file = con, what = what0, sep = "\t", skip = 1, quote = "", na.strings = "NA", quiet = TRUE)
    close(con)
    ## Remove all non CpG rows
    out <- lapply(out, "[", which(out[["mt"]] %in% c(1, 4)))
    ## Create GRanges instance from 'out'
    gr <- GRanges(seqnames = out[["rname"]], ranges = IRanges(start = out[["pos"]], width = 1))
    out[["Cov"]] <- out[["M2"]] + out[["M1"]] + out[["M"]] + out[["C"]] + out[["N"]]
    out[["rname"]] <- out[["pos"]] <- out[["mt"]] <- out[["M1"]] <- out[["M"]] <- out[["C"]] <- out[["N"]] <- NULL
    ## Should leave "M2" and "Cov"
    out <- out[!sapply(out, is.null)]
    df <- DataFrame(out)
    elementMetadata(gr) <- df
    gr
}
