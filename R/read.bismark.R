read.bismark <- function(files, sampleNames, rmZeroCov = FALSE, verbose = TRUE){
    read.bedgraph(files, sampleNames, zeroBased = FALSE, rmZeroCov = FALSE, verbose = TRUE)
}
