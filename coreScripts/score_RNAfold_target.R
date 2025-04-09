##################################################
### compute viral target folding structures

library(optparse)
library(Biostrings, quietly=T)
library(here)

option_list <- list(make_option(c("-i", "--input"), type="character", 
                                default=file.path(here(), "ref_data/NC_045512v2.fa"), 
                                help="genome .fa file name", metavar="character"),
                    make_option(c("-w", "--window"), type="integer", default=20,
                                help="window size", metavar="integer"),
                    make_option(c("-f", "--folding"), type="integer", default=100,
                                help="folding window size", metavar="integer"),
                    make_option(c("-o", "--out"), type="character", default=".", 
                                help="output directory", metavar="character")) 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

RNAfold_path <- "/opt/anaconda3/bin/RNAfold"

cat("\nRNAfold score: viral target structure\n")

# read in genome sequence, process into continuous string
genome_seq <- readLines(opt$input)
genome_seq <- genome_seq[!grepl(">", genome_seq)]
genome_seq <- paste0(genome_seq, collapse="")

# read in folded windows
if(!file.exists(paste0("wuhCor1_windows_", opt$folding, "_RNAfold.txt"))) {
  # break genome into windows
  windows <- sapply(seq.int(nchar(genome_seq)-opt$folding-3), # all windows that have 4 nt following
                    function(x) {
                      substr(genome_seq, start=x, stop=x+opt$folding-1)
                    })
  writeLines(windows, con=file.path(opt$out, paste0("wuhCor1_windows_", opt$folding, ".txt")))
  
  # fold windows
  cat("- computing folding structures\n")
  system(paste0(RNAfold_path, " -i ", file.path(opt$out, paste0("wuhCor1_windows_", opt$folding, ".txt")), 
                " --noPS --outfile=wuhCor1_windows_", opt$folding, "_RNAfold.txt"))
}
RNAfold <- data.frame(matrix(readLines(paste0("wuhCor1_windows_", opt$folding, "_RNAfold.txt")), 
                             ncol=2, byrow=T), stringsAsFactors=F)
colnames(RNAfold) <- c("sequence", "structure")
RNAfold$MFE <- as.numeric(sub(" ", "", sub("\\)", "", sub(".*\\(", "", RNAfold$structure))))
RNAfold$structure <- sub(" .*", "", RNAfold$structure)

# compute propensity to base-pair: nucleotide positions
cat("- calculating base-pairing propensity\n")
propensity_nt <- sapply(seq.int(nchar(genome_seq)),
                        function(x) {
                          which_windows <- seq.int(from=x, to=ifelse(x <= opt$folding, 1, x-opt$folding+1))
                          tmp_windows <- matrix(unlist(strsplit(RNAfold$structure[which_windows], split="")), 
                                                nrow=length(which_windows), byrow=T)
                          tmp_structure <- diag(tmp_windows)
                          return(mean(tmp_structure %in% c("(", ")")))
                        })

# compute propensity to base-pair: windows
if(file.exists(file.path(opt$out, "windows.txt"))) {
  cat("- pulling window positions from windows.txt\n")
  window_starts <- read.table(file.path(opt$out, "windows.txt"), header=T)$start
} else {
  window_starts <- seq.int(from=1, to=nchar(genome_seq)-opt$window+1)
}
propensity_window <- sapply(window_starts,
                            function(x) {
                              mean(propensity_nt[x:(x+opt$window-1)])
                            })

# output target RNAfold scores
write.table(data.frame(start=window_starts,
                       strand=read.table(file.path(opt$out, "windows.txt"), header=T)$strand,
                       target_basepairing_propensity=propensity_window),
            file=file.path(opt$out, "score_RNAfold_target.txt"),
            quote=F, sep="\t", row.names=F)
