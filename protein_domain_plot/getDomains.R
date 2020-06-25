library(ensembldb)
library(EnsDb.Hsapiens.v86)
library("optparse")
 
option_list = list(
    make_option(c("-s", "--start"), type="numeric", default=NULL, 
              help="dataset file name", metavar="character"),
	make_option(c("-w", "--width"), type="numeric", default=NULL, 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-t", "--transcript"), default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "14")

rng_tx <- IRanges(start = c(opt$start), width = c(opt$width),
                  names = c(opt$transcript))
rng_prt <- transcriptToProtein(rng_tx, edbx)

rng_prt
