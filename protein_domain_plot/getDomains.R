library(ensembldb)
library(EnsDb.Hsapiens.v86)
library("optparse")
 
option_list = list(
    make_option(c("-s", "--start"), type="numeric", default=NULL, 
              help="dataset file name", metavar="character"),
	make_option(c("-e", "--end"), type="numeric", default=NULL, 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-a", "--accession"), default=NULL, 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-c", "--chromosome"), default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == opt$chromosome)



uni_rng <- IRanges(start = c(opt$start), end = c(opt$end),
                   names = c(opt$accession))

## We have to specify that the IDs are Uniprot IDs
uni_gnm <- proteinToGenome(uni_rng, edbx, idType = "uniprot_id")

uni_gnm

