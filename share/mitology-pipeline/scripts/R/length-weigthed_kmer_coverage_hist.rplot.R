#!/usr/bin/env Rscript

## Copyright (c) 2014 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.

################################################
##
## @script: kmer_cov_weighted_hist.rplot.R
##
## @author: Joseph Tran (Joseph.Tran@versailles.inra.fr)
##
## @version: dev
##
## @date: 2014-11-09
##
## @description: This R script outputs a kmer coverage node-length-weighted histogram
##
###############################################

version <- "dev"

## execution time: start
T1<-Sys.time()

## functions
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
tryCatch.W.E <- function(expr)
{
	W <- NULL
	w.handler <- function(w)
	{ # warning handler
		W <<- w
		invokeRestart("muffleWarning")
	}
	list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
										warning = w.handler),
										warning = W)
}

## load and install libraries
### getopt
if (is.installed('getopt'))
{
	suppressPackageStartupMessages(require('getopt'));
} else
{
	install.packages('getopt');
    suppressPackageStartupMessages(library('getopt'));
}
### optparse
if (is.installed('optparse'))
{
	suppressPackageStartupMessages(require('optparse'));
} else
{
	install.packages('optparse');
    suppressPackageStartupMessages(library('optparse'));
}
### plotrix
if (is.installed('plotrix'))
{
	suppressPackageStartupMessages(require('plotrix'));
} else
{
	installed.packages('plotrix');
	suppressPackageStartupMessages(library('plotrix'));	
}

## options and usage
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default %default]"),
	make_option(c("-p","--prefix"), action="store", type="character", default="length-weighted_kmer_cov_hist.rplot", metavar="output filename", help="Save the outputs using the given prefix [default %default]"),
	make_option(c("-m","--xlim_min"), action="store", type="integer", default=0, metavar="xlim minimum value", help="The xlim range minimum value [default %default]"),
	make_option(c("-M","--xlim_max"), action="store", type="integer", default=500, metavar="xlim minimum value", help="The xlim range maximum value [default %default]")
)
parser <- OptionParser(usage = "usage:  %prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE, print_help_and_exit=FALSE)
opt <- arguments$options
if (opt$help) {
    cat(paste("author:  Joseph Tran <Joseph.Tran@versailles.inra.fr>\n\nversion: ", version, "\n\n"))
    print_help(parser)
    cat("notes:\n\t1. Copyright (c) 2014 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.\n\n")
    quit()
}

## test input file
if (length(arguments$args) != 1) {
    cat(paste("author:  Joseph Tran <Joseph.Tran@versailles.inra.fr>\n\nversion: ", version, "\n\n"))
    print_help(parser)
    cat("notes:\n\t1. Copyright (c) 2012 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.\n\n")
    stop("Incorrect number of required positional arguments\n\n")
} else {
    file <- arguments$args
}

if ( file.access(file) == -1) {
    stop(sprintf("Specified file ( %s ) does not exist", file))
} else {
    file_text <- readLines(file)
}

## verbosity
if ( opt$verbose ) {
    write("Arguments have been parsed successfully.\n", stderr())
}

## load stats file
# (R) > data = read.table("out-dir/meta-velvetg.Graph2-stats.txt", header=TRUE) 
# (R) > weighted.hist(data$short1_cov, data$lgth, breaks=seq(0, 200, by=1))

if ( opt$verbose ) {
    write("Load Graph2 stats data ... ", stderr())
}

warn_err <- tryCatch.W.E(read.table(file, header=TRUE))

if (is.null(warn_err$warning) && is.null(warn_err$value$message))
{
    data <- warn_err$value
    print(data[1:10,], file=stderr())
    cat("... ", file=stderr())
    write(paste(dim(data),"\n"), stderr())
} else {
    stop(paste(geterrmessage(), str(warn_err)))
}

## weighted histogram

# outdir: ?
outDir <- dirname(file)
pdf(file=file.path(outDir, "length-weighted_kmer_cov_hist.rplot.pdf"),
	   onefile=TRUE,
	   title=paste(strsplit(basename(file),".", fixed=TRUE)[[1]][1], "length-weigthed k-mer coverage histogram", sep=" "),
	   paper="special", height=11.7, width=16.5)

weighted.hist(data$short1_cov, data$lgth, breaks=seq(0, 500, by=1), xlab="K-mer coverage", ylab="Length-weighted Frequency", col=c("red"), main=paste(strsplit(basename(file),".", fixed=TRUE)[[1]][1], "length-weigthed k-mer coverage histogram", sep=" "))

dev.off()

## execution time: end
T2<-Sys.time()
Tdiff= difftime(T2, T1)
cat("Execution time : ", Tdiff,"seconds\n", file=stderr())













