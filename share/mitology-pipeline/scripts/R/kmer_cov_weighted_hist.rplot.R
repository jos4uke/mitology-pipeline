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
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default %default]")
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
if(length(arguments$args) != 1) {
    cat(paste("author:  Joseph Tran <Joseph.Tran@versailles.inra.fr>\n\nversion: ", version, "\n\n"))
    print_help(parser)
    cat("notes:\n\t1. Copyright (c) 2012 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.\n\n")
    stop("Incorrect number of required positional arguments\n\n")
} else {
    file <- arguments$args
}

if( file.access(file) == -1) {
    stop(sprintf("Specified file ( %s ) does not exist", file))
} else {
    file_text <- readLines(file)
}

## verbosity
if ( opt$verbose ) {
    write("Arguments have been parsed successfully.\n", stderr())
}

