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
### ggplot2
if (is.installed('ggplot2'))
{
    suppressPackageStartupMessages(require('ggplot2'));
} else
{
    installed.packages('ggplot2');
    suppressPackageStartupMessages(library('ggplot2'));
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
    cat("notes:\n\t1. Copyright (c) 2014 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.\n\n")
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
#pdf(file=file.path(outDir, "length-weighted_kmer_cov_hist.rplot.pdf"),
#	   onefile=TRUE,
#	   title=paste(strsplit(basename(file),".", fixed=TRUE)[[1]][1], "length-weigthed k-mer coverage histogram", sep=" "),
#	   paper="special", height=11.7, width=16.5)

# plotrix
#weighted.hist(data$short1_cov, data$lgth, breaks=seq(0, 500, by=1), xlab="K-mer coverage", ylab="Length-weighted Frequency", col=c("red"), main=paste(strsplit(basename(file),".", fixed=TRUE)[[1]][1], "length-weigthed k-mer coverage histogram", sep=" "))
# ggplot2
p <- ggplot(data=data, aes(short1_cov, weight = lgth)) 
p <- p + geom_histogram(binwidth=1, aes(fill = ..count..)) +
	scale_fill_gradient("Count", low = "green", high = "red") +
	geom_density() +
	xlim(opt$xlim_min,opt$xlim_max) ### TODO: make it an option xlim_min xlim_max

## get density peaks and (x, y) coordinates
r <- print(p)
peakx=r$data[[2]]$x[which(diff(sign(diff(r$data[[2]]$scaled)))==-2)]
peaks <- data.frame( x=peakx, y=r$data[[2]]$y[which(r$data[[2]]$x %in% peakx)], lgth=rep(0, length(peakx)))
max_peaky <- max(peaks$y)
## plot peaks as dot
#p + geom_point(data=peaks, aes(x=x,y=y, colour="red", size=3)) + guides(colour=FALSE, size=FALSE)
## plot x and y peaks lines and coordinates 
for ( i in 1:length(peaks$x)) {
	#p <- p + geom_point(aes(x=peaks$x[i], y=peaks$y[i], colour="red", size=3)) + guides(colour=FALSE, size=FALSE)
	print(paste("Peak", i, "X:", peaks$x[i], ", Y:", peaks$y[i], sep=" "), file=stderr())
	coordx=c(0,peaks$x[i],peaks$x[i])
	coordy=c(peaks$y[i],peaks$y[i],0)
	p <- p + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=paste("#", i,"\n","x=", round(coordx[2],2),"\n", "y=", round(coordy[2]), sep=''), x=coordx[2], y=ifelse((coordy[2]/max_peaky)<=0.5, max_peaky/2, max_peaky*1.3), size=5, colour="black") + guides(colour=FALSE, size=FALSE)
	p <- p + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=round(coordy[1]), x=ifelse(opt$xlim_max <=100, -2, -10), y=coordy[1]*1.05, size=4, colour="red") + guides(colour=FALSE, size=FALSE)
	p <- p + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=round(coordx[2],2), x=coordx[2], y=ifelse((coordy[2]/max_peaky)<=.5, max_peaky*(-0.1), max_peaky*(-0.05)), size=4, colour="red") + guides(colour=FALSE, size=FALSE)
#})	
	#p <- p + geom_segment(aes(x=coordx[1], y=coordy[1], xend=coordx[2], yend=coordy[2], colour="red")) + guides(colour=FALSE)
	#p <- p + geom_segment(aes(x=coordx[2], y=coordy[2], xend=coordx[3], yend=coordy[3], colour="red")) + guides(colour=FALSE)
}

## render plot	 
#pdf(file=file.path(outDir, paste(opt$prefix, "_xlim_", opt$xlim_min,"_", opt$xlim_max, ".pdf", sep='')),
#       onefile=TRUE,
#       title=paste(opt$prefix, "_xlim_", opt$xlim_min,"_", opt$xlim_max, ".pdf", sep=''),
#	   #title=paste(strsplit(basename(file),".", fixed=TRUE)[[1]][1], "length-weigthed k-mer coverage histogram :", paste(opt$prefix, "_xlim_", opt$xlim_min,"_", opt$xlim_max, ".pdf", sep=''), sep=" "),
#       paper="special", height=11.7, width=16.5)

p + xlab("K-mer coverage") + 
	ylab("Length-weighted Frequency") +
	ggtitle(paste(strsplit(basename(file),".", fixed=TRUE)[[1]][1], "length-weigthed k-mer coverage histogram\n", paste(opt$prefix, "_xlim_", opt$xlim_min,"_", opt$xlim_max, ".pdf", sep=''), sep=" ")) 

ggsave(p, file=paste(opt$prefix, "_xlim_", opt$xlim_min,"_", opt$xlim_max, ".pdf", sep=''),
	   path=outDir,
	   height=11.7, width=16.5, unit="in")

#dev.off()

# save peaks to file
write.table(peaks, file=file.path(outDir, paste(opt$prefix, "_xlim_", opt$xlim_min,"_", opt$xlim_max, "_density_peaks.tab",sep='')), 
				col.names=T, row.names=T, quote=FALSE, sep="\t"
)

## execution time: end
T2<-Sys.time()
Tdiff= difftime(T2, T1)
cat("Execution time : ", Tdiff,"seconds\n", file=stderr())













