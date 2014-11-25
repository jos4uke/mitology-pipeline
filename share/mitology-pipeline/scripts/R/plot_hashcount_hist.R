#!/usr/bin/env Rscript

## Copyright (c) 2014 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.

################################################
##
## @script: plot_hashcount_hist.R
##
## @author: Joseph Tran (Joseph.Tran@versailles.inra.fr)
##
## @version: dev
##
## @date: 2014-11-20
##
## @description: This R script outputs a kmer abundance histogram
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
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print extra output [default %default]"),
	make_option(c("-x","--xlim_min"), action="store", type="integer", default=0, metavar="xlim_minimum_value", help="The xlim range minimum value [default %default]"),
	make_option(c("-X","--xlim_max"), action="store", type="integer", default=NULL, metavar="xlim_minimum_value", help="The xlim range maximum value [default %default]"),
	make_option(c("-y","--ylim_min"), action="store", type="integer", default=0, metavar="ylim_minimum_value", help="The ylim range minimum value [default %default]"),
	make_option(c("-Y","--ylim_max"), action="store", type="integer", default=NULL, metavar="ylim_minimum_value", help="The ylim range maximum value [default %default]"),
	make_option(c("-c","--cutoffestim"), action="store_true", default=FALSE, help="Estimate k-mer coverage cutoff [default %default]"),
	make_option(c("-d","--desc_pctl_min"), action="store", type="integer", default=50, metavar="descending_percentile_minimim_value", help="The descending percentile minimum value [default %default]"),
	 make_option(c("-D","--desc_pctl_max"), action="store", type="integer", default=75, metavar="descending_percentile_maximum_value", help="The descending percentile maximum value [default %default]"),
	 make_option(c("-a","--asc_pctl_min"), action="store", type="integer", default=25, metavar="ascending_percentile_minimum_value", help="The ascending percentile minimum value [default %default]"),
	 make_option(c("-A","--asc_pctl_max"), action="store", type="integer", default=50, metavar="ascending_percentile_maximum_value", help="The ascending percentile maximum value [default %default]")
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

if ( opt$verbose ) {
    write("Load hashcount hist data ... ", stderr())
}

warn_err <- tryCatch.W.E(read.table(file, header=FALSE))

if (is.null(warn_err$warning) && is.null(warn_err$value$message))
{
    data <- warn_err$value
    print(data[1:10,], file=stderr())
    cat("... ", file=stderr())
    write(paste(dim(data),"\n"), stderr())
} else {
    stop(paste(geterrmessage(), str(warn_err)))
}

# header
colnames(data) <- c("kmer_abundance", "kmer_count", "cumulative_count", "fraction_of_total_distinct_kmers") 
write("Update hashcount hist data colnames ... ", stderr())
print(data[1:10,], file=stderr())

# x/y lims max
opt$xlim_max <- ifelse(is.null(opt$xlim_max), max(data$kmer_abundance), opt$xlim_max)
opt$ylim_max <- ifelse(is.null(opt$ylim_max), max(data$kmer_count), opt$ylim_max)
##
print(paste("xlim max: ", opt$xlim_max, sep=''))
print(paste("ylim max: ", opt$ylim_max, sep=''))

# outdir: ?
outDir <- dirname(file)
p0 <- ggplot(data=data, aes(x=kmer_abundance, y=kmer_count))

## render plot
#pdf(file=file.path(outDir, paste(basename(file), ".pdf", sep='')),
#       onefile=TRUE,
#       title=paste(basename(file), ".pdf", sep=''),
#       paper="special", height=11.7, width=16.5)

p1 <- p0 + geom_point(colour="red") + guides(colour= FALSE) +
	geom_smooth(colour="blue") + guides(colour= FALSE) +
	xlim(opt$xlim_min, opt$xlim_max) +
	ylim(opt$ylim_min, opt$ylim_max) +
	xlab("Genomic k-mer abundance") +
    ylab("Genomic k-mers count at that abundance") +
    ggtitle(paste(basename(file), sep=''))

#dev.off()

# render and save plot
ggsave(p1, file=paste(basename(file), "_xlims_", opt$xlim_min, "_", opt$xlim_max, "_ylims_", opt$ylim_min, "_", opt$ylim_max, ".pdf", sep=''),
       path=outDir,
       height=11.7, width=16.5, unit="in")

## estimate cutoff
if ( opt$cutoffestim ) {
	write("Estimate k-mer coverage cut-off ... ", stderr())

	r <- print(p1)
    valleyx=r$data[[2]]$x[which(diff(sign(diff(r$data[[2]]$y)))==2)+1]
	
	# algo 1:
	# 1) search for the first local minimum from 0 on the smoothed data
	# 2) from the local minimum, search for the first local maximum
	# 3) first linear regression
	# 3.1) get the middle point between the point and the local minimum
	# 3.2) get data between the 1st and 3rd quartiles around the middle point
	# 3.3) plot the linear fit of these 2 quartile points using non smoothed data
	# 4) second linear regression: idem but between the local minimum and the local maximum
	# 5) find intersection between the 2 linear fits
	# 6) plot the x coordinate and vline of the intersection point  
	flocmin <- r$data[[1]]$x[which(r$data[[1]]$x == round(valleyx[1]))]
	peakx <- r$data[[2]]$x[which(diff(sign(diff(r$data[[2]]$y)))==-2)]
	flocmax <- r$data[[1]]$x[which(r$data[[1]]$x == round(peakx[1]))]
	descending <- subset(r$data[[1]], x >= 0 & x <= flocmin & !is.na(x), x:y)
	ascending <- subset(r$data[[1]], x >= flocmin & x <= flocmax & !is.na(x), x:y)
	#descending.quantx <- quantile(descending$x, na.rm=T)
	#ascending.quantx <- quantile(ascending$x, na.rm=T)
	#descending.quanty <- quantile(descending$y, na.rm=T)
	#ascending.quanty <- quantile(ascending$y, na.rm=T)
	# fit on y 1st and 3rd quartiles: not as useful as x fit
	#descending.y1quart <- descending.quanty[2]
	#descending.y3quart <- descending.quanty[4]
	#descending.tofit <- subset(descending, y >= descending.y1quart & y <= descending.y3quart)
	#ascending.y1quart <- ascending.quanty[2]
	#ascending.y3quart <- ascending.quanty[4]
	#ascending.tofit <- subset(ascending, y >= ascending.y1quart & y <= ascending.y3quart)
	#descending.yfit <- lm(y ~ x, data=subset(descending, y >= descending.y1quart & y <= descending.y3quart, x:y))
	#ascending.yfit <- lm(y ~ x, data=subset(ascending, y >= ascending.y1quart & y <= ascending.y3quart, x:y)) 
	# fit on x 2nd and 3rd quartiles for descending: best
	#descending.x1quart <- descending.quantx[2]
	#descending.x2quart <- descending.quantx[3]
	#descending.x3quart <- descending.quantx[5]
	descending.pctl.xlims <- quantile(descending$x, c(opt$desc_pctl_min, opt$desc_pctl_max)/100, na.rm=T)
	# fit on x 1st and 2nd quartiles for ascending: best
	#ascending.x1quart <- ascending.quantx[2]	
	#ascending.x2quart <- ascending.quantx[3]
	ascending.pctl.xlims <- quantile(ascending$x, c(opt$asc_pctl_min, opt$asc_pctl_max)/100, na.rm=T)
	## fits
	descending.xfit <- lm(y ~ x, data=subset(descending, x >= descending.pctl.xlims[[1]] & x <= descending.pctl.xlims[[2]], x:y))
	ascending.xfit <- lm(y ~ x, data=subset(ascending, x >= ascending.pctl.xlims[[1]] & x <= ascending.pctl.xlims[[2]], x:y))
	# save fits in a dataframe then plot; search for intercept between fits and finally plot vline
	cm <- rbind(coef(descending.xfit), coef(ascending.xfit))
	fits.intersect <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
	cm.df <- as.data.frame(cm)
	colnames(cm.df)[1] <- "y.intercept"
	p2 <- p1 + geom_abline(data=cm.df, aes(intercept=y.intercept, slope=x, colour=x), linetype=2) +
				 geom_vline(xintercept = fits.intersect[1], colour="dark green", linetype="longdash") +
				 geom_text(aes(fits.intersect[1], 0, label = paste("cut-off=", round(fits.intersect[1]), sep=''), hjust=1.25), colour="dark green")
	# plot peakx
	peaks <- data.frame( x=peakx, y=r$data[[2]]$y[which(r$data[[2]]$x %in% peakx)])
	max_peaky <- max(peaks$y)
	## plot x and y peaks lines and coordinates
	for ( i in 1:length(peaks$x)) {
		print(paste("Peak", i, "X:", peaks$x[i], ", Y:", peaks$y[i], sep=" "), file=stderr())
		xx=c(0,peaks$x[i],peaks$x[i])
		yy=c(peaks$y[i],peaks$y[i],0)
		p4 <- p2 + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=paste("#", i,"\n","x=", round(xx[2],2),"\n", "y=", round(yy[2]), sep=''), x=xx[2], y=ifelse((yy[2]/max_peaky)<=0.5, max_peaky/2, max_peaky*1.3), size=5, colour="black") + guides(colour=FALSE, size=FALSE)
		p4 <- p4 + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=round(yy[1]), x=ifelse(opt$xlim_max <=100, -2, -10), y=yy[1]*1.05, size=4, colour="red") + guides(colour=FALSE, size=FALSE)
		p4 <- p4 + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=round(xx[2],2), x=xx[2], y=ifelse((yy[2]/max_peaky)<=.5, max_peaky*(-0.1), max_peaky*(-0.05)), size=4, colour="red") + guides(colour=FALSE, size=FALSE)
	}

	# render and save plot
    ggsave(p4, file=paste(basename(file), 
			"_xlims_", opt$xlim_min, "_", opt$xlim_max, 
			"_ylims_", opt$ylim_min, "_", opt$ylim_max, 
			"_fitdesclims_", opt$desc_pctl_min, "_", opt$desc_pctl_max,
			"_fitasclims_", opt$asc_pctl_min, "_", opt$asc_pctl_max,
			"_fits_intersect_", round(fits.intersect[1]), ".pdf", sep=''),
           path=outDir,
           height=11.7, width=16.5, unit="in")

	# algo 2:
	# 1) use the smoothed data to find the first local minimum from 0
	# 2) plot the x coordinate and vline of the local minimum
	p3 <- p1 + geom_vline(xintercept = valleyx, colour="blue", linetype = "longdash") +
			 	 geom_text(aes(valleyx, 0, label = paste("cut-off=", round(valleyx,2), sep=''), hjust = -0.1), colour="blue")
	# plot peakx
	## plot x and y peaks lines and coordinates
    for ( i in 1:length(peaks$x)) {
        print(paste("Peak", i, "X:", peaks$x[i], ", Y:", peaks$y[i], sep=" "), file=stderr())
        xx=c(0,peaks$x[i],peaks$x[i])
        yy=c(peaks$y[i],peaks$y[i],0)
        p5 <- p3 + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=paste("#", i,"\n","x=", round(xx[2],2),"\n", "y=", round(yy[2]), sep=''), x=xx[2], y=ifelse((yy[2]/max_peaky)<=0.5, max_peaky/2, max_peaky*1.3), size=5, colour="black") + guides(colour=FALSE, size=FALSE)
        p5 <- p5 + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=round(yy[1]), x=ifelse(opt$xlim_max <=100, -2, -10), y=yy[1]*1.05, size=4, colour="red") + guides(colour=FALSE, size=FALSE)
        p5 <- p5 + geom_text(aes(family="Helvetica", fontface="plain", lineheight=.8), label=round(xx[2],2), x=xx[2], y=ifelse((yy[2]/max_peaky)<=.5, max_peaky*(-0.1), max_peaky*(-0.05)), size=4, colour="red") + guides(colour=FALSE, size=FALSE)
    }

	# render and save plot
	ggsave(p5, file=paste(basename(file), "_xlims_", opt$xlim_min, "_", opt$xlim_max, "_ylims_", opt$ylim_min, "_", opt$ylim_max, "_smoothed_cutoff_", round(valleyx), ".pdf", sep=''),
		   path=outDir,
		   height=11.7, width=16.5, unit="in")

}

## execution time: end
T2<-Sys.time()
Tdiff= difftime(T2, T1)
cat("Execution time : ", Tdiff,"seconds\n", file=stderr())

