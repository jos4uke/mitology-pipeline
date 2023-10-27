#!/usr/bin/env bash

#
# split_files.sh
#

# Copyright (C) 2014 INRAE

# author: Joseph Tran <joseph.tran@inrae.fr>

# This file is part of mitology-pipeline.

# mitology-pipeline is free software: you can redistribute it   
# and/or modify it under the terms of the GNU General Public License as  
# published by the Free Software Foundation, either version 3 of  
# the License, or (at your option) any later version.  
#  
# mitology-pipeline is distributed in the hope that it will be useful,  
# but WITHOUT ANY WARRANTY; without even the implied warranty of  
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  
# GNU General Public License for more details.  
#  
# You should have received a copy of the GNU General Public Licence  
# along with mitology-pipeline. If not, see <https://www.gnu.org/licenses/>  

# Date: 2014-11-18

declare -r VERSION="dev"
PROG_NAME="$(basename ${BASH_SOURCE[0]})"

### LOGGING CONFIGURATION ###

# load log4sh (disabling properties file warning) and clear the default
# configuration
LOG4SH_CONFIGURATION='none' . /usr/local/share/log4sh/log4sh 2>/dev/null
[[ $? != 0 ]] && $(echo "Error loading log4sh lib" >&2; exit 1)
log4sh_resetConfiguration

# set the global logging level
# add and configure a FileAppender that outputs to STDERR
logger_addAppender ${PROG_NAME}.stderr2
appender_setType ${PROG_NAME}.stderr2 FileAppender
appender_file_setFile ${PROG_NAME}.stderr2 STDERR
appender_setLevel ${PROG_NAME}.stderr2 FATAL
appender_setLayout ${PROG_NAME}.stderr2 PatternLayout
appender_setPattern ${PROG_NAME}.stderr2 "%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] ${PROG_NAME} - %m"
appender_activateOptions ${PROG_NAME}.stderr2
# add and configure console appender that outputs to standard output
logger_addAppender ${PROG_NAME}.console2
appender_setType ${PROG_NAME}.console2 ConsoleAppender
appender_setLevel ${PROG_NAME}.console2 INFO
appender_setLayout ${PROG_NAME}.console2 PatternLayout
appender_setPattern ${PROG_NAME}.console2 "%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] ${PROG_NAME} - %m"
appender_activateOptions ${PROG_NAME}.console2
appender_exists ${PROG_NAME}.console2 && logger_debug "Console appender is enabled." || logger_warn "Console appender was not enabled.
 Maybe a log4sh error occured."

### PRINT PROG VERSION AND COMMAND
echo "$PROG_NAME (version: $VERSION)." | tee $ERROR_TMP 2>&1 | logger_info
echo "Executed command: $0 $@" | tee -a $ERROR_TMP 2>&1 | logger_info

### USAGE ###
Usage()
{
printf %s "\
Program: ${PROG_NAME}
Version: $VERSION

Copyright 2014 Joseph Tran <Joseph.Tran@versailles.inra.fr>

Licensed under the CeCILL License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

Usage: ${PROG_NAME} -n|--number CHUNKS -o|--out_dir OUTPUT_DIR INFILE
                   [-d|--debug] [-h|--help]

Mandatory:
	-n|--number CHUNK                       
		generate CHUNKS output files.  See below
                                        
	-o|--out_dir OUTDIR
		the output directory.
	
Options:
	-d|--debug
		Enable debugging mode in the console.

	-h|--help
		Displays this message.

"
}

# Configure opts
### NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately;
CFG_OPTS=$(getopt -o hn:o:d --long help,number:,out_dir:,debug -n 'split_files.sh' -- "$@")

if [[ $? != 0 ]] ; then Usage >&2 ; exit 1 ; fi

# Note the quotes around `$CONFIGURE_OPTS'
eval set -- "$CFG_OPTS"

while true; do
    case "$1" in
        -h | --help ) Usage >&2; exit 1;;
        -n | --number ) CHUNKS="$2"; shift 2 ;;
		-o | --out_dir ) OUTPUT_DIR="$2"; shift 2 ;;
        -d | --debug )
						DEBUG="debug";
						appender_setLevel ${PROG_NAME}.console2 DEBUG
						appender_activateOptions ${PROG_NAME}.console2
						shift ;;
        -- ) shift; break ;;
        * ) break ;;
    esac
done

### VALIDATION ###

# mandatory
if [[ -z $CHUNKS ]]; then
	logger_fatal "Error: the number of chunks must be not null. See Usage using --help option."
	exit 1
else
	re='^[0-9]+$'
	if ! [[ $CHUNKS =~ $re ]] ; then
		logger_fatal "Error: Not a number";  exit 1
	else
		logger_debug "OK: will split file into $CHUNKS chunks."
	fi
fi

# optional
if [[ -z $OUTPUT_DIR ]]; then
    logger_warn "Output directory must be not null. See Usage with --help option.";
    OUTPUT_DIR=$(pwd)
	logger_debug "Will output splitted files into $OUTPUT_DIR "
else
	if ! [[ -d $OUTPUT_DIR ]]; then
		logger_warn "Output directory, ${OUTPUT_DIR}, does not exist.";
		OUTPUT_DIR=$(pwd)
		logger_debug "Will output splitted files into $OUTPUT_DIR "
	else
		logger_debug "OK: will output splitted files into $OUTPUT_DIR "
	fi	
fi

# infile
INFILE=$1
if [[ -s $INFILE ]]; then
	logger_info "OK: infile; ${INFILE}, exists and is not empty."
else
	logger_fatal "Error: infile, ${INFILE}, does not exit or is empty."
	exit 1
fi

###
### SPLIT FILE INTO N CHUNKS
###

# split character
split_char="_"

# get total number of lines
total_nb_lines=$(wc -l <${INFILE})
logger_info "OK: total number of lines = ${total_nb_lines}"

# get number of lines per file 
((lines_per_file = (total_nb_lines + CHUNKS - 1) / CHUNKS))
logger_info "OK: lines per file = ${lines_per_file}"

# split into equivalent N files of X lines
OUTFILES=${OUTPUT_DIR}/$(basename ${INFILE})${split_char}
split --lines=${lines_per_file} -d ${INFILE} ${OUTFILES}
if [[ status=$? -ne 0 ]]; then
	logger_fatal "Error: $PROG_NAME quits without splitting file, ${INFILE}, into $CHUNKS chunks of ${lines_per_file} lines each."
	logger_fatal "exit status code: $status"
	exit 1	
fi

# report
logger_debug "Input file               = ${INFILE}"
logger_debug "Total number of lines    = ${total_nb_lines}"
logger_debug "Number of splitted files = $(ls ${OUTFILES}* | wc -l)"
logger_debug "Number of lines per file = ${lines_per_file}"
logger_debug "$(wc -l ${OUTFILES}*)"













