#!/usr/bin/env bash

#
# MITOLOGY PIPELINE
#

# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>, Delphine Charif <Delphine.Charif@versailles.inra.fr>

# This script provides a pipeline for plant organites (chloroplast and mitochondrion) assembly 
# GRANT: ANR BIOADAPT CYTOPHENO 2013

# This software is governed by the CeCILL license, Version 2.0 (the "License"), under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license, Version 2.0 (the "License"), as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt".

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license, Version 2.0 (the "License"), and that you accept its terms.

# Date: 2014-09-25

declare -r VERSION="dev"

########################
# SECTION CONFIGURATION
#######################

### SESSION VARIABLES ###

NAMESPACE="MITOLOGY"

WORKING_DIR=$(pwd)
DATE=$(date '+%F_%Hh%Mm%Ss')
SESSION_ID=$(date '+%Y%m%d%H%M%S')
EXECUTED_COMMAND="$0 $*"
SESSION_TAG=${NAMESPACE}_${USER}_${SESSION_ID}

LOG_DIR="log"
DEBUGFILE=${SESSION_TAG}.log
ERROR_TMP="/tmp/$(basename ${0%.*})_error_${SESSION_TAG}.log"

[[ $VERSION -eq "dev" ]] && PROG_PATH=$(realpath $(dirname $0));PIPELINE_USER_CONFIG=${PROG_PATH}/../share/mitology-pipeline/etc/mitology-pipeline_user.config || PIPELINE_USER_CONFIG=/usr/local/share/mitology-pipeline/etc/mitology-pipeline_user.config

PIDS_ARR=()
WAITALL_TIMEOUT=259200
WAITALL_INTERVAL=60
WAITALL_DELAY=60

### LOGGING CONFIGURATION ###

# load log4sh (disabling properties file warning) and clear the default
# configuration
LOG4SH_CONFIGURATION='none' . /usr/local/share/log4sh/build/log4sh 2>/dev/null
[[ $? != 0 ]] && $(echo "Error loading log4sh lib" >&2; exit 1)
log4sh_resetConfiguration

# set the global logging level
logger_setLevel DEBUG

# add and configure a FileAppender that outputs to STDERR
logger_addAppender stderr
appender_setType stderr FileAppender
appender_file_setFile stderr STDERR
appender_setLevel stderr FATAL
appender_setLayout stderr PatternLayout
appender_setPattern stderr '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions stderr
appender_exists stderr && logger_debug "Standard error appender is enabled." || logger_warn "Standard error appender was not enabled. Maybe a log4sh error occured."

# add and configure console appender that outputs to standard output
logger_addAppender console
appender_setType console ConsoleAppender
appender_setLevel console INFO
appender_setLayout console PatternLayout
appender_setPattern console '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions console
appender_exists console && logger_debug "Console appender is enabled." || logger_warn "Console appender was not enabled. Maybe a log4sh error occured."

### LOAD LIB ###

# bash-common lib
[[ $VERSION == "dev" ]] && LIB_PATH=$(realpath $(dirname $0))/../../bash-common/share/bash-common/lib/bash-common_lib.inc || LIB_PATH=/usr/local/share/bash-common/lib/bash-common_lib.inc

logger_debug "[Library] Loading $LIB_PATH"
. $LIB_PATH
if [[ $? -ne 0 ]]; then
	logger_fatal "Error loading bash common lib: $LIB_PATH"
	exit 1
fi

### USAGE ###
Usage()
{
printf %s "\
Program: $(basename $0)
Version: $VERSION

Copyright 2014 Joseph Tran <Joseph.Tran@versailles.inra.fr> & Delphine Charif <Delphine.Charif@versailles.inra.fr>

Licensed under the CeCILL License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

Usage: $(basename $0) -c|--configfile CONFIG_FILE -o|--out_dir OUTPUT_DIR [-d|--debug] [-e|--email_address VALID_EMAIL_ADDR]

Options:
-c|--config_file CONFIG_FILE            The user configuration file listing the data samples paths and tetrad analysis parameters.
                                        You can get a copy there: $PIPELINE_USER_CONFIG.
-C|--kmer_abund_cutoff INT				The k-mer abundance cutoff below which k-mers are trimmed with khmer (filter_abund.py) corresponding to errors and contaminants. This value overrides the one given in the CONFIG_FILE. 
-o|--out_dir OUTPUT_DIR                 The output directory.
-d|--debug                              Enable debugging mode in the console.
-e|--email_address VALID_EMAIL_ADDR     An optional but valid email address to send pipeline job/error status notifications
-h|--help                               Displays this message.

"
}

### NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately;
CONFIGURE_OPTS=`getopt -o hc:C:o:e:d --long help,config_file:,kmer_abund_cutoff:,out_dir:,debug,email_address: \
    -n 'mitology-pipeline.sh' -- "$@"`

if [[ $? != 0 ]] ; then Usage >&2 ; exit 1 ; fi

# Note the quotes around `$CONFIGURE_OPTS'
eval set -- "$CONFIGURE_OPTS"

while true; do
    case "$1" in
        -h | --help ) Usage >&2; exit 1;;
        -c | --config_file ) CONFIGFILE="$2"; shift 2 ;;
		-C | --kmer_abund_cutoff ) KMER_ABUND_CUTOFF="$2"; shift 2;;
        -o | --out_dir ) OUTPUT_DIR="$2"; shift 2 ;;
        -d | --debug )
                    appender_setLevel console DEBUG;
                    appender_activateOptions console;
                    shift 1 ;;
        -e | --email_address ) EMAIL="$2"; shift 2 ;;
        -- ) shift; break ;;
        * ) break ;;
    esac
done

### VALIDATION ###
if [[ ! -s $CONFIGFILE ]]; then
    logger_fatal "Config file, $CONFIGFILE, does not exist or is empty. See Usage with --help option.";
    exit 1;
fi

#################
# PIPELINE STEPS
#################

# CREATE OUTPUT DIR
# LOAD CONFIG
# SET GENOME PATH AND INDEXES
# K-MER ABUNDANCE FILTERING
# ASSEMBLY
# DOT PLOT AGAINST REFERENCE GENOME
# INTERNAL CONSISTENCY: MAPPING READS AGAINST REFERENCE GENOME + FILTERING
# CLEANING

#===================
# OUTPUT DIRECTORY
#===================

echo "Start running $NAMESPACE pipeline (version: $VERSION)." | tee $ERROR_TMP 2>&1 | logger_info
echo "Executed command: $0 $*" | tee -a $ERROR_TMP 2>&1 | logger_info

#
# Create a directory named with OUTPUT_DIR value, to save all outputs
#
echo "Creating $OUTPUT_DIR directory ..." | tee -a $ERROR_TMP 2>&1 | logger_info
if [[ -d $OUTPUT_DIR ]]; then
    echo "OK $OUTPUT_DIR directory already exists. Will output all output files in this directory." | tee -a $ERROR_TMP 2>&1  | logger_info
else
    mkdir $OUTPUT_DIR 2>>$ERROR_TMP
    rtrn=$?
    out_dir_failed_msg="[Output directory] Failed. Output directory, $OUTPUT_DIR, was not created."
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
    exit_on_error "$ERROR_TMP" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
    echo "$(date '+%Y-%m-%d %T') [Output directory] OK $OUTPUT_DIR directory was created successfully. Will output all output files in this directory." | tee -a $ERROR_TMP 2>&1 | logger_info
fi

# Create log directory
echo "Creating $LOG_DIR directory ..." | tee -a $ERROR_TMP 2>&1 | logger_info
if [[ -d $OUTPUT_DIR/$LOG_DIR ]]; then
    echo "OK $OUTPUT_DIR/$LOG_DIR directory already exists. Will write log files in this directory." | tee -a $ERROR_TMP 2>&1 | logger_info
else
    mkdir $OUTPUT_DIR/$LOG_DIR 2>>$ERROR_TMP
    rtrn=$?
    log_dir_failed_msg="[Log directory] Failed Log directory, $OUTPUT_DIR/$LOG_DIR, was not created."
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$log_dir_failed_msg"
    exit_on_error "$ERROR_TMP" "$log_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
    echo "$(date '+%Y-%m-%d %T') [Log directory] OK $OUTPUT_DIR/$LOG_DIR directory was created sucessfully. Will write log files in this directory." | tee -a $ERROR_TMP 2>&1 | logger_info
fi

#==================================
# Enable the pipeline debug logger
#==================================

logger_addAppender debuggerF
appender_setType debuggerF FileAppender
appender_file_setFile debuggerF $OUTPUT_DIR/$LOG_DIR/$DEBUGFILE
appender_setLevel debuggerF DEBUG
appender_setLayout debuggerF PatternLayout
appender_setPattern debuggerF '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions debuggerF
appender_exists debuggerF && cat $ERROR_TMP | logger_info
appender_exists debuggerF && logger_info "Debugging infos will be output to $OUTPUT_DIR/$LOG_DIR/$DEBUGFILE file." || logger_warn "The debugger file appender was not enabled. Maybe a log4sh error occured."

#=============
# LOAD CONFIG
#=============

# set backup config file variable
BACKUPED_CONFIG_FILE=$OUTPUT_DIR/$(basename $CONFIGFILE)

# 1. Backup session user config file in session dir if not exist
logger_info "[Check config: session user config file] Backuping session user config file into session directory ..."
cp $CONFIGFILE $OUTPUT_DIR/. 2>$ERROR_TMP
rtrn=$?
cp_user_config_failed_msg="[Check config: session user config file] Failed backuping session user config file into session directory."
[[ "$rtrn" -ne 0 ]] && logger_fatal "$cp_user_config_failed_msg"
exit_on_error "$ERROR_TMP" "$cp_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
logger_info "[Check config: session user config file] Will use backuped session user config file: $BACKUPED_CONFIG_FILE" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Load config parameters from backuped session user config file
logger_info "[Check config: session user config file] Loading session user config parameters from $BACKUPED_CONFIG_FILE file ..."
load_user_config_failed_msg="[Check config: session user config file] Failed loading session user config parameters from $BACKUPED_CONFIG_FILE file."
for cfg in $(get_config_sections $BACKUPED_CONFIG_FILE 2>$ERROR_TMP;); do
    rtrn=$?
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
    exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
    logger_debug "--- Config section [${cfg}] ---"
    unset $(set |awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN {
          cfg = toupper(cfg);
          prefix = toupper(prefix);
          pattern = "\^" prefix "_" cfg "_"
       }
       $0~pattern { print $1 }' 2>$ERROR_TMP ) 2>>$ERROR_TMP
    rtrn=$?
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
    exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
    CONFIG_PARAMS=$(format_config_params $BACKUPED_CONFIG_FILE ${cfg} ${NAMESPACE} 2>$ERROR_TMP)
    eval "${CONFIG_PARAMS}"
    rtrn=$?
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
    exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
    for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
        logger_debug "$params"
    done
    rtrn=$?
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
    exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
done
logger_info "[Check config: session user config file] OK Session user config file, $BACKUPED_CONFIG_FILE, was loaded successfully."

#=============================================
# REFERENCE GENOME SEQUENCE AND INDEXES PATHS
#=============================================

logger_info "[Genome sequences and index paths] set path variables ..."

### set refrence genome sequence base path
declare -r genome_base_path=$(toupper ${NAMESPACE}_paths)_GENOMES_BASE_PATH

if [[ -z ${!genome_base_path} && ! -d ${!genome_base_path} ]]; then
    logger_fatal "An error occured while setting genome base path variable."
    exit 1
fi
logger_debug "[Genome base path] ${genome_base_path}=${!genome_base_path}"

### set refrence genome indexes base path
declare -r genome_index_path=$(toupper ${NAMESPACE}_paths)_INDEXES_BASE_PATH
if [[ -z ${!genome_index_path} ]]; then
    logger_fatal "An error occured while setting genome indexes path variable."
    exit 1
fi
logger_debug "[Genome index path] ${genome_index_path}=${!genome_index_path}"

### SET GENOME SAMTOOLS INDEX PATH RELATIVE TO CURRENT VERSION/TOOL

declare -r genome_samtools_path=$(toupper ${NAMESPACE}_paths)_SAMTOOLS_INDEXES
if [[ -z ${!genome_samtools_path} ]]; then
    logger_fatal "An error occured while setting genome samtools indexes path variable."
    exit 1
fi
logger_debug "[Genome index path] ${genome_samtools_path}=${!genome_samtools_path}"

#### reference
declare -r ga_ref=$(toupper ${NAMESPACE}_genome_alias )_ref
if [[ -z ${!ga_ref} ]]; then
    logger_fatal "An error occured while setting genome alias variable for reference genome."
    exit 1
fi
logger_debug "[Genome alias] ${ga_ref}=${!ga_ref}"

eval "$(toupper ${NAMESPACE}_paths)_ref_samtools_index=${!genome_index_path}/${!genome_samtools_path}/$(get_tool_version samtools)/${!ga_ref}/${!ga_ref}"
declare -r ref_samtools_index_path=$(toupper ${NAMESPACE}_paths)_ref_samtools_index
if [[ ! -s ${!ref_samtools_index_path} ]]; then
    logger_fatal "An error occured while setting genome samtools index path for reference genome."
    exit 1
fi
IDX_FILES=($(ls ${!ref_samtools_index_path}*))
if [[ ${#IDX_FILES[@]} -le 0 ]]; then
    logger_fatal "An error occured while checking genome samtools index files for the reference genome."
    exit 1
fi
logger_info "[Genome samtools index path] ${ref_samtools_index_path}=${!ref_samtools_index_path}"

# call reference directly
#eval echo -e \$"$(toupper ${NAMESPACE}_paths)_ref_fasta"























