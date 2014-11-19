#!/usr/bin/env bash

#
# MITOLOGY PIPELINE
#

# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This script provides a pipeline for plant organites (chloroplast and mitochondrion) assembly 
# GRANT: ANR BIOADAPT CYTOPHENO 2012

# based on the preliminary work of 2 interns (03-07/2014): Alexandrina Brodrug and Myriam Shafie

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
ERROR_TMP_MODEL="/tmp/$(basename ${0%.*})_error_${SESSION_TAG}.XXXXXX"
ERROR_TMP=$(mktemp "$ERROR_TMP_MODEL")
[[ $? -eq 0 ]] && echo "Execute $(basename ${0%.*})">>$ERROR_TMP || exit 1 

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
# add and configure a FileAppender that outputs to STDERR
logger_addAppender stderr
appender_setType stderr FileAppender
appender_file_setFile stderr STDERR
appender_setLevel stderr FATAL
appender_setLayout stderr PatternLayout
appender_setPattern stderr '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions stderr
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

# mitology-pipeline lib
[[ $VERSION == "dev" ]] && LIB_PATH=$(realpath $(dirname $0))/../share/mitology-pipeline/lib/mitology-pipeline_lib.inc || LIB_PATH=/usr/local/share/mitology-pipeline/lib/mitology-pipeline_lib.inc

logger_debug "[Library] Loading $LIB_PATH"
. $LIB_PATH
if [[ $? -ne 0 ]]; then
	logger_fatal "Error loading mitology pipeline lib: $LIB_PATH"
	exit 1
fi

### SCRIPTS ###

# scripts path
[[ $VERSION == "dev" ]] && SCRIPTS_PATH=$(realpath $(dirname $0))/../share/mitology-pipeline/scripts || SCRIPTS_PATH=/usr/local/share/mitology-pipeline/scripts

logger_debug "[Scripts] Setting $SCRIPTS_PATH"


### SUPPORTED ASSEMBLERS/SCAFFOLDERS ###

# supported assemblers/scaffolders list
SUPPORTED_ASSEMBLERS=( meta-velvetg )
SUPPORTED_SCAFFOLDERS=( meta-velvetg )
ASSEMBLER_DEFAULT=meta-velvetg
SCAFFOLDER_DEFAULT=meta-velvetg

# supported assemblers/scaffolders combinations

### USAGE ###
Usage()
{
printf %s "\
Program: $(basename $0)
Version: $VERSION

Copyright 2014 Joseph Tran <Joseph.Tran@versailles.inra.fr>

Licensed under the CeCILL License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

Usage: $(basename $0) -c|--configfile CONFIG_FILE -o|--out_dir OUTPUT_DIR [-A|--assembler ASSEMBLER] [-S|--scaffolder SCAFFOLDER] [-d|--debug] [-C|--kmer_abund_cutoff INT] [-e|--email_address VALID_EMAIL_ADDR]

Mandatory:
-c|--config_file CONFIG_FILE            The user configuration file, CONFIG_FILE, listing the pipeline parameters by section.
                                        You can get a copy there: $PIPELINE_USER_CONFIG.
-o|--out_dir OUTPUT_DIR                 The output directory.

Options:
-C|--kmer_abund_cutoff INT              The k-mer abundance cutoff below which k-mers are trimmed with khmer (filter_abund.py) 
                                        corresponding to errors and contaminants. This value overrides the one given 
                                        in the CONFIG_FILE .
-A|--assembler ASSEMBLER                The assembler, given by ASSEMBLER, to use. [default: $ASSEMBLER_DEFAULT]
                                        Supported list of assemblers: "${SUPPORTED_ASSEMBLERS[@]}"
-S|--scaffolder SCAFFOLDER              The scaffolder, given by SCAFFOLDER, to use. [default: $SCAFFOLDER_DEFAULT]
                                        Supported list of scaffolders: "${SUPPORTED_SCAFFOLDERS[@]}"
-d|--debug                              Enable debugging mode in the console.
-e|--email_address VALID_EMAIL_ADDR     An optional but valid email address to send pipeline job/error status notifications
-h|--help                               Displays this message.

"
}

### NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately;
CONFIGURE_OPTS=`getopt -o hc:A:S:C:o:e:d --long help,config_file:,assembler:,scaffolder:,out_dir:,kmer_abund_cutoff:,debug,email_address: \
    -n 'mitology-pipeline.sh' -- "$@"`

if [[ $? != 0 ]] ; then Usage >&2 ; exit 1 ; fi

# Note the quotes around `$CONFIGURE_OPTS'
eval set -- "$CONFIGURE_OPTS"

while true; do
    case "$1" in
        -h | --help ) Usage >&2; exit 1;;
        -c | --config_file ) CONFIGFILE="$2"; shift 2 ;;
        -o | --out_dir ) OUTPUT_DIR="$2"; shift 2 ;;
		-A | --assembler ) ASSEMBLER_CLI="$2"; shift 2 ;;
		-S | --scaffolder ) SCAFFOLDER_CLI="$2"; shift 2 ;;
		-C | --kmer_abund_cutoff ) KMER_ABUND_CUTOFF="$2"; shift 2;;
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

# mandatory
if [[ ! -s $CONFIGFILE ]]; then
    logger_fatal "Config file, $CONFIGFILE, does not exist or is empty. See Usage with --help option.";
    exit 1;
fi

if [[ -z $OUTPUT_DIR ]]; then
	logger_fatal "Output directory must be not null. See Usage with --help option.";
	exit 1;
fi

# not mandatory
if [[ -z $ASSEMBLER_CLI ]]; then
	logger_warn "Assembler string must be not null. See Usage with --help option."
else
	# check if assembler is supported else use default
	if [[ ! $(elementIn "$ASSEMBLER_CLI" "${SUPPORTED_ASSEMBLERS[@]}") ]]; then 
		logger_warn "Given assembler, $ASSEMBLER_CLI, is not currently supported. See Usage with --help option."
		exit 1
	else
		logger_info "Replace the default assembler, $ASSEMBLER_DEFAULT, by the user provided assembler, $ASSEMBLER_CLI ."
	fi
fi

if [[ -z $SCAFFOLDER_CLI ]]; then
	logger_warn "Scaffolder string must be not null. See Usage with --help option."
else
	# check if scaffolder is supported else use default
	if [[ ! $(elementIn "$SCAFFOLDER_CLI" "${SUPPORTED_SCAFFOLDERS[@]}") ]]; then
		logger_warn "Given scaffolder, $SCAFFOLDER_CLI, is not currently supported. See Usage with --help option."
		exit 1
	else
		logger_info "Replace the default scaffolder, $SCAFFOLDER_DEFAULT, by the user provided scaffolder, $SCAFFOLDER_CLI ."
	fi
fi


#################
# PIPELINE STEPS
#################

# CREATE OUTPUT DIR
# LOAD CONFIG
# OVERRIDE CONFIG
# SET GENOME PATH AND INDEXES
# CHECKING SAMPLE
# K-MER ABUNDANCE FILTERING: FILTERED READS
# ASSEMBLY: 
# - CONTIGING 
# - SCAFFOLDING
# STATS:
# - COMPASS
# - QUAST
# - AMOS/ASTATS
# ALIGNMENTS: 
# - NUCMER (DOT PLOT AGAINST ONE REFERENCE GENOME)
# - BWA (FOR INTERNAL CONSISTENCY: MAPPING FILTERED READS AGAINST CONTIGS AND SCAFFOLDS)
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
appender_file_setFile debuggerF $(realpath $OUTPUT_DIR)/$LOG_DIR/$DEBUGFILE
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
logger_info "[Check config: session user config file] Will use backuped session user config file: $BACKUPED_CONFIG_FILE" 2>&1

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
        logger_debug "export $params"
		export "$params"
    done
    rtrn=$?
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
    exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
done
logger_info "[Check config: session user config file] OK Session user config file, $BACKUPED_CONFIG_FILE, was loaded successfully."

#=================
# OVERRIDE CONFIG
#=================

logger_info "[Override config] checking for options to override loaded config parameters ..."

### khmer filter abund cutoff
declare -r khmer_filter_abund_cutoff=$(toupper ${NAMESPACE}_khmer_filter_abund)_C
if [[ -z ${!khmer_filter_abund_cutoff} ]]; then
	logger_warn "[Override config] Config khmer filter abundance cutoff variable, ${!khmer_filter_abund_cutoff}, is null. Search for option to override ..."
	if [[ -z $KMER_ABUND_CUTOFF ]]; then
		logger_fatal "[Override config] Khmer abundance cutoff option value is null. Please fill in a cutoff value in config file or on the command line, see usage with --help option."
		exit 1
	else
		eval "$(toupper ${NAMESPACE}_khmer_filter_abund)_C=${KMER_ABUND_CUTOFF}"
		logger_info "[Override config] Overrides config khmer filter abundance cutoff value, NULL, by ${!khmer_filter_abund_cutoff}"
	fi
else
	if [[ -n $KMER_ABUND_CUTOFF ]]; then
		cutoff_old=${!khmer_filter_abund_cutoff}
		eval "$(toupper ${NAMESPACE}_khmer_filter_abund)_C=${KMER_ABUND_CUTOFF}"
		logger_info "[Override config] Overrides config khmer filter abundance cutoff value, $cutoff_old, by ${!khmer_filter_abund_cutoff}"
	else
		logger_info "[Override config] Not overriding the config khmer filter abundance cutoff value, ${!khmer_filter_abund_cutoff}."
	fi
fi

### assembler
declare -r CFG_ASSEMBLER=$(toupper ${NAMESPACE}_contig_assembler)_name
if [[ -z ${!CFG_ASSEMBLER} ]]; then
	logger_warn "[Override config] Config contig_assembler variable, ${!CFG_ASSEMBLER}, is null. Search for option to override ..."
	if [[ -z $ASSEMBLER_CLI ]]; then
		logger_fatal "[Override config] User provided assembler option value is null. Please fill in a contig assembler name value in config file or on the command line, see usage with --help option."
		exit 1
	else
		eval "$(toupper ${NAMESPACE}_contig_assembler)_name=$ASSEMBLER_CLI"
		logger_info "[Override config] Overrides config contig_assembler name value, NULL, by ${!CFG_ASSEMBLER}"
	fi
else
	if [[ -n $ASSEMBLER_CLI ]]; then
		assembler_old=${!CFG_ASSEMBLER}
		eval "$(toupper ${NAMESPACE}_contig_assembler)_name=$ASSEMBLER_CLI"
		logger_info "[Override config] Overrides config contig_assembler name value, $assembler_old, by ${!CFG_ASSEMBLER}"
	else
		logger_info "[Override config] Not overriding the config contig_assembler name value, ${!CFG_ASSEMBLER}."
	fi
fi

### scaffolder
declare -r CFG_SCAFFOLDER=$(toupper ${NAMESPACE}_scaffolder)_name
if [[ -z ${!CFG_SCAFFOLDER} ]]; then
	logger_warn "[Override config] Config scaffolder variable, ${!CFG_SCAFFOLDER}, is null. Search for option to override ..."
	if [[ -z $SCAFFOLDER_CLI ]]; then
		logger_fatal "[Override config] User provided scaffolder option value is null. Please fill in a scaffolder name value in config file or on the command line, see usage with --help option."
		exit 1
	else
		eval "$(toupper ${NAMESPACE}_scaffolder)_name=$SCAFFOLDER_CLI"
		logger_info "[Override config] Overrides config scaffolder name value, NULL, by ${!CFG_SCAFFOLDER}"
	fi
else
	if [[ -n $SCAFFOLDER_CLI ]]; then
		scaffolder_old=${!CFG_SCAFFOLDER}
		eval "$(toupper ${NAMESPACE}_scaffolder)_name=$SCAFFOLDER_CLI"
		logger_info "[Override config] Overrides config scaffolder name value, $scaffolder_old, by ${!CFG_SCAFFOLDER}"
	else
		logger_info "[Override config] Not overriding the config scaffolder name value, ${!CFG_SCAFFOLDER}."
	fi
fi


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

#### reference
declare -r ga_ref=$(toupper ${NAMESPACE}_genome_alias)_ref
if [[ -z ${!ga_ref} ]]; then
    logger_fatal "An error occured while setting genome alias variable for reference genome."
    exit 1
fi
logger_debug "[Genome alias] ${ga_ref}=${!ga_ref}"

##### mito
declare -r ga_mito_ref=$(toupper ${NAMESPACE}_genome_alias)_ref_mito
if [[ -z ${!ga_mito_ref} ]]; then
    logger_fatal "An error occured while setting genome alias variable for mitochondria reference genome."
    exit 1
fi
logger_debug "[Genome alias] ${ga_mito_ref}=${!ga_mito_ref}"
ga_mito_ref_path=${!genome_base_path}/${!ga_mito_ref}/${!ga_mito_ref}.fas
if [[ ! -s ${ga_mito_ref_path} ]]; then
    logger_fatal "Reference genome fasta sequence for mitochondria does not exist or is empty."
	logger_fatal "Please check for given sequence path: ${ga_mito_ref_path}"
    exit 1
fi
logger_info "[Genome alias] mito_ref_path=${ga_mito_ref_path}"
ga_mito_gff_ref=$(toupper ${NAMESPACE}_genome_alias)_ref_mito_gff
if [[ -z ${!ga_mito_gff_ref} ]]; then
    logger_fatal "An error occured while setting gff genes alias variable for mitochondria reference genes gff."
    exit 1
fi
logger_info "[Genome alias] ${ga_mito_gff_ref}=${!ga_mito_gff_ref}"
ga_mito_gff_ref_path=${!genome_base_path}/${!ga_mito_ref}/${!ga_mito_gff_ref}
if [[ ! -s ${ga_mito_gff_ref_path} ]]; then
    logger_fatal "Reference gff genes for mitochondria does not exist or is empty."
    logger_fatal "Please check for given gff genes path: ${ga_mito_gff_ref_path}"
    exit 1
fi
logger_info "[Genome alias] mito_gff_ref_path=${ga_mito_gff_ref_path}"

##### chloro
declare -r ga_chloro_ref=$(toupper ${NAMESPACE}_genome_alias)_ref_chloro
if [[ -z ${!ga_chloro_ref} ]]; then
    logger_fatal "An error occured while setting genome alias variable for chloroplast reference genome."
    exit 1
fi
logger_debug "[Genome alias] ${ga_chloro_ref}=${!ga_chloro_ref}"
ga_chloro_ref_path="${!genome_base_path}/${!ga_chloro_ref}/${!ga_chloro_ref}.fas"
if [[ ! -s ${ga_chloro_ref_path} ]]; then
    logger_fatal "Reference genome fasta sequence for chloroplast does not exist or is empty."
    logger_fatal "Please check for given sequence path: ${ga_chloro_ref_path}"
    exit 1
fi
logger_info "[Genome alias] chloro_ref_path=${ga_chloro_ref_path}"
ga_chloro_gff_ref=$(toupper ${NAMESPACE}_genome_alias)_ref_chloro_gff
if [[ -z ${!ga_chloro_gff_ref} ]]; then
    logger_fatal "An error occured while setting gff genes alias variable for chloroplast reference genes gff."
    exit 1
fi
logger_info "[Genome alias] ${ga_chloro_gff_ref}=${!ga_chloro_gff_ref}"
ga_chloro_gff_ref_path=${!genome_base_path}/${!ga_chloro_ref}/${!ga_chloro_gff_ref}
if [[ ! -s ${ga_chloro_gff_ref_path} ]]; then
    logger_fatal "Reference gff genes for chloroplast does not exist or is empty."
    logger_fatal "Please check for given gff genes path: ${ga_chloro_gff_ref_path}"
    exit 1
fi
logger_info "[Genome alias] chloro_gff_ref_path=${ga_chloro_gff_ref_path}"


### SET GENOME SAMTOOLS INDEX PATH RELATIVE TO CURRENT VERSION/TOOL

declare -r genome_samtools_path=$(toupper ${NAMESPACE}_paths)_SAMTOOLS_INDEXES
if [[ -z ${!genome_samtools_path} ]]; then
    logger_fatal "An error occured while setting genome samtools indexes path variable."
    exit 1
fi
logger_debug "[Genome index path] ${genome_samtools_path}=${!genome_samtools_path}"

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

### SET GENOME BWA INDEX PATH RELATIVE TO CURRENT VERSION/TOOL

#### set current tool version index for the reference genome
declare -r genome_bwa_path=$(toupper ${NAMESPACE}_paths)_BWA_INDEXES
if [[ -z ${!genome_bwa_path} ]]; then
    logger_fatal "An error occured while setting genome bwa indexes path variable."
    exit 1
fi
logger_debug "[Genome index path] ${genome_bwa_path}=${!genome_bwa_path}"

eval "$(toupper ${NAMESPACE}_paths)_ref_bwa_index=${!genome_index_path}/${!genome_bwa_path}/$(get_tool_version bwa)/${!ga_ref}/${!ga_ref}"
declare -r ref_bwa_index_path=$(toupper ${NAMESPACE}_paths)_ref_bwa_index
if [[ -z ${!ref_bwa_index_path} ]]; then
    logger_fatal "An error occured while setting genome bwa index path variable for the ref genome."
    exit 1
fi
IDX_FILES=($(ls ${!ref_bwa_index_path}*))
if [[ ${#IDX_FILES[@]} -le 0 ]]; then
    logger_fatal "An error occured while checking genome bwa index files for the ref genome."
    exit 1
fi
logger_info "[Genome index path] ${ref_bwa_index_path}=${!ref_bwa_index_path}"

# call directly
#eval echo -e \$"$(toupper ${NAMESPACE}_paths)_ref_bwa_index"
# test
#eval ls -lh \$"$(toupper ${NAMESPACE}_paths)_ref_bwa_index*"

#=================
# CHECKING SAMPLE
#=================

declare -r current_sample_alias=$(toupper ${NAMESPACE}_sample)_name_alias
declare -r current_sample_seq_dir=$(toupper ${NAMESPACE}_sample)_seqfile_parent_dir

# check if seq dir exists
if [[ ! -d ${!current_sample_seq_dir} ]]; then
	logger_fatal "[Checking sample] Sample directory, ${!current_sample_seq_dir}, does not exist. Please check the sample directory path."
	exit 1
fi
logger_info "[Checking sample] Sample directory, ${!current_sample_seq_dir}, exists."

# check if all seq files exist
declare -r current_sample_seq_R1=$(toupper ${NAMESPACE}_sample)_seqfile_R1
eval "$(toupper ${NAMESPACE}_sample)_seqfile_R1_path=${!current_sample_seq_dir}/${!current_sample_seq_R1}"
declare -r current_sample_seq_R1_path=$(toupper ${NAMESPACE}_sample)_seqfile_R1_path

declare -r current_sample_seq_R2=$(toupper ${NAMESPACE}_sample)_seqfile_R2
eval "$(toupper ${NAMESPACE}_sample)_seqfile_R2_path=${!current_sample_seq_dir}/${!current_sample_seq_R2}"
declare -r current_sample_seq_R2_path=$(toupper ${NAMESPACE}_sample)_seqfile_R2_path

if [[ ! -s "${!current_sample_seq_R1_path}" ]]; then
	logger_fatal "[Checking sample] Sample R1 seq file, ${!current_sample_seq_R1}, does not exist or is empty."
	exit 1
fi	
logger_info "[Checking sample] Sample R1 seq file, ${!current_sample_seq_R1}, exists."

if [[ ! -s "${!current_sample_seq_R2_path}" ]]; then
	logger_fatal "[Checking sample] Sample R2 seq file, ${!current_sample_seq_R2}, does not exist or is empty."
	exit 1
fi
logger_info "[Checking sample] Sample R2 seq file, ${!current_sample_seq_R2}, exists."

#===============================
# 01. K-MER ABUNDANCE FILTERING
#===============================

# STEPS
## CREATE OUTPUT DIR
## INTERLEAVE READS 
## COUNTING KMERS
## FILTER K-MER ABUNDANCE
## EXTRACT AND SPLIT PAIRED READS

#
# Create K-mer abundance filtering output directory
#
KMER_FILTER_ABUND_OUTDIR="01.K-mer_filter_abund"
logger_info "Creating $KMER_FILTER_ABUND_OUTDIR directory ..." 
if [[ -d $OUTPUT_DIR/$KMER_FILTER_ABUND_OUTDIR ]]; then
    logger_debug "OK $KMER_FILTER_ABUND_OUTDIR directory already exists. Will output all k-mer abundance filtering output files in this directory."
else
    mkdir $OUTPUT_DIR/$KMER_FILTER_ABUND_OUTDIR 2>$ERROR_TMP
    rtrn=$?
    out_dir_failed_msg="[$KMER_FILTER_ABUND_OUTDIR] Failed. K-mer abundance filtering output directory, $KMER_FILTER_ABUND_OUTDIR, was not created."
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
    exit_on_error "$ERROR_TMP" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
    logger_debug "[$KMER_FILTER_ABUND_OUTDIR] OK $KMER_FILTER_ABUND_OUTDIR directory was created successfully. Will output all k-mer abundance filtering output files in this directory."
fi

### Enable the k-mer abundance filtering debug logger
KMER_FILTER_ABUND_DEBUGF=${KMER_FILTER_ABUND_OUTDIR}_debug.log
logger_addAppender kmerFiltAbundF
appender_setType kmerFiltAbundF FileAppender
appender_file_setFile kmerFiltAbundF $(realpath $OUTPUT_DIR)/$KMER_FILTER_ABUND_OUTDIR/$KMER_FILTER_ABUND_DEBUGF
appender_setLevel kmerFiltAbundF DEBUG
appender_setLayout kmerFiltAbundF PatternLayout
appender_setPattern kmerFiltAbundF '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions kmerFiltAbundF
appender_exists kmerFiltAbundF && logger_info "[$KMER_FILTER_ABUND_OUTDIR] Debugging infos on k-mer abundance filtering will be output to $OUTPUT_DIR/$KMER_FILTER_ABUND_OUTDIR/$KMER_FILTER_ABUND_DEBUGF file." || logger_warn "The kmerFiltAbundF debugger file appender was not enabled. Maybe a log4sh error occured."
### error handling
KMER_FILTER_ABUND_ERROR=$OUTPUT_DIR/$KMER_FILTER_ABUND_OUTDIR/${KMER_FILTER_ABUND_OUTDIR}.err

### Substep counter
KMER_FILT_ABUND_COUNTER=0

#
# LINK step
#
# link config loading output to interleaving input 
eval "declare -A $(toupper ${NAMESPACE}_interleaving_input)"
interleaving_input=$(toupper ${NAMESPACE}_interleaving_input)
eval "${interleaving_input}=( [R1]=${!current_sample_seq_R1_path} )"
eval "${interleaving_input}+=( [R2]=${!current_sample_seq_R2_path} )"

### check for interleaving input sequences
if [[ "${!current_sample_seq_R1_path}" == "${!current_sample_seq_R2_path}" ]]; then
	logger_fatal "[TEST] Interleaving error R1 == R2"; exit 1
else
	logger_debug "[TEST] Interleaving OK R1 != R2"
fi
R1_input="${interleaving_input}[R1]"
R2_input="${interleaving_input}[R2]"
logger_debug "Interleaving input R1: ${!R1_input}"
logger_debug "Interleaving input R2: ${!R2_input}"

#
# Interleave reads
#

logger_info "[$KMER_FILTER_ABUND_OUTDIR] Interleaving reads ..."

# set interleaving vars
((KMER_FILT_ABUND_COUNTER+=1))
INTERLEAVING_SUBDIR="${KMER_FILT_ABUND_COUNTER}.Interleaving"
INTERLEAVING_SUBDIR_PATH=$OUTPUT_DIR/$KMER_FILTER_ABUND_OUTDIR/$INTERLEAVING_SUBDIR
INTERLEAVING_ERROR=$INTERLEAVING_SUBDIR_PATH/${!current_sample_alias}_interleaved.err
## define interleaving output
eval "$(toupper ${NAMESPACE}_interleaving_output)=(['filename']=${INTERLEAVING_SUBDIR_PATH}/${!current_sample_alias}_interleaved.fastq)"
declare -r interleaving_output=$(toupper ${NAMESPACE}_interleaving_output)
## build cli
declare -r khmer_interleave_reads=$(toupper ${NAMESPACE}_paths)_khmer_interleave_reads
R1_input="${interleaving_input}[R1]"
R2_input="${interleaving_input}[R2]"
interleaving_cli="${!khmer_interleave_reads} ${!R1_input} ${!R2_input} >${!interleaving_output["filename"]}"
#interleaving_cli+=" 2>${INTERLEAVING_ERROR}"

# check if interleaving subdir exists else create
# set interleaving step input vars (no need to check seq file: already done in sample checking block)
# exists: check interleaving step output var => if output file exists else run step
# not exist: run step
logger_info "Creating $INTERLEAVING_SUBDIR_PATH directory ..."
if [[ -d $INTERLEAVING_SUBDIR_PATH ]]; then
	logger_debug "OK $INTERLEAVING_SUBDIR_PATH directory already exists. Will check for already existing output files."
	# check for existing output files
	if [[ -s ${!interleaving_output["filename"]} ]]; then
		# skip step
		logger_info "[$INTERLEAVING_SUBDIR] Output file already exists: ${!interleaving_output["filename"]}."
		logger_info "[$INTERLEAVING_SUBDIR] Skip interleaving step."
	else
		# run step
		logger_info "[$INTERLEAVING_SUBDIR] Will run interleaving step ..."
		run_cli -c "$interleaving_cli" -t "$INTERLEAVING_SUBDIR" -e "$INTERLEAVING_ERROR" -d
	fi
else
	# create subdir
	mkdir $INTERLEAVING_SUBDIR_PATH 2>$KMER_FILTER_ABUND_ERROR
	rtrn=$?
	out_dir_failed_msg="[$KMER_FILTER_ABUND_OUTDIR] Failed. Interleaving output directory, $INTERLEAVING_SUBDIR_PATH, was not created."
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
	exit_on_error "$KMER_FILTER_ABUND_ERROR" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
	logger_debug "[$KMER_FILTER_ABUND_OUTDIR] OK $INTERLEAVING_SUBDIR_PATH directory was created successfully. Will output interleaving output files in this directory."
	# run step
	logger_info "[$INTERLEAVING_SUBDIR] Will run interleaving step ..."
	run_cli -c "$interleaving_cli" -t "$INTERLEAVING_SUBDIR" -e "$INTERLEAVING_ERROR" -d
fi

#
# LINK step
#
# link interleaving output and hashcount input 
eval "$(toupper ${NAMESPACE}_hashcount_input)=(['input_sequence_filename']=${!interleaving_output['filename']})"
declare -r hashcount_input=$(toupper ${NAMESPACE}_hashcount_input)

#
# Counting k-mers
#

logger_info "[$KMER_FILTER_ABUND_OUTDIR] Counting k-mers ..."

# set hashcount vars
((KMER_FILT_ABUND_COUNTER+=1))
HASHCOUNT_SUBDIR="${KMER_FILT_ABUND_COUNTER}.Hashcount"
HASHCOUNT_SUBDIR_PATH="$OUTPUT_DIR/$KMER_FILTER_ABUND_OUTDIR/$HASHCOUNT_SUBDIR"
HASHCOUNT_ERROR=$HASHCOUNT_SUBDIR_PATH/${!current_sample_alias}_counting.err
declare -r khmer_load_into_counting_k=$(toupper ${NAMESPACE}_khmer_load_into_counting)_k
# define hashcount output
eval "$(toupper ${NAMESPACE}_hashcount_output)=(['output_countingtable_filename']=$HASHCOUNT_SUBDIR_PATH/${!current_sample_alias}_k${!khmer_load_into_counting_k}.hashcount)"
declare -r hashcount_output=$(toupper ${NAMESPACE}_hashcount_output)
# build cli options
khmer_load_into_counting_opts=($(buildCommandLineOptions "khmer_load_into_counting" "$NAMESPACE" 2>$KMER_FILTER_ABUND_ERROR))
rtrn=$?
cli_opts_failed_msg="[$KMER_FILTER_ABUND_OUTDIR] An error occured while building the khmer_load_into_counting command line options for current sample ${!current_sample_alias}."
exit_on_error "$KMER_FILTER_ABUND_ERROR" "$cli_opts_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
hashcount_opts="${khmer_load_into_counting_opts[@]}"
logger_debug "[$KMER_FILTER_ABUND_OUTDIR] khmer_load_into_counting options: $hashcount_opts"
# build cli
declare -r khmer_load_into_counting=$(toupper ${NAMESPACE}_paths)_khmer_load_into_counting
hashcount_cli="${!khmer_load_into_counting} $hashcount_opts ${!hashcount_output['output_countingtable_filename']} ${!hashcount_input['input_sequence_filename']}"
#hashcount_cli+=" 2>${HASHCOUNT_ERROR} | logger_debug"

# check if hashcount subdir exists else create
# set hashcount step input vars (see previous LINK step)
# exists: check hashcount step output var => if output file exists skip step
# not exist: run step
logger_info "[$KMER_FILTER_ABUND_OUTDIR] Creating $HASHCOUNT_SUBDIR_PATH directory ..."
if [[ -d $HASHCOUNT_SUBDIR_PATH ]]; then
	logger_debug "[$KMER_FILTER_ABUND_OUTDIR] OK $HASHCOUNT_SUBDIR_PATH directory already exists. Will check for already existing output files."

	# check for existing output files
	if [[ -s ${!hashcount_output['output_countingtable_filename']} ]]; then
		# skip step
		logger_info "[$HASHCOUNT_SUBDIR] Output file already exists: ${!hashcount_output['output_countingtable_filename']}."
		logger_info "[$HASHCOUNT_SUBDIR] Skip hashcount step."
	else
		# run step
		logger_info "[$HASHCOUNT_SUBDIR] Will run hashcount step ... "
		run_cli -c "$hashcount_cli" -t "$HASHCOUNT_SUBDIR" -e "$HASHCOUNT_ERROR" -d
	fi
else
	# create subdir
	mkdir $HASHCOUNT_SUBDIR_PATH 2>$KMER_FILTER_ABUND_ERROR
	rtrn=$?
	out_dir_failed_msg="[$KMER_FILTER_ABUND_OUTDIR] Failed. Hashcount output directory, $HASHCOUNT_SUBDIR_PATH, was not created."
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
	exit_on_error "$KMER_FILTER_ABUND_ERROR" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
	logger_debug "[$KMER_FILTER_ABUND_OUTDIR] OK $HASHCOUNT_SUBDIR_PATH directory was created successfully. Will output hashcount output files in this directory."
	# run step
	logger_info "[$HASHCOUNT_SUBDIR] Will run hashcount step ... "
	run_cli -c "$hashcount_cli" -t "$HASHCOUNT_SUBDIR" -e "$HASHCOUNT_ERROR" -d

	### check for potential errors at start
	# khmer version 1.1 / screed version 0.7
	# load-into-counting.py throws "IOError: InvalidFASTQFileFormat: sequence and quality scores length mismatch"
	# Not very informative and not related to fastq format considering this issue (https://github.com/ged-lab/khmer/issues/249)
	# checked input fastq: ok, readlength == qualitylength
	# awk '{if(NR%4==2) print NR"\t"$0"\t"length($0)}' test/01.K-mer_filter_abund/cvi_interleaved.fastq > cvi_readlength.txt
	# awk '{if(NR%4==0) print NR"\t"$0"\t"length($0)}' test/01.K-mer_filter_abund/cvi_interleaved.fastq > cvi_qualitylength.txt
	# awk 'NR==FNR{a[$3]++;next}!a[$3]' cvi_readlength.txt cvi_qualitylength.txt | wc -l # => result: 0
	# Proposed fix: do not use multi-threading with T>1
	# fix works but that's a pitty!
fi

#
# LINK step
#
# link hashcount output and filterabund input
# link sample inputs and filterabund input
# add prefix for output files
filterabund_prefix_out=$(basename ${!hashcount_output['output_countingtable_filename']%.*})
eval "declare -A $(toupper ${NAMESPACE}_filterabund_input)" 
filterabund_input=$(toupper ${NAMESPACE}_filterabund_input)
eval "${filterabund_input}=( [input_presence_table_filename]=${!hashcount_output['output_countingtable_filename']} )"
eval "${filterabund_input}+=([prefix_out]=$filterabund_prefix_out)"
eval "${filterabund_input}+=([input_sequence_filename_R1]=${!current_sample_seq_R1_path})"
eval "${filterabund_input}+=([input_sequence_filename_R2]=${!current_sample_seq_R2_path})"
eval "${filterabund_input}+=([input_sequence_filename_interleaved]=${!interleaving_output["filename"]})"

## test filterabund_input hash
#echo -e "hash_ori: ${!hashcount_output['output_countingtable_filename']}"
#echo -e "R1_ori: ${!current_sample_seq_R1_path}"
#echo -e "R2_ori: ${!current_sample_seq_R2_path}"
#echo -e "prefix_out: $filterabund_prefix_out"
#echo ""
#echo -e "array: ${filterabund_input}"
#hashcount="${filterabund_input}[input_presence_table_filename]"
#echo -e "hash: ${!hashcount}" # works!
#R1="${filterabund_input}[input_sequence_filename_R1]"
#echo -e "R1: ${!R1}"
#R2="${filterabund_input}[input_sequence_filename_R1]"
#echo -e "R2: ${!R2}"
#prefix="${filterabund_input}[prefix_out]"
#echo -e "prefix_out: ${!prefix}"
#vals="${filterabund_input}[@]"
#keys=$(eval echo '${!'$vals'}') # works!
#echo -e "keys: $keys"
#
## hash
#declare -A arr
#arr=([prefix_out]=$filterabund_prefix_out)
#arr+=([input_presence_table_filename]=${!hashcount_output['output_countingtable_filename']})
#
#echo -e "prefix_out: ${arr[prefix_out]}"
#echo -e "input_presence_table_filename: ${arr[input_presence_table_filename]}"
#echo -e "keys: ${!arr[@]}"
#exit 0
## end test

#
# Filtering k-mer abundance
#

logger_info "[$KMER_FILTER_ABUND_OUTDIR] Filtering k-mer abundance ..."

# set filterabund vars
((KMER_FILT_ABUND_COUNTER+=1))
FILTERABUND_SUBDIR="${KMER_FILT_ABUND_COUNTER}.Filterabund"
FILTERABUND_SUBDIR_PATH="$OUTPUT_DIR/$KMER_FILTER_ABUND_OUTDIR/$FILTERABUND_SUBDIR"
FILTERABUND_ERROR="${FILTERABUND_SUBDIR_PATH}/${!current_sample_alias}_filterabund.err"
# define filterabund output
declare -r khmer_filter_abund_C=$(toupper ${NAMESPACE}_khmer_filter_abund)_C
filterabund_input_prefix="${filterabund_input}[prefix_out]"
eval "$(toupper ${NAMESPACE}_filterabund_output)=(['out']=${FILTERABUND_SUBDIR_PATH}/${!filterabund_input_prefix}_C${!khmer_filter_abund_C}.filterabund)"
declare -r filterabund_output=$(toupper ${NAMESPACE}_filterabund_output)
# build cli options
khmer_filterabund_opts=($(buildCommandLineOptions "khmer_filter_abund" "$NAMESPACE" 2>$KMER_FILTER_ABUND_ERROR))
rtrn=$?
cli_opts_failed_msg="[$KMER_FILTER_ABUND_OUTDIR] An error occured while building the khmer_filter_abund command line options for current sample ${!current_sample_alias}."
exit_on_error "$KMER_FILTER_ABUND_ERROR" "$cli_opts_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
filterabund_opts="${khmer_filterabund_opts[@]}"
logger_debug "[$KMER_FILTER_ABUND_OUTDIR] khmer_filter_abund options: $filterabund_opts"
# build cli
declare -r khmer_filter_abund=$(toupper ${NAMESPACE}_paths)_khmer_filter_abund
filterabund_input_hashcount="${filterabund_input}[input_presence_table_filename]"
filterabund_input_R1="${filterabund_input}[input_sequence_filename_R1]"
filterabund_input_R2="${filterabund_input}[input_sequence_filename_R2]"
filterabund_input_interleaved="${filterabund_input}[input_sequence_filename_interleaved]"
filterabund_cli="${!khmer_filter_abund} $filterabund_opts"
filterabund_cli+=" -o ${!filterabund_output['out']}"
filterabund_cli+=" ${!filterabund_input_hashcount}"
filterabund_cli+=" ${!filterabund_input_interleaved}"
#filterabund_cli+=" 2>${FILTERABUND_ERROR} | logger_debug"

# check if filterabund subdir exists else create 
# set filterabund step input vars (see previous LINK step)
# exists: check filterabund step output var => if output file exists skip step 
# not exist: run step 
logger_info "[$KMER_FILTER_ABUND_OUTDIR] Creating $FILTERABUND_SUBDIR_PATH directory ..."
if [[ -d $FILTERABUND_SUBDIR_PATH ]]; then
	logger_debug "[$KMER_FILTER_ABUND_OUTDIR] OK $FILTERABUND_SUBDIR_PATH directory already exists. Will check for already existing output files."
	
	# check for existing output files
	if [[ -s ${!filterabund_output['out']} ]]; then
		# skip step
		logger_info "[$FILTERABUND_SUBDIR] Output file already exists: ${!filterabund_output['out']}."
		logger_info "[$FILTERABUND_SUBDIR] Skip filterabund step."
	else
		# run step
		logger_info "[$FILTERABUND_SUBDIR] Will run filterabund step ... "
		run_cli -c "$filterabund_cli" -t "$FILTERABUND_SUBDIR" -e "$FILTERABUND_ERROR" -d
	fi
else
	# create subdir
	mkdir $FILTERABUND_SUBDIR_PATH 2>$KMER_FILTER_ABUND_ERROR
	rtrn=$?
	out_dir_failed_msg="[$KMER_FILTER_ABUND_OUTDIR] Failed. Filterabund output directory, $FILTERABUND_SUBDIR_PATH, was not created."
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
	exit_on_error "$KMER_FILTER_ABUND_ERROR" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
	# run step
	logger_info "[$FILTERABUND_SUBDIR] Will run filterabund step ... "
	run_cli -c "$filterabund_cli" -t "$FILTERABUND_SUBDIR" -e "$FILTERABUND_ERROR" -d
fi

#
# LINK step
#
# link filterabund output and extractpairs input
eval "declare -A $(toupper ${NAMESPACE}_extractpairs_input)"
declare -r extractpairs_input=$(toupper ${NAMESPACE}_extractpairs_input)
eval "${extractpairs_input}=( [infile]=${!filterabund_output['out']} )"
extractpairs_prefix_out=$(basename ${!filterabund_output['out']})
eval "${extractpairs_input}+=( [prefix_out]=$extractpairs_prefix_out )"

#
# Extracting paired-end reads
#
logger_info "[$KMER_FILTER_ABUND_OUTDIR] Extracting paired-end and orphans reads ... "

# set extractpairs vars
EXTRACTPAIRS_SUBDIR=$FILTERABUND_SUBDIR
EXTRACTPAIRS_SUBDIR_PATH=$FILTERABUND_SUBDIR_PATH
EXTRACTPAIRS_ERROR="$EXTRACTPAIRS_SUBDIR_PATH/${!current_sample_alias}_extractpairs.err"
# define extractpairs output
eval "declare -A $(toupper ${NAMESPACE}_extractpairs_output)"
declare -r extractpairs_output=$(toupper ${NAMESPACE}_extractpairs_output)
eval "${extractpairs_output}=( [pe]=${!filterabund_output['out']}.pe )"
eval "${extractpairs_output}+=( [se]=${!filterabund_output['out']}.se )"

# build cli
declare -r khmer_extract_paired_reads=$(toupper ${NAMESPACE}_paths)_khmer_extract_paired_reads
extractpairs_input_filterabund="${extractpairs_input}[infile]"
extractpairs_cli="cd $EXTRACTPAIRS_SUBDIR_PATH; ${!khmer_extract_paired_reads} $(basename ${!extractpairs_input_filterabund})"
extractpairs_cli+=" 2>$(basename ${EXTRACTPAIRS_ERROR}) | logger_debug"

# no need to check for extractpairs output dir, nor to create it
# set extractpairs step input vars (see previous LINK step)
# check for extractpairs step output vars => if already exist skip step else run step
logger_info "[$KMER_FILTER_ABUND_OUTDIR] Checking for $EXTRACTPAIRS_SUBDIR_PATH directory ..."
if [[ -d $EXTRACTPAIRS_SUBDIR_PATH ]]; then
	logger_debug "[$KMER_FILTER_ABUND_OUTDIR] OK $EXTRACTPAIRS_SUBDIR_PATH directory already exists. Will check for already existing output files."

	# check for existing output files
	filterabund="${extractpairs_input}[infile]"
	filterabund_pe="${!filterabund}.pe"
	filterabund_se="${!filterabund}.se"
	extractpairs_skip=false
	## pe : mandatory
	if [[ -s ${filterabund_pe} ]]; then
		extractpairs_skip=true
		logger_info "[$EXTRACTPAIRS_SUBDIR] Output file already exists: ${filterabund_pe}."
	else
		extractpairs_skip=false
		logger_warn "[$EXTRACTPAIRS_SUBDIR] Output file, ${filterabund_pe}, does not exist. Should run extractpairs step."
	fi
	## se : not mandatory
	if [[ -e ${filterabund_se} ]]; then
		#extractpairs_skip=true
		logger_info "[$EXTRACTPAIRS_SUBDIR] Output file already exists: ${filterabund_se}."
	else
		#extractpairs_skip=false
		logger_warn "[$EXTRACTPAIRS_SUBDIR] Output file, ${filterabund_se}, does not exist. Maybe, should run extractpairs step."
	fi
	
	# run or skip step
	case $extractpairs_skip in
		# skip
		(true) 
			logger_info "[$EXTRACTPAIRS_SUBDIR] Skip extractpairs step."	
			;;
		# run
		(false)
			logger_info "[$EXTRACTPAIRS_SUBDIR] Will run extractpairs step ... "
			run_cli -c "$extractpairs_cli" -t "$EXTRACTPAIRS_SUBDIR" -e "$EXTRACTPAIRS_ERROR" -d 
			cd "$WORKING_DIR"
			;;
	esac
else
	logger_fatal "[$KMER_FILTER_ABUND_OUTDIR] The $EXTRACTPAIRS_SUBDIR_PATH directory does not exist."
	logger_fatal "[$KMER_FILTER_ABUND_OUTDIR] Cannot run the extractpairs step."

	exit_on_error "$KMER_FILTER_ABUND_ERROR" "$EXTRACTPAIRS_SUBDIR_PATH directory does not exist." 1 "" $SESSION_TAG $EMAIL
fi

#
# LINK step
#
# link extractpairs output and splitpairs input
eval "declare -A $(toupper ${NAMESPACE}_splitpairs_input)"
declare -r splitpairs_input=$(toupper ${NAMESPACE}_splitpairs_input)
extractpairs_pe="${extractpairs_output}[pe]"
eval "${splitpairs_input}=( [infile]=${!extractpairs_pe} )"
splitpairs_prefix_out=$(basename ${!extractpairs_pe})
eval "${splitpairs_input}+=( [prefix_out]=$splitpairs_prefix_out )"

#
# Splitting paired reads
#
logger_info "[$KMER_FILTER_ABUND_OUTDIR] Splitting paired reads ... "
logger_info "[$KMER_FILTER_ABUND_OUTDIR] current dir: $(pwd)"

# set splitpairs vars
SPLITPAIRS_SUBDIR=$FILTERABUND_SUBDIR
SPLITPAIRS_SUBDIR_PATH=$FILTERABUND_SUBDIR_PATH
SPLITPAIRS_ERROR="$SPLITPAIRS_SUBDIR_PATH/${!current_sample_alias}_splitpairs.err"
# define splitpairs output
eval "declare -A $(toupper ${NAMESPACE}_splitpairs_output)"
declare -r splitpairs_output=$(toupper ${NAMESPACE}_splitpairs_output)
eval "${splitpairs_output}=( [pe.1]=${!extractpairs_pe}.1 )"
eval "${splitpairs_output}+=( [pe.2]=${!extractpairs_pe}.2 )"

# build cli
declare -r khmer_split_paired_reads=$(toupper ${NAMESPACE}_paths)_khmer_split_paired_reads
splitpairs_input_filterabund_pe="${splitpairs_input}[infile]"
splitpairs_cli="cd $SPLITPAIRS_SUBDIR_PATH; ${!khmer_split_paired_reads} $(basename ${!splitpairs_input_filterabund_pe})"
#splitpairs_cli+=" 2>$(basename $SPLITPAIRS_ERROR) | logger_debug"

# no need to check for splitpairs output dir, nor to create it
# set splitpairs step input vars (see previous LINK step)
# check for splitpairs step output vars => if already exist skip step else run step
logger_info "[$KMER_FILTER_ABUND_OUTDIR] Checking for $SPLITPAIRS_SUBDIR_PATH directory ..."
if [[ -d $SPLITPAIRS_SUBDIR_PATH ]]; then
	logger_debug "[$KMER_FILTER_ABUND_OUTDIR] OK $SPLITPAIRS_SUBDIR_PATH directory already exists. Will check for already existing output files."

	# check for existing output files
	filterabund_pe="${splitpairs_input}[infile]"
	filterabund_pe_1="${!filterabund_pe}.1"
	filterabund_pe_2="${!filterabund_pe}.2"
	splitpairs_skip=false
	## pe.1 && pe.2 : mandatory
	if [[ -s ${filterabund_pe_1} && -s ${filterabund_pe_2} ]]; then
		splitpairs_skip=true
		logger_info "[$SPLITPAIRS_SUBDIR] Output file already exists: ${filterabund_pe_1}"
		logger_info "[$SPLITPAIRS_SUBDIR] Output file already exists: ${filterabund_pe_2}"
	else
		splitpairs_skip=false
		[[ ! -s ${filterabund_pe_1} ]] && logger_warn "[$SPLITPAIRS_SUBDIR] Output file, ${filterabund_pe_1}, does not exist. Should run splitpairs step."
		[[ ! -s ${filterabund_pe_2} ]] && logger_warn "[$SPLITPAIRS_SUBDIR] Output file, ${filterabund_pe_2}, does not exist. Should run splitpairs step."
	fi

	# run or skip step
	case $splitpairs_skip in 
		# skip
		(true)
			logger_info "[$SPLITPAIRS_SUBDIR] Skip splitpairs step."
			;;
		(false)
			logger_info "[$SPLITPAIRS_SUBDIR] Will run splitpairs step ... "
			run_cli -c "$splitpairs_cli" -t "$SPLITPAIRS_SUBDIR" -e "$SPLITPAIRS_ERROR" -d
			cd "$WORKING_DIR"
			;;
	esac
else
	# trouble
	logger_fatal "[$KMER_FILTER_ABUND_OUTDIR] The $SPLITPAIRS_SUBDIR_PATH directory does not exist."
	logger_fatal "[$KMER_FILTER_ABUND_OUTDIR] Cannot run the splitpairs step."

	exit_on_error "$KMER_FILTER_ABUND_ERROR" "$SPLITPAIRS_SUBDIR_PATH does not exist." 1 "" $SESSION_TAG $EMAIL
fi

### close appender
appender_exists kmerFiltAbundF && appender_close kmerFiltAbundF


#========================================
# 02.ASSEMBLY: CONTIGING AND SCAFFOLDING 
#========================================

# STEPS
## CONTIGING
## SCAFFOLDING

#
# Create Assembly output directory
#
ASSEMBLY_OUTDIR="02.Assemblies"
logger_info "Creating $ASSEMBLY_OUTDIR directory ..."
if [[ -d $OUTPUT_DIR/$ASSEMBLY_OUTDIR ]]; then
	logger_debug "OK $ASSEMBLY_OUTDIR directory already exists. Will output all assembly output files in this directory."
else
	mkdir $OUTPUT_DIR/$ASSEMBLY_OUTDIR 2>$ERROR_TMP
	rtrn=$?
	out_dir_failed_msg="[$ASSEMBLY_OUTDIR] Failed. Assembly output directory, $ASSEMBLY_OUTDIR, was not created."
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
	exit_on_error "$ERROR_TMP" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
	logger_debug "[$ASSEMBLY_OUTDIR] OK $ASSEMBLY_OUTDIR directory was created successfully. Will output all output files in this directory."
fi
	
### Enable the assembly debug logger
ASSEMBLY_DEBUGF=${ASSEMBLY_OUTDIR}_debug.log
logger_addAppender assemblyDebugF
appender_setType assemblyDebugF FileAppender
appender_file_setFile assemblyDebugF $(realpath $OUTPUT_DIR)/$ASSEMBLY_OUTDIR/$ASSEMBLY_DEBUGF
appender_setLevel assemblyDebugF DEBUG
appender_setLayout assemblyDebugF PatternLayout
appender_setPattern assemblyDebugF '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions assemblyDebugF
appender_exists assemblyDebugF && logger_info "[$ASSEMBLY_OUTDIR] Debugging infos on assembly will be output to $OUTPUT_DIR/$ASSEMBLY_OUTDIR/$ASSEMBLY_DEBUGF file." || logger_warn "The assemblyDebugF debugger file appender was not enabled. Maybe a log4sh error occured."
### error handling
ASSEMBLY_ERROR=$OUTPUT_DIR/$ASSEMBLY_OUTDIR/${ASSEMBLY_DEBUGF}.err

#
# LINKS step
#
# link kmer filter abund output to assembly input
eval "declare -A $(toupper ${NAMESPACE}_assembly_input)"
assembly_input=$(toupper ${NAMESPACE}_assembly_input)
filterabund_pe="${splitpairs_input}[infile]"
eval "${assembly_input}=( [pe]=${!filterabund_pe} )"
filterabund_pe_1="${splitpairs_output}[pe.1]"
eval "${assembly_input}+=( [pe_1]=${!filterabund_pe_1} )"
filterabund_pe_2="${splitpairs_output}[pe.2]"
eval "${assembly_input}+=( [pe_2]=${!filterabund_pe_2} )"
filterabund_se="${extractpairs_output}[se]"
eval "${assembly_input}+=( [se]=${!filterabund_se} )"
declare -r assembly_k=$(toupper ${NAMESPACE}_contig_assembler)_hash_length
eval "${assembly_input}+=( [k]=${!assembly_k} )"

# 
# ASSEMBLY
#
logger_info "[$ASSEMBLY_OUTDIR] De novo assembly on filtered reads ... "

# set assembly vars
filterabund_pe="${assembly_input}[pe]"
#sample_id="$(basename ${!filterabund_pe%.pe})"
filterabund_se="${assembly_input}[se]"
[[ -s ${!filterabund_se} ]] && sample_id="$(basename ${!filterabund_pe}).se" || sample_id="$(basename ${!filterabund_pe})"
sample_assemblies_outdir="${sample_id}.assemblies"
sample_id+=".assembly_meta-velvetg_k${!assembly_k}"
# define assembly
eval "declare -A $(toupper ${NAMESPACE}_assembly_output)" 
assembly_output=$(toupper ${NAMESPACE}_assembly_output)

#
# CONTIGING
#

declare -r ASSEMBLER=$(toupper ${NAMESPACE}_contig_assembler)_name
logger_info "[$ASSEMBLY_OUTDIR] Use ${!ASSEMBLER} as the contig assembler."

# contiging output directory
contigs_by="contiging_by_${!ASSEMBLER}"

# run assembler
case "${!ASSEMBLER}" in
	"meta-velvetg")
		# set assembler vars
		klen="${assembly_input}[k]"
		## build contiging options
		# use meta_velvetg as written in the config file not meta-velevetg, the actual assembler name, because config parser does not allow hyphen in section name
		contiging_opts=($(buildCommandLineOptions "meta_velvetg" "$NAMESPACE" 2>$ASSEMBLY_ERROR))
		contiging_opts_sorted=($(shortenAndSortOptions "${contiging_opts[@]}" 2>$ASSEMBLY_ERROR))
		contiging_opts_sorted_cat="opts"$(echo "${contiging_opts_sorted[@]}" | sed -e 's/[ =]/_/g')
		## contiging outdir
		CONTIGING_OUTDIR=${OUTPUT_DIR}/${ASSEMBLY_OUTDIR}/${sample_assemblies_outdir}/${contigs_by}/k${!klen}/${contiging_opts_sorted_cat}
		### no meta-velvetg exp_covs defined
		if [[ ! "$contiging_opts_sorted_cat" =~ /_-exp_covs__/ ]]; then
			contiging_opts_sorted_cat_no_mv_exp_covs=$(echo "$contiging_opts_sorted_cat" | sed 's/_-exp_covs_[0-9_]*_/_-exp_covs__/')
			CONTIGING_NO_MV_EXP_COVS_OUTDIR=${OUTPUT_DIR}/${ASSEMBLY_OUTDIR}/${sample_assemblies_outdir}/${contigs_by}/k${!klen}/${contiging_opts_sorted_cat_no_mv_exp_covs}
		fi

		rename_pre_assembly_outdir()
		{
			if [[ -d $CONTIGING_NO_MV_EXP_COVS_OUTDIR ]]; then
				logger_info "[$contigs_by] A previous pre-assembly output directory exists: $CONTIGING_NO_MV_EXP_COVS_OUTDIR"
				logger_info "[$contigs_by] Rename directory: $CONTIGING_NO_MV_EXP_COVS_OUTDIR to $CONTIGING_OUTDIR"
				mv $CONTIGING_NO_MV_EXP_COVS_OUTDIR $CONTIGING_OUTDIR 2>$ERROR_TMP
				rtrn=$?
				failed_rename_dir_msg="Cannot rename directory: $CONTIGING_NO_MV_EXP_COVS_OUTDIR to $CONTIGING_OUTDIR"
				[[ "$rtrn" -ne 0 ]] && logger_fatal "$failed_rename_dir_msg"
				exit_on_error "$ERROR_TMP" "$failed_rename_dir_msg" $rtrn "" $SESSION_TAG $EMAIL
				logger_info "[$contigs_by] Done"
			fi	
		}

		# define expected assembler output
		eval "declare -A $(toupper ${NAMESPACE}_assembler_output)"
		declare -r assembler_output=$(toupper ${NAMESPACE}_assembler_output)
		assembler_contigs=$CONTIGING_OUTDIR/meta-velvetg.contigs.fa
		eval "${assembler_output}=( [contigs]=$assembler_contigs )"

		# build cli
		filterabund_se="${assembly_input}[se]"
		scaffold=no
		run_meta_velvetg_cli="$SCRIPTS_PATH/run_meta-velvetg.sh -o $CONTIGING_OUTDIR -N $NAMESPACE -P ${!filterabund_pe} -S ${!filterabund_se} --scaffold $scaffold --skip_config -d"

		# check if contiging outdir exists 
		# exists: check for expected assembler contigs output => exists skip contiging else run contiging
		# not exist: run contiging
		logger_info "[$contigs_by] Checking for $CONTIGING_OUTDIR directory ..."
		if [[ -d $CONTIGING_OUTDIR ]]; then
			logger_debug "[$contigs_by] OK $CONTIGING_OUTDIR directory already exists. Will check for already existing output files."
	
			# check for existing output files
			contiging_skip=false
			## assembler_contigs : mandatory
			if [[ -s $assembler_contigs ]]; then
				contiging_skip=true
				logger_info "[$contigs_by] Output file already exists: $assembler_contigs"
			else
				contiging_skip=false
				logger_warn "[$contigs_by] Output file, ${assembler_contigs}, does not exist. Should run contiging step."
			fi

			# run or skip step
			case $contiging_skip in
				(true)
					logger_info "[$contigs_by] Skip contiging step."
					;;
				(false)
					logger_info "[$contigs_by] Will run contiging step ..."
					
					# if MV no exp covs dir exist then rename it to $CONTIGING_OUTDIR
					rename_pre_assembly_outdir
					
					#logger_debug "[$contigs_by] meta-velvetg cli: $run_meta_velvetg_cli"
					#run_one_cli -c "$run_meta_velvetg_cli" -t "$contigs_by" -e "$ERROR_TMP"
					$run_meta_velvetg_cli 2>$ERROR_TMP &
					pid=$!
					wait $pid
					rtrn=$?
					pid_file="run_meta-velvetg.sh.pid"
            		[[ -s $CONTIGING_OUTDIR/$pid_file ]] && pid=$(cat $CONTIGING_OUTDIR/$pid_file)
            		status_file="${pid_file%.pid}.exit-status"
            		echo $rtrn >$CONTIGING_OUTDIR/$status_file
            		logger_debug "[$contigs_by] pid, ${pid}, exit status: $rtrn"
					run_MV_failed_msg="[$contigs_by] Meta-velvetg script returns a non-zero status exit code. See $MV_ERROR file for more details."
					exit_on_error "$ERROR_TMP" "$run_MV_failed_msg" "$rtrn" "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
					;;
			esac
		else
			# run contiging
			logger_info "[$contigs_by] Will run contiging step ..."
			
			# if MV no exp covs dir exist then rename it to $CONTIGING_OUTDIR
			rename_pre_assembly_outdir
			
			#logger_debug "[$contigs_by] meta-velvetg cli: $run_meta_velvetg_cli"
			#run_one_cli -c "${run_meta_velvetg_cli}" -t "$contigs_by" -e "$ERROR_TMP"
			$run_meta_velvetg_cli 2>$ERROR_TMP &
			pid=$!
			wait $pid
			rtrn=$?
			pid_file="run_meta-velvetg.sh.pid"
			[[ -s $CONTIGING_OUTDIR/$pid_file ]] && pid=$(cat $CONTIGING_OUTDIR/$pid_file)
			status_file="${pid_file%.pid}.exit-status"
			echo $rtrn >$CONTIGING_OUTDIR/$status_file
			logger_debug "[$contigs_by] pid, ${pid}, exit status: $rtrn"
			run_MV_failed_msg="[$contigs_by] Meta-velvetg script returns a non-zero status exit code. See $OUTPUT_DIR/$LOG_DIR/$DEBUGFILE file for more details."
			exit_on_error "$ERROR_TMP" "$run_MV_failed_msg" "$rtrn" "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
		fi
		;;
	# put here other assemblers
esac


#
# LINKS step
#
# case assembler/scaffolder combinations
# add contiging output to assembly output hash
#eval "declare -A $(toupper ${NAMESPACE}_assembly_output)"
#assembly_output=$(toupper ${NAMESPACE}_assembly_output)
assembler_contigs="${assembler_output}[contigs]"
eval "${assembly_output}=( [assembler_contigs]=${!assembler_contigs} )"
assembly_assembler_contigs="${assembly_output}[assembler_contigs]"
logger_debug "[$contigs_by] Add contiging output to assembly output hash"
logger_debug "[$contigs_by] - ${!assembly_assembler_contigs} => ${assembly_output} hash"


# 
# SCAFFOLDING
#

declare -r SCAFFOLDER=$(toupper ${NAMESPACE}_scaffolder)_name
logger_info "[$ASSEMBLY_OUTDIR] Use ${!SCAFFOLDER} as the scaffolder."

# scaffolding output directory
scaffolds_by="scaffolding_by_${!ASSEMBLER}"

# run scaffolder
case "${!SCAFFOLDER}" in
    "meta-velvetg")
        # set scaffolder vars
        klen="${assembly_input}[k]"
        ## build scaffolding options
        # use meta_velvetg as written in the config file not meta-velevetg, the actual assembler name, because config parser does not allow hyphen in section name
        scaffolding_opts=($(buildCommandLineOptions "meta_velvetg" "$NAMESPACE" 2>$ASSEMBLY_ERROR))
        scaffolding_opts_sorted=($(shortenAndSortOptions "${scaffolding_opts[@]}" 2>$ASSEMBLY_ERROR))
        scaffolding_opts_sorted_cat="opts"$(echo "${scaffolding_opts_sorted[@]}" | sed -e 's/[ =]/_/g')
        ## scaffolding outdir
        SCAFFOLDING_OUTDIR=${OUTPUT_DIR}/${ASSEMBLY_OUTDIR}/${sample_assemblies_outdir}/${scaffolds_by}/k${!klen}/${scaffolding_opts_sorted_cat}
 
		# define expected scaffolder output
        eval "declare -A $(toupper ${NAMESPACE}_scaffolder_output)"
        declare -r scaffolder_output=$(toupper ${NAMESPACE}_scaffolder_output)
        scaffolder_contigs=$SCAFFOLDING_OUTDIR/meta-velvetg.contigs.fa
        eval "${scaffolder_output}=( [contigs]=$scaffolder_contigs )"

        # build cli
        filterabund_se="${assembly_input}[se]"
        scaffold=yes
        run_meta_velvetg_cli="$SCRIPTS_PATH/run_meta-velvetg.sh -o $SCAFFOLDING_OUTDIR -N $NAMESPACE -P ${!filterabund_pe} -S ${!filterabund_se} --scaffold $scaffold --skip_config -d --pre_assembly_dir $CONTIGING_OUTDIR"

		# check if scaffolding outdir exists
        # exists: check for expected scaffolder contigs output => exists skip scaffolding else run scaffolding
        # not exist: run scaffolding 
        logger_info "[$scaffolds_by] Checking for $SCAFFOLDING_OUTDIR directory ..."
        if [[ -d $SCAFFOLDING_OUTDIR ]]; then
            logger_debug "[$scaffolds_by] OK $SCAFFOLDING_OUTDIR directory already exists. Will check for already existing output files."

            # check for existing output files
            scaffolding_skip=false
            ## scaffolder_contigs : mandatory
            if [[ -s $scaffolder_contigs ]]; then
                scaffolding_skip=true
                logger_info "[$scaffolds_by] Output file already exists: $scaffolder_contigs"
            else
                scaffolding_skip=false
                logger_warn "[$scaffolds_by] Output file, ${scaffolder_contigs}, does not exist. Should run scafflding step."
            fi

			# run or skip step
            case $scaffolding_skip in
                (true)
                    logger_info "[$scaffolds_by] Skip scaffolding step."
                    ;;
                (false)
                    logger_info "[$scaffolds_by] Will run scaffolding step ..."

                    $run_meta_velvetg_cli 2>$ERROR_TMP &
                    pid=$!
                    wait $pid
                    rtrn=$?
                    pid_file="run_meta-velvetg.sh.pid"
                    [[ -s $SCAFFOLDING_OUTDIR/$pid_file ]] && pid=$(cat $SCAFFOLDING_OUTDIR/$pid_file)
                    status_file="${pid_file%.pid}.exit-status"
                    echo $rtrn >$SCAFFOLDING_OUTDIR/$status_file
                    logger_debug "[$scaffolds_by] pid, ${pid}, exit status: $rtrn"
                    run_MV_failed_msg="[$scaffolds_by] Meta-velvetg script returns a non-zero status exit code. See $MV_ERROR file for more details."
                    exit_on_error "$ERROR_TMP" "$run_MV_failed_msg" "$rtrn" "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
                    ;;
            esac
		else
			logger_info "[$scaffolds_by] Will run scaffolding step ..."

            $run_meta_velvetg_cli 2>$ERROR_TMP &
            pid=$!
            wait $pid
            rtrn=$?
            pid_file="run_meta-velvetg.sh.pid"
            [[ -s $SCAFFOLDING_OUTDIR/$pid_file ]] && pid=$(cat $SCAFFOLDING_OUTDIR/$pid_file)
            status_file="${pid_file%.pid}.exit-status"
            echo $rtrn >$SCAFFOLDING_OUTDIR/$status_file
            logger_debug "[$scaffolds_by] pid, ${pid}, exit status: $rtrn"
            run_MV_failed_msg="[$scaffolds_by] Meta-velvetg script returns a non-zero status exit code. See $MV_ERROR file for more details."
            exit_on_error "$ERROR_TMP" "$run_MV_failed_msg" "$rtrn" "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
		fi
		;;
esac


#
# LINKS step
#
# add scaffolding output to assembly output hash
scaffolder_contigs="${scaffolder_output}[contigs]"
eval "${assembly_output}+=( [scaffolder_contigs]=${!scaffolder_contigs} )"
assembly_scaffolding_contigs="${assembly_output}[scaffolder_contigs]"
logger_debug "[$scaffolds_by] Add scaffolding output to assembly output hash"
logger_debug "[$scaffolds_by] - ${!assembly_scaffolding_contigs} => ${assembly_output} hash"


### close assembly debug logger
appender_exists assemblyDebugF && appender_close assemblyDebugF


#
# 03.ALIGNMENTS: NUCMER and bwa
#

# STEPS
## NUCMER
## BWA

#
# Create Alignments output directory
#
ALIGNMENTS_OUTDIR="03.Alignments"
logger_info "Creating $ALIGNMENTS_OUTDIR directory ..."
#if [[ -d $OUTPUT_DIR/$ALIGNMENTS_OUTDIR ]]; then
#    logger_debug "OK $ALIGNMENTS_OUTDIR directory already exists. Will output all alignments output files in this directory."
#else
#    mkdir $OUTPUT_DIR/$ALIGNMENTS_OUTDIR 2>$ERROR_TMP
#    rtrn=$?
#    out_dir_failed_msg="[$ALIGNMENTS_OUTDIR] Failed. Alignments output directory, $ALIGNMENTS_OUTDIR, was not created."
#    [[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
#    exit_on_error "$ERROR_TMP" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
#    logger_debug "[$ALIGNMENTS_OUTDIR] OK $ALIGNMENTS_OUTDIR directory was created successfully. Will output all output files in this directory."
#fi
createDir -n "$ALIGNMENTS_OUTDIR" -t "$ALIGNMENTS_OUTDIR" -e "$ERROR_TMP" -d

### Enable the alignments debug logger
ALIGNMENTS_DEBUGF=${ALIGNMENTS_OUTDIR}_debug.log
logger_addAppender alignDebugF
appender_setType alignDebugF FileAppender
appender_file_setFile alignDebugF $(realpath $OUTPUT_DIR)/$ALIGNMENTS_OUTDIR/$ALIGNMENTS_DEBUGF
appender_setLevel alignDebugF DEBUG
appender_setLayout alignDebugF PatternLayout
appender_setPattern alignDebugF '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions alignDebugF
appender_exists alignDebugF && logger_info "[$ALIGNMENTS_OUTDIR] Debugging infos on alignments will be output to $OUTPUT_DIR/$ALIGNMENTS_OUTDIR/$ALIGNMENTS_DEBUGF file." || logger_warn "The alignDebugF debugger file appender was not enabled. Maybe a log4sh error occured."
### error handling
ALIGNMENTS_ERROR=$OUTPUT_DIR/$ALIGNMENTS_OUTDIR/${ALIGNMENTS_DEBUGF}.err

#
# LINKS step
#
# cf assembly_output hash for contigs
# cf sample_id for naming files

#
# NUCMER
#
NUCMER="nucmer"
logger_info "[$ALIGNMENTS_OUTDIR] Align contigs against reference genomes using $NUCMER ... "
NUCMER_OUTDIR=$OUTPUT_DIR/$ALIGNMENTS_OUTDIR/$NUCMER
logger_info "[$ALIGNMENTS_OUTDIR] Will output $NUCMER alignment results into directory $NUCMER_OUTDIR ... "
createDir -n "$NUCMER_OUTDIR" -t "$NUCMER" -e "$ERROR_TMP" -d

function run_nucmer() {
    local REF=$1
    local SEQ=$2
    local OUTDIR=$3
    local TITLE=$4
    local CWD=$(pwd)
    echo "curDir: $CWD"
	cd $(realpath $OUTDIR)
    nucmer --prefix=$TITLE $REF $SEQ
    show-coords ${TITLE}.delta > ${TITLE}.coords
    mummerplot -l ${TITLE}.delta --small --postscript --title="$TITLE"
    show-tiling -c ${TITLE}.delta > ${TITLE}.tiling
    mv "out.ps" $TITLE.ps
    ps2pdf $TITLE.ps $TITLE.pdf
    cd $(realpath $CWD)
}

# vs CONTIGING
assembler_contigs="${assembly_output}[assembler_contigs]"
sample_id+="_scaffolding_no"
sample_dir=$(dirname ${!assembler_contigs##$OUTPUT_DIR/$ASSEMBLY_OUTDIR/})
logger_info "$[${NUCMER}_contiging] run $NUCMER on ${!assembler_contigs} ... "
NUCMER_CONTIGING_OUTDIR=$OUTPUT_DIR/$ALIGNMENTS_OUTDIR/$NUCMER/$sample_dir
createDir -n "$NUCMER_CONTIGING_OUTDIR" -t "${NUCMER}_contiging" -e "$ERROR_TMP" -d

for organite in mito chloro; do
	case $organite in
		mito)
			REF_GENOME=${ga_mito_ref_path}
			;;
		chloro)
			REF_GENOME=${ga_chloro_ref_path}
			;;
	esac
	NUCMER_CONTIGING_ORGANITE=$NUCMER_CONTIGING_OUTDIR/$organite
	NUCMER_CONTIGING_ORGANITE_ERROR=$NUCMER_CONTIGING_ORGANITE/${NUCMER}_${organite}.err
	createDir -n "$NUCMER_CONTIGING_ORGANITE" -t "${NUCMER}_$organite" -e "$ERROR_TMP" -d
	run_nucmer "$REF_GENOME" "$(realpath ${!assembler_contigs})" "$(realpath $NUCMER_CONTIGING_ORGANITE)" "${sample_id}_${organite}" 2>$NUCMER_CONTIGING_ORGANITE_ERROR
	rtrn=$?
	nucmer_failed_msg="[${NUCMER}_$organite] Failed. Errors occured while running $NUCMER on ${!assembler_contigs} against $organite genome"
	exit_on_error "$NUCMER_CONTIGING_ORGANITE_ERROR" "$nucmer_failed_msg" "$rtrn" "$OUTPUT_DIR/$ALIGNMENTS_OUTDIR/$ALIGNMENTS_DEBUGF" $SESSION_TAG $EMAIL
	logger_debug "[${NUCMER}_$organite] PWD: $(pwd)"
	cd $WORKING_DIR
done

## mito
#run_nucmer "$REF_MT" "$(realpath ${!assembler_contigs})" "$(realpath $NUCMER_OUTDIR)" "${sample_id}_mito" 
#echo "PWD: $(pwd)"
#cd $WORKING_DIR

## chloro
#run_nucmer "$REF_PT" "$(realpath ${!assembler_contigs})" "$(realpath $NUCMER_OUTDIR)" "${sample_id}_chloro"
#echo "PWD: $(pwd)"
#cd $WORKING_DIR

# vs SCAFFOLDING
scaffolder_contigs="${assembly_output}[scaffolder_contigs]"
sample_id=${sample_id/scaffolding_no/scaffolding_yes}
sample_dir=$(dirname ${!scaffolder_contigs##$OUTPUT_DIR/$ASSEMBLY_OUTDIR/})
logger_info "$[${NUCMER}_scaffolding] run $NUCMER on ${!scaffolder_contigs} ... "
NUCMER_SCAFFOLDING_OUTDIR=$OUTPUT_DIR/$ALIGNMENTS_OUTDIR/$NUCMER/$sample_dir
createDir -n "$NUCMER_SCAFFOLDING_OUTDIR" -t "$NUCMER_SCAFFOLDING_OUTDIR" -e "$ERROR_TMP" -d

for organite in mito chloro; do
	case $organite in
        mito)
            REF_GENOME=${ga_mito_ref_path}
            ;;
        chloro)
            REF_GENOME=${ga_chloro_ref_path}
            ;;
    esac
	NUCMER_SCAFFOLDING_ORGANITE=$NUCMER_SCAFFOLDING_OUTDIR/$organite
	NUCMER_SCAFFOLDING_ORGANITE_ERROR=$NUCMER_SCAFFOLDING_ORGANITE/${NUCMER}_${organite}.err
	createDir -n "$NUCMER_SCAFFOLDING_ORGANITE" -t "${NUCMER}_$organite" -e "$ERROR_TMP" -d
    run_nucmer "$REF_GENOME" "$(realpath ${!scaffolder_contigs})" "$(realpath $NUCMER_SCAFFOLDING_ORGANITE)" "${sample_id}_${organite}" 2>$NUCMER_SCAFFOLDING_ORGANITE_ERROR
	rtrn=$?
    nucmer_failed_msg="[${NUCMER}_$organite] Failed. Errors occured while running $NUCMER on ${!scaffolder_contigs} against $organite genome"
    exit_on_error "$NUCMER_SCAFFOLDING_ORGANITE_ERROR" "$nucmer_failed_msg" "$rtrn" "$OUTPUT_DIR/$ALIGNMENTS_OUTDIR/$ALIGNMENTS_DEBUGF" $SESSION_TAG $EMAIL
    logger_debug "[${NUCMER}_$organite] PWD: $(pwd)"
    cd $WORKING_DIR	
done

## mito
#run_nucmer "$REF_MT" "$(realpath ${!scaffolder_contigs})" "$(realpath $NUCMER_OUTDIR)" "${sample_id}_mito"
#echo "PWD: $(pwd)"
#cd $WORKING_DIR

## chloro
#run_nucmer "$REF_PT" "$(realpath ${!scaffolder_contigs})" "$(realpath $NUCMER_OUTDIR)" "${sample_id}_chloro"
#echo "PWD: $(pwd)"
#cd $WORKING_DIR

### close alignments debug logger
appender_exists alignDebugF && appender_close alignDebugF


#===============
# 04.Statistics
#===============

# STEPS
## QUAST
## AMOS/ASTATS ### TODO ###
## COMPASS ### TODO ###

#
# Create Statistics output directory
# 
STATS_OUTDIR="04.Statistics"
logger_info "Creating $STATS_OUTDIR directory ..."
createDir -n "$STATS_OUTDIR" -t "$STATS_OUTDIR" -e "$ERROR_TMP" -d

### Enable the alignments debug logger
STATS_DEBUGF=${STATS_OUTDIR}_debug.log
logger_addAppender statsDebugF
appender_setType statsDebugF FileAppender
appender_file_setFile statsDebugF $(realpath $OUTPUT_DIR)/$STATS_OUTDIR/$STATS_DEBUGF
appender_setLevel statsDebugF DEBUG
appender_setLayout statsDebugF PatternLayout
appender_setPattern statsDebugF '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions statsDebugF
appender_exists statsDebugF && logger_info "[$STATS_OUTDIR] Debugging infos on alignments will be output to $OUTPUT_DIR/$STATS_OUTDIR/$STATS_DEBUGF file." || logger_warn "The statsDebugF debugger file appender was not enabled. Maybe a log4sh error occured."
### error handling
STATS_ERROR=$OUTPUT_DIR/$STATS_OUTDIR/${STATS_DEBUGF}.err

#
# QUAST
#
QUAST="quast"
logger_info "[$STATS_OUTDIR] Will compute assemblies statistics using $QUAST ... "

## build quast cli options 
quast_path=$(toupper ${NAMESPACE}_paths)_quast
quast_opts=($(buildCommandLineOptions "quast" "$NAMESPACE" "remove_equal" 2>$ERROR_TMP))
quast_opts_sorted=($(shortenAndSortOptions "${quast_opts[@]}" 2>$ERROR_TMP))
quast_opts_sorted_cat="opts"$(echo "${quast_opts_sorted[@]}" | sed -e 's/[ ]/_/g')
QUAST_OUTDIR=$OUTPUT_DIR/04.Statistics/$QUAST

## create quast output directory
createDir -n "$QUAST_OUTDIR" -t "$QUAST_OUTDIR" -e "$ERROR_TMP" -d 

### CONTIGING
sample_dir=$(dirname ${!assembler_contigs##$OUTPUT_DIR/$ASSEMBLY_OUTDIR/})
logger_info "$[${QUAST}_contiging] run $QUAST on ${!assembler_contigs} ... "
for organite in mito chloro; do
	QUAST_CONTIGING_OUTDIR=$QUAST_OUTDIR/$sample_dir/$organite/${quast_opts_sorted_cat}
	QUAST_ERROR=$QUAST_CONTIGING_OUTDIR/${QUAST}_${organite}.err
	case $organite in
		chloro)
			REF_GENOME=${ga_chloro_ref_path}
			REF_GFF=${ga_chloro_gff_ref_path}
			;;
		mito)
			REF_GENOME=${ga_mito_ref_path}
			REF_GFF=${ga_mito_gff_ref_path}
			;;
	esac	
	createDir -n "$QUAST_CONTIGING_OUTDIR" -t "${QUAST}_$organite" -e "$ERROR_TMP" -d

	quast_cli="${!quast_path} -o $QUAST_CONTIGING_OUTDIR -R $REF_GENOME -G $REF_GFF ${quast_opts[@]} $(realpath ${!assembler_contigs})"
	run_cli -c "$quast_cli" -t "${QUAST}_$organite" -e "$QUAST_ERROR" -d
done

### SCAFFOLDING
sample_dir=$(dirname ${!scaffolder_contigs##$OUTPUT_DIR/$ASSEMBLY_OUTDIR/})
logger_info "$[${QUAST}_scaffolding] run $QUAST on ${!scaffolder_contigs} ... "
for organite in mito chloro; do
    QUAST_SCAFFOLDING_OUTDIR=$QUAST_OUTDIR/$sample_dir/$organite/${quast_opts_sorted_cat}
    QUAST_ERROR=$QUAST_SCAFFOLDING_OUTDIR/${QUAST}_${organite}.err
    case $organite in
        chloro)
            REF_GENOME=${ga_chloro_ref_path}
            REF_GFF=${ga_chloro_gff_ref_path}
            ;;
        mito)
            REF_GENOME=${ga_mito_ref_path}
            REF_GFF=${ga_mito_gff_ref_path}
            ;;
    esac
    createDir -n "$QUAST_SCAFFOLDING_OUTDIR" -t "${QUAST}_$organite" -e "$ERROR_TMP" -d

    quast_cli="${!quast_path} -o $QUAST_SCAFFOLDING_OUTDIR -R $REF_GENOME -G $REF_GFF ${quast_opts[@]} $(realpath ${!scaffolder_contigs})"
    run_cli -c "$quast_cli" -t "${QUAST}_$organite" -e "$QUAST_ERROR" -d
done

### close stats debugger file
appender_exists statsDebugF && appender_close statsDebugF










#=====
# END
#=====

logger_info "[End] Run successfully the pipeline."
logger_info "[End] Will exit now."

# close all appenders
appender_exists stderr && appender_close stderr
appender_exists console && appender_close console
appender_exists debuggerF && appender_close debuggerF

exit 0



























