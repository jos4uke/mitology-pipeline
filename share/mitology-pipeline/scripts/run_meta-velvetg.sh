#!/usr/bin/env bash

#
# run_meta-velvetg.sh
#

# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This script provides a wrapper around meta-velvetg assembler pipeline
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

# Date: 2014-11-02

declare -r VERSION="dev"

########################
# SECTION CONFIGURATION
#######################

### SESSION VARIABLES ###

# DEFAULTS
NAMESPACE_DEFAULT="MITOLOGY"
SCAFFOLD_DEFAULT="no"
CONFIG_SECTIONS=("paths" "contig_assembler" "scaffolder" "velveth" "velvetg" "meta_velvetg")

WORKING_DIR=$(pwd)
DATE=$(date '+%F_%Hh%Mm%Ss')
SESSION_ID=$(date '+%Y%m%d%H%M%S')
EXECUTED_COMMAND="${BASH_SOURCE[0]} $@"
PROG_NAME="$(basename ${BASH_SOURCE[0]})"
SESSION_TAG=${NAMESPACE_DEFAULT}_${USER}_${SESSION_ID}

DEBUGFILE=${SESSION_TAG}.log
ERROR_TMP_MODEL="/tmp/${PROG_NAME%.*}_error_${SESSION_TAG}.XXXXXX"
ERROR_TMP=$(mktemp "$ERROR_TMP_MODEL")
[[ $? -eq 0 ]] && echo -e "Execute ${PROG_NAME%.*} pipeline. Version: $VERSION"| tee -a $ERROR_TMP || exit 1

[[ $VERSION -eq "dev" ]] && PROG_PATH=$(realpath $(dirname ${BASH_SOURCE[0]}));PIPELINE_USER_CONFIG=$(find ${PROG_PATH} -iname "mitology-pipeline_user.config") || PIPELINE_USER_CONFIG=/usr/local/share/mitology-pipeline/etc/mitology-pipeline_user.config

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
appender_exists ${PROG_NAME}.console2 && logger_debug "Console appender is enabled." || logger_warn "Console appender was not enabled. Maybe a log4sh error occured."

### PRINT PROG VERSION AND COMMAND
echo "$PROG_NAME pipeline (version: $VERSION)." | tee $ERROR_TMP 2>&1 | logger_info
echo "Executed command: ${EXECUTED_COMMAND}" | tee -a $ERROR_TMP 2>&1 | logger_info

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

Usage: ${PROG_NAME} -c|--configfile CONFIG_FILE -o|--out_dir OUTPUT_DIR -P|--paired_end PAIRED_END 
                   [-S|--singletons SINGLETONS] [--scaffold no|yes] [-N|--namespace NAMESPACE] 
                   [--pre_assembly_dir PA_DIRECTORY] [--skip_config] 
                   [-d|--debug] [-e|--email_address VALID_EMAIL_ADDR] 
                   [-h|--help]

Mandatory:
-c|--config_file CONFIG_FILE            The user configuration file, CONFIG_FILE, listing the pipeline parameters by section.
                                        You can get a copy there: $PIPELINE_USER_CONFIG.
                                        Only the sections (and following parameters) listed here are mandatory:
                                        $(echo "${CONFIG_SECTIONS[@]}" | sed -e 's/ /, /g')
-o|--out_dir OUTPUT_DIR                 The output directory.
-P|--paired_end PAIRED_END              The interleaved paired end sequences.

Options:
-N|--namespace NAMESPACE                The namespace to use for the pipeline section/parameters.
-S|--singletons SINGLETONS              The singletons sequences file.
--scaffold no|yes                       Disable or enable the scaffolding process.
--pre_assembly_dir PA_DIRECTORY         The pre-assembly directory, PA_DIRECTORY, containing the velveth/velvetg output files 
                                        mandatory to meta-velvetg. Useful to avoid redundancy when run a 2 rounds assembly
                                        first without scaffolding then with scaffolding. 
--skip_config                           Will skip loading the mandatory config file.
                                        Useful when the config was pre-loaded by a caller script.
-d|--debug                              Enable debugging mode in the console.
-e|--email_address VALID_EMAIL_ADDR     An optional but valid email address to send pipeline job/error status notifications
-h|--help                               Displays this message.

"
}

# Configure opts
### NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately;
CFG_OPTS=$(getopt -o hc:N:P:S:o:e:d --long help,config_file:,namespace:,paired_end:,singletons:,out_dir:,email_address:,debug,scaffold:,skip_config,pre_assembly_dir: -n 'run_meta-velvetg.sh' -- "$@")

if [[ $? != 0 ]] ; then Usage >&2 ; exit 1 ; fi

# Note the quotes around `$CONFIGURE_OPTS'
eval set -- "$CFG_OPTS"

while true; do
    case "$1" in
        -h | --help ) Usage >&2; exit 1;;
        -c | --config_file ) CONFIGFILE="$2"; shift 2 ;;
        -o | --out_dir ) OUTPUT_DIR="$2"; shift 2 ;;
        -N | --namespace ) NAMESPACE="$2"; shift 2 ;;
        -P | --paired_end ) PAIRED_END="$2"; shift 2 ;;
        -S | --singletons ) SINGLETONS="$2"; shift 2;;
        -d | --debug )
					DEBUG="debug";
					appender_setLevel ${PROG_NAME}.console2 DEBUG 
					appender_activateOptions ${PROG_NAME}.console2
					shift ;;
        --scaffold ) SCAFFOLD="$2"; shift 2;;
        --pre_assembly_dir ) PA_DIR="$2"; shift 2;;
		--skip_config ) SKIP_CONFIG=true; shift;;
		-e | --email_address ) EMAIL="$2"; shift 2 ;;
        -- ) shift; break ;;
        * ) break ;;
    esac
done

### VALIDATION ###

# mandatory
if [[ ! -s $CONFIGFILE ]]; then
    case $SKIP_CONFIG in
		(true)
			logger_info "Set skipping loading config file, $([[ -n $CONFIGFILE ]] && echo $CONFIGFILE)."
			;;
		(false)	
			logger_fatal "Config file, $CONFIGFILE, does not exist or is empty. See Usage with --help option.";
    		exit 1;
			;;
	esac
else
	case $SKIP_CONFIG in
		(true)
			logger_info "Set skipping loading config file, $CONFIGFILE."
			;;
	esac
fi

if [[ -z $OUTPUT_DIR ]]; then
    logger_fatal "Output directory must be not null. See Usage with --help option.";
    exit 1;
fi

if [[ ! -s $PAIRED_END ]]; then
	logger_fatal "Interleaved paired end file, $PAIRED_END, does not exist or is empty. See Usage with --help option."
	exit 1;
fi

# not mandatory
if [[ -z $NAMESPACE ]]; then
	NAMESPACE=$NAMESPACE_DEFAULT
	logger_info "Use default namespace: $NAMESPACE"
fi

if [[ ! -s $SINGLETONS ]]; then
	logger_warn "Singletons file, $SINGLETONS, does not exist or is empty. See Usage with --help option."
fi

if [[ "$SCAFFOLD" -eq "no" || "$SCAFFOLD" -eq "yes" ]]; then
	logger_info "The user scaffold value is now, ${SCAFFOLD}."
else
	logger_warn "The user scaffold value, $SCAFFOLD, is not recognized. The authorized values are no or yes."
	SCAFFOLD=$SCAFFOLD_DEFAULT
	logger_warn "Will use the default scaffold value, ${SCAFFOLD}."
fi

if [[ ! -s "$PA_DIR" ]]; then
	if [[ -z "$PA_DIR" ]]; then
		logger_info "Set the pre-assembly directory to the output directory: $OUTPUT_DIR"	
	else
		logger_warn "The user pre-assembly directory, ${PA_DIR}, does not exist."
		logger_info "Set the pre-assembly directory to the output directory: $OUTPUT_DIR"
	fi
	PA_DIR=$OUTPUT_DIR
else
	logger_debug "Set the pre-assembly directory to ${PA_DIR}"
	logger_debug "Velveth/velvetg output files are supposed to be in this directory."
fi
 
### LOAD LIB ###

# bash-common lib
[[ $VERSION == "dev" ]] && LIB_PATH=$(realpath $(dirname ${BASH_SOURCE[0]}))/../../../../bash-common/share/bash-common/lib/bash-common_lib.inc || LIB_PATH=/usr/local/share/bash-common/lib/bash-common_lib.inc

logger_debug "[Library] Loading $LIB_PATH"
. $LIB_PATH
if [[ $? -ne 0 ]]; then
    logger_fatal "Error loading bash common lib: $LIB_PATH"
    exit 1
fi

# mitology-pipeline lib
[[ $VERSION == "dev" ]] && LIB_PATH=$(realpath $(dirname $0))/../lib/mitology-pipeline_lib.inc || LIB_PATH=/usr/local/share/mitology-pipeline/lib/mitology-pipeline_lib.inc

logger_debug "[Library] Loading $LIB_PATH"
. $LIB_PATH
if [[ $? -ne 0 ]]; then
    logger_fatal "Error loading mitology pipeline lib: $LIB_PATH"
    exit 1
fi


#################
# PIPELINE STEPS
#################

# CREATE OUTPUT DIR
# LOAD CONFIG
# ASSEMBLY:
# - CONTIGING
# - SCAFFOLDING

#===================
# OUTPUT DIRECTORY
#===================

echo "$PROG_NAME pipeline (version: $VERSION)." | tee $ERROR_TMP 2>&1 | logger_info
echo "Executed command: ${EXECUTED_COMMAND}" | tee -a $ERROR_TMP 2>&1 | logger_info

#
# Create a directory named with OUTPUT_DIR value, to save all outputs
#
echo "Creating $OUTPUT_DIR directory ..." | tee -a $ERROR_TMP 2>&1 | logger_info
if [[ -d $OUTPUT_DIR ]]; then
    echo "OK $OUTPUT_DIR directory already exists. Will output all assembly output files in this directory." | tee -a $ERROR_TMP 2>&1  | logger_info
else
    mkdir -p $OUTPUT_DIR 2>>$ERROR_TMP
    rtrn=$?
    out_dir_failed_msg="[Output directory] Failed. Output directory, $OUTPUT_DIR, was not created."
    [[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
    exit_on_error "$ERROR_TMP" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
    echo "$(date '+%Y-%m-%d %T') [Output directory] OK $OUTPUT_DIR directory was created successfully. Will output all output files in this directory." | tee -a $ERROR_TMP 2>&1 | logger_info
fi

#==================================
# Enable the pipeline debug logger
#==================================

DEBUGFILE="${PROG_NAME}_$(date '+%F_%Hh%Mm%Ss').log"
DEBUGFILE_PATH=$(realpath $OUTPUT_DIR)/$DEBUGFILE
logger_addAppender ${PROG_NAME}.debuggerF2
appender_setType ${PROG_NAME}.debuggerF2 FileAppender
appender_file_setFile ${PROG_NAME}.debuggerF2 $DEBUGFILE_PATH
appender_setLevel ${PROG_NAME}.debuggerF2 DEBUG
appender_setLayout ${PROG_NAME}.debuggerF2 PatternLayout
appender_setPattern ${PROG_NAME}.debuggerF2 '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions ${PROG_NAME}.debuggerF2
appender_exists ${PROG_NAME}.debuggerF2 && cat $ERROR_TMP | logger_info
appender_exists ${PROG_NAME}.debuggerF2 && logger_info "Debugging infos will be output to $DEBUGFILE_PATH file." || logger_warn "The debugger file appender was not enabled. Maybe a log4sh error occured."

#=============
# LOAD CONFIG
#=============

case $SKIP_CONFIG in
	(false)
		# set backup config file variable
		BACKUPED_CONFIG_FILE=$OUTPUT_DIR/$(basename $CONFIGFILE)

		# 1. Backup session user config file in session dir if not exist
		logger_info "[Check config: session user config file] Backuping session user config file into session directory ..."
		cp $CONFIGFILE $OUTPUT_DIR/. 2>$ERROR_TMP
		rtrn=$?
		cp_user_config_failed_msg="[Check config: session user config file] Failed backuping session user config file into session directory."
		[[ "$rtrn" -ne 0 ]] && logger_fatal "$cp_user_config_failed_msg"
		exit_on_error "$ERROR_TMP" "$cp_user_config_failed_msg" $rtrn "$DEBUGFILE_PATH" $SESSION_TAG $EMAIL
		logger_info "[Check config: session user config file] Will use backuped session user config file: $BACKUPED_CONFIG_FILE" 2>&1

		# 2. Load config parameters from backuped session user config file
		logger_info "[Check config: session user config file] Loading session user config parameters from $BACKUPED_CONFIG_FILE file ..."
		load_user_config_failed_msg="[Check config: session user config file] Failed loading session user config parameters from $BACKUPED_CONFIG_FILE file."
		for cfg in $(get_config_sections $BACKUPED_CONFIG_FILE 2>$ERROR_TMP;); do
			rtrn=$?
			[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
			exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$DEBUGFILE_PATH" $SESSION_TAG $EMAIL
			logger_debug "--- Config section [${cfg}] ---"
			unset $(set |awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN {
				  cfg = toupper(cfg);
				  prefix = toupper(prefix);
				  pattern = "\^" prefix "_" cfg "_"
			   }
			   $0~pattern { print $1 }' 2>$ERROR_TMP ) 2>>$ERROR_TMP
			rtrn=$?
			[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
			exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$DEBUGFILE_PATH" $SESSION_TAG $EMAIL
			CONFIG_PARAMS=$(format_config_params $BACKUPED_CONFIG_FILE ${cfg} ${NAMESPACE} 2>$ERROR_TMP)
			eval "${CONFIG_PARAMS}"
			rtrn=$?
			[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
			exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$DEBUGFILE_PATH" $SESSION_TAG $EMAIL
			for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
				logger_debug "$params"
			done
			rtrn=$?
			[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
			exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$DEBUGFILE_PATH" $SESSION_TAG $EMAIL
		done
		logger_info "[Check config: session user config file] OK Session user config file, $BACKUPED_CONFIG_FILE, was loaded successfully."
		;;
	(true)
		logger_info "Skip loading config file. Will use given namespace, $NAMESPACE, to recover pipeline parameters."
		# print pipeline parameters
		params_counter=0
		for cfg in "${CONFIG_SECTIONS[@]}"; do
			logger_debug "--- Config section [${cfg}] ---"
			for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
				logger_debug "$params"
				((params_counter+=1))
			done
		done
		if [[ $params_counter == 0 ]]; then
			logger_warn "[Check config: inherit env variables] No pipeline parameters were inherited from the current session."
		else
			 logger_info "[Check config: inherit env variables] OK $params_counter pipeline parameters were inherited from the current session."
		fi
		;;
esac		


#=======================
# CONTIGING/SCAFFOLDING
#=======================

### STEPS ###
# PRE-ASSEMBLY
## VELVETH
## VELVETG
## RPLOT
# ASSEMBLY
## META-VELVETG
### SCAFFOLD NO|YES

### expected output ###
# if not passed as options to the cli, search in the current output dir
# PRE-ASSEMBLY
## VELVETH
# velveth/Roadmaps
# velveth/Sequences 

## VELVETG
# velvetg/Graph2
# velvetg/stats.txt

## RPLOT
# *.rplot.pdf

# ASSEMBLY
## META-VELVETG
### SCAFFOLD NO|YES
# meta-velvetg.contigs.fa
# meta-velvetg.asm.afg

# assembly vars 
ASSEMBLY_STEP="assembly"

# expected output hashes
declare -A velveth_output
velveth_output=( [Roadmaps]="Roadmaps" )
velveth_output+=( [Sequences]="Sequences" )
declare -A velvetg_output
velvetg_output=( [Graph2]="Graph2" )
velvetg_output+=( [stats]="stats.txt" )
velvetg_output+=( [weighted_hist]="weighted_hist.rplot.pdf" )
declare -A metavelvetg_output 
metavelvetg_output=( [contigs]="meta-velvetg.contigs.fa" )
metavelvetg_output+=( [afg]="meta-velvetg.asm.afg" )

eval "declare -A $(toupper ${NAMESPACE}_pre_assembly_output)"
pre_assembly_output=$(toupper ${NAMESPACE}_pre_assembly_output)
eval "${pre_assembly_output}=( [velveth_Roadmaps]=${PA_DIR}/${velveth_output[Roadmaps]} )"
eval "${pre_assembly_output}+=( [velveth_Sequences]=${PA_DIR}/${velveth_output[Sequences]} )"
eval "${pre_assembly_output}+=( [velvetg_Graph2]=${PA_DIR}/${velvetg_output[Graph2]} )"
eval "${pre_assembly_output}+=( [velvetg_stats]=${PA_DIR}/${velvetg_output[stats]} )"
eval "${pre_assembly_output}+=( [rplot]=${PA_DIR}/${velvetg_output[weighted_hist]} )" 

eval "declare -A $(toupper ${NAMESPACE}_assembly_output)"
assembly_output=$(toupper ${NAMESPACE}_assembly_output)
eval "${assembly_output}=( [metavelvetg_contigs]=${OUTPUT_DIR}/${metavelvetg_output[contigs]} )"
eval "${assembly_output}+=( [metavelvetg_afg]=${OUTPUT_DIR}/${metavelvetg_output[afg]} )"

# check for existing pre-assembly output files
## velveth
VH_Roadmaps="${pre_assembly_output}[velveth_Roadmaps]"
VH_Sequences="${pre_assembly_output}[velveth_Sequences]"
[[ ! -s ${!VH_Roadmaps} ]] && SKIP_VELVETH=false || SKIP_VELVETH=true
[[ ! -s ${!VH_Sequences} ]] && SKIP_VELVETH=false || SKIP_VELVETH=true
## velvetg
VG_Graph2="${pre_assembly_output}[velvetg_Graph2]"
VG_stats="${pre_assembly_output}[velvetg_stats]"
[[ ! -s ${!VG_Graph2} ]] && SKIP_VELVETG=false || SKIP_VELVETG=true
[[ ! -s ${!VG_stats} ]] && SKIP_VELVETG=false || SKIP_VELVETG=true
## rplot
VG_rplot="${pre_assembly_output}[rplot]"
[[ ! -s ${!VG_rplot} ]] && SKIP_RPLOT=false || SKIP_RPLOT=true
## pre-assembly
logger_debug "[$ASSEMBLY_STEP] SKIP_VELVETH: $SKIP_VELVETH"
logger_debug "[$ASSEMBLY_STEP] SKIP_VELVETG: $SKIP_VELVETG"
logger_debug "[$ASSEMBLY_STEP] SKIP_RPLOT: $SKIP_RPLOT"
[[ $SKIP_VELVETH == true && $SKIP_VELVETG == true && $SKIP_RPLOT == true ]] && SKIP_PA=true || SKIP_PA=false
logger_debug "[$ASSEMBLY_STEP] SKIP_PRE_ASSEMBLY: $SKIP_PA"

# check for existing assembly output files
## metavelvetg
MV_contigs="${assembly_output}[metavelvetg_contigs]"
MV_afg="${assembly_output}[metavelvetg_afg]"
[[ ! -s ${!MV_contigs} ]] && SKIP_MV=false || SKIP_MV=true
amos=$(toupper ${NAMESPACE}_meta_velvetg)_amos_file
[[ ${!amos} -eq "yes" && ! -s ${!MV_afg} ]] && SKIP_MV=false || SKIP_MV=true
logger_debug "[$ASSEMBLY_STEP] SKIP_METAVELVETG: $SKIP_MV"














#=====
# END
#=====

logger_info "[End] Run successfully the $PROG_NAME pipeline."
logger_info "[End] Will exit now."

# close all appenders
appender_exists ${PROG_NAME}.stderr2 && appender_close ${PROG_NAME}.stderr2
appender_exists ${PROG_NAME}.console2 && appender_close ${PROG_NAME}.console2
appender_exists ${PROG_NAME}.debuggerF2 && appender_close ${PROG_NAME}.debuggerF2

exit 0

