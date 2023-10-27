#!/usr/bin/env bash

#
# run_meta-velvetg.sh
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


#=====
# PID
#=====
trap "exit 1" TERM
pipeline_pid="$$"
echo -e "$pipeline_pid" | tee $OUTPUT_DIR/$(basename ${BASH_SOURCE[0]}).pid | logger_debug


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

#========
# sample
#========
assembly_k="$(toupper ${NAMESPACE}_contig_assembler)_hash_length"
SAMPLE_ID=$(basename ${PAIRED_END})
[[ -s ${SINGLETONS} ]] && SAMPLE_ID+=".se"
SAMPLE_ID+=".assembly_meta-velvetg_k${!assembly_k}"
logger_info "Set the sample id to: $SAMPLE_ID"

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
velvetg_output+=( [weighted_hist]="length-weighted_kmer_cov_hist.rplot_xlim_0_500.pdf" )
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

#
# PRE-ASSEMBLY
#
PA="Pre-assembly"
PA_ERROR=$OUTPUT_DIR/${PA}.err
case $SKIP_PA in
	(true)
		logger_info "[$PA] Skip pre-assembly step ... "
		
		# velveth
		logger_debug "[$PA] Expected velveth output files already exist: "
		logger_debug "[$PA] - ${!VH_Roadmaps}"
		logger_debug "[$PA] - ${!VH_Sequences}"
        if [[ $(realpath $OUTPUT_DIR) != $(realpath $PA_DIR) ]]; then
            logger_debug "[$PA] Will link files to the current output directory"
            ln -s $(realpath ${!VH_Roadmaps}) ${OUTPUT_DIR}/${velveth_output[Roadmaps]}
            ln -s $(realpath ${!VH_Sequences}) ${OUTPUT_DIR}/${velveth_output[Sequences]}
            [[ -e ${OUTPUT_DIR}/${velveth_output[Roadmaps]} ]] && logger_debug "[$PA] ${velveth_output[Roadmaps]} symlinked" || logger_warn "[$PA] ${velveth_output[Roadmaps]} not symlinked"
            [[ -e ${OUTPUT_DIR}/${velveth_output[Sequences]} ]] && logger_debug "[$PA] ${velveth_output[Sequences]} symlinked" || logger_warn "[$PA] ${velveth_output[Sequences]} not symlinked"
        fi

		# velvetg
		logger_debug "[$PA] Expected velvetg output files already exist: "
        logger_debug "[$PA] - ${!VG_Graph2}"
        logger_debug "[$PA] - ${!VG_stats}"
        if [[ $(realpath $OUTPUT_DIR) != $(realpath $PA_DIR) ]]; then
            logger_debug "[$PA] Will link files to the current output directory"
            ln -s $(realpath ${!VG_Graph2}) ${OUTPUT_DIR}/${velvetg_output[Graph2]}
            ln -s $(realpath ${!VG_stats}) ${OUTPUT_DIR}/${velvetg_output[stats]}
            [[ -e ${OUTPUT_DIR}/${velvetg_output[Graph2]} ]] && logger_debug "[$PA] ${velvetg_output[Graph2]} symlinked" || logger_warn "[$PA] ${velvetg_output[Graph2]} not symlinked"
            [[ -e ${OUTPUT_DIR}/${velvetg_output[stats]} ]] && logger_debug "[$PA] ${velvetgh_output[Graph2]} symlinked" || logger_warn "[$PA] ${velvetg_output[stats]} not symlinked"
        fi
		;;
	(false)
		logger_info "[$PA] Running the pre-assembly step ... "
		
		# velveth
		velveth_path="$(toupper ${NAMESPACE}_paths)_velveth"
		VH="${!velveth_path##*/}"
		VH_ERROR=$OUTPUT_DIR/${VH}.err
		case $SKIP_VELVETH in
			(true)
				logger_debug "[$VH] Expected velveth output files already exist: "
				logger_debug "[$VH] - ${!VH_Roadmaps}"
				logger_debug "[$VH] - ${!VH_Sequences}"
				if [[ $(realpath $OUTPUT_DIR) != $(realpath $PA_DIR) ]]; then
					logger_debug "[$VH] Will link files to the current output directory"
					ln -s $(realpath ${!VH_Roadmaps}) ${OUTPUT_DIR}/${velveth_output[Roadmaps]}
					ln -s $(realpath ${!VH_Sequences}) ${OUTPUT_DIR}/${velveth_output[Sequences]}
					[[ -e ${OUTPUT_DIR}/${velveth_output[Roadmaps]} ]] && logger_debug "[$VH] ${velveth_output[Roadmaps]} symlinked" || logger_warn "[$VH] ${velveth_output[Roadmaps]} not symlinked" 
					[[ -e ${OUTPUT_DIR}/${velveth_output[Sequences]} ]] && logger_debug "[$VH] ${velveth_output[Sequences]} symlinked" || logger_warn "[$VH] ${velveth_output[Sequences]} not symlinked"
				fi

				logger_info "[$VH] Skip velveth pre-assembly step ... "
				;;
			(false)
				logger_info "[$VH] Run velveth step ... "
		
				# build cli options
				velveth_opts=($(buildCommandLineOptions "velveth" "$NAMESPACE" "remove_equal" 2>${PA_ERROR}))
				velveth_opts_sorted=($(shortenAndSortOptions "${velveth_opts[@]}" 2>${PA_ERROR}))
				rtrn=$?
				cli_opts_failed_msg="[$VH] An error occured while building $VH command line options for current sample ${SAMPLE_ID}."
				exit_on_error "${PA_ERROR}" "$cli_opts_failed_msg" "$rtrn" "$OUTPUT_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
				logger_debug "[$VH] $VH options: ${velveth_opts_sorted[@]}"
				# build cli
				velveth_format_type="-fastq" 
				velveth_read_type_pe="-shortPaired" 
				velveth_read_type_se="-short"
				velveth_cli="${!velveth_path} $OUTPUT_DIR ${!assembly_k} ${velveth_format_type} ${velveth_read_type_pe} ${PAIRED_END}"
				[[ -s $SINGLETONS ]] && velveth_cli+=" ${velveth_format_type} ${velveth_read_type_se} ${SINGLETONS}"
				velveth_cli+=" ${velveth_opts_sorted[@]}"
				#velveth_cli+=" 2>$VH_ERROR | logger_debug"

				# run cli
				run_cli -c "$velveth_cli" -t "$VH" -e "$VH_ERROR" -d					
				;;	
		esac	

		# velvetg
		velvetg_path="$(toupper ${NAMESPACE}_paths)_velvetg"
		VG="${!velvetg_path##*/}"
		VG_ERROR=$OUTPUT_DIR/${VG}.err
		case $SKIP_VELVETG in
			(true)
				logger_debug "[$VG] Expected velvetg output files already exist: "
				logger_debug "[$VG] - ${!VG_Graph2}"
				logger_debug "[$VG] - ${!VG_stats}"
				if [[ $(realpath $OUTPUT_DIR) != $(realpath $PA_DIR) ]]; then
                    logger_debug "[$VG] Will link files to the current output directory"
                    ln -s $(realpath ${!VG_Graph2}) ${OUTPUT_DIR}/${velvetg_output[Graph2]}
                    ln -s $(realpath ${!VG_stats}) ${OUTPUT_DIR}/${velvetg_output[stats]}
                    [[ -e ${OUTPUT_DIR}/${velvetg_output[Graph2]} ]] && logger_debug "[$VG] ${velvetg_output[Graph2]} symlinked" || logger_warn "[$VG] ${velvetg_output[Graph2]} not symlinked"
                    [[ -e ${OUTPUT_DIR}/${velvetg_output[stats]} ]] && logger_debug "[$VG] ${velvetgh_output[Graph2]} symlinked" || logger_warn "[$VG] ${velvetg_output[stats]} not symlinked"
                fi

				logger_info "[$VG] Skip velvetg pre-assembly step ... "
				;;
			(false)
				logger_debug "[$VG] Run velvetg step ... "

				# build cli options
				velvetg_opts=($(buildCommandLineOptions "velvetg" "$NAMESPACE" "remove_equal" 2>${PA_ERROR}))
				velvetg_opts_sorted=($(shortenAndSortOptions "${velvetg_opts[@]}" 2>${PA_ERROR}))
				rtrn=$?
				cli_opts_failed_msg="[$VG] An error occured while building $VG  command line options for current sample ${SAMPLE_ID}."
				exit_on_error "${PA_ERROR}" "$cli_opts_failed_msg" "$rtrn" "$OUTPUT_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
				logger_debug "[$VG] $VG options: ${velvetg_opts_sorted[@]}"
				# build cli
				velvetg_cli="${!velvetg_path} $OUTPUT_DIR ${velvetg_opts_sorted[@]}"
				#velvetg_cli+=" 2>${VG_ERROR} | logger_debug"
				
				# run cli
				run_cli -c "$velvetg_cli" -t "$VG" -e "$VG_ERROR" -d
				;;
		esac

		# rplot
		RP="rplot"
		RP_ERROR=$OUTPUT_DIR/${RP}.err
		case $SKIP_RPLOT in
			(true)
				logger_debug "[$RP] Expected rplot output file already exist: "
				logger_debug "[$RP] - ${!VG_rplot}"
				logger_info "[$RP] Skip rplot pre-assembly step ... "
				;;
			(false)
				logger_info "[$RP] Run rplot step ... "

				# build cli options
				Rplot_script="length-weigthed_kmer_coverage_hist.rplot.R"
				Rplot_path="$(realpath $(dirname $BASH_SOURCE[0]))/R/$Rplot_script"
				Rplot_stats="${pre_assembly_output}[velvetg_stats]"
				Rplot_cli="Rscript --vanilla $Rplot_path $OUTPUT_DIR/$(basename ${!Rplot_stats})"
				#Rplot_cli+=" 2>${RP_ERROR} | logger_debug"
				#Rplot_cli="$Rplot_path $OUTPUT_DIR/$(basename ${!Rplot_stats}) 2>${RP_ERROR}"
				# run cli
				run_cli -c "$Rplot_cli" -t "$RP" -e "$RP_ERROR" -d
				;;
		esac

		# end pre-assembly
		logger_info "[$PA] Pre-assembly completed."
		;;
esac

#
# ASSEMBLY
#
ASBL="Assembly"
ASBL_ERROR=$OUTPUT_DIR/${ASBL}.err

# if MV_exp_covs not null
# warn the user to check the current exp_covs values with the expected coverage peaks in the rplot output
# else
# exit with the same warning message and add to fill in the exp_covs values in the config file
MV_exp_covs="$(toupper ${NAMESPACE}_meta_velvetg)_exp_covs"
MV_exp_covs_warning="[$ASBL] Current meta-velvetg expected coverage peaks are, $([[ -z ${!MV_exp_covs} ]] && echo "NULL" || echo ${!MV_exp_covs}). \
Please check those values to fit those in the length-weighted kmer coverage histogram, $OUTPUT_DIR/${velvetg_output[weighted_hist]}."
if [[ -n ${!MV_exp_covs} ]]; then
	logger_warn "$MV_exp_covs_warning"
else
	logger_warn "$MV_exp_covs_warning"
	MV_exp_covs_recom="[$ASBL] Then fill in the config file with the expected coverage peaks values, and run again the pipeline to complete the assembly."
	logger_warn "$MV_exp_covs_recom"
	echo -e "$MV_exp_covs_warning" 1>&2
	echo -e "$MV_exp_covs_recom" 1>&2
	#echo 1 >$OUTPUT_DIR/${PA}.exit-status
	#kill -s TERM $pipeline_pid
	exit 1
fi

case $SKIP_MV in
	(true)
		logger_info "[$ASBL] Skip assembly step ... "
		logger_debug "[$ASBL] Expected meta-velvetg output files already exist: "
		logger_debug "[$ASBL] - ${!MV_contigs}"
		[[ ${!amos} -eq "yes" && -s ${!MV_afg} ]] && logger_debug "[$ASBL] - ${!MV_afg}"
		logger_info "[$ASBL] Skip meta-velvetg assembly step ... "
		;;
	(false)
		logger_info "[$ASBL] Running the assembly step ... "
		logger_info "[$ASBL] Run meta-velvetg step ... "

		# meta-velvetg
		meta_velvetg_path="$(toupper ${NAMESPACE}_paths)_metavelvetg"
		MV="${!meta_velvetg_path##*/}"
		MV_ERROR=$OUTPUT_DIR/${MV}.err

		## build cli options
		meta_velvetg_opts=($(buildCommandLineOptions "meta_velvetg" "$NAMESPACE" "remove_equal" 2>${MV_ERROR}))
		meta_velvetg_opts_sorted=($(shortenAndSortOptions "${meta_velvetg_opts[@]}" 2>${MV_ERROR}))
		[[ $SCAFFOLD == "yes" ]] && meta_velvetg_opts_sorted_str=$(echo "${meta_velvetg_opts_sorted[@]}" | sed -e 's/scaffolding no/scaffolding yes/') || meta_velvetg_opts_sorted_str=$(echo "${meta_velvetg_opts_sorted[@]}" | sed -e 's/scaffolding yes/scaffolding no/')
		rtrn=$?
		cli_opts_failed_msg="[$MV] An error occured while building $MV command line options for current sample ${SAMPLE_ID}."
		exit_on_error "${MV_ERROR}" "$cli_opts_failed_msg" "$rtrn" "$OUTPUT_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
		logger_debug "[$MV] $MV options: ${meta_velvetg_opts_sorted_str}"
		## build cli
		meta_velvetg_cli="${!meta_velvetg_path} $OUTPUT_DIR $meta_velvetg_opts_sorted_str"
		#meta_velvetg_cli+=" 2>${MV_ERROR} | logger_debug"
		
		## run cli
		run_cli -c "$meta_velvetg_cli" -t "$MV" -e "$MV_ERROR" -d

		## symlink
		[[ -s ${!MV_contigs} ]] && ln -s $(basename ${!MV_contigs}) $OUTPUT_DIR/${SAMPLE_ID}.contigs.fa || logger_warn "[$MV] Meta-velvetg contigs fasta file, ${!MV_contigs}, does not exist. Cannot symlink ${!MV_contigs} to $OUTPUT_DIR/${SAMPLE_ID}.contigs.fa"
		;;
esac


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

