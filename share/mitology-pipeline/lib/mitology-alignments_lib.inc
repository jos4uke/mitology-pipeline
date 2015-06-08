#!/usr/bin/env bash

#
# Alignments 
#

#-----------------------------------
# Function: run_nucmer 
# Purpose: run nucmer command
function run_nucmer() 
{
	# usage
	run_cli_usage() {
		printf %s "\
		Usage: run_nucmer [-h|--help] -r|--reference_genome REF_GENOME -t|--tag TAG -e|--cmd_err COMMAND_ERR_FILE [-d|--debug]
	"
	}

    local REF=$1
    local SEQ=$2
    local OUTDIR=$3
    local TITLE=$4
    local PWD=$(pwd)
    echo "curDir: $PWD"
    cd $(realpath $OUTDIR)
    nucmer --prefix=$TITLE $REF $SEQ
    show-coords ${TITLE}.delta > ${TITLE}.coords
    mummerplot -l ${TITLE}.delta --small --postscript --title="$TITLE"
    show-tiling -c ${TITLE}.delta > ${TITLE}.tiling
    mv "out.ps" $TITLE.ps
    ps2pdf $TITLE.ps $TITLE.pdf
    cd $(realpath $PWD)
}
