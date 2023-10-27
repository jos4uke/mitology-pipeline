#!/usr/bin/env bash
#
# MITOLOGY PIPELINE FUNCTIONS LIB

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

# date: 2014-10-13

# dependencies: nucmer, ps2pdf 

# VERSION="RC1"

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
