#!/usr/bin/env bash

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

declare -r VERSION="RC1"

# vars
hashcount=$1
## interleaved reads
reads=$2

# paths
[[ $VERSION == "dev" ]] && PREFIX=$HOME/dev/unix/bash/projects/mitology || PREFIX=/usr/local/src/mitology
SCRIPTS=$PREFIX/share/mitology-pipeline/scripts
R_SCRIPTS=$SCRIPTS/R
KHMER_ABUND_DIST_PATH=/usr/local/bin/abundance-dist.py

# run
${KHMER_ABUND_DIST_PATH} $hashcount $reads ${hashcount}.hist

# plot with R
$R_SCRIPTS/plot_hashcount_hist.R -x 0 -X 1000 -y 0 -Y 3000 -c -d 30 -D 55 -a 25 -A 50 ${hashcount}.hist
