#!/usr/bin/env bash

# author: Joseph Tran <Joseph.Tran@versailles.inra.fr>
# license: GPL >= v3 (see mitology-pipeline/LICENSE)
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
