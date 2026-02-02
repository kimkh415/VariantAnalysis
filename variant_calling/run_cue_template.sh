#! /bin/bash

#$ -cwd

#$ -q broad
#$ -P regevlab
#$ -l h_vmem=8g
#$ -l h_rt=4:00:00
#$ -l os=RedHat7
#$ -pe smp 16
#$ -binding linear:16
#$ -R y

source /broad/software/scripts/useuse
source ~/kwanho/git_repos/cue2/env/bin/activate
export PYTHONPATH=${PYTHONPATH}:~/kwanho/git_repos/cue2

use UGER
use .openssl-1.1.1w
use .gcc-5.2.0
use .python-3.9.2


cue call --config config.yaml

