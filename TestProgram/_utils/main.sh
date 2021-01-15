#!/usr/bin/env bash
HIREP_REPO=$1
MASTER_TEST_DIR=$2
HEADER=$3

_UTILS=$(dirname ${BASH_SOURCE[0]})

set -eu

source $_UTILS/create_runs.sh

header_check(){
    if [ -z $HEADER ]
    then
        HEADER=$_UTILS/script_header_sunbird.sh
    fi
    if [ ! -f $HEADER ]
    then
        echo "ERROR: $HEADER does not exist."
    fi
}

header_check

NTESTS=$(create_run_scripts $MASTER_TEST_DIR $HEADER) &&
    create_numbered_hirep_repo_clones $HIREP_REPO $MASTER_TEST_DIR $NTESTS





