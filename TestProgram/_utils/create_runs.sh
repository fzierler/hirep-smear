#!/usr/bin/bash


_UTILS=$(dirname ${BASH_SOURCE[0]})
source $_UTILS/test_build_lists.sh

create_run_scripts(){
    # Generates a script for each 
    MASTER_TEST_DIR=$1
    HEADER=$2 # e.g., script_header_sunbird.sh
    I=0
    mkdir -p $MASTER_TEST_DIR
    all_run_option_list > .runs
    while IFS='' read -r LINE
    do
        I=$((I+1))
        (
            cat $HEADER
            echo cd $(pwd)/$MASTER_TEST_DIR/$I/TestProgram
            echo echo run_test.sh arguments: $LINE
            echo echo test script: \$0
            echo ./run_tests.sh $LINE
        ) > $MASTER_TEST_DIR/script$I.sh
    done < .runs
    N=$(wc -l < .runs)
    echo "Created $N test scripts." 1>&2
    echo $N # This can be fed to create_numbered_hirep_repo_clones
}

create_numbered_hirep_repo_clones(){
    # Copies s 
    HIREP_REPO=$1
    MASTER_TEST_DIR=$2
    N=$3
    mkdir -p $MASTER_TEST_DIR
    for I in $(seq $N)
    do 
        DIR=$MASTER_TEST_DIR/$I
        git clone --depth=1 $HIREP_REPO $DIR
    done
}



