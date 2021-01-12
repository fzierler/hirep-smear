#!/usr/bin/env bash
set -euo pipefail

run_test_suite(){
    inputs_ccache=$1
    inputs_mpi=$2
    inputs_ecsw=$3
    inputs_nc=$4
    inputs_repr=$5
    inputs_gaugegroup=$6
    
    ./TestProgram/run_tests.sh \
        ${inputs_ccache} \
        ${inputs_mpi} \
        ${inputs_ecsw} \
        -n \
        ${inputs_nc} \
        -r \
        ${inputs_repr} \
        --gauge \
        ${inputs_gaugegroup} \

}

run_sun_config(){
    # creates list of tests to execute
    for mpi in -no-mpi -mpi
    do 
        for matrix_nc in 2 3
        do
            for matrix_repr in FUND ADJ
            do
                for matrix_ecsw in expclover -no-expclover
                do
                    echo \
                        -no-ccache\
                        ${mpi} \
                        ${matrix_ecsw} \
                        ${matrix_nc} \
                        ${matrix_repr} \
                        SUN
                done
            done
        done
    done
}

run_spn_config(){
    # creates list of tests to execute
    for mpi in -no-mpi -mpi
    do 
        for matrix_nc in 4 6
        do
            for matrix_repr in FUND 2A 
            do
                for matrix_ecsw in -no-clover
                do
                    echo \
                        -no-ccache\
                        ${mpi} \
                        ${matrix_ecsw} \
                        ${matrix_nc} \
                        ${matrix_repr} \
                        SPN
                done
            done
        done
    done
}

all_runs(){
    run_sun_config
    run_spn_config
}


create_clones(){
    HIREP_BASE=$1
    TEST_DIR_BASE=$2
    mkdir -p $TEST_DIR_BASE
    I=0
    all_runs | while IFS='' read -r LINE
    do
        DIR=$TEST_DIR_BASE/$I
        git clone --depth=1 $HIREP_BASE $DIR
        I=$((I+1))
    done
}

create_run_scripts(){
    TEST_DIR_BASE=$1
    HEADER=$2 # run_tests_header.sh
    I=0
    all_runs | while IFS='' read -r LINE
    do
        (
            cat $HEADER
            echo cd $(pwd)/$TEST_DIR_BASE/$I/TestProgram
            echo ./run_tests.sh $LINE
        ) > $TEST_DIR_BASE/script$I.sh
        I=$((I+1))
    done
}


