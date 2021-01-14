#!/usr/bin/env bash

# These functions create list of tests to execute,
# to be passed to TestProgram/run_tests.sh
# One group of options per line
# 
sun_run_option_list(){
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

spn_run_option_list(){
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

all_run_option_list(){
    sun_run_option_list
    spn_run_option_list
}
