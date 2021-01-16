#!/usr/bin/env bash

get_options(){
    local FILENAME=$1
    head -n 1 $FILENAME | cut -d':' -f2
}

check_ok(){
    local FILENAME=$1
    if grep -i fail $FILENAME &> /dev/null
    then
        echo FAIL
        grep -i fail $FILENAME
    else
        echo OK
    fi
}

analyse_file(){
    local FILENAME=$1
    printf "%s %s\n" "$(get_options $FILENAME )" "$(check_ok $FILENAME)"
}

for FILENAME in out*.txt
do
    analyse_file $FILENAME
done 
