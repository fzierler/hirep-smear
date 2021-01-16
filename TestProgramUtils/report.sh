#!/usr/bin/env bash

BOLD="\e[0;1m"
NORM="\e[0;0m"
BOLDGREEN="\e[1;32m"
BOLDRED="\e[1;31m"

get_options(){
    local FILENAME=$1
    echo -e $BOLD $(head -n 1 $FILENAME | cut -d':' -f2) $NORM
}


check_ok(){
    local FILENAME=$1
    if grep -i fail $FILENAME &> /dev/null
    then
        echo -e $BOLDRED FAIL $NORM
        echo "=> What failed:"
        grep -i fail $FILENAME | grep '\.\.\.'
    else
        echo -e $BOLDGREEN OK $NORM
    fi
}

analyse_file(){
    local FILENAME=$1
    printf "## %s %s\n" "$(get_options $FILENAME )" "$(check_ok $FILENAME)"
}

for FILENAME in out*.txt
do
    analyse_file $FILENAME
done 
