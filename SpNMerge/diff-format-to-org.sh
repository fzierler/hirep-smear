#!/bin/bash
# Problem: git diff is not syntax aware,
#          so trivial formatting changes clutter its output.
# Solution: use clang-format on all the files before diffing.
# This script clones the repository in /tmp in multiple copies
# checking out different revisions in each copy.
# It then uses the system diff, but feeds it with clang-formatted files.
# Usage example:
# bash ./diff-format-to-org.sh .. > notes-merge-diffs-9e944b8-0af1e75.org
# where:
# 0af1e75 is the commit on master
# 9e944b8 is the merge-base commit 
BASE=$1

setup_dir(){
    local BASE=$1
    local REVISION=$2
    local TMPDIR=/tmp/$REVISION
    if [ ! -d $TMPDIR ]
    then
        git clone $BASE $TMPDIR
    else
        echo "Already cloned $BASE to $TMPDIR" 1>&2
    fi
    (cd $TMPDIR
    git checkout $REVISION 1>&2
    )
}

setup_dirs(){
    local BASE=$1
    local MASTER=$2
    local LASTMERGED=$3
    for revision in $MASTER $LASTMERGED
    do
        setup_dir $BASE $revision 
    done
}

diff2_formatted(){
    local FNAME=$1
    local REVISION1=$2
    local REVISION2=$3
    local R1=/tmp/$REVISION1/$FNAME
    local R2=/tmp/$REVISION2/$FNAME
    if [[ -f $R2 ]] && [[ -f $R1 ]]
    then
        if [[ $FNAME =~ .*\.[ch].* ]]
        then
        echo "$FNAME: formatted" 1>&2
        clang-format -i $R1
        clang-format -i $R2

        else
        echo "$FNAME: not formatted" 1>&2
        fi
        diff -w $R1 $R2
    else
        echo "$FNAME not in both branches:"
        ls $R1
        ls $R2
    fi
}


MERGE_BASE=$(git merge-base spn master)
MASTER=$(git rev-parse master)
echo "Merge Base: $MERGE_BASE" 1>&2
# MAIN
setup_dirs $BASE $MASTER $MERGE_BASE
(echo "* master ($MASTER) to $MERGE_BASE"
for fname in $(git diff --name-only $MASTER $MERGE_BASE \
    | grep -vE '.*.nb' # removing Mathematica notebooks 
)
do 
    echo '** '$fname ; 
    diff2_formatted $fname $MERGE_BASE $MASTER
done) 
