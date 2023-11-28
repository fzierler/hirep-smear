#!/bin/bash

TESTDIR=$(pwd)
TOPDIR=$(dirname $(dirname $TESTDIR))
MAKEDIR=$TOPDIR/Make
MAKEUTILS=$MAKEDIR/Utils
AUTOSUNDIR=$MAKEDIR/Utils/autosun


checkN(){
    N=$1
    #Create symplectic headers
    $MAKEUTILS/write_suN_headers.pl $N REPR_FUNDAMENTAL 0 GAUGE_SPN
    mv suN.h SP.h
    mv suN_types.h SP_types.h
    $MAKEUTILS/write_suN_headers.pl $N REPR_FUNDAMENTAL 0 SUN 
    mv suN.h SU$N.h
    mv suN_types.h SU$N\_types.h
 
    #replace the basename suN with sp for compatibility
    sed -i 's/suN/SP/g' SP.h  SP_types.h
    #replace the basename suN with sp for compatibility
    sed -i 's/SUN_H/SP_H/' SP.h  SP_types.h
    sed -i 's/SUN_TYPES_H/SP_TYPES_H/' SP.h  SP_types.h

    mv SP.h SP$N.h
    mv SP_types.h SP$N\_types.h
}


checkN 4 && checkN 6 && checkN 8  && checkN 10 && checkN 12 && checkN 16 && checkN 20 && checkN 22 

