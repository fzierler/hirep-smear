#!/bin/bash

#Check macros in suN.h. There are several different cases, so do an exhaustive check.
# This test takes ~25 seconds to run on a i5 from 2012 

TESTDIR=$(pwd)
TOPDIR=$(dirname $(dirname $TESTDIR))
MAKEDIR=$TOPDIR/Make
MAKEUTILS=$MAKEDIR/Utils
AUTOSUNDIR=$MAKEDIR/Utils/autosun

create_spn_headers(){
    N=$1
    rep=$2
    echo "Creating SP$N headers ($rep)"
    #Create symplectic headers
    $MAKEUTILS/write_suN_headers.pl $N $rep 0 GAUGE_SPN  || { echo Problem writing symplectic headers ; exit ; }
    mv suN.h SP.h
    mv suN_types.h SP_types.h 
    #replace the basename suN with sp for compatibility
    sed -i 's/suN/SP/g' SP.h  SP_types.h
    sed -i 's/SUN_H/SP_H/' SP.h  SP_types.h
    sed -i 's/SUN_TYPES_H/SP_TYPES_H/' SP.h  SP_types.h
    sed -i 's/NF/SPNF/' SP_types.h

}
    
check_headers(){
    N=$1
    rep=$2
    echo "Checking headers for N=$N ($rep)"
    create_spn_headers $N $rep

    #Write the normal SU(N) header to test against
    $MAKEUTILS/write_suN_headers.pl $N $rep || { echo Problem writing suN headers ; exit ; }
    #compile the test
    gcc -o test_headers -O0 -g -I $TOPDIR/Include/ test_headers.c  || { echo Problem compiling the test ; exit ; }
    # run the test
    ./test_headers || {  echo Found a bug ; exit ; }
}

#Check the representation 

cp $MAKEDIR/MkFlags $MAKEDIR/MkFlags.bu

restore(){
    echo 'Reverting MkFlags...'
    mv $MAKEDIR/MkFlags.bu $MAKEDIR/MkFlags
    touch  $MAKEDIR/MkFlags
}

trap restore INT EXIT TERM

run_spnalgtest(){
    cd $TESTDIR/autosun
    make spnalgtest   || { echo Problem compiling spnalgtest ; exit ; }
    for N in 4 6 8
    do
        echo "Checking algebra for N=$N"
        ./spnalgtest $N   || { echo Problem running spnalgtest ; exit ; }
    done
    cd $TESTDIR
}


check_reps(){

    N=$1
    rep=$2
    echo "Checking representation for N=$N ($rep)"

    sed 's/REPRESENTATION/'${rep}'/' testflags  > $MAKEDIR/MkFlags

    echo Testing $rep at N=$N

    create_spn_headers $N $rep  

    #Write the normal SU(N) header to test against
    $MAKEUTILS/write_suN_headers.pl $N $rep  || { echo Problem writing suN headers ; exit ; }

    #Write the representation specific headers
    cd $TESTDIR/autosun
    make write_spn_sun_algconv   || { echo Problem compiling write_spn_sun_algconv; exit ; }
    ./write_spn_sun_algconv $N   || { echo Problem running write_spn_sun_algconv; exit ; }
    mv spn_sun_algconv.h $TESTDIR
    echo "Building WriteREPR - GAUGE_SUN"
    sed -i 's/GAUGETYPE/GAUGE_SUN/'  $MAKEDIR/MkFlags || exit
    cd $AUTOSUNDIR
    make 
    ./writeREPR $N $TESTDIR/suN_repr_func.h.tmpl > $TESTDIR/suN_func.h || { exit ; } 
    ./writeREPR $N $TESTDIR/suN_exp.c.tmpl > $TESTDIR/suN_exp.h || { exit ; } 
    echo "Building WriteREPR - GAUGE_SPN"
    sed -i 's/GAUGE_SUN/GAUGE_SPN/'  $MAKEDIR/MkFlags
    make 
    ./writeREPR $N $TESTDIR/SP_repr_func.h.tmpl > $TESTDIR/spN_func.h || { exit ; }
    ./writeREPR $N $TESTDIR/SP_exp.c.tmpl> $TESTDIR/spN_exp.h || { exit ; }
    cd $TESTDIR

    sed -i 's/suN/SP/g' spN_exp.h
    #compile and run the test
    gcc -o test_reps -O0 -g -I $TOPDIR/Include/ test_reps.c -D$rep || { echo Problem compiling the representation test ; exit ; }
    echo "Running test"
    ./test_reps || { echo Found a bug ; exit ; }

}

for N in 4 6 8 10 # 12 14 16 18 20 22 24 26 # cut for brevity
do
    check_headers $N REPR_FUNDAMENTAL || exit
done

run_spnalgtest || exit

for rep in REPR_FUNDAMENTAL REPR_ADJOINT REPR_ANTISYMMETRIC
do
    for N in 4 #6 8  
    do
       check_reps $N $rep || exit
    done #rep
done #N


