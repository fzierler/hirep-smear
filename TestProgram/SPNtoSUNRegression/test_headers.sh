#!/bin/bash

#Check macros in suN.h. There are several different cases, so do an exhaustive check.
TESTDIR=$(pwd)
TOPDIR=$(dirname $(dirname $TESTDIR))
MAKEDIR=$TOPDIR/Make
MAKEUTILS=$MAKEDIR/Utils
AUTOSUNDIR=$MAKEDIR/Utils/autosun

for N in 2 4 6 8 10 12 14 16 18 20 22 24 26 ; do

  echo Testing N=$N
  
  #Create symplectic headers
  $MAKEUTILS/write_suN_headers.pl $N REPR_FUNDAMENTAL 0 GAUGE_SPN  || { echo Problem writing symplectic headers ; exit ; }
  mv suN.h SP.h
  mv suN_types.h SP_types.h 
  #replace the basename suN with sp for compatibility
  sed -i 's/suN/SP/g' SP.h  SP_types.h
  #replace the basename suN with sp for compatibility
  sed -i 's/SUN_H/SP_H/' SP.h  SP_types.h
  sed -i 's/SUN_TYPES_H/SP_TYPES_H/' SP.h  SP_types.h

  #Write the normal SU(N) header to test against
  $MAKEUTILS/write_suN_headers.pl $N REPR_FUNDAMENTAL || { echo Problem writing suN headers ; exit ; }

  #compile and run the test
  gcc -o test_headers -O0 -g -I $TOPDIR/Include/ test_headers.c  || { echo Problem compiling the test ; exit ; }
  ./test_headers || {  echo Found a bug ; exit ; }

done


#Check the representation 

cp $MAKEDIR/MkFlags $MAKEDIR/MkFlags.bu

restore(){
    echo 'Reverting MkFlags...'
    mv $MAKEDIR/MkFlags.bu $MAKEDIR/MkFlags
    touch  $MAKEDIR/MkFlags
}

trap restore INT EXIT TERM

for N in 2 4 6 8  ; do
  cd $TESTDIR/autosun
    make spnalgtest   || { echo Problem compiling spnalgtest ; exit ; }
    ./spnalgtest $N   || { echo Problem running spnalgtest ; exit ; }
  cd $TESTDIR
  for rep in REPR_FUNDAMENTAL REPR_ADJOINT; do
  
    sed 's/REPRESENTATION/'${rep}'/' testflags  > $MAKEDIR/MkFlags
  
    echo 
    echo Testing $rep at N=$N
    echo 
    
    #Create symplectic headers
    $MAKEUTILS/write_suN_headers.pl $N  $rep 0  GAUGE_SPN  || { echo Problem writing symplectic headers ; exit ; }
    mv suN.h SP.h
    mv suN_types.h SP_types.h 
    #replace the basename suN with sp for compatibility
    sed -i 's/suN/SP/g' SP.h  SP_types.h
    sed -i 's/SUN_H/SP_H/' SP.h  SP_types.h
    sed -i 's/SUN_TYPES_H/SP_TYPES_H/' SP.h  SP_types.h
    sed -i 's/NF/SPNF/' SP_types.h
    
    #Write the normal SU(N) header to test against
    $MAKEUTILS/write_suN_headers.pl $N $rep  || { echo Problem writing suN headers ; exit ; }
    
    #Write the representation specific headers
    cd $TESTDIR/autosun
      make spnalgtest   || { echo Problem compiling spnalgtest ; exit ; }
      ./spnalgtest $N   || { echo Problem running spnalgtest ; exit ; }
      mv spn_sun_algconv.h $TESTDIR
      sed -i 's/GAUGETYPE/GAUGE_SUN/'  $MAKEDIR/MkFlags || exit
      cd $AUTOSUNDIR
        echo "Building WriteREPR"
        make 
        ./writeREPR $N $TESTDIR/suN_repr_func.h.tmpl > $TESTDIR/suN_func.h || { exit ; } 
        ./writeREPR $N $TESTDIR/suN_exp.c.tmpl > $TESTDIR/suN_exp.h || { exit ; } 
        sed -i 's/GAUGE_SUN/GAUGE_SPN/'  $MAKEDIR/MkFlags
        make 
        ./writeREPR $N $TESTDIR/SP_repr_func.h.tmpl > $TESTDIR/spN_func.h || { exit ; }
        ./writeREPR $N $TESTDIR/SP_exp.c.tmpl> $TESTDIR/spN_exp.h || { exit ; }
    cd $TESTDIR
  
    sed -i 's/suN/SP/g' spN_exp.h
    #compile and run the test
    gcc -o test_reps -O0 -g -I $TOPDIR/Include/ test_reps.c -D$rep || { echo Problem compiling the representation test ; exit ; }
    ./test_reps || { echo Found a bug ; exit ; }
  
  done #rep
done #N






