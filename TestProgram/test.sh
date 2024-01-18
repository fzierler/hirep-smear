


MKFILE=$(pwd)/../Make/MkFlags
cp $MKFILE $MKFILE.bk

sed -e 's/#.*$//;s/^\s*//' $MKFILE.bk | grep -vE '^\s*$|^REPR|^GAUGE_GROUP|^NG|^CFLAGS' > $MKFILE

cat $MKFILE

CLEANUP(){

    cp $MKFILE.bk $MKFILE

}

trap CLEANUP SIGHUP SIGINT SIGQUIT SIGILL SIGTRAP SIGABRT SIGKILL SIGSEGV SIGTERM 

TEST_ALGEBRA_1(){
  GROUP=$1
  shift;
  for REPR in REPR_FUNDAMENTAL REPR_SYMMETRIC REPR_ANTISYMMETRIC REPR_ADJOINT
  do
    for NG in $@
    do
    if [ "$NG" -ne 2 ] || [ "$REPR" != "REPR_ANTISYMMETRIC" ]; then
      MAKEVARS="REPR=$REPR GAUGE_GROUP=$GROUP NG=$NG"
      echo $MAKEVARS
      touch $MKFILE
      cd Algebra
      echo "Test check_algebra_1"
      echo $MAKEVARS > $GROUP\_$REPR\_$NG\_check_algebra_1.out
      make $MAKEVARS CFLAGS="-Wall -std=c99 -g -O2" &> /dev/null && ./check_algebra_1 >> $GROUP\_$REPR\_$NG\_check_algebra_1.out
      if [ ! $? -eq 0 ];then CLEANUP; exit 1; fi
      if [ "$NG" -eq 3 ] && [ "$REPR" == "REPR_ANTISYMMETRIC" ] && [ "$GROUP" == "GROUP_SUN" ]; then
          ./check_algebra_2 >> $GROUP\_$REPR\_$NG\_check_algebra_2.out
          if [ ! $? -eq 0 ];then CLEANUP; exit 1; fi
      fi


      cd ..
      cd PureGauge
      echo "Test check_puregauge_1"
      echo $MAKEVARS > $GROUP\_$REPR\_$NG\_check_puregauge_1.out
      make $MAKEVARS CFLAGS="-Wall -std=c99 -g -O2" &> /dev/null && ./check_puregauge_1 >> $GROUP\_$REPR\_$NG\_check_puregauge_1.out
      if [  ! $? -eq 0 ];then CLEANUP; exit 1; fi
      if [ "$NG" -eq 2 ] ;then
          echo "Test check_puregauge_2"
          echo $MAKEVARS > $GROUP\_$REPR\_$NG\_check_puregauge_2.out
        ./check_puregauge_2
        cat out_0 >> $GROUP\_$REPR\_$NG\_check_puregauge_2.out
      fi
      if [  ! $? -eq 0 ];then CLEANUP; exit 1; fi
      echo "Test check_puregauge_3"
      echo $MAKEVARS > $GROUP\_$REPR\_$NG\_check_puregauge_3.out
      ./check_puregauge_3 >> $GROUP\_$REPR\_$NG\_check_puregauge_3.out
      if [  ! $? -eq 0 ];then CLEANUP; exit 1; fi

      cd ..
    fi
    done
  done
  
}

TEST_ALGEBRA_1 GAUGE_SUN 2 3 4 
TEST_ALGEBRA_1 GAUGE_SPN 2 4 6 8

