#!/bin/bash

TMPFILENAME=cinfo.tmp.$BASHPID
FILENAME=cinfo.c.$BASHPID
FILENAME_TRUE=cinfo.c
MKDIR=$1
TOPDIR=$2
PWD=`pwd`

shift; shift
MACROS=$@


echo "// created automatically by $0 on $(date)" > ${FILENAME}


echo -n ${MACROS} | tr '"' '@' > $TMPFILENAME
len=$((`cat $TMPFILENAME | wc -c`+1))
echo -n "static char MACROS[${len}] = \"" >> ${FILENAME}
cat $TMPFILENAME >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm $TMPFILENAME

awk '{printf "%s\\n",$0}' ${MKDIR}/MkFlags > $TMPFILENAME
len=$((`cat $TMPFILENAME | wc -c`+1))
echo -n "static char CI_mkflags[${len}] = \"" >> ${FILENAME}
cat $TMPFILENAME >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm $TMPFILENAME

if [[ -a /proc/cpuinfo ]] 
then
 awk '{printf "%s\\n",$0}' /proc/cpuinfo > $TMPFILENAME
else
 echo -n "No CPU info\n" > $TMPFILENAME
fi
len=$((`cat $TMPFILENAME | wc -c`+1))
echo -n "static char CI_cpuinfo[${len}] = \"" >> ${FILENAME}
cat $TMPFILENAME >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm $TMPFILENAME

if [[ -a /proc/version ]] 
then
 awk '{printf "%s\\n",$0}' /proc/version > $TMPFILENAME
else
 echo -n "No VERSION info\n" > $TMPFILENAME
fi
len=`cat $TMPFILENAME | wc -c`+1
echo -n "static char CI_linux[${len}] = \"" >> ${FILENAME}
cat $TMPFILENAME >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm $TMPFILENAME

gcc -v 2>&1 | awk '{printf "%s\\n",$0}' > $TMPFILENAME
len=`cat $TMPFILENAME | wc -c`+1
echo -n "static char CI_gcc[${len}] = \"" >> ${FILENAME}
cat $TMPFILENAME >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm $TMPFILENAME

#svn --version >&- ; ret=$?
#if [ "${ret}" -eq "0" ] ; then
#  svn info ${TOPDIR} | awk '{printf "%s\\n",$0}' > $TMPFILENAME
#  len=`cat $TMPFILENAME | wc -c`+1
#  echo -n "static char CI_svninfo[${len}] = \"" >> ${FILENAME}
#  cat $TMPFILENAME >> ${FILENAME}
#  echo "\";" >> ${FILENAME}
#  echo "" >> ${FILENAME}
#  rm $TMPFILENAME

#  svn st -q ${TOPDIR} | awk '{printf "%s\\n",$0}' > $TMPFILENAME
#  len=`cat $TMPFILENAME | wc -c`+1
#  echo -n "static char CI_svnstatus[${len}] = \"" >> ${FILENAME}
#  cat $TMPFILENAME >> ${FILENAME}
#  echo "\";" >> ${FILENAME}
#  echo "" >> ${FILENAME}
#  rm $TMPFILENAME
#else
  echo -n "static char CI_svninfo[1] = \"\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
  echo -n "static char CI_svnstatus[1] = \"\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
#fi

#REV=$(svn info | grep Revision | awk '{ print $2 }')
#if [ -z "$REV" ]; then
  REV="0"
#fi
echo "static int CI_svnrevision = ${REV};" >> ${FILENAME}
echo "" >> ${FILENAME}

cat ${MKDIR}/Utils/${FILENAME_TRUE}.tmpl >> ${FILENAME}

mv $FILENAME $FILENAME_TRUE
