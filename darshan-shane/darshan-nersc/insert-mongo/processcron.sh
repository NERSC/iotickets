#!/bin/bash

#usage: processcron.sh `date +%s` hostname

SCRIPTDIR=`dirname $0`

datenow=$1
NERSC_HOST=$2
date +%Y/%-m/%-d -d@$1
#echo -n `date` ":"

if [ $NERSC_HOST == "cori" ]
  then
    logpath=/global/cscratch1/sd/darshanlogs/`date +%Y/%-m/%-d -d@$datenow`
elif [ $NERSC_HOST == "edison" ]
  then
    logpath=/scratch1/scratchdirs/darshanlogs/`date +%Y/%-m/%-d -d@$datenow`
fi
#logpath=/scratch1/scratchdirs/darshanlogs/darshanmongo-das
#echo $logpath

mkdir -p $SCRIPTDIR/nosvn_log_$NERSC_HOST/

oldlog=$SCRIPTDIR/nosvn_log_$NERSC_HOST/daylog_`date +%Y_%-m_%-d -d@$datenow`.old
newlog=$SCRIPTDIR/nosvn_log_$NERSC_HOST/daylog_`date +%Y_%-m_%-d -d@$datenow`.new
errlog=$SCRIPTDIR/nosvn_log_$NERSC_HOST/daylog_`date +%Y_%-m_%-d -d@$datenow`.stderr
oufname=$SCRIPTDIR/nosvn_log_$NERSC_HOST/daylog_`date +%Y_%-m_%-d -d@$datenow`.latest.js

#echo $oldlog
#echo $newlog

if [ ! -e $oldlog ]
  then
    touch $oldlog
fi

find $logpath -name "*.darshan.gz" -o -name "*.darshan" > $newlog

#echo `which python`

diff $newlog $oldlog| grep "^<" |cut -b 3-|xargs -I '{}' -n1 python $SCRIPTDIR/../log2mongoinsert.py slog $NERSC_HOST {} 2>> $errlog > $oufname 2>> $errlog

rm -f $oldlog
mv $newlog $oldlog

#Insert into MongoDB
echo `wc $oufname|awk '{print $1}'` "new records"
$SCRIPTDIR/connectrwmongo.sh $oufname


