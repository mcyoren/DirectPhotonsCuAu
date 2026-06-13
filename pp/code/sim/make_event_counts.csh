#!/bin/csh

source /opt/phenix/core/bin/phenix_setup.csh -n
set LIST = (`ls -lhtr /phenix/plhf/mlario/taxi/pp/Run15pp200CAnoVTXMBPro104/17265/data/*.root | awk '{print $9}'`)

set OUT = event_counts.txt
rm -f $OUT

set NUM = 0

foreach file ($LIST)

  set NAME = `basename $file`
  set RUNNO = `echo $NAME | awk -F . '{print $1}'`

  set TMP = count_tmp_${NUM}.txt
  rm -f $TMP

root -l -b << EOF > $TMP
.L count_one.C
count_one("$file");
.q
EOF

  set NEVT = `grep NEVENTS $TMP | tail -1 | awk '{print $2}'`

  if ("$NEVT" == "") then
    set NEVT = -1
  endif

  echo "$NUM $file $RUNNO $NEVT" | tee -a $OUT

  rm -f $TMP

  @ NUM = $NUM + 1

end

echo "Wrote $OUT"