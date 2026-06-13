#!/bin/csh
source /opt/phenix/core/bin/phenix_setup.csh -n

if ( $#argv < 2 ) then
  echo "Usage: ./get_Ntag_uncorr_balanced.csh JOBID NJOBS"
  echo "Example: ./get_Ntag_uncorr_balanced.csh 0 300"
  exit 1
endif

set JOBID = $1
set NJOBS = $2
set SYSTEM = 2
set MODE = 0

setenv MYINSTALL /direct/phenix+u/tongzhouguo/install
setenv LD_LIBRARY_PATH ${MYINSTALL}/lib:${LD_LIBRARY_PATH}

# ============================================================
# Compile-only mode.
# Run this once before submitting many Condor jobs.
# ============================================================
if ( $JOBID == -1 ) then

  echo "Compile-only mode"
  echo "Compiling get_Ntag_uncorr.C in current/original folder: $PWD"

  set ROOTMACRO = compile_get_Ntag_uncorr.C

  cat > $ROOTMACRO << EOF
{
  gSystem->Load("/phenix/plhf/roli/install/lib/libPhotonConversionAnalysisEvent");
  gSystem->Load("/phenix/plhf/roli/install/lib/libPhotonConversionAnalysisReco");
  gROOT->LoadMacro("get_Ntag_uncorr.C+");
}
EOF

  root -l -b -q $ROOTMACRO

  rm -f $ROOTMACRO

  echo "Compile step finished."
  echo "Expected library: get_Ntag_uncorr_C.so"
  exit 0
endif

# ============================================================
# Normal running mode.
# Uses precompiled get_Ntag_uncorr_C.so.
# ============================================================

set COUNTFILE = event_counts.txt
if ( ! -e $COUNTFILE ) then
  echo "Missing $COUNTFILE. First run: ./make_event_counts.csh"
  exit 1
endif

set TOTEVENTS = `awk '{if($4>0) s+=$4} END{printf("%d",s)}' $COUNTFILE`
if ( $TOTEVENTS <= 0 ) then
  echo "Bad total event count: $TOTEVENTS"
  exit 1
endif

set GSTART = `awk -v j=$JOBID -v n=$NJOBS -v tot=$TOTEVENTS 'BEGIN{printf("%d", int(tot*j/n))}'`
set GEND   = `awk -v j=$JOBID -v n=$NJOBS -v tot=$TOTEVENTS 'BEGIN{printf("%d", int(tot*(j+1)/n))}'`

if ( $GEND <= $GSTART ) then
  echo "Empty job range for JOBID=$JOBID NJOBS=$NJOBS"
  exit 0
endif

echo "JOBID=$JOBID NJOBS=$NJOBS total=$TOTEVENTS global_range=[$GSTART,$GEND)"

set PLAN = job_${JOBID}_plan.txt
set AWKSCRIPT = make_plan_${JOBID}.awk

cat > $AWKSCRIPT << EOF
BEGIN {
  cum = 0;
  part = 0;
}
{
  idx = \$1;
  file = \$2;
  runno = \$3;
  nev = \$4;

  if (nev <= 0) next;

  fstart = cum;
  fend = cum + nev;

  s = (gs > fstart ? gs : fstart);
  e = (ge < fend ? ge : fend);

  if (s < e) {
    local_s = s - fstart;
    local_e = e - fstart;
    printf("%d|%s|%s|%d|%d|%d\\n", part, file, runno, local_s, local_e, idx);
    part++;
  }

  cum = fend;
}
EOF

awk -v gs=$GSTART -v ge=$GEND -f $AWKSCRIPT $COUNTFILE > $PLAN
rm -f $AWKSCRIPT

if ( ! -s $PLAN ) then
  echo "No file segments for this job."
  rm -f $PLAN
  exit 0
endif

pushd output_new

ln -sf $PWD/../lookup_3D_one_phi.root .
#ln -sf $PWD/../get_Ntag_uncorr_C.so .
ln -sf $PWD/../dc_dead_map.h .

foreach line (`cat ../$PLAN`)

  set PART  = `echo "$line" | awk -F'|' '{print $1}'`
  set file  = `echo "$line" | awk -F'|' '{print $2}'`
  set RUNNO = `echo "$line" | awk -F'|' '{print $3}'`
  set START = `echo "$line" | awk -F'|' '{print $4}'`
  set END   = `echo "$line" | awk -F'|' '{print $5}'`
  set FIDX  = `echo "$line" | awk -F'|' '{print $6}'`

  set output = ${PWD}/first_TTree_${RUNNO}_job${JOBID}_part${PART}_ev${START}_${END}.root

  echo "Running file_index=$FIDX run=$RUNNO local_range=[$START,$END) output=$output"

  if ( -e $output ) then
    echo "File already exists: $output"
  else

    set ROOTMACRO = run_job_${JOBID}_${PART}.Cget_Ntag_uncorr_C

    cat > $ROOTMACRO << EOF
{
  gSystem->Load("/phenix/plhf/roli/install/lib/libPhotonConversionAnalysisEvent");
  gSystem->Load("/phenix/plhf/roli/install/lib/libPhotonConversionAnalysisReco");
  gSystem->Load("/phenix/plhf/mitran/Analysis/Roli/real/get_Ntag_uncorr_C.so");
  get_Ntag_uncorr("$file", "$output", "$RUNNO", $SYSTEM, $MODE, $START, $END);
}
EOF

    root -l -b -q $ROOTMACRO

    if ( -e $output ) then
      echo "Setting permissions for $output"
      chmod 777 $output
    endif

    rm -f $ROOTMACRO

  endif

end

popd

rm -f $PLAN