#!/bin/csh
source /opt/phenix/core/bin/phenix_setup.csh -n

set LIST = (`ls -lhtr /phenix/plhf/tongzhouguo/taxi/Run12CuAu200MinBias/19325/data/*.root | awk '{printf("%s\n",$9)}'`)

set NUM = 0
set SYSTEM = 2
set MODE = 0

setenv MYINSTALL /direct/phenix+u/tongzhouguo/install
setenv LD_LIBRARY_PATH ${MYINSTALL}/lib:${LD_LIBRARY_PATH}

pushd output

ln -sf $PWD/../lookup_3D_one_phi.root .
ln -sf $PWD/../get_Ntag_uncorr_C.so .
#ln -sf $PWD/../Run12CuAuEMCal_deadmap.txt .

foreach file ($LIST)

  if ( $NUM == $1 ) then

    set NAME = `echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    set RUNNO = `echo $NAME | awk -F \. '{printf("%s\n",$1)}'`

    echo $NAME $RUNNO

    set output = ${PWD}/first_TTree_${RUNNO}.root

    if ( -e $output ) then
      echo "File already exists"
    else

root -l -b << EOF
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisReco");
.L get_Ntag_uncorr.C+
gSystem->Load("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Analysis/Tongzhou/real/get_Ntag_uncorr_C.so");
get_Ntag_uncorr("$file","$output","$RUNNO",$SYSTEM,$MODE)
.q
EOF

    endif
  endif

  @ NUM = $NUM + 1

end

popd
