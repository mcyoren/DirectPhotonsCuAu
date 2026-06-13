#!/bin/sh

LIST=`ls -lhtr /phenix/plhf/mlario/taxi/pp/Run15pp200CAnoVTXMBPro104/17265/data/*.root | awk '{printf("%s\n",$9)}'`
NUM=0
SYSTEM=2
MODE=0

export MYINSTALL=/direct/phenix+u/tongzhouguo/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH

pushd output

ln -s $PWD/../lookup_3D_one_phi.root .
ln -s $PWD/../get_Ntag_uncorr_C.so .
#ln -s $PWD/../Run12CuAuEMCal_deadmap.txt .

for file in $LIST
do
  if (( $NUM == $1 ))
  then
    NAME=`echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    RUNNO=`echo $NAME | awk -F \. '{printf("%s\n",$1)}'`
    echo $NAME $RUNNO 
    output=$PWD/first_TTree_$RUNNO.root
    if [ -a $PWD/first_TTree_$RUNNO.root ]
    then
      echo "File already exists"
    else
      root -l -b << EOF
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisReco");
.L get_Ntag_uncorr.C+
gSystem->Load("get_Ntag_uncorr_C.so");
get_Ntag_uncorr("$file","$output","$RUNNO",$SYSTEM,$MODE)
EOF
    fi
  fi
  NUM=$(( $NUM + 1 ))
done

popd

