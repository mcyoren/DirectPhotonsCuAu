#!/bin/sh


export INSTALL=/direct/phenix+u/tongzhouguo/install
export LD_LIBRARY_PATH=$INSTALL/lib:$LD_LIBRARY_PATH

echo "==============================================="
echo "============= START SIMULATION  ==============="
echo "==============================================="


root -l -b << EOF
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisReco");
gSystem->Load("/direct/phenix+u/roli/scratch/HELIOS/WriteEvent_C.so");
.L get_Ntag_uncorr_sim_fastmask.C+
gSystem->Load("get_Ntag_uncorr_sim_fastmask_C.so");
get_Ntag_uncorr_sim_fastmask($1, $2)
EOF

