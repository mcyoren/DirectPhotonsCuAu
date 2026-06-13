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
EOF



echo "==============================================="
echo "================ BUILD IS MADE ================"
echo "==============================================="

./get_Ntag_uncorr_sim.sh $1 1 &
./get_Ntag_uncorr_sim.sh $1 2 &
./get_Ntag_uncorr_sim.sh $1 3 &
./get_Ntag_uncorr_sim.sh $1 4 &
./get_Ntag_uncorr_sim.sh $1 5 &
./get_Ntag_uncorr_sim.sh $1 6 &
./get_Ntag_uncorr_sim.sh $1 7 &
./get_Ntag_uncorr_sim.sh $1 8 &
./get_Ntag_uncorr_sim.sh $1 9 &
