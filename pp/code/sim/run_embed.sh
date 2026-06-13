#!/bin/sh


export INSTALL=/direct/phenix+u/tongzhouguo/install
export LD_LIBRARY_PATH=$INSTALL/lib:$LD_LIBRARY_PATH


echo "==============================================="
echo "============= START SIMULATION  ==============="
echo "==============================================="


root -l -b <<EOF
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonEvent.so");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libsvxcentana.so");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonConvReco.so");
.L get_pi0_from_gamma_data.C+
.q
EOF


echo "==============================================="
echo "================ BUILD IS MADE ================"
echo "==============================================="

./get_pi0_from_gamma_data_mt.sh 0 &
./get_pi0_from_gamma_data_mt.sh 1 &
./get_pi0_from_gamma_data_mt.sh 2 &
./get_pi0_from_gamma_data_mt.sh 3 &
./get_pi0_from_gamma_data_mt.sh 4 &
