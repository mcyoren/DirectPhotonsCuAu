#!/bin/bash

cd ./
INPUT=$1
root -l -b <<'EOF'
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonEvent.so");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libsvxcentana.so");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libRun14AuAuLeptonConvReco.so");
.x get_pi0_from_gamma_data_cand.C+
.q
EOF
