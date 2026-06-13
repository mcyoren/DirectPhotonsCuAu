#!/bin/csh
# source the alma9 setup script to get our alma9 environment
source /opt/phenix/core/bin/phenix_alma9_setup.csh -n
# run our SL7 script in an SL7 singularity container
sl7cmd ./get_Ntag_uncorr.csh $*
echo wrapper done
