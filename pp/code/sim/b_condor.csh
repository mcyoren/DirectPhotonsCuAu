#!/bin/csh
#source /opt/phenix/core/bin/phenix_setup.csh -n
if (-f /opt/phenix/core/bin/phenix_setup.csh) then
  source /opt/phenix/core/bin/phenix_setup.csh -n
endif

if (-f /opt/phenix/core/bin/phenix_alma9_setup.csh) then
  source /opt/phenix/core/bin/phenix_alma9_setup.csh -n
endif

setenv MYLOG /gpfs/mnt/gpfs02/phenix/plhf/plhf1/$USER/log
setenv MYINSTALL /direct/phenix+u/$USER/install
setenv MYCVS /direct/phenix+u/$USER/cvsdir

setenv PATH .:$HOME/bin:$MYINSTALL/bin:$PATH
setenv LD_LIBRARY_PATH .:$MYINSTALL/lib:$LD_LIBRARY_PATH
setenv TSEARCHPATH /direct/phenix+u/$USER/install
sh get_Ntag_uncorr.sh $*
