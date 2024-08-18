#! /bin/tcsh -f
                        
setenv BELLE_LEVEL b20090127_0910
setenv BASF_USER_IF basfsh.so
setenv BASF_NPROCESS 0
setenv BELLE_DEBUG opt
setenv USE_GRAND_REPROCESS_DATA 1
setenv BASF_USER_INIT geant_init

source /sw/belle/local/etc/cshrc_general

setenv RECO_HOME /home/belle/varghese/IPHC/B_Kspipigamma_belle/reconstruction

#Adjust stream number format according to exp-no
if ($3 < 29) then
setenv stream 1$5
else
setenv stream $5
endif



setenv FILE exp$3_run$1_$2_MC_st$5_$4

cd $RECO_HOME/src


basf <<EOF >! $RECO_HOME/log/log_GMC$5/log_$FILE.out



path create main
path create anal

module register fix_mdst RecoDecays
path add_module main fix_mdst
path add_module anal RecoDecays


path add_condition main >:0:anal
path add_condition main =<:0:KILL
path add_condition anal <:0:KILL




initialize
table save belle_event
table save mdst_all
table save mid_all
table save ext_all
table save evtcls_all
table save evtvtx_all


histogram define $RECO_HOME/output/hbook_GMC$5/$FILE.hbk  


 
process_url  http://bweb3/montecarlo.php?ex=$3&rs=$1&re=$2&ty=evtgen-$4&dt=on_resonance&bl=caseB&st=$stream
 

terminate 

EOF






