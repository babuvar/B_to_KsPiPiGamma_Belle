#! /bin/tcsh -f
                        
setenv BELLE_LEVEL b20090127_0910
setenv BASF_USER_IF basfsh.so
setenv BASF_NPROCESS 0
setenv BELLE_DEBUG opt
setenv USE_GRAND_REPROCESS_DATA 1
setenv BASF_USER_INIT geant_init

source /sw/belle/local/etc/cshrc_general

setenv RECO_HOME /home/belle/varghese/IPHC/B_Kspipigamma_belle/reconstruction
setenv INPUT_DIR /group/belle/bdata_b/mcprod/specialmc2010/mixedrare

setenv FILE MCrare_$2

cd $RECO_HOME/src


basf <<EOF >! $RECO_HOME/log/log_MCrare/log_$FILE.out



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


histogram define $RECO_HOME/output/hbook_MCrare/$FILE.hbk  


process_dir $INPUT_DIR/$1/$2
 

terminate 

EOF






