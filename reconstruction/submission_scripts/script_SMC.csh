#! /bin/tcsh -f

setenv RECO_HOME /home/belle/varghese/IPHC/B_Kspipigamma_belle/reconstruction
#setenv SMC_Mdst /home/belle/varghese/IPHC/B_Kspipigamma_belle/mcproduzh/gsim/mdst
setenv SMC_Mdst /group/belle/users/varghese/mcproduzh/gsim/mdst


cd $RECO_HOME/src


source  /home/belle/varghese/.cshrc
setenv BASF_USER_INIT geant_init

         
setenv FILE $1_$2

basf <<EOF >! $RECO_HOME/log/log_SMC/log_$FILE.out


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


histogram define $RECO_HOME/output/hbook_SMC/$FILE.hbk  


process_dir $SMC_Mdst/$1/dir$2/



terminate 

EOF






