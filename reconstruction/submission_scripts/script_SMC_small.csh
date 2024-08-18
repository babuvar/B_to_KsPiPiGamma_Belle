#! /bin/tcsh -f

setenv RECO_HOME /home/belle/varghese/IPHC/B_Kspipigamma_belle/reconstruction 
setenv SMC_Mdst /group/belle/users/varghese/mcproduzh/gsim/mdst 


cd $RECO_HOME/src

source /home/belle/varghese/.cshrc
setenv BASF_USER_INIT geant_init



foreach FILE (b_kspipigam_S0p0_1M)


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


histogram define $RECO_HOME/output/hbook_SMC/hbook_$FILE.hbk  

process_dir $SMC_Mdst/$FILE/dir1/

terminate 

EOF


end



