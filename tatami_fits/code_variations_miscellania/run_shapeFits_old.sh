source tatami_fit_setenv.sh

shape_rootdir="/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles/for_shape"
rootdir="/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles"

#Continuum
./exe/cont_shape -f $shape_rootdir/continuum_5streams.root > shape_logs/cont_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/cont_log.txt > results/fix_shape/Correlation_matrices/cont_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/cont_log.txt > results/fix_shape/Correlation_matrices/cont_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/cont_log.txt > results/fix_shape/Correlation_matrices/cont_dT.txt 
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/cont_log.txt > results/fix_shape/Covariance_matrices/cont_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/cont_log.txt > results/fix_shape/Covariance_matrices/cont_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/cont_log.txt > results/fix_shape/Covariance_matrices/cont_dT.txt


#B-Bbar
./exe/bb_shape -f $shape_rootdir/bb_5streams.root > shape_logs/bb_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/bb_log.txt > results/fix_shape/Correlation_matrices/bb_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/bb_log.txt > results/fix_shape/Correlation_matrices/bb_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/bb_log.txt > results/fix_shape/Correlation_matrices/bb_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/bb_log.txt > results/fix_shape/Covariance_matrices/bb_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/bb_log.txt > results/fix_shape/Covariance_matrices/bb_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/bb_log.txt > results/fix_shape/Covariance_matrices/bb_dT.txt


#Rare-bkg
./exe/rarebkg_shape -f $shape_rootdir/merged_MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches_shuffled_49PartsIn50_bkg.root > shape_logs/rare_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/rare_log.txt > results/fix_shape/Correlation_matrices/rare_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/rare_log.txt > results/fix_shape/Correlation_matrices/rare_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/rare_log.txt > results/fix_shape/Correlation_matrices/rare_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/rare_log.txt > results/fix_shape/Covariance_matrices/rare_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/rare_log.txt > results/fix_shape/Covariance_matrices/rare_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/rare_log.txt > results/fix_shape/Covariance_matrices/rare_dT.txt


#SCF
./exe/scf_shape -f $shape_rootdir/b_kspipigam_sig_500k_wRMVA_wB2MVA_sel_BCS_branches_SCF.root > shape_logs/scf_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/scf_log.txt > results/fix_shape/Correlation_matrices/scf_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/scf_log.txt > results/fix_shape/Correlation_matrices/scf_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/scf_log.txt > results/fix_shape/Correlation_matrices/scf_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/scf_log.txt > results/fix_shape/Covariance_matrices/scf_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/scf_log.txt > results/fix_shape/Covariance_matrices/scf_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/scf_log.txt > results/fix_shape/Covariance_matrices/scf_dT.txt


#Signal
./exe/signal_shape -f $shape_rootdir/b_kspipigam_sig_500k_wRMVA_wB2MVA_sel_BCS_branches_Puresignal.root > shape_logs/signal_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/signal_log.txt > results/fix_shape/Correlation_matrices/signal_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/signal_log.txt > results/fix_shape/Correlation_matrices/signal_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/signal_log.txt > results/fix_shape/Correlation_matrices/signal_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/signal_log.txt > results/fix_shape/Covariance_matrices/signal_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/signal_log.txt > results/fix_shape/Covariance_matrices/signal_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/signal_log.txt > results/fix_shape/Covariance_matrices/signal_dT.txt


#Full fit
#./exe/fullFit -f $rootdir/mc_cocktail_4.root






