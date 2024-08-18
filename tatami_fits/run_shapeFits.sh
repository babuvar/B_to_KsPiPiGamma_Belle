source tatami_fit_setenv.sh


#Continuum
./exe/cont_random_shape  > shape_logs/cont_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/cont_log.txt > results/fix_shape_new/Correlation_matrices/cont_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/cont_log.txt > results/fix_shape_new/Correlation_matrices/cont_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/cont_log.txt > results/fix_shape_new/Correlation_matrices/cont_dT.txt 
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/cont_log.txt > results/fix_shape_new/Covariance_matrices/cont_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/cont_log.txt > results/fix_shape_new/Covariance_matrices/cont_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/cont_log.txt > results/fix_shape_new/Covariance_matrices/cont_dT.txt


#B-Bbar missing FSP
./exe/bb_realBmother_shape > shape_logs/bb_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/bb_log.txt > results/fix_shape_new/Correlation_matrices/bb_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/bb_log.txt > results/fix_shape_new/Correlation_matrices/bb_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/bb_log.txt > results/fix_shape_new/Correlation_matrices/bb_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/bb_log.txt > results/fix_shape_new/Covariance_matrices/bb_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/bb_log.txt > results/fix_shape_new/Covariance_matrices/bb_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/bb_log.txt > results/fix_shape_new/Covariance_matrices/bb_dT.txt


#B-Bbar random
./exe/bb_random_shape > shape_logs/bb_random_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/bb_random_log.txt > results/fix_shape_new/Correlation_matrices/bb_random_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/bb_random_log.txt > results/fix_shape_new/Correlation_matrices/bb_random_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/bb_random_log.txt > results/fix_shape_new/Correlation_matrices/bb_random_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/bb_random_log.txt > results/fix_shape_new/Covariance_matrices/bb_random_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/bb_random_log.txt > results/fix_shape_new/Covariance_matrices/bb_random_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/bb_random_log.txt > results/fix_shape_new/Covariance_matrices/bb_random_dT.txt


#SCF
./exe/scf_random_shape > shape_logs/scf_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/scf_log.txt > results/fix_shape_new/Correlation_matrices/scf_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/scf_log.txt > results/fix_shape_new/Correlation_matrices/scf_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/scf_log.txt > results/fix_shape_new/Correlation_matrices/scf_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/scf_log.txt > results/fix_shape_new/Covariance_matrices/scf_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/scf_log.txt > results/fix_shape_new/Covariance_matrices/scf_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/scf_log.txt > results/fix_shape_new/Covariance_matrices/scf_dT.txt


#Signal
./exe/signal_pure_shape > shape_logs/signal_log.txt
sed '/correlationMBC/,/correlationMBC/!d; //d' shape_logs/signal_log.txt > results/fix_shape_new/Correlation_matrices/signal_MBC.txt
sed '/correlationDE/,/correlationDE/!d; //d' shape_logs/signal_log.txt > results/fix_shape_new/Correlation_matrices/signal_dE.txt
sed '/correlationDT/,/correlationDT/!d; //d' shape_logs/signal_log.txt > results/fix_shape_new/Correlation_matrices/signal_dT.txt
sed '/covarianceMBC/,/covarianceMBC/!d; //d' shape_logs/signal_log.txt > results/fix_shape_new/Covariance_matrices/signal_MBC.txt
sed '/covarianceDE/,/covarianceDE/!d; //d' shape_logs/signal_log.txt > results/fix_shape_new/Covariance_matrices/signal_dE.txt
sed '/covarianceDT/,/covarianceDT/!d; //d' shape_logs/signal_log.txt > results/fix_shape_new/Covariance_matrices/signal_dT.txt


#Full fit
#./exe/fullFit -f $rootdir/mc_cocktail_4.root






