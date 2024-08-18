source /home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/run_linearityFit_v2.sh 


echo " cleaning toy-plots in $plot_dir "
#rm -f $plot_dir/*.pdf
rm -rf $plot_dir
mkdir $plot_dir

echo " cleaning results in $result_dir "
#rm -f $result_dir/*.txt
rm -rf $result_dir 
mkdir $result_dir

echo " cleaning logs in $log_dir "
#rm -f $log_dir/*.txt
rm -rf $log_dir
mkdir $log_dir

echo " cleaning consolidated results in $consolidated_results_dir "
#rm -f $consolidated_results_dir/*.txt
rm -rf $consolidated_results_dir
mkdir $consolidated_results_dir


