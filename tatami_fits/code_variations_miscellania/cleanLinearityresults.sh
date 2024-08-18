source /home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/run_linearityFit.sh

echo " cleaning toy-plots in $plot_dir "
rm -f $plot_dir/*.pdf

echo " cleaning results in $result_dir "
rm -f $result_dir/*.txt

echo " cleaning logs in $log_dir "
rm -f $log_dir/*.txt

echo " cleaning consolidated results in $linearity_dir"
rm -f $linearity_dir/con*/*.txt $linearity_dir/con*/*.pdf

