cd /gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way

export toyindex=$1
export smc_sample=$2 


export BELLE_LEVEL=b20090127_0910
export BELLE_DEBUG=opt
export BELLE_MSG_MAX_SHLVL=1
export USE_I386_EXP5X=
export MALLOC_CHECK_=0
if [ -f /sw/belle/local/etc/bashrc_general ]; then
    . /sw/belle/local/etc/bashrc_general
fi
export ROOTSYS=/home/belle/varghese/IPHC/old-root/root-tatami
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:.:$LD_LIBRARY_PATH
export MANPATH=$MANPATH:$ROOTSYS/man
export BELLE_FLC=/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way
export LD_LIBRARY_PATH=$BELLE_FLC/lib:$LD_LIBRARY_PATH



export input_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/bootstrapped_datasets
export plot_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/plots/linearity_fits_v2
export result_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/results/linearity_fits_v2
export log_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/linearity_v2_logs
export consolidated_results_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/consolidated_results


#execute command
if [ "$0" = "$BASH_SOURCE" ] ; then
./exe/linearityToyFit_v2_correctMCMistag -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

#./exe/linearityToyFit_v2_binBybin -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

#./exe/linearity_sigOnly_reco1D_correctedMistag -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

#./exe/linearity_sigOnly_reco1D_genFlav -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

#./exe/linearity_sigOnly_gen -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

#./exe/linearity_sigOnly_reco1D -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

#./exe/linearityTest_v2 -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

#./exe/linearityToyFit_v2_2plus1D -f $input_dir/toy_dataset_${smc_sample}_${toyindex}.root -p $plot_dir/toy_fit_${smc_sample}_${toyindex}.pdf -r $result_dir/result_toy_${smc_sample}_${toyindex}.txt -i $toyindex -c $consolidated_results_dir -t $smc_sample   >> $log_dir/log_toy_${smc_sample}_${toyindex}.txt

fi



