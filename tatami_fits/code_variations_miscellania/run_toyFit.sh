cd /gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way



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



export input_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part1/bootstrapped_datasets
export plot_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/plots/toy_fits
export result_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/results/toy_fits
export log_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/toy_logs
export consolidated_results_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part1/consolidated_results


#execute command
if [ "$0" = "$BASH_SOURCE" ] ; then

#./exe/toyFit -f $input_dir/toy_dataset_$1.root -p $plot_dir/toy_fit$1.eps -r $result_dir/result_toy$1.txt -i $1 -c $consolidated_results_dir   >> $log_dir/log_toy$1.txt

./exe/toyFit -f $input_dir/toy_dataset_$1.root -p $plot_dir/toy_fit$1.pdf -r $result_dir/result_toy$1.txt -i $1 -c $consolidated_results_dir   >> $log_dir/log_toy$1.txt

fi



