cd /gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way

export toyindex=$1
export Acp_A=$2
export Acp_S=$3

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
export plot_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/plots/linearity_fits
export result_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/results/linearity_fits
export log_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/linearity_logs
export consolidated_results_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2/consolidated_results_a${Acp_A}_s${Acp_S}
export linearity_dir=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2 #for cleaning script

#execute command
if [ "$0" = "$BASH_SOURCE" ] ; then


./exe/linearityTest -f $input_dir/toy_dataset_${toyindex}.root -p $plot_dir/linearity_fit${toyindex}_a${Acp_A}_s${Acp_S}.pdf -r $result_dir/result_linearity${toyindex}_a${Acp_A}_s${Acp_S}.txt -i $toyindex -c $consolidated_results_dir -a $Acp_A -s $Acp_S  >> $log_dir/log_linearity${toyindex}_a${Acp_A}_s${Acp_S}.txt

fi



