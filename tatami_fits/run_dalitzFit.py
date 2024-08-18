#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess

#General config
os.system('source /home/belle/varghese/.bashrc')
work_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way'
os.chdir(work_dir)

os.environ['BELLE_LEVEL'] = 'b20090127_0910'
os.environ['BELLE_DEBUG'] = 'opt'
os.environ['BELLE_MSG_MAX_SHLVL'] = '1'
os.environ['USE_I386_EXP5X'] = ''
os.environ['MALLOC_CHECK_'] = '0'
if os.path.exists('/sw/belle/local/etc/bashrc_general'):
	subprocess.call('. /sw/belle/local/etc/bashrc_general', shell=True)
os.environ['ROOTSYS'] = '/home/belle/varghese/IPHC/old-root/root-tatami'
os.system('export PATH=$ROOTSYS/bin:$PATH')
os.system('export LD_LIBRARY_PATH=$ROOTSYS/lib:.:$LD_LIBRARY_PATH')
os.system('export MANPATH=$MANPATH:$ROOTSYS/man')
os.system('export BELLE_FLC=/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way')
os.system('export LD_LIBRARY_PATH=$BELLE_FLC/lib:$LD_LIBRARY_PATH')


#Specific config
input_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/bootstrapped_datasets'
plot_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/plots/dalitz_fits'
result_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/results/dalitz_fits'
log_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/dalitz_logs'
consolidated_results_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/consolidated_results'

numFits = int(sys.argv[1])
jobindex = int(sys.argv[2])
smc_sample = sys.argv[3]
fitindex_low = jobindex * numFits
fitindex_high = (jobindex + 1) * numFits
fitindices = [*range(fitindex_low, fitindex_high)]
executable = sys.argv[4]


if __name__ == '__main__':
	
	for fitindex in fitindices:

		os.system('%s/exe/%s'%(work_dir, executable) + ' -f %s/toy_dataset_%s_%s.root'%(input_dir, smc_sample, fitindex) + ' -p %s/toy_fit_%s_%s.pdf'%(plot_dir, smc_sample, fitindex) + ' -r %s/result_toy_%s_%s.txt'%(result_dir, smc_sample, fitindex) + ' -i %s -c %s'%(fitindex, consolidated_results_dir) + ' -t %s   > %s/log_toy_%s_%s.txt'%(smc_sample, log_dir, smc_sample, fitindex))




