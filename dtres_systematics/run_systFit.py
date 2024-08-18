#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess

#General config
os.system('source /home/belle/varghese/.bashrc')
work_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/dtres_systematics'
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
os.system('export BELLE_FLC=/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/dtres_systematics')
os.system('export LD_LIBRARY_PATH=$BELLE_FLC/lib:$LD_LIBRARY_PATH')


#Specific config
input_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/bootstrapped_datasets'
plot_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/dtres_systematics/plots/systematic_fits'
result_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/dtres_systematics/results/systematic_fits'
log_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/dtres_systematics/systematic_logs'
consolidated_results_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/systematic_results'

numFits = int(sys.argv[1])
jobindex = int(sys.argv[2])
component = sys.argv[3]
fitindex_low = jobindex * numFits
fitindex_high = (jobindex + 1) * numFits
fitindices = [*range(fitindex_low, fitindex_high)]
executable = sys.argv[4]


if __name__ == '__main__':
	
	for fitindex in fitindices:

		#save plot in pdf file
		os.system('%s/exe/%s'%(work_dir, executable) + ' -f %s/toy_dataset_S0p4_501.root'%input_dir + ' -p %s/toy_fit_%s_%s.pdf'%(plot_dir, component, fitindex) + ' -r %s/result_toy_%s_%s.txt'%(result_dir, component, fitindex) + ' -i %s -c %s'%(fitindex, consolidated_results_dir) + ' -t %s   > %s/log_toy_%s_%s.txt'%(component, log_dir, component, fitindex))
	
		#print('%s/exe/%s'%(work_dir, executable) + ' -f %s/toy_dataset_S0p4_501.root'%input_dir + ' -p %s/toy_fit_%s_%s.pdf'%(plot_dir, component, fitindex) + ' -r %s/result_toy_%s_%s.txt'%(result_dir, component, fitindex) + ' -i %s -c %s'%(fitindex, consolidated_results_dir) + ' -t %s   > %s/log_toy_%s_%s.txt'%(component, log_dir, component, fitindex))



