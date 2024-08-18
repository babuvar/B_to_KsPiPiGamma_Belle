#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

#Choose the executable
# CP FITS
#executable = 'sigOnly_genDt_simFit'
#executable = 'linearityToyFit_v2_3DKDE_correctMCMistag'
#executable = 'linearityToyFit_v2_correctMCMistag' #latest & greatest
#executable = 'mbckDE_2DFit_yieldBiasCheck'
#executable = 'linearityToyFit_v2_binBybin'
#executable = 'linearity_sigOnly_reco1D_correctedMistag'
#executable = 'linearity_sigOnly_reco1D_genFlav'
#executable = 'linearity_sigOnly_gen'
#executable = 'linearity_sigOnly_reco1D'
#executable = 'linearityTest_v2'
#executable = 'linearityToyFit_v2_2plus1D'
#executable = 'fullFit_explicitPDFdefinitions'
executable = 'dalitzFit_explicitPDFdefinitions'

#numFits_per_job = 2; numJobs = 500
#numFits_per_job = 5; numJobs = 200
numFits_per_job = 10; numJobs = 100
#numFits_per_job = 20; numJobs = 50
#numFits_per_job = 50; numJobs = 20
samples = ['S0p2', 'S0p4', 'S0p6', 'Sm0p2', 'Sm0p4', 'Sm0p6', 'S0p0']

#numFits_per_job = 1; numJobs = 5; samples = ['S0p6']
#numFits_per_job = 1; numJobs = 5; samples = ['S0p2']
numFits_per_job = 1; numJobs = 5; samples = ['S0p4']



#Normal jobs
for job in range(numJobs):
	for smc_sample in samples:

		#os.system('bsub -q s python3 run_linearityFit_v2.py %s %s %s %s'%(numFits_per_job, job, smc_sample, executable))
		print('bsub -q s python3 run_linearityFit_v2.py %s %s %s %s'%(numFits_per_job, job, smc_sample, executable))		


'''
#Dalitz jobs
executable = 'dalitzFit'
for job in range(numJobs):
	for smc_sample in samples:

		os.system('bsub -q s python3 run_dalitzFit.py %s %s %s %s'%(numFits_per_job, job, smc_sample, executable))
'''

'''
#bin by bin jobs
executable = 'linearityToyFit_v2_correctMCMistag_byBin'
for job in range(numJobs):
	for smc_sample in samples:
		for rBin in range(1,7):
			os.system('bsub -q s python3 run_linearityFit_v2_bin.py %s %s %s %s %s'%(numFits_per_job, job, smc_sample, executable, rBin))
'''			






