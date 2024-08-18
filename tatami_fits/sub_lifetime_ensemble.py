#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

#Choose the executable
# LIFETIMEFITS
executable = 'lifetime_fit'


#numFits_per_job = 2; numJobs = 500
#numFits_per_job = 5; numJobs = 200
#numFits_per_job = 10; numJobs = 100
numFits_per_job = 20; numJobs = 50
#numFits_per_job = 50; numJobs = 20
#samples = ['S0p2', 'S0p4', 'S0p6', 'Sm0p2', 'Sm0p4', 'Sm0p6', 'S0p0']
samples = ['S0p0']

#numFits_per_job = 1; numJobs = 5; samples = ['S0p2']
#numFits_per_job = 1; numJobs = 5; samples = ['S0p2']


#Normal jobs
for job in range(numJobs):
	for smc_sample in samples:

		os.system('bsub -q s python3 run_lifetimeFit.py %s %s %s %s'%(numFits_per_job, job, smc_sample, executable))








