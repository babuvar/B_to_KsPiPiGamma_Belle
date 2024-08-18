#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

#executable = 'dtres_syst'
executable = 'dtres_syst_dalitz'

#numFits_per_job = 2; numJobs = 500
numFits_per_job = 5; numJobs = 200
#numFits_per_job = 10; numJobs = 100
#numFits_per_job = 20; numJobs = 50
#numFits_per_job = 50; numJobs = 20

components = ['all']#for dtres systematics


#numFits_per_job = 2; numJobs = 2


#Normal jobs
for job in range(numJobs):
	for component  in components:

		os.system('bsub -q s python3 run_systFit.py %s %s %s %s'%(numFits_per_job, job, component, executable))

		#print('bsub -q s python3 run_systFit.py %s %s %s %s'%(numFits_per_job, job, component, executable))







