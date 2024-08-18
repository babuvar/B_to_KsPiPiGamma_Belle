#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

#Choose the executable
#executable = 'fixedShape_syst'
#executable = 'wtag_syst'
#executable = 'deltaMtauB_syst'
#executable = 'fixedShape_syst_dalitz'
#executable = 'wtag_syst_dalitz'
executable = 'deltaMtauB_syst_dalitz'

#numFits_per_job = 2; numJobs = 500
numFits_per_job = 5; numJobs = 200
#numFits_per_job = 10; numJobs = 100
#numFits_per_job = 20; numJobs = 50
#numFits_per_job = 50; numJobs = 20

#components = ['bb', 'cont', 'scf', 'rare', 'signal', 'all']
#components = ['bbrand', 'cont', 'scf', 'bbmissfsp', 'signal', 'all']
#components = ['dummy']#for wtag systematics
components = ['all']

#numFits_per_job = 1; numJobs = 3; components = ['all']
#numFits_per_job = 1; numJobs = 3



#Normal jobs
for job in range(numJobs):
	for component  in components:

		os.system('bsub -q s python3 run_systFit.py %s %s %s %s'%(numFits_per_job, job, component, executable))

		#print('bsub -q s python3 run_systFit.py %s %s %s %s'%(numFits_per_job, job, component, executable))







