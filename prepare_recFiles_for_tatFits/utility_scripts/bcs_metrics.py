#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TH1F, TCanvas, gROOT, gStyle
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import json
import subprocess
import shlex
import uproot
import sys

gROOT.SetBatch(True)
gStyle.SetOptStat(0)

def shuffle(df):
	df = df.reindex(np.random.permutation(df.index))
	df.reset_index(drop=True, inplace=True)
	return df

def split_prop(var):
	if isinstance(var, str):
		return (var, {})
	else:
		return var


'''
gmc_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/genericMC_x6_wRMVA_wB2MVA_sel_BCS_branches.root'
rare_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches.root'
'''

gmc_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/genericMC_x6_wRMVA_wB2MVA_sel.root'
rare_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/MCrare_x50_wRMVA_wB2MVA_sel.root'


variables = ['nexp', 'nrun', 'eventid', 'vtxcl', 'issignal', 'issgnevt']

#load dataframes
df_gmc = uproot.open(gmc_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
#df_gmc = shuffle(df_gmc)
df_rare = uproot.open(rare_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
#print(df_rare)
#df_rare = shuffle(df_rare)

# estimate true - yields
true_yields = {'gmc_yield' : len(df_gmc)/6.0,
               'rare_yield' : len(df_rare)/50.0}

#print(true_yields)

#print('gmc_yield = ', len(df_gmc))
#print('rare_yield = ', len( df_rare))

dfg_gmc = df_gmc.groupby(['nexp', 'nrun', 'eventid'])
dfg_rare = df_rare.groupby(['nexp', 'nrun', 'eventid'])
'''
#print('df_rare post BCS yield = ', dfg_rare.ngroups)
BCS_yield = dfg_rare.ngroups
dfg_rare_counts = dfg_rare.size().reset_index(name='counts')
dfg_gmc_counts = dfg_gmc.size().reset_index(name='counts')
multiplicity_vals = [*range(8)]
multiplicity_rare = np.zeros(8)
multiplicity_gmc = np.zeros(8)
multiplicity_combined = np.zeros(8)
n_multicands = 0
n_bcsyield_nominal = 0

for mult_val in multiplicity_vals:

	multiplicity_rare[mult_val] = np.count_nonzero(dfg_rare_counts['counts'] == mult_val)
	multiplicity_gmc[mult_val] = np.count_nonzero(dfg_gmc_counts['counts'] == mult_val)
	multiplicity_combined[mult_val] = (multiplicity_rare[mult_val] / 50) + (multiplicity_gmc[mult_val] / 6)

	#count number of events with multiple candidates
	n_bcsyield_nominal = n_bcsyield_nominal + multiplicity_combined[mult_val]
	if mult_val > 1:
		n_multicands = n_multicands + multiplicity_combined[mult_val]

print('number of events with multiple candidates = ', n_multicands)
print('number of events post bcs = ', n_bcsyield_nominal)
print('% of events with multiple candidates = ', float((n_multicands*100)/n_bcsyield_nominal))
#print('rare-multiplicity = ', multiplicity_rare)
#print('gmc-multiplicity = ', multiplicity_gmc)

#Multiplicity plot
hmult = TH1F("multiplicity","Candidate multiplicity", 5, 0.5, 5.5)
for i in range(4):
	j=i+1
	for val in range(int(multiplicity_combined[j])):
		hmult.Fill(j)

hmult.SetXTitle("Candidate multiplicity")

can = TCanvas("can", "can", 900, 600)
can.cd(); can.SetLogy() 
hmult.Draw("E1")
can.SaveAs('multiplicity.png')
'''

#success rate of BCS
num_evts_withSignal = 0
num_evts_correctBCS = 0
for name, group in dfg_rare:
	grp_idx = -1; containsSig_flag = 0
	if len(group) > 1:
		for index, row in group.iterrows():
			grp_idx = grp_idx + 1
			if row['issignal'] == 1:
				containsSig_flag = 1
				#Is bcs-candidate a true-signal		?
				if grp_idx == 0:
					num_evts_correctBCS = num_evts_correctBCS + 1
		if containsSig_flag == 1:
			num_evts_withSignal = num_evts_withSignal + 1


print('num_evts_withSignal = ', num_evts_withSignal)
print('num_evts_correctBCS = ', num_evts_correctBCS)



