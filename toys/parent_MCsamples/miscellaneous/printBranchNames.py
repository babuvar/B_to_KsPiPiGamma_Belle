import os
import sys
from ROOT import TFile


f = TFile('/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/genericMC_x6_wRMVA_wB2MVA_sel_BCS_branches.root')
tree = f.Get('h1')
branchlist = tree.GetListOfBranches().Clone()
namelist = []
for branch in branchlist:
	#print(branch.GetName())
	namelist.append(branch.GetName())

print(namelist)
