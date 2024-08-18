#Need to setup the environment for basf2 before running this program

from ROOT import TFile
from array import array
import sys
import os
import glob
sys.path.append('../csmva/basf2_version/')
sys.path.append('../csmva/rootTMVA_version/')


from trainTest import add_addMVA_variable as add_rMVA 
from mva_training import add_addMVA_variable as add_b2MVA

rootdir = 'recFiles'
outdir = 'outFiles'
#files = ['merged_GMC1_charged.root', 'merged_GMC1_charm.root', 'merged_GMC1_mixed.root', 'merged_GMC1_uds.root', 'merged_GMC2_charged.root', 'merged_GMC2_charm.root', 'merged_GMC2_mixed.root', 'merged_GMC2_uds.root', 'merged_GMC3_charged.root', 'merged_GMC3_charm.root', 'merged_GMC3_mixed.root', 'merged_GMC3_uds.root', 'merged_GMC4_charged.root', 'merged_GMC4_charm.root', 'merged_GMC4_mixed.root', 'merged_GMC4_uds.root', 'merged_GMC5_charged.root', 'merged_GMC5_charm.root', 'merged_GMC5_mixed.root', 'merged_GMC5_uds.root', 'merged_GMC_st0_charged.root', 'merged_GMC_st0_charm.root', 'merged_GMC_st0_mixed.root', 'merged_GMC_st0_uds.root', 'b_kspipigam_sig_500k.root', 'merged_data_sideband.root']

#files = ['MCrare_shuffled_1PartsIn50.root', 'MCrare_shuffled_49PartsIn50.root']
#files = ['merged_MCrare_x50.root']
#files = ['merged_SMC_S0p2_100k.root', 'merged_SMC_S0p4_100k.root',  'merged_SMC_S0p6_100k.root', 'merged_SMC_Sm0p2_100k.root', 'merged_SMC_Sm0p4_100k.root',  'merged_SMC_Sm0p6_100k.root']

#files = ['SMC_S0p2_1M.root', 'SMC_S0p6_1M.root', 'SMC_Sm0p2_1M.root', 'SMC_Sm0p6_1M.root', 'SMC_S0p4_1M.root', 'SMC_Sm0p0_1M.root', 'SMC_Sm0p4_1M.root']

#files = ['Sm0p6_1M_correct.root', 'Sm0p4_1M_correct.root', 'Sm0p2_1M_correct.root', 'S0p0_1M_correct.root', 'S0p2_1M_correct.root', 'S0p4_1M_correct.root', 'S0p6_1M_correct.root']

#files = ['S0p0_1M_correct.root']

#Continuum-only
#files = ['merged_GMC1_charm.root',  'merged_GMC3_charm.root',  'merged_GMC5_charm.root', 'merged_GMC2_charm.root', 'merged_GMC4_charm.root', 'merged_GMC_st0_charm.root', 'merged_GMC1_uds.root',  'merged_GMC3_uds.root', 'merged_GMC5_uds.root', 'merged_GMC2_uds.root', 'merged_GMC4_uds.root',  'merged_GMC_st0_uds.root']

#files = ['S0p0_1M_correct_new.root']

#Better way
#files = ['merged_data.root']
#files = ['merged_MCrare_x50.root']
#files = [os.path.basename(x) for x in glob.glob('recFiles/S*.root')]
#files = [os.path.basename(x) for x in glob.glob('recFiles/merged_GMC*.root')]
files = [os.path.basename(x) for x in glob.glob('recFiles/M*.root')]

files2 = []; files3 = []; files4 = []; files5 = []; files6 = []; filesExp = []
for f in files:
	f2 = f[:len(f) - 5] + '_wRMVA.root'; files2.append(f2)
	f3 = f[:len(f) - 5] + '_wRMVA_wB2MVA.root'; files3.append(f3)
	f4 = f[:len(f) - 5] + '_wRMVA_wB2MVA_sel.root'; files4.append(f4)
	f5 = f[:len(f) - 5] + '_wRMVA_wB2MVA_sel_BCS.root'; files5.append(f5)
	f6 = f[:len(f) - 5] + '_wRMVA_wB2MVA_sel_BCS_branches.root'; files6.append(f6)
	fexp = f[:len(f) - 5] + '_wRMVA_expert.root'; filesExp.append(fexp)
file_dic = {'f1' : files, 'f2' : files2, 'f3' : files3, 'f4' : files4, 'f5' : files5, 'f6' : files6, 'fexp' : filesExp}


variables = ['coscms', 'costh', 'k0et', 'k0hso00', 'k0hso02', 'k0hso04', 'k0hso10', 'k0hso12', 'k0hso20', 'k0hso22', 'k0hso24']
b2_fBDT_xml = 'mva_XMLfiles/MyContinuumFBDT.xml'
#b2_fBDT_xml = 'mva_XMLfiles/Sviat_Tristan_b2FBDT.xml'
rTMVA_xml = 'mva_XMLfiles/BDT_BDT.weights.xml'

#Selections
#cuts = ''
#cuts = 'b2MVA > 0.0' #if you want to cut on the basf2-fBDT variable
#cuts = 'b2MVA > 0.05' #if you want to cut on the basf2-fBDT variable
cuts = 'rMVA > 0.5' #if you want to cut on the root-TMVA variable

cuts += ' && pi0veto < 0.2 && etaveto < 0.25'
cuts += ' && g_e9e25 > 0.95 && g_cos > -0.65 && g_cos < 0.86'
cuts += ' && tagcl > 0 && vtxcl > 0'
cuts += ' && m23 > 0.6 && m23 < 0.9'
cuts += ' && x_p0_kid < 0.25 && x_p1_kid < 0.25'
cuts += ' && x_p0_eid < 0.25 && x_p1_eid < 0.25'
cuts += ' && x_m < 1.8'
cuts += ' && ( tagh < 50 || tagntrk == 1) && vtxh < 50'


def apply_selections(input_files, cuts, rootdir, outdir, key):

	for input_file in input_files:
		
		inputfilename = '%s/%s'%(rootdir, input_file)
		file_truncated = input_file[:len(input_file) - 5]
		outputfilename = '%s/%s_sel.root'%(outdir, file_truncated)
		f1 = TFile(inputfilename, 'read')
		initialtree = f1.Get(key)
		f2 = TFile(outputfilename, 'recreate')
		tree = initialtree.CopyTree(cuts)
		f2.Write(); f2.Close()
	print('Selections have been applied to the given files')


def apply_BCS(input_files, rootdir, outdir, key):

	for input_file in input_files:

		inputfilename = '%s/%s'%(rootdir, input_file)
		file_truncated = input_file[:len(input_file) - 5]
		outputfilename = '%s/%s_BCS.root'%(outdir, file_truncated)
		oldfile = TFile(inputfilename); oldtree = oldfile.Get(key)
		newfile = TFile(outputfilename, 'recreate')
		newtree = oldtree.CloneTree(0)
		lasteventid = 0
		#loop over candidates
		for i in range(oldtree.GetEntries()):
			oldtree.GetEntry(i);
			if lasteventid == oldtree.eventid:
				continue
			lasteventid = oldtree.eventid
			newtree.Fill()
		newfile.Write(); newfile.Close()

	print('Finished BCS for the given files')		


def add_branches(input_files, rootdir, outdir, key):

	for input_file in input_files:

		inputfilename = '%s/%s'%(rootdir, input_file)
		file_truncated = input_file[:len(input_file) - 5]
		outputfilename = '%s/%s_branches.root'%(outdir, file_truncated)
		oldfile = TFile(inputfilename); oldtree = oldfile.Get(key)

		mc_type = array('f',[0]); bmc_type = oldtree.Branch('mc_type', mc_type, 'mc_type/F');
		inexp = array('i',[0]); binexp = oldtree.Branch('inexp', inexp, 'inexp/I');
		irbin = array('i',[0]); birbin = oldtree.Branch('irbin', irbin, 'irbin/I');
		iflavor = array('i',[0]); biflavor = oldtree.Branch('iflavor', iflavor, 'iflavor/I');
		dalitzCategory = array('i',[0]); bDC = oldtree.Branch('dalitzCategory', dalitzCategory, 'dalitzCategory/I');

		newfile = TFile(outputfilename, 'recreate')
		newtree = oldtree.CloneTree(0)

		#determine branch MC-type
		mc_type[0] = 50 #default for data and signal-MC
		if 'rare' in input_file:
			mc_type[0] = 0
		if 'charged' in input_file:
			mc_type[0] = 10
		if 'charm' in input_file:
			mc_type[0] = 20
		if 'mixed' in input_file:
			mc_type[0] = 30
		if 'uds' in input_file:
			mc_type[0] = 40

		#loop over candidates
		for i in range(oldtree.GetEntries()):
			oldtree.GetEntry(i);

			# nexp
			if oldtree.nexp < 31:
				inexp[0] = 7
			else:
				inexp[0] = 31
			# rbin
			if oldtree.rbin >= 0.0 and oldtree.rbin <= 0.1: 
				irbin[0] = 0
			if oldtree.rbin > 0.1 and oldtree.rbin <= 0.25:
				irbin[0] = 1
			if oldtree.rbin > 0.25 and oldtree.rbin <= 0.5:
				irbin[0] = 2
			if oldtree.rbin > 0.5 and oldtree.rbin <= 0.625:
				irbin[0] = 3
			if oldtree.rbin > 0.625 and oldtree.rbin <= 0.75:
				irbin[0] = 4
			if oldtree.rbin > 0.75 and oldtree.rbin <= 0.875:
				irbin[0] = 5
			if oldtree.rbin > 0.875 and oldtree.rbin <= 1.0:
				irbin[0] = 6
			# flavor
			iflavor[0] = int(oldtree.flavor)
			# Dalitz Category
			if oldtree.m12 > oldtree.m13:
				dalitzCategory[0] = 1
			else:
				dalitzCategory[0] = -1

			bmc_type.Fill()
			binexp.Fill()
			birbin.Fill()
			biflavor.Fill()
			bDC.Fill()
			newtree.Fill()
		newfile.Write(); newfile.Close()
		del newfile; del newtree; del oldfile; del oldtree;
		del bmc_type; del binexp; del birbin; del biflavor; del bDC 

	print('Finished adding extra branches')

		
def clean_output(filename_dic, outdir = outdir):

	keys = ['f2', 'f3', 'f4', 'f5', 'fexp']

	for idx, outfile in enumerate(filename_dic['f1']):
		#rename files from last step so as to a sane length
		oldname = '%s/%s'%(outdir, filename_dic['f6'][idx])
		newname = '%s/%s'%(outdir, filename_dic['f1'][idx])
		newname = newname[:len(newname) - 5] + '_toFit.root'
		if os.path.isfile(oldname):
			os.rename(oldname, newname)
		#remove files from intermediate steps
		for key in keys:
			intermediate_file = '%s/%s'%(outdir, filename_dic[key][idx])
			if os.path.isfile(intermediate_file):
				os.remove(intermediate_file)

	print('Finished cleaning intermediate files')


if __name__ == '__main__':

	#add root-TMVA variable
	add_rMVA(do_addMVAvar = True, variables = variables, addMVA_inputFiles = files , xmlFile = rTMVA_xml, rootdir = rootdir, outdir = outdir, key = 'h1')

	#add basf2-fBDT variable
	add_b2MVA(do_addMVAvar = True, addMVA_inputFiles = files2, xmlFile = b2_fBDT_xml, rootdir = outdir, outdir = outdir, key = 'h1')

	#Apply selection criteria
	apply_selections(input_files = files3, cuts = cuts, rootdir = outdir, outdir = outdir, key = 'h1')

	#Apply best candidate selection for each event
	apply_BCS(input_files = files4, rootdir = outdir, outdir = outdir, key = 'h1')

	#add extra branches to help with analysis
	add_branches(input_files = files5, rootdir = outdir, outdir = outdir, key = 'h1')

	#optional : remove files from intermediate steps
	#clean_output(filename_dic = file_dic, outdir = outdir)


