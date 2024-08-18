#Need to setup the environment for basf2 before running this program

from ROOT import TChain, TFile
from array import array
import sys
import os
sys.path.append('../csmva/basf2_version/')


from mva_training import read_shuffle_store

#filenames = ['merged_MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches.root']
#filenames = ['S0p0_1M_correct_wRMVA_wB2MVA_sel_BCS_branches.root']
#filenames = ['b_kspipigam_sig_500k_wRMVA_wB2MVA_sel_BCS_branches_Puresignal.root']
filenames = ['MCrare_full_bkg.root']
rootdir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/prepare_recFiles_for_tatFits/splitfiles'
#filenames = ['merged_MCrare_x50_wRMVA_wB2MVA.root']
#filenames = ['merged_MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches.root']
#rootdir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/prepare_recFiles_for_tatFits/outFiles'
#rootdir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/prepare_recFiles_for_tatFits/outFiles_with_myb2XMLFile'
InputFiles = []
for filename in filenames:
	InputFile = '%s/%s'%(rootdir, filename)
	InputFiles.append(InputFile)
#MCrarefilename = 'MCrare_wRMVA_wB2MVA_shuffled.root'
#MCrarefilename = 'merged_MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches_shuffled.root'
#MCrarefilename = 'MCrare_shuffled_wRMVA_wB2MVA_sel_BCS_branches.root'
#MCrarefilename = 'S0p0_1M_correct_wRMVA_wB2MVA_sel_BCS_branches_shuffled.root'
#MCrarefilename = 'b_kspipigam_sig_500k_wRMVA_wB2MVA_sel_BCS_branches_Puresignal_shuffled.root'
MCrarefilename = 'MCrare_full_bkg_shuffled.root'
splitdir = 'splitfiles'
MCrarefile = '%s/%s'%(splitdir, MCrarefilename)
key = 'h1'


def split(filename, key, part1, part2):

	chain = TChain(key)
	chain.Add(filename)
	
	num_entries = chain.GetEntries()
	print('number of events = ', num_entries)

	nentries1 = num_entries * (part1/(part1+part2)); nentries1 = int(nentries1)
	nentries2 = num_entries * (part2/(part1+part2)); nentries2 = int(nentries2)


	print('Entries split as %s, %s'%(nentries1, nentries2))

	#Create 'splitfiles' directory if it does not exist
	if not os.path.exists('splitfiles'):
		os.makedirs('splitfiles')

	#Create new files and trees
	filename_truncated = filename[:len(filename) - 5]
	#part-1
	newfilename1 = '%s_%sPartsIn%s.root'%(filename_truncated, part1, part1+part2)
	newfile1 = TFile(newfilename1, 'recreate'); newtree1 = chain.CloneTree(0)
	for i in range(nentries1):
		chain.GetEntry(i)
		newtree1.Fill()
	newfile1.Write(); newfile1.Close()
	#part-2
	newfilename2 = '%s_%sPartsIn%s.root'%(filename_truncated, part2, part1+part2)
	newfile2 = TFile(newfilename2, 'recreate'); newtree2 = chain.CloneTree(0)
	for i in range(nentries1, num_entries):
		chain.GetEntry(i)
		newtree2.Fill()
	newfile2.Write(); newfile2.Close()
	print('Finished splitting files')


def clean_output(MCrarefile):

	if os.path.isfile(MCrarefile):
		os.remove(MCrarefile)
	print('Finished cleaning intermediate file')


if __name__ == '__main__':

	#Randomize rare-MC input before splitting it
	read_shuffle_store(do_shuffle = True, inputFilenames = InputFiles, key = key, shuffledFilename = MCrarefile)

	#Split MC-rare into two parts : 49 + 1
	split(MCrarefile, key, 49, 1)

	#Clean intermediate shuffled file
	#clean_output(MCrarefile)





