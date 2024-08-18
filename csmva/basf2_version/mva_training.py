from root_pandas import read_root
from ROOT import gROOT, TFile, TTree, TCanvas
import os.path
import sys
from basf2 import *
import sysconfig
import os
gROOT.ProcessLine(".include " + sysconfig.get_path("include"))
import basf2_mva
import basf2_mva_util
from array import array

##################################################################################################
# Configuration area
##################################################################################################

#Flags to call functions
do_shuffle = False
do_split = False
do_trainMVA = False
do_addMVAvar = True

#For read_shuffle_store
key = 'h1'
rootdir = 'rootfiles_copy_forMVA'
outdir = 'rootfiles_copy_forMVA' # same as rootdir by default
inputFilenames = ['%s/merged_GMC_st0_charm.root'%rootdir, '%s/merged_GMC_st0_uds.root'%rootdir, '%s/puresignal_b_spipigam_500k.root'%rootdir]
shuffledFilename = '%s/allCSMVA_shuffled.root'%rootdir

#For split_TestTrain
ratio = 0.5
TestingFilename = '%s/mva_Testing.root'%rootdir
TrainingFilename = '%s/mva_Training.root'%rootdir
trainVars = ['coscms', 'costh', 'k0et', 'k0hso00', 'k0hso02', 'k0hso04', 'k0hso10', 'k0hso12', 'k0hso20', 'k0hso22', 'k0hso24']

#For train_MVA
mva_outdir = 'mva_output'
xmlname = 'MyContinuumFBDT'
pdfname = 'evalTest'

#For add_addMVA_variable
addMVA_inputFiles = ['merged_SMC_500K.root', 'merged_GMC_st0_uds.root', 'merged_GMC_st0_charm.root', 'merged_GMC_st0_charged.root', 'merged_GMC_st0_mixed.root'] #Note that these may be independent of testing/training samples
#addMVA_inputFiles = ['merged_GMC_st0_mixed.root']
xmlFile = '%s/%s.xml'%(mva_outdir, xmlname)

##################################################################################################
# Function definitions
##################################################################################################

def read_shuffle_store(do_shuffle, inputFilenames, key, shuffledFilename):

	if do_shuffle == False:
		return

	print("Starting reading and shuffling input root files")
	dataframe = read_root(inputFilenames, key=key)
	dataframe = dataframe.sample(frac=1).reset_index(drop=True)
	dataframe.to_root(shuffledFilename, key = key, store_index=False)	
	print("Finished shuffling. Output stored at %s"%shuffledFilename)
	return 

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def split_TestTrain(do_split, shuffledFilename, key, ratio, TestingFilename, TrainingFilename, trainVars): 

	if do_split == False: 
		return 

	print("Starting to split into testing and training samples")
	nancuttemp = '!TMath::IsNaN(%s)'
	andtemp = ' && '
	nancut = ''
	for i in range(len(trainVars)):
		nancut += nancuttemp%trainVars[i]
		if(i < len(trainVars)-1):
			nancut += andtemp
	fshuffled = TFile(shuffledFilename, "READ")
	shuffled_tree = fshuffled.Get(key)
	ftrain = TFile(TrainingFilename, "RECREATE")
	training_tree = shuffled_tree.CopyTree(nancut, "", nentries = int(shuffled_tree.GetEntries()*ratio))
	ftrain.Write(); ftrain.Close()
	ftest = TFile(TestingFilename, "RECREATE")
	testing_tree = shuffled_tree.CopyTree(nancut, "", nentries = int(shuffled_tree.GetEntries()*(1.0-ratio)), firstentry = int(shuffled_tree.GetEntries()*ratio))
	ftest.Write(); ftest.Close()
	print("MC-data has been split into a testing sample:%s and a training sample:%s"%(TestingFilename, TrainingFilename)) 
	return 

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def train_MVA(do_trainMVA, TestingFilename, TrainingFilename, trainVars, mva_outdir, xmlname, pdfname):
	
	if do_trainMVA ==False:
		return

	print("Starting MVA training and testing")
	# Global/General options
	go = basf2_mva.GeneralOptions()
	go.m_datafiles = basf2_mva.vector(TrainingFilename)  # training sample
	go.m_treename = key # ntuple tree name
	go.m_identifier = "%s/%s.xml"%(mva_outdir, xmlname) # name of the file with the training info

	#go.m_identifier = "MVADatabaseIdentifier"
	go.m_variables = basf2_mva.vector(*trainVars) # input variables
	go.m_target_variable = "issignal" # target for training
	go.m_weight_variable = ""

	# Specific options
	sp = basf2_mva.FastBDTOptions() # here we use the FastBDT method
	sp.transform2probability = False;
	#sp.m_nTrees = 100 # number of trees in the FastBDT forest
	#sp.m_nCuts = 8 # number of cuts for each tree
	#spbdt.m_shrinkage = 0.2 # i dont know...
	#sp.m_nLevels = 4 # number of levels of a specific tree
	#sp.m_nTrees = 100 # number of trees in the FastBDT forest
	#sp.m_nCuts =  # number of cuts for each tree
	#sp.m_nLevels = 4 # number of levels of a specific tree

	tmva_bdt_options = basf2_mva.TMVAOptionsClassification()
	#tmva_bdt_options.m_method = "Fisher"
	tmva_bdt_options.m_method = "BDT"
	#tmva_bdt_options.m_type = "Fisher"
	tmva_bdt_options.m_type = "BDT"
	#tmva_bdt_options.m_config = ("H:V:CreateMVAPdfs:NTrees=1000:nCuts=16:MaxDepth=10:VarTransform=D,G,G")
	#tmva_bdt_options.m_config = ("H:V:CreateMVAPdfs:VarTransform=P:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10")
	tmva_bdt_options.m_config = ("H:V:CreateMVAPdfs:VarTransform=D,G,P:NTrees=1000:nCuts=16:MaxDepth=10")
	tmva_bdt_options.transform2probability = False;

	basf2_mva.teacher(go,tmva_bdt_options)
	#basf2_mva.teacher(go,sp)
	outputPdf = '%s/%s.pdf'%(mva_outdir, pdfname)
	#os.system("basf2 mva_helper.py -id %s -tree %s -train %s -data %s -out %s"%(go.m_identifier + ' mva-addition/_old_MyContinuumFBDT.xml', go.m_treename, TrainingFilename, TestingFilename, outputPdf))
	os.system("basf2 mva_helper.py -id %s -tree %s -train %s -data %s -out %s"%(go.m_identifier, go.m_treename, TrainingFilename, TestingFilename, outputPdf))
	print("Finished MVA training and testing")
	return

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def add_addMVA_variable(do_addMVAvar, addMVA_inputFiles, xmlFile, rootdir, outdir, key):

	if do_addMVAvar == False:
		return

	print("Starting to add MVA variable to specified rootfiles")

	for addMVA_inputFile in addMVA_inputFiles:
		print("Creating expert-root-file for %s"%addMVA_inputFile)
		inputRootFile = '%s/%s'%(rootdir, addMVA_inputFile); 
		iRF_truncated = addMVA_inputFile[:len(addMVA_inputFile) - 5]
		expertRootFile = '%s/%s_expert.root'%(outdir, iRF_truncated)
		outputRootFile = '%s/%s_wB2MVA.root'%(outdir, iRF_truncated)

		basf2_mva.expert(basf2_mva.vector(xmlFile), basf2_mva.vector(inputRootFile), key, expertRootFile)
		fexp = TFile(expertRootFile,"READ"); texp = fexp.Get("variables")
		finp = TFile(inputRootFile,"READ"); tinp = finp.Get('h1')
		fout = TFile(outputRootFile, "RECREATE"); tout = tinp.CopyTree("", "")
	
		branch = texp.GetListOfBranches().At(0).GetName(); mva=array('f',[0]); texp.SetBranchAddress(branch, mva)
		newvar=array('f',[0]); newbranch = tout.Branch('b2MVA', newvar, 'b2MVA/F')
		print("Adding branches...")
		for i in range( texp.GetEntries() ):
			texp.GetEntry(i); newvar[0] = mva[0]; newbranch.Fill()
		#tout.Write(); 
		fout.Write(); 
		finp.Close(); fexp.Close(); fout.Close()

	print("Finished adding MVA variable to specified rootfiles")
	
##################################################################################################
# Main program
##################################################################################################

if __name__ == '__main__':
	
	#read input files and make a combined and shuffled output file
	read_shuffle_store(do_shuffle, inputFilenames, key, shuffledFilename)

	#split into testing and training samples with a given ratio
	split_TestTrain(do_split, shuffledFilename, key, ratio, TestingFilename, TrainingFilename, trainVars)

	#train and test MVA
	train_MVA(do_trainMVA, TestingFilename, TrainingFilename, trainVars, mva_outdir, xmlname, pdfname)

	#read list of root files, output rootfiles with the MVA variable added (in the same directory) 
	add_addMVA_variable(do_addMVAvar, addMVA_inputFiles, xmlFile, rootdir, outdir, key)


##################################################################################################








