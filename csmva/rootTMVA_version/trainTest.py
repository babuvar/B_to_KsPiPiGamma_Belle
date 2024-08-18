from ROOT import TMVA, TFile, TChain, gROOT, gInterpreter, Float_t, TString
from array import array

##################################################################################################
# Configuration area
##################################################################################################

#Flags to call functions
do_traintest = False
do_addMVAvar = True

#common
rootdir = '/home/varghese/Research_KEK/IPHC/Analysis/RootFiles/MyRootFiles'
outdir = '/home/varghese/Research_KEK/IPHC/Analysis/RootFiles/MyRootFiles' #same as rootdir by default
key = 'h1'
variables = ['coscms', 'costh', 'k0et', 'k0hso00', 'k0hso02', 'k0hso04', 'k0hso10', 'k0hso12', 'k0hso20', 'k0hso22', 'k0hso24']


#for train_test_MVA
sigFile_list = ['%s/puresignal_b_spipigam_500k.root'%rootdir]
#sigFile_list = ['%s/full_b_spipigam_500k.root'%rootdir, '%s/merged_GMC_st0_charged.root'%rootdir, '%s/merged_GMC_st0_mixed.root'%rootdir]
bkgFile_list = ['%s/merged_GMC_st0_charm.root'%rootdir, '%s/merged_GMC_st0_uds.root'%rootdir]
performance_rootfile = 'performance.root'
w_s = 1 # signal weight
w_b = 1 # background weight


#for add_addMVA_variable
addMVA_inputFiles = ['merged_SMC_500K.root', 'merged_GMC_st0_uds.root', 'merged_GMC_st0_charm.root', 'merged_GMC_st0_charged.root', 'merged_GMC_st0_mixed.root'] #Note that these may be independent of testing/training samples
xmlFile = 'dataset/weights/BDT_BDT.weights.xml'


##################################################################################################
# Function definitions
##################################################################################################


def train_test_MVA(do_traintest, variables, sigFile_list, bkgFile_list, key, performance_rootfile, w_s, w_b):

	if do_traintest == False:
		return

	TMVA.Tools.Instance();

	print('Starting MVA training-testing')
	fo = TFile(performance_rootfile,'recreate')

	factory = TMVA.Factory('BDT', fo,'!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D')
	dataloader = TMVA.DataLoader('dataset')


	#Variables for training/testing (same set of 11 vars used by S.Bilokin)
	for variable in variables:
		dataloader.AddVariable(variable)


	#Signal and background events for training/testing
	#Signal
	print('reading signal files')
	sig_chain = TChain(key)
	for sigFile in sigFile_list:
		sig_chain.Add(sigFile)
	sig_chain.Merge('all_sig.root')
	ts = sig_chain.GetTree()
  
	#Background
	bkg_chain = TChain(key)
	print('reading background files')
	for bkgFile in bkgFile_list:
		bkg_chain.Add(bkgFile)
	bkg_chain.Merge('all_bkg.root')
	tb = bkg_chain.GetTree()

 
	#the 1's are event weights (setting both to 1 ignores that)
	dataloader.AddSignalTree(ts, w_s)
	dataloader.AddBackgroundTree(tb, w_b)


	dataloader.PrepareTrainingAndTestTree('','','nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V')

	#AdaBoost BDT
	#factory.BookMethod(dataloader,TMVA::Types::kBDT,'BDT','!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20')

	#Gradient boost BDT
	factory.BookMethod(dataloader,TMVA.Types.kBDT,'BDT','!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2')

	#MLP neural net
	#factory->BookMethod(dataloader,TMVA::Types::kMLP,'MLP','H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator')

	#train
	factory.TrainAllMethods()
	factory.TestAllMethods()
	factory.EvaluateAllMethods()

	#clean up
	fo.Close()
	del factory
	del dataloader

	TMVA.TMVAGui(performance_rootfile) # Meant to use when running in desktop/laptop. If used in a remote (ssh -Y ...) session, it may be very slow.
	print('Finished MVA training-testing')
	
	return

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def add_addMVA_variable(do_addMVAvar, variables, addMVA_inputFiles, xmlFile, rootdir, outdir, key):

	if do_addMVAvar == False:
		return

	print('Starting to add MVA variable to specified rootfiles')
	
	for addMVA_inputFile in addMVA_inputFiles :
	
		fullname = '%s/%s'%(rootdir, addMVA_inputFile)
		name_truncated = addMVA_inputFile[:len(addMVA_inputFile) - 5]
		outname = '%s/%s_wRMVA.root'%(outdir, name_truncated)
		finp = TFile(fullname,"READ"); tinp = finp.Get(key)
		fout = TFile(outname, "RECREATE"); tout = tinp.CopyTree("", "")
	
		
		reader = TMVA.Reader("!Color:!Silent");
	
		branches = {}
		for branchName in variables:
			branches[branchName] = array('f', [-999])
			reader.AddVariable(branchName, branches[branchName])
			tinp.SetBranchAddress(branchName, branches[branchName])
		
		newvar=array('f',[0]); newbranch = tout.Branch('rMVA', newvar, 'rMVA/F')	
		reader.BookMVA('BDT', TString(xmlFile))
		
		for i in range(tinp.GetEntries()):
			tinp.GetEntry(i); 
			newvar[0] = reader.EvaluateMVA('BDT')
			newbranch.Fill()

		fout.Write();
		finp.Close(); fout.Close()
		
	print('Finished adding MVA variable to specified rootfiles')
	return

##################################################################################################
# Main program
##################################################################################################

if __name__ == '__main__':
	
	#read signal & background files and implement the ROOT-TMVA training and testing
	train_test_MVA(do_traintest, variables, sigFile_list, bkgFile_list, key, performance_rootfile, w_s, w_b)

	#read list of root files and output rootfiles with the MVA variable added (in the same directory) 
	add_addMVA_variable(do_addMVAvar, variables, addMVA_inputFiles, xmlFile, rootdir, outdir, key)



##################################################################################################

















