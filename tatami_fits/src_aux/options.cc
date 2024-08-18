/**
 * @file	options.cc
 * @date	Aug 15, 2011
 * @author	mprim
 * @brief	brief description
 *
 * long description
 */

#include <iostream>

#include "options.h"
// boost library
#include <boost/program_options.hpp>
namespace po = boost::program_options;
// Utilities
#include "utility.h"


bool Options::GetDoLifetimeFit() const {
    return m_dltf;
}

std::vector<std::string> Options::GetFilenames() const {
    return m_filenames;
}

unsigned int Options::GetNumcpu() const {
    return m_numcpu;
}


std::string Options::GetPlotname() const {
    return m_plotname;
}

std::string Options::GetResultname() const {
    return m_resultname;
}

std::string Options::GetConsolidatedResultsDir() const {
    return m_consolidated_results_dir;
}

int Options::GetToyIndex() const{
    return m_toyindex;
}

float Options::GetAcp_Ainput() const{
    return m_Ainput;
}

float Options::GetAcp_Sinput() const{
    return m_Sinput;
}

std::string Options::GetSMCSample() const{
    return m_SMCSample;
}

int Options::GetrBin() const{
    return m_rBin;
}

bool Options::ParseOptions(int argc, char *argv[]) {
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Produce help message)")
		("Plotname,p", po::value<std::string>(&m_plotname), "Name of plot file")
		("Resultname,r", po::value<std::string>(&m_resultname), "Name of result file")
		("filenames,f", po::value<std::vector<std::string> >()->multitoken(), "Input filename(s) or regular expression for TChain")
		("lifetime,l", po::value<bool>(&m_dltf)->default_value(false), "Do lifetime")
		("numcpu", po::value<unsigned int>(&m_numcpu)->default_value(10), "Number of CPUs used during minimization")
		("Consolidatedresultsdir,c", po::value<std::string>(&m_consolidated_results_dir), "Name of consolidated result directory")
		("toyindex,i", po::value<int>(&m_toyindex), "Index of toy-fit")
		("Acp_Ainput,a", po::value<float>(&m_Ainput), "Input Acp-A value")
		("Acp_Sinput,s", po::value<float>(&m_Sinput), "Input Acp-S value")
		("smc_sample,t", po::value<std::string>(&m_SMCSample), "SMC sample used for toy set")
		("rBin,b", po::value<int>(&m_rBin), "Input wr-bin")
		;
			

	// create variables map and run parser
	bool parser_succeeded = true;
	po::variables_map vm;
	try {
		po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
		po::notify(vm);
	} catch (po::error e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		parser_succeeded = false;
	}

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		parser_succeeded = false;
	}

	if (vm.count("filenames")) {
		m_filenames = vm["filenames"].as<std::vector<std::string> >();
	} else {
		if(m_file_required) {
			std::cout << "ERROR: At least one file required!" << std::endl;
			parser_succeeded = false;
		}
	}

	if(m_numcpu > sysconf( _SC_NPROCESSORS_ONLN )) {
		std::cout << "WARNING: NumCPU set is larger than available CPUs, value was reset to " << sysconf( _SC_NPROCESSORS_ONLN ) << std::endl;
		m_numcpu = sysconf( _SC_NPROCESSORS_ONLN );
	}


	return parser_succeeded;
}

void Options::PrintOptions() {
	std::cout << "Number of CPUs used for minimization: " << m_numcpu << std::endl;
	std::cout << util::vector_to_string(m_filenames, true) << std::endl;
}
