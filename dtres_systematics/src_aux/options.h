#pragma once
/**
 * @file	options.h
 * @date	Aug 15, 2011
 * @author	mprim
 * @brief	brief description
 *
 * long description
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>
#include <vector>

class Options {
public:
	Options(bool file_required = true) : m_file_required(file_required) { };
	~Options() { };

	/**
	 * @brief Parse the program options and set all members
	 * @param argc
	 * @param argv
	 * @return true if success, false if error occured
	 */
	bool ParseOptions(int argc, char *argv[]);

	/**
	 * @brief Prints current options
	 */
	void PrintOptions();

    std::string GetSMCSample() const;
    std::string GetConsolidatedResultsDir() const;
    std::string GetPlotname() const;
    std::string GetResultname() const;
    std::vector<std::string> GetFilenames() const;
    bool GetDoLifetimeFit() const;
    unsigned int GetNumcpu() const;
    int GetToyIndex() const;
    float GetAcp_Ainput() const;
    float GetAcp_Sinput() const;
    int GetrBin() const;

private:

	std::string m_SMCSample;
	std::string m_consolidated_results_dir;
	std::string m_plotname;
	std::string m_resultname;
	bool m_file_required;
	std::vector<std::string> m_filenames;
	bool m_dltf;
	unsigned int m_numcpu;
	int m_toyindex;
	float m_Ainput;
	float m_Sinput;
	int m_rBin;
};

#endif /* OPTIONS_H_ */
