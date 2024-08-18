#pragma once
/**
 * @brief Header file with some useful utilities
 *
 * Can be included into any file/project, no extra implementation file necessary
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <stdlib.h>

namespace util {

/**
 * @brief Returns a string of t if t overloads << operator
 *
 * @param t will be converted to a string
 * @return string of t
 */
template <class T>
inline std::string to_string (const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}

/**
 * @brief Returns string representation of var
 *
 * @param var boolean to be converted
 * @return "true" or "false"
 */
inline std::string bool_string(bool var) {
    if(var)
        return "true";
    else
        return "false";
}

/**
 * @brief Creates a std::string representation of std::vector<T> if T overloads << operator
 *
 * @param vec vector to convert
 * @param newline true if newline should be used for each value
 * @return string of vec
 */
template <class T>
inline std::string vector_to_string(const std::vector<T> &vec, bool newline = true) {
    std::stringstream ss;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        ss << "Vector[" << i << "] " << vec.at(i);
        if(newline)
            ss << std::endl;
        else
            ss << " ";
    }
    return ss.str();
}

/**
 * @brief Creates a std::string representation of std::pair<T1,T2> if T1 and T2 overload << operator
 *
 * @param pair pair to convert
 * @return string of pair
 */
template <class T1, class T2>
inline std::string pair_to_string(const std::pair<T1, T2> &pair) {
    std::stringstream ss;
    ss << "Pair[" << pair.first << "," << pair.second << "]";
    return ss.str();
}

/**
 * @brief Calculate mean of a given vector
 * @param vec values
 * @return mean
 */
template <class T>
inline double calc_mean(const std::vector<T> &vec){
	double sum = 0;
	for(unsigned int j = 0; j < vec.size();j++){
		sum += vec[j];
	}
	sum /= vec.size();
	return sum;
}

/**
 * @brief Calculate sigma of a given vector
 * @param vec values
 * @return sigma
 */
template <class T>
inline double calc_sigma(const std::vector<T> &vec){
	double sum = 0;
	double mean = util::calc_mean(vec);

	for(unsigned int j = 0; j < vec.size();j++){
		sum += (vec[j] - mean)*(vec[j] - mean);
	}
	sum /= vec.size()-1;
	return std::sqrt(sum);
}

/**
 * @brief Calculate weighted mean of a given vector
 * @param vec values
 * @param err weights
 * @return weighted mead
 */
template <class T, class S>
inline double calc_weighted_mean(const std::vector<T> &vec, const std::vector<S> &err){
	double sum = 0;
	double weight;
	double sumweight = 0;

	for(unsigned int j = 0; j < vec.size();j++){
		weight = 1.0/(err[j]*err[j]);
		sum += weight*vec[j];
		sumweight += weight;
	}
	sum /= sumweight;
	return sum;
}

/**
 * @brief Executes command cmd and in case of error, prints command and exists with EXIT_FAILURE code
 * @param cmd command to execute
 */
inline void execute_command(const std::string &cmd) {
	int err = system(cmd.c_str());
	if(err != 0) {
		std::cerr << "ERROR: Failed to execute command: " << cmd << std::endl;
		exit(EXIT_FAILURE);
	}
}

} // namespace util

#endif /* UTILITIES_H_ */
