#ifndef __INPUT_OPTIONS__
#define __INPUT_OPTIONS__
/**
 *  @file   input_options.h
 *  @brief  Utility functions for input options processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include "datatypes.h"
#include "utilities.h"
#include <getopt.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <regex>

using namespace std;


/**
 * @brief Reads the options to the fastqprocess program
 *
 * @param argc  no of arguments to the main function
 * @param argv arguments array to the main function
 * @param options the structure for holding the options for getopt
*/
void read_options_fastqprocess(int, char **, INPUT_OPTIONS_FASTQPROCESS &);


/**
 * @brief Reads the options to the tagsort program
 *
 * @param argc  no of arguments to the main function
 * @param argv arguments array to the main function
 * @param options the structure for holding the options for getopt
*/
void read_options_tagsort(int, char **, INPUT_OPTIONS_TAGSORT &);

#endif
