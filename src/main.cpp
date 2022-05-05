/**
 * @file main.cpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief A program for preparing vcf files for benchmarking 
 * @brief This program supports ETCHING, DELLY, LUMPY, Manta, SvABA, novoBreak, and GRIDSS.
 * @version 0.2.11
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

/**
 * @mainpage ETCHING_BENCH: Benchmarking tool of ETCHING project
 * 
 * @section Usage
 * @subsection 
 * 
 */

#include "prepr.hpp"
#include "bench.hpp"
#include "merge.hpp"

/**
 * @brief main function
 * 
 */
int main (int argc, char ** argv){   
  /** @brief Set input arguments */
  args CONF(argc, argv);

  /** @brief Preprocessing */
  prepr(CONF);

  /** @brief Merge preprocessed vcf files */
  merge(CONF);

  /** @brief Calculate performance */
  benchmark(CONF);

  return 0;
}
