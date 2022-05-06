/**
 * @file prepr.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.12
 * @date 2022-05-06
 * 
 * @copyright Copyright (c) 2022 Bioinformatic and Genomics Lab. Hanyang University
 * 
 */

#ifndef PREPR
#define PREPR

#include "DEBUG_BLOCKS.hpp"

#include "prepr_vcf.hpp"
#include "prepr_typer.hpp"
#include "args.hpp"

void prepr (args CONF);

void prepr (args CONF, int i);

void apply_cutoff (args CONF, int i);

#endif
