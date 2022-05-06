/**
 * @file prepr_typer.cpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.12
 * @date 2022-05-06
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef PREPR_TYPER_CPP
#define PREPR_TYPER_CPP

#include "prepr_typer.hpp"

int typer ( std::string infile , std::string outfile ){

  VCF container;
  VCF container_SV;

  // Read VCF
  container.read_vcf_file ( infile );

  // add missed meta information
  container.add_missed_header_lines();

  // Typing SVs
  container_SV = container;
  container_SV.typing_SV();

  // print result
  container_SV.fwrite(outfile);

  return 0;
}

#endif
