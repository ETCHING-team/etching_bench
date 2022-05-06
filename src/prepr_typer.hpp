/**
 * @file prepr_typer.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.12
 * @date 2022-05-06
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef PREPR_TYPER_HPP
#define PREPR_TYPER_HPP

#include "prepr_vcf.hpp"

/**
 * @brief Typing SVs into DEL, DUP, INV, or TRA
 * 
 * @param [in] infile
 * @param [in,out] outfile
 * @return int
 * 
 * @details Print outfile after sv typing
 */
int typer ( std::string infile , std::string outfile );

#endif
