/**
 * @file merge.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.11
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef MERGE_HPP
#define MERGE_HPP

#include "args.hpp"

void merge (args CONF);


/**
 * @brief Merge for silver standard without truth set
 * 
 * @note Skip this, if truth set (or golden standard) is given by user.
 * 
 * @param [in,out] CONF
 * 
 * @details Hidden parameters in CONF:
 * @details merge_window: The SVs overlapped in this window is treated as the same SV.
 * @details truth_cutoff (default: CONSENSUS_CUTOFF)
 * @details min_sv_size (default: MIN_SV_SIZE)
 * @details initial_merge.outfile: output file
 * @details initial_merge.list: list of input files
 * @details initial_merge.cutoff = truth_cutoff
 * 
 * @details Steps in this function:
 * @details 1. Input: ${prefix}.initial.list: a list file of preprocessed vcf files
 * @details 2. Output: ${prefix}.truth.vcf: Merged SV calls with comine_calls_svs of SURVIVOR merge
 * 
 * @return ${outfile_prefix}.truth.vcf
 */
void merge_for_silver_standard ( args CONF ) ;


/**
 * @brief 
 * 
 */
void apply_truth_cutoff(STR outfile, int truth_cutoff);

/**
 * @brief Add calls to ID
 * 
 * @param [in,out] &CONF
 * 
 * @details Steps of this function:
 * @details 1. Read truth set: ${outfile_prefix}.truth.vcf for silver standard, or given by user for golden standard.
 * @details 2. Rename ID to TRUTH${ID}
 * @details 3. Remove mates
 * @details 4. Make file name (${prefix}.benchmark.annotated.vcf)
 * @details 5. Annotate truth set
 * 
 * @return ${outfile_prefix}.truth.annotated.vcf: Annotated truth file 
 */
void add_binary_calls_to_id ( args CONF ) ;

void add_binary_calls_to_id ( STR infile, STR outfile );



/**
 * @brief Merge for golden standard with given truth set
 * 
 * @param CONF 
 */
void merge_with_truth_set ( args CONF ) ;





/** 
 * @brief Modified combine_calles_svs in SURVIVOR/src/merge_vcf/combine_svs.cpp
 * @detauks Equivalent to SURVIVOR merge (except for print_log)
 * @details Merging SV module of SURVIVOR merge
 * @details print_log: true if you want to print log
 * @param [in] files : list of vcf files
 * @param [in] max_dist : maximum distance between two SVs (resolution)
 * @param [in] min_support : minimum number of supporting reads
 * @param [in] type_save : type of SVs to save
 * @param [in] strand_save : Strand aware
 * @param [in] dynamic_size : dynamic size (1 or 0)
 * @param [in] min_svs : minimum size of SVs 
 * @param [in,out] output : output file
 * @param [in] print_log : 1 if you want to print log (introduced)
 */
void combine_calls_svs(std::string files, double max_dist, int min_support, int type_save, int strand_save, int dynamic_size, int min_svs, std::string output, bool print_log);


#endif
