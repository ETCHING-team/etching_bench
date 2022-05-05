/**
 * @file DEBUG_BLOCKS.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.11
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef DEBUG_BOCKS
#define DEBUG_BOCKS

#include <stdio.h>
#include <iostream>

#define DEBUG_MAIN 
#define DEBUG_PREPR_VCF
#define DEBUG_PREPR
#define DEBUG_PREPR_DELLY
#define DEBUG_TYPING_SV
#define DEBUG_REMOVE_MATE
#define DEBUG_REARRANGE_BP
#define DEBUG_FWRITE_SHORT
#define DEBUG_MERGE
#define DEBUG_MERGE_FOR_BENCH
#define DEBUG_MAKE_TRUTH_FILE
#define DEBUG_SET_ANNOTATED_TRUTH
#define DEBUG_COMBINE_CALLS_SVS
#define DEBUG_GET_FORMAT_SIZE
#define DEBUG_ADD_BINARY_CALLS_TO_ID

#define DEBUG_MERGE
#define DEBUG_BENCH

#define CHECK printf("CHECK...%s(%d)::%s\n", __FILE__, __LINE__, __FUNCTION__);
#define DEBUG_COUT(STR)    std::cout << __FILE__ << "(" << __LINE__ << ")::"<< __FUNCTION__ << "::" << STR << std::endl 
#endif
