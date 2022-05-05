/**
 * @file args.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.11
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef ARGS_HPP
#define ARGS_HPP

#include "DEBUG_BLOCKS.hpp"
#include "define.hpp"
#include <getopt.h>
#include <math.h>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

/**
 * @brief Struct for a vcf file
 */
struct source_info{
  STR tool;
  STR infile;
  double cutoff;
  bool noimprecise;
  STR outfile;
  STR cutoff_file;
  STR list;
};


/**
 * @class args args.hpp "args.hpp"
 * @brief A class for configuration of the program
 * 
 * @details All options, input vcf files, output prefix, etc. are stored in this class.
 */
class args{
private:
  /** 
   * @var STR default_source
   * @brief Default settings of supporting tools
   * @details Updated in intializaion
   */
  std::map<STR,source_info > default_source;

  /** @brief A set of help strings */
  std::set< STR > help_set;

  
  /** 
   * @var int remove_own_prediction
   * @brief 1 to remove own prediction, 0 for else
   */  
  int remove_own_prediction;

  /** 
   * @var int keep_own_prediction
   * @brief 1 to keep own prediction, 0 for else
   * @note If remove_own_prediction == 1 , keep_own_prediction will be ignored.
   */  
  int keep_own_prediction;

  /**
   * @var bool any_help
   * @brief Check if any help is given
   * 
   * @param argc 
   * @param argv 
   * @return true for any help is given, false otherwise
   */
  bool any_help(int argc, char ** argv);

  /**
   * @brief Reorder source vector
   * 
   * @details Following order:
   * 
   * 0. truth (after update or if provided)
   * 1. etching
   * 2. delly
   * 3. lumpy
   * 4. manta
   * 5. svaba
   * 6. novobreak
   * 7. gridss
   * 
   * Any tool not in the list will be ignored.
   */
  void reorder();

  /**
   * @brief convert uppercase to lowercase
   * @param str tool name
   */
  STR to_lower_case(STR str);

  /**
   * @brief Preprocessed file name of input vcf file ( )
   * @ref source_info and @ref source
   * 
   * @details add .prep.vcf to the end of the file name
   * 
   * @param [in] source.infile
   * @return ${source.infile_prefix}.prep.vcf
   */
  STR preprocessed_file_name (STR infile );


  /**
   * @brief VCF file name after cutoff applied
   * @ref source_info and @ref source
   * 
   * @details add .cutoff.vcf to the end of the file name
   * 
   * @param [in] source.infile
   * @return ${source.infile_prefix}.cutoff.vcf
   */
  STR cutoff_file_name (STR file);


  /** 
   * @brief Initialize
   * @ref config_file 
   */
  void init();

  /**
   * @brief Read options from argv
   * 
   * @param argc 
   * @param argv 
   */
  void read_options(int argc, char ** argv);

  /**
   * @brief Check if the required options are set
   * 
   * @details Required options are:
   * @details -c, --config for config file
   * @details -o, --outfile for outfile prefix
   */
  void check_required();

  /**
   * @brief Check if the tool is supported
   * 
   * @param [in] tool name
   * @details Exit if the tool is not supported.
   */
  void is_supporting(STR &tool);


  /**
   * @brief Update merge information
   * 
   * @details Merge information:
   * @details source_info initial_merge
   * @details source_info truth
   * @details source_info benchmark
   */
  void update_merge_info();

  /**
   * @brief Check if input string is number
   *
   * @param [in] std::string input string
   * @return bool 1 for number, 0 for else
   */
  bool isNumber(const std::string s);

public:
  /** @brief Constructor */
  args();
    
  /** @brief Destructor */
  ~args();

  /**
   * @brief Constructor with arguments
   * @param [in] "int argc" number of arguments
   * @param [in] "char *argv[]" arguments
   */
  args ( int argc , char ** argv );

  /**
   * @var STR config_file
   * @brief configuration file for the source array
   * 
   * @details Format:
   * @details tool_name    file_name.vcf    [cutoff]
   * 
   * @note cutoff is optional. 
   */
  STR config_file;

  /** @brief Minimum SV size */
  int min_sv_size;
  int max_sv_size;
  int intra_SV;
  int inter_SV;

  /** @brief Output file prefix */
  STR prefix;

  /** @brief Merge window */
  int merge_window;

  /** @brief Consensus cutoff */
  double truth_cutoff;

  /** @brief Truth file */
  STR truth_file;

  /**
   * @ brief 1 for maximizing f1 score, else 0 [default 0]
   *
   */
  bool maximizing_f1;
  
  /**
   * @var std::vector<source_info> source
   * @brief Source files
   * 
   * @details Format example:
   * @details source[1].tool = "etching"
   * @details source[1].infile = "etching_input.vcf"
   * @details source[1].outfile = "prefix.etching.prep.vcf"
   * @details source[1].cutoff_file = "prefix.etching.cutoff.vcf" (This is for initial merge)
   * @details source[1].cutoff = 0.5
   */
  std::vector<source_info> source;

  /**
   * @var source initial_merge
   * @brief Initial merge for silver standard
   * 
   * @details Format example:
   * @details initial_merge.tool = "survivor"
   * @details initial_merge.outfile = "${prefix}.merge.vcf"
   * @details initial_merge.cutoff = truth_cutoff
   * @details initial_merge.list = "${prefix}.merge.list"
   * 
   * @note ${prefix}.merge.list is a list of source.outfile
   * @note ${prefix}.merge.vcf is a SURVIVOR merge file with min. support
   */
  source_info initial_merge;

  /**
   * @var source_info truth
   * @brief truth (annotated) set information
   * 
   * @details truth.tool = "truth"
   * @details truth.infile = ${prefix}.merge.vcf (initial merge)
   * @details truth.outfile = ${prefix}.truth.vcf (annotated truth)
   */
  source_info truth;   

  /**
   * @var source_info benchmark
   * @brief 
   * 
   * @details benchmark.tool = "benchmark"
   * @details benchmark.infile = ${prefix}.benchmark.vcf
   * @details benchmark.outfile = ${prefix}.benchmark.annotated.vcf
   * @details benchmark.list = ${prefix}.benchmark.list
   * 
   * @note benchmark.infile is a list of truth.outfile and source.outfile
   */
  source_info benchmark;

  /**
   * @brief Read arguments
   * 
   * @param [in] argc Number of arguments
   * @param [in] argv Arguments
   * 
   * @details 
   * 1. Check if any help string was given.
   * 2. Read options
   * 3. Check if required options were given.
   * 4. Read config file given by -c option.
   * 5. Update all file names
   * 
   * @return Update args objects
   */
  void read_args ( int argc, char ** argv );

  /**
   * @var int is_silver
   * @brief 1 for silver standard, 0 for a given true set (golden standard, or simulation)
   */
  int is_silver;


  /**
   * @var int precise_only
   * @brief 1 to remove imprecise (precise only), or 0 to keep imprecise
   *
   */
  int precise_only;

  /**
   * @var int remove_self_bias
   * @brief 1 to remove the bias from own prediction
   *
   */
  int remove_self_bias;

  /**
   * @var int gen_TRA_mate
   * @brief 1 to generate a mate for TRA, or 0 for else
   *
   */
  std::string reuse_anno_vcf;
    
  /**
   * @brief read the file of input vcf files, and store in class args
   * 
   * @details
   * config file format (tab- or space-delimited):
   * tool_name     input_file     [cutoff  [noimprecise] ]
   * 
   * Available tool names: etching, delly, lumpy, manta, svaba, novobreak, and gridss.
   * If cutoff was not given, it will be set to the default value of the tool.
   * 
   * @note Other tools are not supported, yet.
   */
  void read_conf_file ();

  /** @brief Print USAGE and exit(1) */
  void usage();

  /** @brief Print source information mainly for debugging */
  void print_source();

  /** @brief Print settings, mainly for debugging */ 
  void print_settings();

  /** @brief return begin of source */
  std::vector<source_info>::iterator begin();

  /** @brief return end of source */
  std::vector<source_info>::iterator end();
    
  /** @brief Operator for iterating over source */
  source_info & operator [](unsigned int i);

  /** @brief find "tool" over source */
  std::vector<source_info>::iterator find(STR tool);

  /** @brief iterator */
  using iterator=std::vector<source_info>::iterator;

  /** @brief const iterator */
  using const_iterator=std::vector<source_info>::const_iterator;

  /** @brief Number of files to be merged */
  std::size_t size();

  /** @brief Make initial merge list file */
  void make_initial_merge_list_file();

  /** @brief Make benchmark list file */
  void make_benchmark_list_file();
};

#endif
