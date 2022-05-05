/**
 * @file prepr.cpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.11
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022 Bioinformatic and Genomics Lab. Hanyang University
 * 
 */

#ifndef PREPR_CPP
#define PREPR_CPP

#include <chrono>
#include "prepr.hpp"
#include "prepr_vcf.hpp"

void prepr ( args CONF ){
  if ( CONF.reuse_anno_vcf.empty() ){
    CERR << "[Preprocessing]" << std::endl;
    for ( int i = 0; i < CONF.source.size(); i++ ){
      prepr(CONF,i);
    }
  }

}


/**
 * @brief For a vcf file
 * 
 */

// void prepr (args CONF, int i){
//     VCF VCF(CONF[i]);
void prepr (args CONF, int i){
  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  STR TIME = std::ctime ( & end_time );
  TIME.pop_back();
  CERR << "Date: " << TIME << "\n";

  VCF VCF(CONF,i);
  // add missed meta information
  VCF.add_missed_header_lines();
  // Preprocessing for Delly, Lumpy, Manta, and novoBreak
  // Add score to id for all tools
  VCF.prepr_tools ( );
  // typing
  VCF.typing_SV ();
  // remove mate SV
  VCF.remove_mate ( );

  // remove short SVs (<min_sv_size)
  VCF.remove_short( CONF.min_sv_size );

  // remove large SVs (>=max_sv_size)
  VCF.remove_large( CONF.max_sv_size );

  // remove large SVs (>=max_sv_size)
  VCF.remove_inter_chromosome( CONF.intra_SV );

  // remove large SVs (>=max_sv_size)
  VCF.remove_intra_chromosome( CONF.inter_SV );

  // rearrange BP
  VCF.rearrange_bp();
  
  VCF.add_score_to_id();

  int benchmarking_size = VCF.size();

  // write VCF file without FORMATs
  VCF.fwrite_short( );
  VCF.remove_low_qual_records();
    
  VCF.fwrite_short(VCF.cutoff_file);

  int above_cutoff_size = VCF.above_cutoff_size();

  std::string prepr_note = "prepared for benchmarking: SV-typing, and TRA mate removed";
  // if ( ! VCF.gen_TRA_mate ){
  //     prepr_note += ", remove TRA mate";
  // }
  CERR << VCF.tool_name << "\t" << VCF.outfile << "\t" << benchmarking_size << "\t" << prepr_note << "\n";
  CERR << VCF.tool_name << "\t" << VCF.cutoff_file << "\t" << above_cutoff_size << "\tprepared for silver standard: remove low-quality SVs with cutoff " << VCF.cutoff << "\n";
  CERR << "\n";
} 

#endif
