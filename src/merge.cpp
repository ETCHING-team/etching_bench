/**
 * @file merge.cpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief
 * @version 0.2.11
 * @date 2022-04-29
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef MERGE_CPP
#define MERGE_CPP

#include <chrono>
#include "merge.hpp"
#include "define.hpp"
#include "prepr_vcf.hpp"
#include "SURVIVOR/src/merge_vcf/combine_svs.cpp"

void merge (args CONF){
  if ( CONF.reuse_anno_vcf.empty() ){
    CERR << "[Merging]" << std::endl;
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    STR TIME = std::ctime ( & end_time );
    TIME.pop_back();
    CERR << "Date: " << TIME << "\n";
    
    if (!CONF.truth_file.size()) {
      merge_for_silver_standard(CONF);
      CONF.truth_file=CONF.initial_merge.outfile;
    }
    else{
      add_binary_calls_to_id(CONF.truth_file,CONF.truth.outfile);
      merge_with_truth_set(CONF);
    }
  }
}

  
/**
 * @brief [Merge Step 1] Initial merge for silver standard
 * 
 * @param CONF 
 */
void merge_for_silver_standard ( args CONF ){
  CONF.make_initial_merge_list_file();

  std::cerr << "merge_for_silver_standard::Equivalent to: SURVIVOR merge "
	    << CONF.initial_merge.list << " " << CONF.merge_window << " 1 1 1 0 " << CONF.min_sv_size << " " << CONF.initial_merge.outfile << "\n";

  bool log_switch = false;
  combine_calls_svs(CONF.initial_merge.list, CONF.merge_window, 1, 1, 1, 0, CONF.min_sv_size, CONF.initial_merge.outfile,log_switch);

  VCF benchmarking_vcf(CONF.initial_merge.outfile);
  VCF truth_vcf;
  truth_vcf=benchmarking_vcf;
  truth_vcf.clear();
  benchmarking_vcf.extract_id_from_SURVIVOR_merge_format();
  benchmarking_vcf.update_format_name_for_silver_standard();
  int truth_id_count=1;
  for ( auto &i : benchmarking_vcf ){
    for ( auto &j : i.second ){
      VCF_LINE vcf_line=j.second;
      vcf_line.format.push_back("");
      int count=0;
      for ( std::size_t k = vcf_line.format.size()-1 ; k > 1 ; k--){
	vcf_line.format[k] = vcf_line.format[k-1];
      }
      STR truth_tag="COMBINED"+std::to_string(truth_id_count++);
      for ( std::size_t k = 2 ; k < vcf_line.format.size() ; k++){
	STS ss(vcf_line.format[k]);
	STR tmp;
	double score;
	GETLINE ( ss, tmp, ':' );
	GETLINE ( ss, tmp, ':' );
	score = std::stod(tmp);
	if ( score >= CONF.source[k-2].cutoff ){
	  count++;
	  truth_tag += "_1";
	}
	else{
	  truth_tag += "_0";
	}
      }
      if ( count >= CONF.truth_cutoff ){
	truth_tag += ":1";
	truth_vcf[i.first][j.first]=vcf_line;
      }
      else{
	truth_tag += ":0";
      }
      vcf_line.format[1]=truth_tag;
      vcf_line.format[0]="ID:SCORE";
      j.second=vcf_line;
    }
  }
  // update header
  STS header_ss(benchmarking_vcf.header);
  VEC<STR> header_tokens;
  STR header_token;
  while ( header_ss >> header_token ){
    header_tokens.push_back(header_token);
  }
  STR new_header;
  for ( std::size_t i = 0 ; i < 9 ; i++ ){
    new_header += header_tokens[i] + "\t";
  }
  new_header += "TRUTH\t";
  for ( std::size_t i = 9 ; i < header_tokens.size() ; i++ ){
    new_header += header_tokens[i] + "\t";
  }
  new_header.pop_back();
  benchmarking_vcf.header=new_header;

  truth_vcf.fwrite(CONF.truth.outfile);
  benchmarking_vcf.fwrite(CONF.benchmark.outfile);
  CERR << "Benchmarking vcf file: " << CONF.benchmark.outfile << "\n";
  CERR << "\n";
}


void add_binary_calls_to_id ( STR infile, STR outfile ){
  VCF annotated_truth(infile); 
#ifdef DEBUG_ADD_BINARY_CALLS_TO_ID
  DEBUG_COUT ( "add_binary_calls_to_id: " << infile << "\t" << outfile << "\n");
#endif
  annotated_truth.rename_id("TRUTH"); 
  annotated_truth.typing_SV(); 
  annotated_truth.remove_mate(); 
  annotated_truth.rearrange_bp();
  annotated_truth.add_binary_calls_to_id(outfile); 
  annotated_truth.fwrite_short(outfile);
}


void merge_with_truth_set ( args CONF ){
  /// 1. Make CONF.benchmark.list: annotated_truth_file + preprocessed vcfs
  CONF.make_benchmark_list_file();

  /// 2. Merge to ${outfile_prefix}.benchmark.vcf using SURVIVOR merge 
  std::cerr << "merge_with_truth_set::Equivalent to: SURVIVOR merge "
	    << CONF.benchmark.list << " " << CONF.merge_window << " " << 1 << " 1 1 0 " << CONF.min_sv_size << " " << CONF.benchmark.infile << "\n";
  bool log_switch = true;
  combine_calls_svs(CONF.benchmark.list, CONF.merge_window, 1, 1, 1, 0, CONF.min_sv_size, CONF.benchmark.infile,log_switch);
    
  /// 3. Read benchmark.infile
  VCF final_set (CONF.benchmark.infile); 
    
  final_set.extract_id_from_SURVIVOR_merge_format();

  final_set.update_format_name();

  final_set.rename_id("COMBINED");

  final_set.fwrite(CONF.benchmark.outfile);

  CERR << "Benchmarking vcf file: " << CONF.benchmark.outfile << "\n";
  CERR << "\n";
}


void combine_calls_svs(std::string files, double max_dist, int min_support, int type_save, int strand_save, int dynamic_size, int min_svs, std::string output, bool log_switch) {
  std::vector<std::string> names = parse_filename(files);

  Parameter::Instance()->max_caller = names.size();
  Parameter::Instance()->max_dist = max_dist;
  Parameter::Instance()->use_type = (type_save == 1);

  Parameter::Instance()->use_strand = (strand_save == 1);
  Parameter::Instance()->min_length = min_svs;
  Parameter::Instance()->dynamic_size = false;   //(dynamic_size==1);
  Parameter::Instance()->min_support = min_support;

  IntervallTree bst;
  TNode *root = NULL;
  std::map<std::string, int> chrs;
  for (size_t id = 0; id < names.size(); id++) {
    parse_vcf_header(chrs, names[id]);
    std::vector<strvcfentry> entries = parse_vcf(names[id], min_svs);
    if ( log_switch ) std::cerr << "merging entries: " << entries.size() << std::endl;
    for (size_t j = 0; j < entries.size(); j++) {
      breakpoint_str start = convert_position(entries[j].start);
      breakpoint_str stop = convert_position(entries[j].stop);
      if (entries[j].type == 4) { 
	stop.position += entries[j].sv_len;
      }
      meta_data_str tmp;
      tmp.caller_id = (int) id;
      tmp.genotype = entries[j].genotype;
      tmp.QV = entries[j].quality;
      tmp.num_reads = entries[j].num_reads;
      tmp.sv_len = entries[j].sv_len;
      tmp.pre_supp_vec = entries[j].prev_support_vec;
      tmp.vcf_ID = entries[j].sv_id;
      tmp.allleles = entries[j].alleles;

      if(start.position==151164){
	std::cout<<"FOUND "<<std::endl;
      }
      bst.insert(start, stop, entries[j].type, entries[j].strands, tmp, root);
    }
    entries.clear();
  }

  std::map<std::string, std::vector<SVS_Node *> > union_set;
  bst.get_breakpoints(root, union_set);

  FILE * file;
  file = fopen(output.c_str(), "w");
  print_header(file, names, chrs);

  int id = 0;

  std::vector<std::string> keys;
  for (std::map<std::string, std::vector<SVS_Node *> >::iterator i = union_set.begin(); i != union_set.end(); i++) {
    keys.push_back((*i).first);
  }
	
  std::sort(keys.begin(), keys.end(), numeric_string_compare);
	
  for (size_t i = 0; i < keys.size(); i++) {
    std::vector<SVS_Node *> points = union_set[keys[i]];
    for (std::vector<SVS_Node *>::reverse_iterator i = points.rbegin(); i != points.rend(); i++) {
      int support = get_support((*i)->caller_info);
      int len = 100000;
      if ((*i)->type != 3) {
	len = get_avglen((*i)->caller_info);
      }
      short type = (*i)->type;
      if ((*i)->type == -1) {
	type = 5;
      }

      if (support >= min_support && len > min_svs) {
	print_entry_overlap(file, (*i), id);
      }

      id++;
    }
  }
  fclose(file);
}

#endif
