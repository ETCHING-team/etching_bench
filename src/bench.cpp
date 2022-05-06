/**
 * @file bench.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.12
 * @date 2022-05-06
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef BENCH_CPP
#define BENCH_CPP

#include "bench.hpp"
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime>

void benchmark (args CONF){
  CERR << "[Benchmarking]" << std::endl;

  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  STR TIME = std::ctime ( & end_time );
  TIME.pop_back();
  CERR << "Date: " << TIME << "\n";

  benchmark_t BENCH(CONF);

  BENCH.calc_performance();
  BENCH.write_performance();

}

benchmark_t::benchmark_t(){}
benchmark_t::benchmark_t(args CONF){
  get_args(CONF);
}
benchmark_t::~benchmark_t(){}

void benchmark_t::get_args(args CONF){
  prefix = CONF.prefix;
  is_silver = CONF.is_silver;
  truth_cutoff = CONF.truth_cutoff;
  remove_self_bias = CONF.remove_self_bias;
  maximizing_f1 = CONF.maximizing_f1;
  if ( ! CONF.reuse_anno_vcf.empty() ) benchmark_annotated_vcf = CONF.reuse_anno_vcf;
  else benchmark_annotated_vcf = CONF.benchmark.outfile;

  CERR << "Annotated VCF file for benchmarking: " << benchmark_annotated_vcf << "\n";
  
  if ( is_silver ) {    
    CERR << "Silver standard" << "\t" << "True if predicted by >=" << truth_cutoff << " callers" << "\n";
    if ( remove_self_bias ) {
      CERR << "Exclude own caller from silver standard" << "\n";
    }
    else{
      CERR << "Include own caller in silver standard" << "\n";
    }
  }
    
  tool_score_vec.resize(CONF.size());
  tool.resize(CONF.size());
  for(int i=0;i<CONF.size();i++){
    tool[i]=CONF.source[i].tool;
    cutoff.push_back(CONF.source[i].cutoff);
  }
  CHECK_FIN(benchmark_annotated_vcf);
  IFS fin(benchmark_annotated_vcf);
  STR line;
  STR svtype;
  while(GETLINE(fin,line)){
    if(! line.empty() && line[0] != '#') {
      VEC<STR> line_vec=line_to_vec(line);
      if ( line_vec.size() != 10 + CONF.size() ){
	std::cerr << "Error: " << line << std::endl;
	std::cerr << "Error: " << line_vec.size() << " " << CONF.size() << std::endl;
	exit(1);
      }
      svtype_vec.push_back(extract_svtype(line_vec[7]));
      truth_vec.push_back(truth_value(line_vec[9]));
      for(int i=0;i<CONF.size();i++){
	tool_score_vec[i].push_back(second_element(line_vec[10+i]));
      }
    }
  }
  CLOSE (fin);
}

VEC<STR> benchmark_t::line_to_vec(STR line){
  VEC<STR> vec;
  STS ss(line);
  STR token;
  while ( ss >> token ){
    vec.push_back(token);
  }
  return vec;
}


STR benchmark_t::extract_svtype(STR info){
  STS ss(info);
  STR token;
  STR svtype="NONE";
  while ( GETLINE ( ss , token , ';' ) ){
    STS ss2(token);
    STR key;
    GETLINE ( ss2 , key , '=' );
    if ( key == "SVTYPE" ){
      GETLINE ( ss2 , svtype , '=' );
      break;
    }
  }
  return svtype;
}


double benchmark_t::truth_value(STR str){
  if ( is_silver ){
    return prediction_number(str);
  }
  return second_element(str);
}

double benchmark_t::prediction_number(STR str){
  STS ss(str);
  STR::size_type sz;
  STR token;
  GETLINE(ss,token,':');

  int output = 0;
  STS ss2(token);
  STR token2;
  GETLINE(ss2,token2,'_');
  while ( GETLINE(ss2,token2,'_') ){
    output += atoi ( token2.c_str() );
  }
  return (double) output;
}

double benchmark_t::second_element(STR str){
  STS ss(str);
  STR::size_type sz;
  STR token;
  GETLINE(ss,token,':');
  GETLINE(ss,token,':');
  return std::stod(token,&sz);
}


void benchmark_t::calc_performance(){
  for(std::size_t tool_idx = 0 ; tool_idx < tool.size() ; tool_idx ++){
    double tool_cutoff=0;

    if ( ! maximizing_f1 ){
      tool_cutoff = cutoff[tool_idx];
    }

    std::size_t Size = truth_vec.size();
    VEC<double> score_vec = tool_score_vec[tool_idx];
    std::set < double > score_set;
    for ( auto score : score_vec ){
      score_set.insert(score);
    }

    if ( maximizing_f1 ){
      double F1=0;
      truth_table TT_for_max_f1 ( truth_vec, score_vec, is_silver, truth_cutoff, remove_self_bias );
      bool is_first = 1;
      for (auto score : score_set){
	if ( score != -1 ){
	  if ( F1 < TT_for_max_f1.F1 ) {
	    F1 = TT_for_max_f1.F1;
	    tool_cutoff = score ;
	  }
	  TT_for_max_f1.calc(is_first);
	  TT_for_max_f1.update(score);
	  is_first = 0;
	}
      }
    }


    for ( std::size_t svtype_idx = 0 ; svtype_idx < 5 ; svtype_idx ++ ){
      STR svtype = svtype_list[svtype_idx];

      VEC<double> truth_vec_for_svtype;
      VEC<double> score_vec_for_svtype;
      std::set < double > score_set_for_svtype;

      if ( svtype_idx > 0 ) {
	for ( std::size_t j = 0 ; j < Size ; j ++ ){
	  if ( svtype_vec[j] == svtype ){
	    truth_vec_for_svtype.push_back(truth_vec[j]);
	    score_vec_for_svtype.push_back(score_vec[j]);
	    score_set_for_svtype.insert(score_vec[j]);
	  }
	}
      }
      else{
	truth_vec_for_svtype = truth_vec;
	score_vec_for_svtype = score_vec;
	score_set_for_svtype = score_set;
      }
      
      truth_table TT_for_cutoff_performance ( truth_vec_for_svtype, score_vec_for_svtype, is_silver, truth_cutoff, remove_self_bias );
      STR cutoff_performance ;
      
      cutoff_performance = TT_for_cutoff_performance.to_string("0");
      bool is_first = 1;
      for (auto score : score_set_for_svtype){
	if ( score != -1 ){
	  TT_for_cutoff_performance.calc(is_first);
	  cutoff_performance = TT_for_cutoff_performance.to_string ( std::to_string ( tool_cutoff ) ) ;
	  if ( score >= tool_cutoff){
	    break;
	  }
	  TT_for_cutoff_performance.update(score);
	  is_first = 0 ;
	}
      }
      cutoff_performance_map[tool[tool_idx]][svtype] = cutoff_performance ;
      
      truth_table TT_for_PR_curve ( truth_vec_for_svtype, score_vec_for_svtype, is_silver, truth_cutoff, remove_self_bias );
      STR performance_table_in_string ="#Score\tTP\tTN\tFP\tFN\tRec\tPre\tSpe\tF1\tbAcc\tnMCC\n";
      is_first = 1;
      for (auto score : score_set_for_svtype){
	if ( score != -1 ){
	  TT_for_PR_curve.calc(is_first);
	  performance_table_in_string += TT_for_PR_curve.to_string(std::to_string(score)) ;
	  TT_for_PR_curve.update(score);
	  is_first = 0 ;
	}
      }
      TT_for_PR_curve.last_update();
      
      
      performance_table_in_string += TT_for_PR_curve.to_string_aupr_auroc() ;
      performance_table_map[tool[tool_idx]][svtype] = performance_table_in_string; 
    }
  }
}

void benchmark_t::write_performance(){
  STR outfile = prefix + ".benchmark.txt";
  CHECK_FOUT(outfile);
  OFS fout(outfile);

  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  STR TIME = std::ctime ( & end_time );
  TIME.pop_back();
  COUT << "[" << TIME << "]\n";
  COUT << "\n";
  COUT << "Result summary:\t" << outfile << "\n";
    
  fout << "#Tool\tSV_type\tCutoff\tTP\tTN\tFP\tFN\tRec\tPre\tSpe\tF1\tbAcc\tnMCC\n";
  for ( std::size_t i = 0 ; i < tool.size() ; i ++ ){
    for ( std::size_t j = 0 ; j < 5 ; j ++ ){
      fout << tool[i] << "\t" << svtype_list[j] << "\t" << cutoff_performance_map[tool[i]][svtype_list[j]];	    
    }
  }
  CLOSE(fout);

  for(std::size_t i = 0 ; i < tool.size() ; i ++){
    COUT << "\n";
    for ( std::size_t svtype_idx = 0 ; svtype_idx < 5 ; svtype_idx ++ ){
      STR svtype = svtype_list[svtype_idx];
      STR outfile ;
      outfile = prefix + "." + tool[i] + ".performance." + svtype + ".txt";
      CHECK_FOUT(outfile);
      COUT << tool[i] << " performance table:\t" << outfile << "\n";
      OFS fout_tool_table(outfile);
      fout_tool_table << performance_table_map[tool[i]][svtype];
      CLOSE(fout_tool_table);
    }
  }
}

#endif
