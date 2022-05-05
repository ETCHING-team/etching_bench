/**
 * @file args.cpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief
 * @version 0.2.11
 * @date 2022-04-29
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef ARGS_CPP
#define ARGS_CPP

#include "args.hpp"
#include "prepr.hpp"
#include "define.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

args::args (){
  init();
}

args::~args(){
  init();
};
args::args( int argc , char ** argv ){
  init();
  read_args(argc, argv);
#undef DEBUG_ARGS
#ifdef DEBUG_ARGS
  print_settings();
#endif
}


void USAGE(){
  std::cerr << PROGRAM_NAME << "-" << PROGRAM_VERSION << "\n"
	    << "\n"
	    << "Usage: " << PROGRAM_NAME << " [options]\n"
	    << "\n"
	    << "Required:\n"
	    << "\t-c FILE    Config file (required)\n"
	    << "\t           Format of config file:\n"
	    << "\t               caller_name  file_name.vcf\n"
	    << "\t           Example)\n"
	    << "\t               etching ETCHING.vcf\n"
	    << "\t               delly   Delly.vcf\n"
	    << "\t               LUMPY   /path/to/sample_1/lumpy/colon_cancer.vcf\n"
	    << "\t               gridss  gridss.vcf\n"
	    << "\t               MANta   manta_tumor.vcf\n"
	    << "\t               sVaba   somatic_sv.vcf\n"
	    << "\t               nOvobreak   novobreak.pass.flt.vcf\n"
	    << "\t           You can ignore the letter cases in the caller_name.\n"
	    << "\t           (etching, ETCHING, and even etchiNg)\n"
	    << "\n"
	    << "Options:\n"
	    << "\t-o STR     Outfile prefix [etching_bench]\n"
	    << "\t-s INT     Consensus cutoff for silver standard [" << CONSENSUS_CUTOFF << "]\n"
	    << "\t           TRUE if detected by >=3 (in default) callers.\n"
	    << "\t-t FILE    Truth set\n"
	    << "\t-w INT     Merge window size [" << MERGE_WINDOW << "]\n"
	    << "\t-m INT     Minimun SV size [" << MIN_SV_SIZE <<"]\n"
	    << "\t-M INT     Maximum SV size\n"
	    << "\t-I         Remove IMPRECISE SVs\n"
	    << "\t-X         Exclude own prediction in a silver standard set [default]\n"
	    << "\t-K         Keep own prediction in a silver standard set\n"
	    << "\t-F         Calculate performance at F1-score-maximizing cutoff\n"
	    << "\t-R FILE    Re-use a given.benchmark.annotated.vcf (skip to calculation)\n"
	    << "\t--intra    Only intra-chromosomal SVs (DEL, DUP, or INV)\n"
	    << "\t--inter    Only inter-chromosomal SVs (TRA)\n"
	    << "\t--version  Print version\n"
	    << "\t--version  Print version\n"
	    << "\t-h         Print this message\n"
	    << "\n"
	    << "Note: Not for general use yet.\n"
	    << "      We guarantee only for the programs listed below:\n"
	    << "      ETCHING, DELLY, LUMPY, Manta, SvABA, novoBreak, and GRIDSS\n";
  exit(1);
}


void args::init(){
  default_source[ETCHING_TOOL] = {ETCHING_TOOL, "", ETCHING_CUTOFF, 0, ""};
  default_source[DELLY_TOOL] = {DELLY_TOOL, "", DELLY_CUTOFF, 0, ""};
  default_source[LUMPY_TOOL] = {LUMPY_TOOL, "", LUMPY_CUTOFF, 0, ""};
  default_source[MANTA_TOOL] = {MANTA_TOOL, "", MANTA_CUTOFF, 0, ""};
  default_source[SVABA_TOOL] = {SVABA_TOOL, "", SVABA_CUTOFF, 0, ""};
  default_source[NOVOBREAK_TOOL] = {NOVOBREAK_TOOL, "", NOVOBREAK_CUTOFF, 0, ""};
  default_source[GRIDSS_TOOL] = {GRIDSS_TOOL, "", GRIDSS_CUTOFF, 0, ""};

  help_set = {"h", "-h", "--h", "help", "-help", "--help", "?", "--?", "-?", 
	      "--H", "-H", "H", "--HELP", "-HELP", "HELP", "Help", "-Help", "--Help"};

  config_file=CONFIG_FILE;
  prefix=DEFAULT_PREFIX;
  truth_file=TRUTH_FILE;
  truth_cutoff=CONSENSUS_CUTOFF;
  merge_window=MERGE_WINDOW;
  min_sv_size=MIN_SV_SIZE;
  max_sv_size=MAX_SV_SIZE;
  intra_SV=0;
  inter_SV=0;

  reuse_anno_vcf=""; // R
  
  is_silver=1;
  precise_only=0; // I
  remove_self_bias=1; // X or K
  remove_own_prediction=0; // X
  keep_own_prediction=0; // K
  maximizing_f1=0; // F
}


void args::read_args ( int argc, char ** argv ){
  CHECK_HELP (argc, argv);
  read_options(argc, argv);
  check_required();
  read_conf_file();
  update_merge_info();
}


/**
 * @brief Update merge information
 * 
 * @details Merge information:
 * @details source_info initial_merge
 * @details source_info truth
 * @details source_info benchmark
 */
void args::update_merge_info(){
  initial_merge.tool = "survivor";
  initial_merge.outfile = prefix+".merge.vcf";
  initial_merge.cutoff = truth_cutoff;
  initial_merge.list = prefix+".merge.list";

  truth.tool = "truth";
  truth.infile = prefix + ".merge.vcf";
  truth.outfile = prefix + ".truth.vcf";

  benchmark.tool = "benchmark";
  benchmark.list = prefix+".benchmark.list";
  benchmark.infile = prefix+".benchmark.vcf";
  benchmark.outfile = prefix+".benchmark.annotated.vcf";
}


void args::read_options(int argc, char ** argv){
  int opt=0;
  double cutoff=truth_cutoff;
  int check_true_set=0;
  int check_silver=0;
  int check_selfbias=0;
  while (1){
    static struct option long_options[]=
      {
       {"config",           required_argument, 0, 'c'},
       {"outfile",          required_argument, 0, 'o'},
       {"truth-set",        optional_argument, 0, 't'},
       {"consensus-cutoff", optional_argument, 0, 's'},
       {"merge-window",     optional_argument, 0, 'w'},
       {"min-sv-size",      optional_argument, 0, 'm'},
       {"max-sv-size",      optional_argument, 0, 'M'},
       {"reuse-anno-vcf",   optional_argument, 0, 'R'},
       {"remove-imprecise",       no_argument, 0, 'I'},
       {"exclude-own-prediction", no_argument, 0, 'X'},
       {"keep-own-prediction",    no_argument, 0, 'K'},
       {"f1-maximizing",          no_argument, 0, 'F'},
       {"intra",                  no_argument, 0, '1'},
       {"inter",                  no_argument, 0, '2'},
       {"version",                no_argument, 0, 'v'},
       {"help",                   no_argument, 0, 'h'},
       {0, 0, 0, 0}
      };
    int opt_ind=0;
    // opt=getopt_long(argc,argv,"c:o:t:s:w:m:IXKRFvh",long_options,&opt_ind);
    opt=getopt_long(argc,argv,"c:o:t:s:w:m:M:R:IXKF12vh",long_options,&opt_ind);
    if (opt==-1) break;
#undef DEBUG_ARGS_READ_OPTIONS
#ifdef DEBUG_ARGS_READ_OPTIONS
    DEBUG_COUT("opt:"<<(char)opt<<" opt_ind:"<<opt_ind<<" argv[opt_ind]:"<<argv[opt_ind]);
#endif
    switch(opt){
    case 'c': config_file = optarg; break;
    case 'o': prefix = optarg; break;
    case 't': truth_file = optarg; truth_cutoff = 0 ; break;
    case 's': truth_cutoff = atoi(optarg); check_silver = 1 ; break;
    case 'w': merge_window = atoi(optarg); break;
    case 'm': min_sv_size = atoi(optarg); break;
    case 'M': max_sv_size = atoi(optarg); break;
    case 'R': reuse_anno_vcf = optarg; break; // reuse_anno_vcf
    case 'I': precise_only = 1 ; break; // remove IMPRECISE
    case 'X': remove_own_prediction = 1; break; // Exclude own prediction
    case 'K': keep_own_prediction = 1; break; // Keep own prediction
    case 'F': maximizing_f1 = 1 ; break; // F1-score maximizing
    case '1': intra_SV = 1 ; break; // intra only
    case '2': inter_SV = 1 ; break; // inter only
    case 'v': COUT << PROGRAM_VERSION << "\n" ; exit(0);
    case 'h': USAGE();
    default:
      CERR_ERROR("Unknown option: "<<(char)opt);
      USAGE();
    }
  }
  
  if ( truth_file.size() > 0 ){
    if ( check_silver ){
      CERR << "Warning!!! You used -t and -s together. The option \"-s " << truth_cutoff << "\" will be ignored.\n";
    }
    is_silver = 0;
    truth_cutoff = 0;
    remove_self_bias = 0;
  }

  if ( remove_own_prediction == 1 ){
    if( keep_own_prediction == 1 ){
      std::cerr << "Warning!!! -K option will be ignored because of -X option.\n";
    }
  }
  else{
    if( keep_own_prediction == 1 ){
      remove_self_bias = 0 ;
    }    
  }

  // if ( precise_only == 1 ){
  //   // do in read_conf_file();
  // }
}


void args::check_required(){
  if ( config_file == "" ){
    std::cerr << "ERROR: The config file is not set" << std::endl;
    exit(1);
  }
  if ( prefix == "" ){
    std::cerr << "ERROR: The output prefix (-o) is not set" << std::endl;
    exit(1);
  }
}


bool args::isNumber(const std::string s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

void args::read_conf_file (){
  CHECK_FIN(config_file);
  IFS fin(config_file);
  STR line;
  STR tool;
  STR infile;
  double cutoff(-1);
  while(std::getline(fin, line)){
    bool noimprecise = precise_only;
    if ( line.empty() || line[0] == '#' ) continue;
    STS ss(line);
    ss >> tool ;
    is_supporting ( tool );
    
    if ( ss >> infile ){
      if ( reuse_anno_vcf.empty() ) CHECK_FIN(infile);
      STR tmp;
      cutoff = default_source[tool].cutoff;
      if ( ss >> tmp ){
	if ( isNumber(tmp) ){
	  std::size_t sz;
	  cutoff = std::stod ( tmp , &sz);
	}
	else if (tmp == "noimprecise") {
	  noimprecise = 1;
	}
	if ( ss >> tmp ) {
	  if ( isNumber(tmp) ){
	    std::size_t sz;
	    cutoff = std::stod ( tmp , &sz);
	  }
	  else if (tmp == "noimprecise") {
	    noimprecise = 1;
	  }
	}
      }
    }
    else{
      std::cerr << "Error: " << tool << " has no file name.\n";
      CHECK_CONFIG;
    }
    source.push_back({tool, infile, cutoff, noimprecise, preprocessed_file_name (tool), cutoff_file_name (tool) } );
  }
  CLOSE(fin);
}


void args::is_supporting(STR &tool){
  STR lower_case_name = to_lower_case(tool);
  if ( default_source.find(lower_case_name) == default_source.end() ){
    std::cerr << "ERROR: " << tool << " is not supported now.\n";
    std::cerr << "This version supports: ";
    for ( auto &it : default_source )
      std::cerr << it.first << " ";
    std::cerr << "\n";
    CHECK_CONFIG;
  }
  tool = lower_case_name;
}


STR args::to_lower_case(STR str){
  STR new_str = "";
  for(int i = 0; i < str.length(); i++){
    if(str[i] >= 'A' && str[i] <= 'Z'){
      new_str += str[i] + 32;
    }
    else{
      new_str += str[i];
    }
  }
  return new_str;
}


bool args::any_help(int argc, char ** argv) {
  for(int i=0;i<argc;i++) if(help_set.find(argv[i]) != help_set.end()) return true; 
  return false; 
}


STR args::preprocessed_file_name (STR tool){
  // if ( file.substr(file.length()-4) == ".vcf" )
  //   file = file.substr(0, file.length()-4);
  // return file + ".prep.vcf";
  return prefix + "." + tool + ".prep.vcf";
}

STR args::cutoff_file_name (STR tool){
  // if ( file.substr(file.length()-4) == ".vcf" )
  //   file = file.substr(0, file.length()-4);
  // return file + ".cutoff.vcf";
  return prefix + "." + tool + ".cutoff.vcf";
}


void args::reorder(){
  std::vector < STR > tv = {"etching", "delly", "lumpy", "manta", "svaba", "novobreak", "gridss"};
  std::vector < source_info > new_source;
  for(int j = 0; j < tv.size(); j++){
    for(int i = 0; i < source.size(); i++){
      if(source[i].tool == tv[j]){
	new_source.push_back(source[i]);
      }
    }
  }
  source.clear();
  source=new_source;
  new_source.clear();
  tv.clear();
}

void args::print_source(){
  for (auto it = source.begin(); it != source.end(); ++it){
    std::cout << it->tool << " " << it->infile << " " << it->cutoff << " " << it->outfile << std::endl;
  }
}


void args::print_settings(){
  std::cerr << "config file: " << config_file << std::endl;
  std::cerr << "Outfile prefix: " << prefix << std::endl;
  if ( truth_file.size() != 0 )
    std::cerr << "Truth set: " << truth_file << std::endl;
  if ( truth_cutoff > 0){
    std::cerr << "Truth cutoff: " << truth_cutoff << std::endl;
  }
  std::cerr << "Merge window size: " << merge_window << std::endl;
  std::cerr << "Minimun SV size: " << min_sv_size << std::endl;
  for (auto it = source.begin(); it != source.end(); ++it){
    std::cerr << it->tool << ": Input file: " << it->infile << std::endl;
    std::cerr << it->tool << ": PASS cutoff: " << it->cutoff << std::endl;
    std::cerr << it->tool << ": Processed file: " << it->outfile << std::endl;
  }
}


source_info & args::operator [](unsigned int i){
  return source[i];
}


std::vector<source_info>::iterator args::begin(){
  return source.begin();
}


std::vector<source_info>::iterator args::end(){
  return source.end();
}


std::size_t args::size(){
  return source.size();
}

void args::make_initial_merge_list_file(){
  CHECK_FOUT(initial_merge.list);
  OFS fout(initial_merge.list);
  for ( std::size_t i = 0 ; i < source.size() ; i ++ ){
    // fout << source[i].cutoff_file << std::endl;
    fout << source[i].outfile << std::endl;
  }
  CLOSE(fout);
}



void args::make_benchmark_list_file(){
  CHECK_FOUT(benchmark.list);
  OFS fout(benchmark.list);
  fout << truth.outfile << std::endl;
  for ( std::size_t i = 0 ; i < source.size() ; i ++ ){
    fout << source[i].outfile << std::endl;
  }
  CLOSE(fout);
}


#endif
