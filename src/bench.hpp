/**
 * @file bench.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.11
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef BENCH_HPP
#define BENCH_HPP

#include "args.hpp"
#include "define.hpp"

#include <sstream>

/**
 * @brief Calculate performance
 * 
 * @param [in] CONF 
 * @return void
 * 
 * @details Step 1: Read vcf files
 * @details Step 2: Calculate performance
 */
void benchmark (args CONF);


/**
 * @struct truth_table bench.hpp "bench.hpp"
 * @brief A struct for truth table
 * 
 * @details TP: True Positive
 * @details FP: False Positive
 * @details FN: False Negative
 * @details TN: True Negative
 * 
 * @details Rec: Recall
 * @details Pre: Precision
 * @details Spe: Specificity
 * 
 * @details F1: F1 score
 * @details bAcc: balanced Accuracy
 * @details nMcc: normalized MCC
 * 
 * @details Rec0: Recall at 0
 * @details Pre0: Precision at 0
 * @details Rec1: Recall at 1
 * @details Pre1: Precision at 1
 * 
 * @details fdr0: FDR at 0
 * @details sen0: Sensitivity at 0
 * @details fdr1: FDR at 1
 * @details sen1: Sensitivity at 1
 * 
 * @details auroc: Area under ROC curve
 * @details aupr: Area under PR curve
 * 
 * @details truth_table(): Constructor
 * @details truth_table(VEC<double>,VEC<double>): Constructor
 * @details init_map(VEC<double>,VEC<double>): Initialize map
 * @details calc(): Calculate performance
 * @details update(double): Update for next iteration
 * @details to_string_aupr_auroc(): Return auPR and auROC in string
 */
struct truth_table{

  std::map < double, double > T_map;
  std::map < double, double > F_map;

  double TP;
  double FP;
  double FN;
  double TN;

  double Rec;
  double Pre;
  double Spe;

  double F1;
  double bAcc;
  double nMcc;

  double Rec0;
  double Pre0;
  double Rec1;
  double Pre1;

  double fdr0;
  double sen0;
  double fdr1;
  double sen1;

  double auroc;
  double aupr;

  /**
   * @brief Construct a new truth table object
   * 
   */
  void init(){
    TP=0;
    FP=0;
    FN=0;
    TN=0;

    Rec=0;
    Pre=0;
    Spe=0;

    F1=0;
    bAcc=0;
    nMcc=0;

    Rec0=1;
    Pre0=0;
    Rec1=1;
    Pre1=0;

    fdr0=1;
    sen0=1;
    fdr1=1;
    sen1=1;

    auroc=0;
    aupr=0;
  }

  /**
   * @brief Construct a new truth table object with given truth table
   * 
   * @param [in] truth_vec 
   * @param [in] score_vec of a tool (of tool_idx)
   * 
   * @details truth_vec[i] is the prediction number (in silver standard) or 1/0 of the i'th SV.
   * @details score_vec[i] is the the score of the i'th SV predicted by the caller (-1 if not predicted).
   */
  truth_table(VEC<double> truth_vec, VEC<double> score_vec,int is_silver,double truth_cutoff, int remove_self_bias){
    init();
    init_map ( truth_vec, score_vec, is_silver, truth_cutoff, remove_self_bias );
    init_update();
  }
  void init_update(){
    fdr1=1-Spe;
    sen1=Rec;
    
    Rec1=Rec;
    Pre1=Pre;
    
    auroc = 0;
    aupr  = 0;
    
    fdr0=fdr1;
    sen0=sen1;
    
    Rec0=Rec1;
    Pre0=Pre1;
    
    if ( T_map.find(-1) != T_map.end()){
      TP -= T_map[-1];
      FN += T_map[-1];
    }
    if ( F_map.find(-1) != F_map.end()){
      FP -= F_map[-1];
      TN += F_map[-1];
    }
  }
  
  // init map
  void init_map(VEC<double> truth_vec, VEC<double> score_vec,int is_silver,double truth_cutoff,int remove_self_bias){
    if ( is_silver == 1 ){
      for ( std::size_t i = 0 ; i < truth_vec.size() ; i ++){
	double prediction_number = truth_vec[i] ;
	if ( remove_self_bias && score_vec[i] >= 0 ){
	  prediction_number --;
	}
	if( prediction_number >= truth_cutoff ){
	  TP++;
	  T_map[score_vec[i]]++;
	}
	else {
	  FP++;
	  F_map[score_vec[i]]++;
	};
      }
    }
    else{
      for ( std::size_t i = 0 ; i < truth_vec.size() ; i ++){
	if(truth_vec[i] == 1) {
	  TP++;
	  T_map[score_vec[i]]++;
	}
	else {
	  FP++;
	  F_map[score_vec[i]]++;
	};
      }
    }
  }
  
  // Calculate performance
  void calc(bool is_first){
    Rec = TP / ( TP + FN ) ;
    Pre = TP / ( TP + FP ) ;
    Spe = TN / ( TN + FP ) ;
    F1 = 2 * Rec * Pre / ( Rec + Pre ) ;
    if ( Rec * Pre == 0 ) F1 = 0 ;
    bAcc = ( Rec + Spe ) / 2 ;
    nMcc = ( ( TP * TN ) - ( FP * FN ) ) / sqrt ((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) ;
    if ( ( TP * TN ) - ( FP * FN ) == 0 ) nMcc = 0 ;
    nMcc = ( nMcc + 1 ) / 2;

    fdr1=1-Spe;
    sen1=Rec;
    
    Rec1=Rec;
    Pre1=Pre;
    
    // auroc+=(fdr0-fdr1)*(sen0+sen1)/2;
    // aupr +=(Rec0-Rec1)*(Pre0+Pre1)/2;
    auroc+=(fdr0-fdr1)*sen1;
    is_first ? aupr +=(Rec0-Rec1)*Pre0 : aupr +=(Rec0-Rec1)*(Pre0+Pre1)/2 ;
    
    fdr0=fdr1;
    sen0=sen1;
    
    Rec0=Rec1;
    Pre0=Pre1;
    
  }
  
  // update
  void update(double score){
    TP -= T_map[score];
    FN += T_map[score];
    
    FP -= F_map[score];
    TN += F_map[score];
  }

  void last_update(){
    aupr +=Rec0*Pre0;
  }
  
  STR to_string(STR score){
    // calc();
    STS ss;
    ss << score << "\t" << TP << "\t" << TN << "\t" << FP << "\t" << FN << "\t" << Rec << "\t" << Pre << "\t" << Spe << "\t" << F1 << "\t" << bAcc << "\t" << nMcc << "\n";
    return ss.str();
  }
  
  STR to_string(){
    // calc();
    STS ss;
    ss << TP << "\t" << TN << "\t" << FP << "\t" << FN << "\t" << Rec << "\t" << Pre << "\t" << Spe << "\t" << F1 << "\t" << bAcc << "\t" << nMcc << "\n";
    return ss.str();
  }
  
  STR to_string_aupr_auroc(){
    STS ss;
    ss << "#auROC\t" << auroc << "\n" << "#auPR\t" << aupr << "\n";
    return ss.str();
  }
};



class benchmark_t{
private:
  /** 
   * @var STR truth
   * @brief truth vcf file 
   */
  STR truth;

  /** @brief true or false */
  VEC<double> truth_vec;

  /** 
   * @var VEC<STR> tool
   * @brief vcf files for benchmarking 
   */
  VEC<STR> tool;

  /** @brief Cutoff for tool */
  VEC<double> cutoff;

  /** 
   * @brief A vector of vector to store SV score
   *
   * @details tool_score_vec[2][66825] is the 2nd tool's score of 66825'th SV in a merged vcf.
   *
   */
  VEC < VEC<double> > tool_score_vec;

  /**
   * @brief A vector to store SV type
   *
   * @details If i'th SV is TRA, then svtype_vec [ i ] = "TRA".
   */
  VEC < STR > svtype_vec;
    
  /** @brief Line to vector */
  VEC<STR> line_to_vec(STR line);

  double truth_value(STR str);
  double prediction_number(STR str);
  
  /** @brief second element of a string delimited by ':' */
  double second_element(STR str);

  /** @brief truth table map */
  VEC<std::map<double,truth_table> > truth_table_map;

  /** 
   * @brief tool's performance map 
   * 
   * @details performance_table_map["tool_name"]["svtype"]
   */
  std::map < STR, std::map < STR, STR > > performance_table_map;
  std::map < STR, std::map < STR, STR > > cutoff_performance_map;

  /**
   * @brief Return SV type
   *
   * @details Find SV type and return it
   * @param [in] std::string INFO (the 8'th column of a VCF file)
   * @return std::string
   */
  STR extract_svtype(STR info);

  /**
   * @brief Calculate performance for each tool
   * 
   * @param [in] index
   */
  void performance(int tool_idx);

  /**
   * @var STR prefix
   * @brief 
   * 
   */
  STR prefix;


  const VEC<STR> svtype_list={"ALL","DEL","DUP","INV","TRA"};
    
public:
  benchmark_t();
  benchmark_t(args CONF);
  ~benchmark_t();

  /**
   * @brief Read vcf files
   * 
   * @param [in] CONF 
   * @return void
   * 
   */
  void get_args (args CONF);

  /**
   * @brief Calculate performance for all tools
   * 
   */
  void calc_performance();

  /**
   * @brief Write performance to file
   * 
   * @return void
   */
  void write_performance();

  /**
   * @var int is_silver
   * @brief 1 for silver standard, 0 for a given true SV set (golden standard or simulation)
   *
   */
  int is_silver;

  /**
   * @var double truth_cutoff
   * @brief true positive if prediction number >= cutoff
   *
   */  
  double truth_cutoff;

  /**
   * @var std:string benchmark_annotated_vcf
   *
   * @brief Annotated VCF file for benchmarking
   */
  STR benchmark_annotated_vcf;

  
  /** use this to remove the bias from self-prediction in benchmarking */
  int remove_self_bias;

  /**
   * @var bool maximizing_f1
   *
   * @brief Indicator for searching the F1-score maximizing cutoff
   * @details 1 for searching the F1-score maximizing cutoff, 0 for else (default)
   */
  bool maximizing_f1;

};



#endif
