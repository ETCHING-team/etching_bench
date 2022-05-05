/**
 * @file prepr_tools.cpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief
 * @version 0.2.11
 * @date 2022-04-29
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "prepr_vcf.hpp"
#include "define.hpp"
#include <sstream>
#include <string>

#ifndef IS_PRECISE
#define IS_PRECISE
int is_precise ( std::string input ){
  STS ss ( input );  
  std::string tmp;
  while ( std::getline ( ss , tmp , ';' ) ){
    if ( tmp == "IMPRECISE" ){
      return 0;
    }
  }
  return 1;
}

void VCF::remove_imprecise(){
  VCF_MAP replace;
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      if ( is_precise( j.second.info ) ){
	replace[i.first][j.first]=j.second;
      }
    }
  }
  int n0=size();
  vcf_map = replace;
  int n1=size();
  CERR << tool_name << "\t" << n0 - n1 << " IMPRECISE SVs removed\n";
}

#endif


#ifndef PREPR_DELLY
#define PREPR_DELLY



/**
 * @brief return delly score
 * 
 * @param [in] input in format "GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV"
 * @return int (DV+RV)
 */
int delly_score(std::string input){
  STS ss ( input );
  std::vector < std::string > tmp(12);
  for ( std::size_t i = 0 ; i < 12 ; i ++ ){
    GETLINE ( ss , tmp[i] , ':' );
  }
  return atoi(tmp[9].c_str()) + atoi(tmp[11].c_str());
}

void VCF::prepr_delly(){
  VCF_MAP replace;
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      int score_1 = delly_score(vcf_line.format[1]); 
      // int score_2 = delly_score(vcf_line.format[2]);
      vcf_line.qual = std::to_string ( score_1 ) ;
      replace[i.first][j.first]=vcf_line;      
      // int score = score_1>score_2?score_1:score_2;
      // vcf_line.qual = std::to_string ( score ) ;
      // replace[i.first][j.first]=vcf_line;      
      // int score_1 = delly_score(vcf_line.format[1]); 
      // int score_2 = delly_score(vcf_line.format[2]);
      // if ( score_2 == 0 ){
      //  	vcf_line.qual = std::to_string ( score_1 ) ;
      //  	replace[i.first][j.first]=vcf_line;
      // }
    }
  }
  // int n0=size();
  vcf_map = replace;
  // int n1=size();
  // CERR << tool_name << "\t" << n0 - n1 << " germline SVs removed\n";
}




#endif


#ifndef PREPR_LUMPY
#define PREPR_LUMPY

/**
 * @brief return SU in the form "GT:SU:PE:SR"
 * @param [in] input string in the form "GT:SU:PE:SR"
 * @return int SU in the form "GT:SU:PE:SR"
 */
int calc_val_lumpy ( std::string input ){
  STS ss ( input );  
  std::string tmp;
  std::getline ( ss , tmp , ':' );
  std::getline ( ss , tmp , ':' );
  return atoi ( tmp.c_str() );
}

void VCF::prepr_lumpy(){
  // FORMAT[0] = "GT:SU:PE:SR"
  // FORMAT[1] SU is somatic score
  // FORMAT[2] is germline score
  // Remove germline SV FORMAT[2].SU > 0
  // Add FORMAT[1].SU to QUAL

  VCF_MAP replace;

  // iterate over all elements of vcf_map in VCF class
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      if ( vcf_line.format.size() == 2 ){
	int score = calc_val_lumpy ( vcf_line.format[1] );
	vcf_line.qual = std::to_string ( score );
	replace[i.first][j.first]=vcf_line;
      }
      if ( vcf_line.format.size() == 3 ){
	int somatic_score = calc_val_lumpy ( vcf_line.format[1] );
	int germline_score = calc_val_lumpy ( vcf_line.format[2] );
	if ( germline_score == 0 ){
	  vcf_line.qual = std::to_string ( somatic_score );
	  replace[i.first][j.first]=vcf_line;
	}
      }
    }
  }
  int n0=size();
  vcf_map = replace;
  int n1=size();
  CERR << tool_name << "\t" << n0 - n1 << " germline SVs removed\n";
}



#endif


#ifndef PREPR_MANTA
#define PREPR_MANTA

/**
 * @brief Return SOMATICSCORE=XXX; from a VCF record
 * @param [in] info 
 * @details 1. Split the INFO field into tokens by ';'
 * @details 2. store the tokens into a vector
 * @details 3. find the SOMATICSCORE=XXX token
 * @details 4. return XXX
 */
std::string manta_score(std::string info) {
  STS ss(info);
  std::string item;
  std::vector<std::string> tokens;
  while (std::getline(ss, item, ';')) {
    tokens.push_back(item);
  }
  for (int i = 0; i < tokens.size(); i++) {
    if (tokens[i].find("SOMATICSCORE") != std::string::npos) {
      return tokens[i].substr(tokens[i].find("=") + 1);
    }
  }
  return "";
}

/**
 * @brief Add manta SOMATICSCORE of info field to qual field for all recods in vcf file
 * 
 */
void VCF::prepr_manta(){
  VCF_MAP replace;
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      vcf_line.qual = manta_score ( vcf_line.info );
      if (vcf_line.qual == "" ){
	vcf_line.qual = std::to_string( manta_PR_plus_SR( vcf_line.format ) );
      }
      replace[i.first][j.first]=vcf_line;
    }
  }
  vcf_map = replace;
}


int VCF::manta_PR_plus_SR ( std::vector < std::string > format_vec ){
  int output;
  for ( auto i : format_vec ){
    std::stringstream ss ( i );
    std::string token;
    while ( std::getline ( ss , token , ':' ) ){
      std::stringstream ss1 ( token );
      std::string token1;
      std::getline ( ss1 , token1, ',' );
      std::getline ( ss1 , token1, ',' );
      output += atoi ( token1.c_str() );
    }
  }

  return output;
}



#endif


#ifndef PREPR_NOVOBREAK
#define PREPR_NOVOBREAK

/**
 * @brief Add id for novoBreak
 * @param [in] id_prefix
 */

void VCF::prepr_novoBreak(){
  rename_id ( "novobreak" );
}


#endif





#ifndef PREPR_GRIDSS
#define PREPR_GRIDSS

/**
 * @brief Add id for GRIDSS
 * @param [in] id_prefix
 */

void VCF::prepr_gridss(){
  // select only "PASS" and "qual" in FILTER field
  VCF_MAP replace;
  
  // iterate over all elements of vcf_map in VCF class
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      if ( vcf_line.filter == "PASS" || vcf_line.filter == "qual" ){
	replace[i.first][j.first]=vcf_line;
      }
    }
  }
  int n0=size();
  vcf_map = replace;
  int n1=size();
  CERR << tool_name << "\t" << n0 - n1 << " low-quality-tagged SVs removed\n";
}


#endif


#ifndef PREPR_TOOLS
#define PREPR_TOOLS

void VCF::prepr_tools(){
  if ( noimprecise ) remove_imprecise();
  if (tool_name == "novobreak") prepr_novoBreak();
  else if ( tool_name == "lumpy") prepr_lumpy();
  else if ( tool_name == "delly") prepr_delly();
  else if ( tool_name == "manta") prepr_manta();
  // else if ( tool_name == "gridss") prepr_gridss();
}




#endif
