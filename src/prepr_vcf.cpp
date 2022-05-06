/**
 * @file prepr_vcf.cpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.12
 * @date 2022-05-06
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef PREPR_VCF_CPP
#define PREPR_VCF_CPP

#undef DEBUG_PREPR_VCF

#include "prepr_vcf.hpp" 
#include "define.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <set>

#define UNKNOWN_STR "."

VCF_LINE::VCF_LINE()
{
  chr1=UNKNOWN_STR;
  chr2=UNKNOWN_STR;
}

VCF_LINE::VCF_LINE(STR line)
{
  chr1=UNKNOWN_STR;
  chr2=UNKNOWN_STR;
  read_vcf_line(line);
}

VCF_LINE::~VCF_LINE(){}



void VCF_LINE::read_vcf_line(STR line)
{

  chr1=UNKNOWN_STR; 
  chr2=UNKNOWN_STR;
  pos2=-1;

  STR tmp;
  STR key;
  // std::size_t sz;
  STS ss ( line );
  
  format.clear();

  std::getline ( ss,   chr1, '\t' );
  std::getline ( ss,    tmp, '\t' ); pos1 = atoi ( tmp.c_str() );
  std::getline ( ss,    id, '\t'  );
  std::getline ( ss,    ref, '\t' );
  std::getline ( ss,    alt, '\t' );
  std::getline ( ss,    tmp, '\t' ); if ( isdigit(tmp[0]) ) qual = tmp ; else qual = UNKNOWN_STR ;
  std::getline ( ss, filter, '\t' );
  std::getline ( ss,   info, '\t' );
  while(std::getline(ss,tmp,'\t') ) {
    format.push_back(tmp);
  }



  STS ss_info(info);
  std::map < STR , STR > info_map;
  while ( std::getline ( ss_info , tmp , ';') ){
    STS atom(tmp);
    STR a,b;
    std::getline ( atom , a , '=' );
    std::getline ( atom , b , '=' );
    info_map[a]=b;
  }

  svtype=info_map["SVTYPE"];
  if ( svtype == "TRA" ) svtype = "BND";
  

  STR chr2_key="CHR2=";
  STR pos2_key=";END=";

  tmp = alt;
  STR tmp1;

  key = "[";
  tmp1 = cut_str(tmp, key);
  if ( tmp.size() > 0 ) {
    chr2 = cut_str(tmp,":");
    tmp = cut_str(tmp,"[");
    pos2 = std::stol(tmp);
  }

  tmp = alt;
  key = "]";
  tmp1 = cut_str(tmp, key);
  if ( tmp.size() > 0 ) { 
    chr2 = cut_str(tmp, ":");
    tmp = cut_str(tmp,key);
    pos2 = std::stol(tmp);
  }  


  if ( info_map.find("CHR2") != info_map.end() ) chr2 = info_map["CHR2"];
  if ( info_map.find("END") != info_map.end() ) pos2 = atoi(info_map["END"].c_str());
    

  if ( chr2 == UNKNOWN_STR ){
    tmp=svtype.substr(0,3);
    if ( tmp=="DUP"||tmp=="DEL"||tmp=="INV"||tmp=="INS"||tmp=="CNV" ){
      chr2 = chr1;
    }
  }


  strand = "";
  
  key = "SVTYPE=BND";
  if ( info.find(key) != STR::npos ) {
    if ( alt[alt.size()-1] == '[' ) strand = "FR"; 
    else if ( alt[alt.size()-1] == ']' ) strand = "FF"; 
    else if ( alt[0] == '[' ) strand = "RR"; 
    else if ( alt[0] == ']' ) strand ="RF"; 
  }
  
  if ( alt[alt.size()-1] == '.' ) {
    chr2=UNKNOWN_STR;
    pos2=-1;
    strand ="F";
  }
  else if ( alt[0] == '.' ) {
    chr2=UNKNOWN_STR;
    pos2=-1;
    strand ="R";
  }

  if (strand.size() == 0) {
    key="STRANDS=";
    tmp=info;
    cut_str ( tmp , key );
    if ( tmp.size() != 0 ) {
      tmp = tmp.substr(0,2);
      if ( tmp == "+-" ) strand = "FR"; 
      else if ( tmp == "++" ) strand = "FF"; 
      else if ( tmp == "--" ) strand = "RR"; 
      else if ( tmp == "-+" ) strand = "RF"; 
    }
    else{
      key="CT=";
      tmp=info;
      cut_str(tmp,key);
      tmp=tmp.substr(0,4);
      if ( tmp == "5to3" ) {
	strand = "RF"; 
      }
      else if ( tmp == "5to5" ) { 
	strand = "RR"; 
      }
      else if ( tmp == "3to3" ) {
	strand = "FF"; 
      }
      else if ( tmp == "3to5" ) {
	strand = "FR"; 
      }
    }
  }

  dir1=0;
  dir2=0;
  if ( strand=="FF") {
    dir1=1;
    dir2=1;
  }
  else if ( strand=="FR") {
    dir1=1;
    dir2=-1;
  }
  else if ( strand=="RF") {
    dir1=-1;
    dir2=1;
  }
  else if ( strand=="RR") {
    dir1=-1;
    dir2=-1;
  }
  else if ( strand=="F") {
    dir1=1;
    dir2=0;
  }
  else if ( strand=="R") {
    dir1=-1;
    dir2=0;
  }


  tmp = info ;
  key = "MATEID=";
  cut_str(tmp,key);
  mate_id = cut_str(tmp,";");

}


STR VCF_LINE::cut_str(STR & input, const STR key)
{
  const std::size_t found (input.find(key));
  STR output;
  if ( found != STR::npos ) {
    output = input.substr(0,found);
    input = input.substr(found+key.size());
  }
  else{
    output = input;
    input.clear();
  }
  return output;
}

  



void VCF_LINE::update_svtype_in_info()
{
  std::size_t start,end;
  const STR svt("SVTYPE=");
  STR front;
  STR back;
  if ( svtype=="BND") {
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "BND" + back;
  }
  else if ( svtype=="DEL" ) {
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "DEL" + back;
  }
  else if ( svtype=="DUP" ) {
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "DUP" + back;
  }
}


void VCF_LINE::update_info()
{
  info.clear();
  if ( svtype.size() != 0 ) info += "SVTYPE=" + svtype + ";";
  if ( chr2.size() == 0 ) chr2 = UNKNOWN_STR;
  info += "CHR2=" + chr2 + ";";
  if ( pos2 > 0 ) info += "END=" + std::to_string(pos2) + ";";
  // if ( strand.size() !=0 ) info += "STRAND=" + strand + ";";
  if ( mate_id.size() !=0 ) info += "MATEID=" + mate_id + ";";
  // if ( strand == "F" ) info+="REPATH="+chr1+":"+std::to_string(pos1)+"(+)-UNKNOWN"+";";
  // else if ( strand == "R" ) info+="REPATH="+chr1+":"+std::to_string(pos1)+"(-)-UNKNOWN"+";";
  // else if ( strand == "FR") info+="REPATH="+chr1+":"+std::to_string(pos1)+"(+)-"+chr2+":"+std::to_string(pos2)+"(+)"+";";
  // else if ( strand == "FF") info+="REPATH="+chr1+":"+std::to_string(pos1)+"(+)-"+chr2+":"+std::to_string(pos2)+"(-)"+";";
  // else if ( strand == "RR") info+="REPATH="+chr1+":"+std::to_string(pos1)+"(-)-"+chr2+":"+std::to_string(pos2)+"(+)"+";";
  // else if ( strand == "RF") info+="REPATH="+chr1+":"+std::to_string(pos1)+"(-)-"+chr2+":"+std::to_string(pos2)+"(-)"+";";

  if ( chr1 == chr2 ) {
    svlen = pos1 - pos2  ;
    if (svlen<0) svlen = -svlen;
    info += "SVLEN=" + std::to_string(svlen) + ";";
  }
  if ( chr1 != chr2 && svtype != "SND" ) info += "INTERCHR;";
  if ( info[info.size()-1] == ';' ) info.pop_back();
}




STR VCF_LINE::to_string()
{
  STR tmp=chr1+"\t"+std::to_string(pos1)+"\t"+id+"\t"+ref+"\t"+alt+"\t"+qual+"\t"+filter+"\t"+info;
  for ( auto & i : format ) {
    tmp += "\t" + i ;
  }
  return tmp;
}

STR VCF_LINE::to_string_short()
{
#undef DEBUG_PREPR_VCF
#ifdef DEBUG_PREPR_VCF
  DEBUG_COUT ( "info: " << info << std::endl );
#endif
  return chr1+"\t"+std::to_string(pos1)+"\t"+id+"\t"+ref+"\t"+alt+"\t"+qual+"\t"+filter+"\t"+info;
}



VCF_LINE return_mate(VCF_LINE input)
{
  VCF_LINE output = input;

  output.chr1 = input.chr2;
  output.pos1 = input.pos2;

  output.chr2 = input.chr1;
  output.pos2 = input.pos1;

  if ( input.strand == "FR" ) output.strand = "RF";
  else if ( input.strand == "RF" ) output.strand = "FR";

  output.alt = mate_alt ( input );
  output.info = mate_info ( input );

  return output;
}


STR mate_alt ( VCF_LINE input ){
  STR alt1 = input.alt;
  STR chr1 = input.chr1;
  STR pos1 = std::to_string(input.pos1);
  STR alt2 = alt1;

  std::size_t found1 = alt1.find('[');
  std::size_t found2 = alt1.find('[',found1+1);
  std::size_t found3 = alt1.find(']');
  std::size_t found4 = alt1.find(']',found3+1);

  if ( found3 == STR::npos && found4 == STR::npos ) {
    if ( found1 != STR::npos && found2 != STR::npos ) {
      if ( found1 != 0 ) {
	alt2 = "]" + chr1 + ":" + pos1 + "]N";
      }
      else {
	alt2 = "[" + chr1 + ":" + pos1 + "[N";
      }
    }
  }
  else {
    if ( found3 != STR::npos && found4 != STR::npos ) {
      if ( found3 != 0 ) {
	alt2 = "N]" + chr1 + ":" + pos1 + "]";
      }
      else {
	alt2 = "N[" + chr1 + ":" + pos1 + "[";
      }
    }
  }
  
  return alt2;
}
  

STR mate_info ( VCF_LINE input ){
  STR info = input.info;
  STR chr1 = input.chr1;
  STR pos1 = std::to_string(input.pos1);
  STR output;
  
  STS ss ( info );
  std::vector<STR> info_vec;
  STR tmp;
  while ( std::getline ( ss , tmp , ';' ) ){
    info_vec.push_back(tmp);
  }
  
  for ( auto & i : info_vec ){
    STS ss1(i);
    STR key;
    std::getline ( ss1 , key , '=' );
    if ( key == "CHR2" ){
      i="CHR2="+chr1;
      break;
    }
  }

  for ( auto & i : info_vec ){
    STS ss1(i);
    STR key;
    std::getline ( ss1 , key , '=' );
    if ( key == "END" ){
      i="END="+pos1;
      break;
    }
  }
  
  for ( auto & i : info_vec ){
    output += i + ";";
  }
  output.pop_back();
  return output;
}


STR VCF::get_first ( STR input )
{
  STS ss(input);
  STR output ;
  ss >> output;
  return output;
}



void VCF::get_genome(STR infile){
  IFS fin(infile.c_str());
  STR id;
  STR seq;
  STR tmp;
  getline ( fin , tmp );
  id = get_first ( tmp.substr(1) ) ;
  while ( getline ( fin , tmp ) ) {
    if ( tmp[0] == '>' ) {
      genome[id]=seq;
      genome_info.push_back(std::make_pair(id,seq.size()));
      id = get_first ( tmp.substr(1) );
      seq.clear();
    }
    else{
      seq += tmp ;
    }
  }
  genome[id]=seq;
  genome_info.push_back(std::make_pair(id,seq.size()));
  fin.close();
}




void VCF::add_metainfo ( STR key , STR input , bool force)
{
  if ( force ) metainfo += input + "\n";
  else if ( metainfo.find ( key ) == STR::npos ) metainfo += input + "\n";
}

void VCF::init()
{
  id_ref_map[UNKNOWN_STR]=-1;
  ref_id_map[-1]=UNKNOWN_STR;
  // gen_TRA_mate=1;
  short_header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
}
VCF::VCF() {
  init();
}
VCF::VCF(const VCF & vcf)
{
  init();*this=vcf;
}
VCF::VCF(const STR infile)
{
  init();
  read_vcf_file(infile);
}
VCF::VCF(const source_info source_info){
  init();
  vcf_file=source_info.infile;
  tool_name=source_info.tool;
  outfile=source_info.outfile;
  read_vcf_file();
}
VCF::VCF(const args CONF, const int i){
  args conf=CONF;
  source_info source_info = conf[i];
  init();
  vcf_file=source_info.infile;
  tool_name=source_info.tool;
  outfile=source_info.outfile;
  read_vcf_file();
  noimprecise=source_info.noimprecise;
  // gen_TRA_mate=conf.gen_TRA_mate;
  cutoff_file=conf.source[i].cutoff_file;
  cutoff=conf.source[i].cutoff;
  CERR << tool_name << "\t" << vcf_file << "\t" << size() << "\tinitially\n";
}


VCF::~VCF()
{
  clear();
}

void VCF::get_header(const STR infile)
{
  metainfo.clear();
  header.clear();
  STR tmp;
  IFS fin ( vcf_file.c_str() );
  while ( std::getline ( fin , tmp ) ) {
    if ( tmp[1] == '#' ) {
      metainfo += tmp + "\n";
    }
    else{
      header = tmp ;
      break;
    }
  }
  fin.close();
}

void VCF::read_vcf_file(const STR infile)
{
  vcf_file = infile;
  read_vcf_file();
}


void VCF::read_vcf_file()
{
  STR tmp;
  VCF_LINE vcf_line;

  get_header ( vcf_file );

  IFS fin ( vcf_file.c_str() );
  while ( std::getline ( fin , tmp ) ) {
    if ( tmp[1] != '#' ) {
      break;
    }
  }

  build_id_ref_map(vcf_file);

  while ( std::getline ( fin , tmp ) ) {
    vcf_line.read_vcf_line(tmp);
    insert(vcf_line);
  }
  fin.close();
}


void VCF::build_id_ref_map(const STR infile)
{
  IFS fin ( infile.c_str() );
  std::size_t found;
  STR tmp;
  STR id;
  int count = 0;

  const STR key="##contig=<ID=";

  while ( std::getline ( fin , tmp ) ) {
    if ( tmp[1] != '#' ) break;
    found = tmp.find(key);
    if ( found != STR::npos ) {
      tmp = tmp.substr(found+key.size());
      found = tmp.find(",");
      tmp = tmp.substr(0,found);
      id_ref_map[tmp]=count;
      ref_id_map[count]=tmp;
      count ++;
    }
  }
  
  fin.close();
}


void VCF::clear()
{
  vcf_map.clear();
  genome_info.clear();
  vcf_file.clear();
  header.clear();
  metainfo.clear();

  id_ref_map.clear();
  ref_id_map.clear();
}

void VCF::update_id_ref_map(STR chr)
{
  if ( id_ref_map.find(chr) == id_ref_map.end()) {
    std::map < int , STR > :: iterator it = ref_id_map.end();
    it -- ;
    int NEW_REF_NUM ( it->first + 1 );
    id_ref_map[chr] = NEW_REF_NUM;
    ref_id_map[NEW_REF_NUM] = chr;
  }
}


char VCF::GetNucl(const STR & refChr, const std::size_t & refPos)
{
  return genome[refChr][refPos-1];
}


void VCF::insert(VCF_LINE vcf_line)
{
  update_id_ref_map(vcf_line.chr1);
  update_id_ref_map(vcf_line.chr2);

  BP bp1;
  BP bp2;
  bp1.insert ( id_ref_map[vcf_line.chr1] , vcf_line.pos1 , vcf_line.dir1 );
  bp2.insert ( id_ref_map[vcf_line.chr2] , vcf_line.pos2 , vcf_line.dir2 );

  vcf_map[bp1][bp2]=vcf_line;
}





BND VCF::find(BP bp1, BP bp2)
{
  BND output;
  if (vcf_map.find(bp1) != vcf_map.end()) {
    if ( vcf_map[bp1].find(bp2) != vcf_map[bp1].end() )  {
      output.insert(bp1,bp2);      
    }
  }
  return output;
}


BND VCF::find(BND bnd)
{
  BND output;
  if (vcf_map.find(bnd.first) != vcf_map.end()) {
    if ( vcf_map[bnd.first].find(bnd.second) != vcf_map[bnd.first].end() )  {
      output=bnd;
    }
  }
  return output;
}



void VCF::write()
{
  std::cout << metainfo << header << "\n";
  write_no_header();
}


void VCF::write_no_header()
{ 
  for ( auto & i : vcf_map ) {
    for ( auto & j : i.second ) {
      VCF_LINE * vcf_line ( & j.second );
      std::cout << vcf_line->to_string() << "\n";
    }
  }
}

void VCF::fwrite(STR OUTFILE)
{
  OFS fout (OUTFILE.c_str());
  fout << metainfo << header << "\n";
  for ( auto & i : vcf_map ) for ( auto & j : i.second ) {
      VCF_LINE * vcf_line ( & j.second );
      fout << vcf_line->to_string() << "\n";
    }
  fout.close();
}

void VCF::fwrite()
{
  OFS fout (outfile.c_str());
  fout << metainfo << header << "\n";
  for ( auto & i : vcf_map ) for ( auto & j : i.second ) {
      VCF_LINE * vcf_line ( & j.second );
      fout << vcf_line->to_string() << "\n";
    }
  fout.close();
}

void VCF::write_short()
{
  std::cout << metainfo << short_header << "\n";
  for ( auto & i : vcf_map ) for ( auto & j : i.second ) {
      VCF_LINE * vcf_line ( & j.second );
      std::cout << vcf_line->to_string_short() << "\n";
    }
}


void VCF::fwrite_short(STR OUTFILE)
{
  OFS fout (OUTFILE.c_str());
  fout << metainfo << short_header << "\n";
  for ( auto & i : vcf_map ) for ( auto & j : i.second ) {
      VCF_LINE * vcf_line ( & j.second );
      fout << vcf_line->to_string_short() << "\n";
    }
  fout.close();
}


void VCF::fwrite_short()
{
  OFS fout (outfile.c_str());
  fout << metainfo << short_header << "\n";
  for ( auto & i : vcf_map ) for ( auto & j : i.second ) {
      VCF_LINE * vcf_line ( & j.second );
      fout << vcf_line->to_string_short() << "\n";
    }
  fout.close();
}

void VCF::rename_id ( STR prefix )
{
  int count(1);
  for ( auto & i : vcf_map ) for ( auto & j : i.second ) {
      VCF_LINE * vcf_line ( & j.second );
      vcf_line->id = prefix + std::to_string(count);
      vcf_line->mate_id.clear();
      count ++;
    }
}


std::size_t VCF::size()
{
  std::size_t SIZE(0);
  for ( auto i : vcf_map ) for ( auto j : i.second ) SIZE++;
  return SIZE;
}


std::size_t VCF::above_cutoff_size()
{
  std::size_t SIZE(0);
  std::size_t sz;
  for ( auto i : vcf_map ) for ( auto j : i.second ) if ( std::stod(j.second.qual.c_str(),&sz) > cutoff ) SIZE++;
  return SIZE;
}



void VCF::typing_SV()
{
  VCF_MAP replace;

  int ins_count = 0;
  int del_count = 0;
  int dup_count = 0;
  int inv_count = 0;

  std::set < STR > eliminated;

  for ( auto & i : vcf_map ) {
    for ( auto & j : i.second ) {
      VCF_LINE vcf_line = j.second;
      VCF_LINE typed_vcf_line = vcf_line;
      
      if ( eliminated.find(vcf_line.id) == eliminated.end() ) {
	
	if ( vcf_line.svtype == "BND") { // intra-Chr. BNDs are classified into DEL, DUP, or INV.
	  if ( vcf_line.chr1 == vcf_line.chr2 ) { // intra-Chr. BNDs are classified into DEL, DUP, or INV.
	    if ( vcf_line.pos1 == vcf_line.pos2 ) {
	      if ( vcf_line.dir1 == 1 && vcf_line.dir2 ==1 ) { // INS
		ins_count ++ ;
		
		typed_vcf_line.id = "INS" + std::to_string(ins_count);
		typed_vcf_line.alt = "<INS>";
		typed_vcf_line.svtype = "INS";
		typed_vcf_line.mate_id = vcf_line.id + "," + vcf_line.mate_id;
		
		eliminated.insert(vcf_line.id);
		eliminated.insert(vcf_line.mate_id);
	      }
	    }
	    else{
	      if ( vcf_line.pos1 < vcf_line.pos2 ) {
		if ( vcf_line.dir1 == 1 && vcf_line.dir2 == -1 ) { // DEL
		  
		  del_count ++ ;
		  
		  typed_vcf_line.id = "DEL" + std::to_string(del_count);
		  typed_vcf_line.alt = "<DEL>";
		  typed_vcf_line.svtype = "DEL";
		  typed_vcf_line.mate_id = vcf_line.id + "," + vcf_line.mate_id;
		  
		  typed_vcf_line.id = "DEL" + std::to_string(del_count);
		  
		  eliminated.insert(vcf_line.id);
		  eliminated.insert(vcf_line.mate_id);
		}
		else if ( vcf_line.dir1 == -1 && vcf_line.dir2 == 1 ) { // DUP
		  dup_count ++ ;
	      
		  typed_vcf_line = vcf_line;
		  typed_vcf_line.id = "DUP" + std::to_string(dup_count);
		  typed_vcf_line.alt = "<DUP>";
		  typed_vcf_line.svtype = "DUP";
		  typed_vcf_line.mate_id = vcf_line.id + "," + vcf_line.mate_id;
	      
		  eliminated.insert(vcf_line.id);
		  eliminated.insert(vcf_line.mate_id);
		}
		else { // INV : if (dir1,dir2) is (1,1) or (-1,-1)
		  inv_count ++ ;
	      
		  typed_vcf_line = vcf_line;
		  typed_vcf_line.id = "INV" + std::to_string(inv_count);
		  typed_vcf_line.alt = "<INV>";
		  typed_vcf_line.svtype = "INV";
		  typed_vcf_line.mate_id = vcf_line.id + "," + vcf_line.mate_id;
	      
		  eliminated.insert(vcf_line.id);
		  eliminated.insert(vcf_line.mate_id);
		}
	      }
	      typed_vcf_line.update_svtype_in_info();
	    }
	  }
	  else {
	    if ( vcf_line.chr2.size() == 0 || vcf_line.chr2=="." ){
	      typed_vcf_line.alt = "<SND>";
	      typed_vcf_line.svtype = "SND";
	      
	      eliminated.insert(vcf_line.id);
	    }
	  }
	}
	typed_vcf_line.update_info();

	BP bp1(id_ref_map[vcf_line.chr1], vcf_line.pos1, vcf_line.dir1);
	BP bp2(id_ref_map[vcf_line.chr2], vcf_line.pos2, vcf_line.dir2);

	replace[bp1][bp2]=typed_vcf_line;
      }
    }
  }
  vcf_map = replace;
  // CERR << tool_name << "\t" << vcf_file << "\t" << size() << "\tafter_typing\n";
}


void VCF::remove_mate()
{
  VCF_MAP replace;

  std::set < STR > eliminated;

  for ( auto & i : vcf_map ) {
    for ( auto & j : i.second ) {
      VCF_LINE vcf_line = j.second;
      if ( eliminated.find(vcf_line.id) == eliminated.end() ){
	
	eliminated.insert(vcf_line.id);
	eliminated.insert(vcf_line.mate_id);

	BP bp1(id_ref_map[vcf_line.chr1], vcf_line.pos1, vcf_line.dir1);
	BP bp2(id_ref_map[vcf_line.chr2], vcf_line.pos2, vcf_line.dir2);

	replace[bp1][bp2]=vcf_line;
      }
    }
  }
  vcf_map = replace;
  // CERR << tool_name << "\t" << vcf_file << "\t" << size() << "\tafter_remove_mate\n";
}

void VCF::rearrange_bp()
{
  VCF_MAP replace;

  for ( auto & i : vcf_map ) {
    for ( auto & j : i.second ) {
      VCF_LINE vcf_line = j.second;
      BP bp1(id_ref_map[vcf_line.chr1], vcf_line.pos1, vcf_line.dir1);
      BP bp2(id_ref_map[vcf_line.chr2], vcf_line.pos2, vcf_line.dir2);
      if ( bp1 > bp2 ) vcf_line = return_mate ( vcf_line );
      replace[bp1][bp2]=vcf_line;
    }
  }
  vcf_map = replace;
}


void VCF::make_mate_for_TRA()
{
  VCF replace;

  for ( auto & i : vcf_map ) {
    for ( auto & j : i.second ) {
      VCF_LINE vcf_line = j.second;
      BP bp1(id_ref_map[vcf_line.chr1], vcf_line.pos1, vcf_line.dir1);
      BP bp2(id_ref_map[vcf_line.chr2], vcf_line.pos2, vcf_line.dir2);

      VCF_LINE mate = return_mate ( vcf_line );
      BP mate_bp1(id_ref_map[mate.chr1], mate.pos1, mate.dir1);
      BP mate_bp2(id_ref_map[mate.chr2], mate.pos2, mate.dir2);
      BND bnd=replace.find(mate_bp1,mate_bp2);
      if ( vcf_line.svtype == "BND" && bnd.is_empty()){
	  vcf_line.id+="_a";
	  mate.id+="_b";
	  replace[bp1][bp2]=vcf_line;
	  replace[mate_bp1][mate_bp2]=mate;
      }
      else{
	  replace[bp1][bp2]=vcf_line;
      }
    }
  }
  vcf_map = replace.vcf_map;
}



double return_depdif(Position Pos, std::vector < double > & dep_vec, pos_type Size, double window)
{
  pos_type start (Pos.second - window - 1);
  pos_type center (Pos.second - 1);
  pos_type end (Pos.second + window - 1);
  pos_type max (dep_vec.size());
  double left (0) ;
  double right (0);

  if ( start < 0  ) start = 0    ;
  if ( end > Size ) end   = Size ;

  if ( center > max ) center = max;
  if ( end > max ) end = max;

  for ( pos_type pos = start ; pos < center ; pos ++ ) left += dep_vec[pos];
  for ( pos_type pos = center + 1 ; pos < end ; pos ++ ) right += dep_vec[pos];

  double output ( ( left - right ) / window );
  return output < 0 ? output * -1 : output;
}



BP_MAP & VCF::operator [](BP bp)
{
  return vcf_map[bp];
}


VCF_MAP::iterator VCF::begin()
{
  return vcf_map.begin();
}

VCF_MAP::iterator VCF::end()
{
  return vcf_map.end();
}


STR VCF::current_date()
{
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct); // YYYY-MM-DD
  return buf;
}


STR VCF::current_date_time()
{
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct); // YYYY-MM-DD.HH:mm:ss
  return buf;
}

void VCF::add_missed_header_lines()
{
  add_metainfo("##fileformat=" , "##fileformat=VCFv4.2", 0);
  add_metainfo("##fileDate=", "##fileDate="+ current_date(), 0);
  for ( auto & i : genome_info ) add_metainfo( "##contig=<ID="+i.first+",", "##contig=<ID="+i.first+",length="+std::to_string(i.second)+">", 0);
  add_metainfo("##INFO=<ID=SVTYPE," , "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant.\">", 0);
  add_metainfo("##INFO=<ID=CHR2,", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation.\">", 0);
  add_metainfo("##INFO=<ID=END,", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant.\">", 0);
  add_metainfo("##INFO=<ID=STRAND,", "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand(s) at break-points, F (5') or R (3').\">", 0);
  add_metainfo("##INFO=<ID=MATEID,", "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate BND, or original IDs of BNDs before SV typing.\">", 0);
  add_metainfo("##INFO=<ID=SVLEN,", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV.\">", 0);
  add_metainfo("##INFO=<ID=REPATH,", "##INFO=<ID=REPATH,Number=1,Type=String,Description=\"Path of rearrangement.\">", 0);
  add_metainfo("##INFO=<ID=INTERCHR,", "##INFO=<ID=INTERCHR,Number=0,Type=Flag,Description=\"Inter-chromosomal variation.\">", 0);
  add_metainfo("##INFO=<ID=SVSCORE," , "##INFO=<ID=SVSCORE,Number=1,Type=Float,Description=\"Quality score of the SV (same with QUAL).\">", 0);
  add_metainfo("##INFO=<ID=SCOREMETHOD," , "##INFO=<ID=SCOREMETHOD,Number=1,Type=String,Description=\"Scoring methods.\">", 0);
  add_metainfo("##ALT=<ID=DEL,", "##ALT=<ID=DEL,Description=\"Deleltion\">", 0);
  add_metainfo("##ALT=<ID=DUP,", "##ALT=<ID=DUP,Description=\"Duplication\">", 0);
  add_metainfo("##ALT=<ID=INV,", "##ALT=<ID=INV,Description=\"Inversion\">", 0);
  add_metainfo("##ALT=<ID=TRA,", "##ALT=<ID=TRA,Description=\"Breakend, which may be translocation or unclssified variation yet\">", 0);
  add_metainfo("##ALT=<ID=SND,", "##ALT=<ID=SND,Description=\"Single breakend, which may be translocation to unknown contig/scaffold or unclssified variation\">", 0);
  add_metainfo("##FILTER=<ID=PASS,", "##FILTER=<ID=PASS,Description=\"High confidential variation greater than or equal to cut-off.\">", 0);
  add_metainfo("##FILTER=<ID=LOWQUAL,", "##FILTER=<ID=LOWQUAL,Description=\"Low quality variation lower than cut-off.\">", 0);
}



void VCF::add_score_to_id(){
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      vcf_line.id = vcf_line.id + "_" + vcf_line.qual;
      vcf_map[i.first][j.first]=vcf_line;
    }
  }
}


void VCF::remove_short (int min_sv_size){
  VCF_MAP replace;

  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      VCF_LINE vcf = j.second;
      if ( vcf.chr1 == vcf.chr2 ){
	if ( vcf.svlen >= min_sv_size ){
	  replace [i.first][j.first] = vcf;
	}
      }
      else {
	if ( ! vcf.is_empty_second() ) replace [i.first][j.first] = vcf;
      }
    }
  }
  int n0 = size();
  vcf_map = replace;
  int n1 = size();
  CERR << tool_name << "\t" << n0 - n1 << " small SVs removed (<" << min_sv_size << ")\n";
}


void VCF::remove_large (int max_sv_size){
  if ( max_sv_size > 0 ){
    VCF_MAP replace;

    for ( auto & i : vcf_map ){
      for ( auto & j : i.second ){
	VCF_LINE vcf = j.second;
	if ( vcf.chr1 == vcf.chr2 ){
	  if ( vcf.svlen < max_sv_size ){
	    replace [i.first][j.first] = vcf;
	  }
	}
	else {
	  if ( ! vcf.is_empty_second() ) replace [i.first][j.first] = vcf;
	}
      }
    }
    int n0 = size();
    vcf_map = replace;
    int n1 = size();
    CERR << tool_name << "\t" << n0 - n1 << " large SVs removed (<" << max_sv_size << ")\n";
  }
}




void VCF::remove_intra_chromosome (int inter_SV){
  if ( inter_SV ){
    VCF_MAP replace;

    for ( auto & i : vcf_map ){
      for ( auto & j : i.second ){
	VCF_LINE vcf = j.second;
	if ( vcf.chr1 != vcf.chr2 && ! vcf.is_empty_second() ){
	  replace [i.first][j.first] = vcf;
	}
      }
    }
    int n0 = size();
    vcf_map = replace;
    int n1 = size();
    CERR << tool_name << "\t" << n0 - n1 << " intra-chromosomal SVs removed\n";
  }
}



void VCF::remove_inter_chromosome (int intra_SV){
  if ( intra_SV ){
    VCF_MAP replace;

    for ( auto & i : vcf_map ){
      for ( auto & j : i.second ){
	VCF_LINE vcf = j.second;
	if ( vcf.chr1 == vcf.chr2 && ! vcf.is_empty_second() ){
	  replace [i.first][j.first] = vcf;
	}
      }
    }
    int n0 = size();
    vcf_map = replace;
    int n1 = size();
    CERR << tool_name << "\t" << n0 - n1 << " inter-chromosomal SVs removed\n";
  }
}




bool VCF_LINE::is_empty_second()
{
  return chr2==UNKNOWN_STR;
}

int VCF::get_format_size(){
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
#undef DEBUG_GET_FORMAT_SIZE
#ifdef DEBUG_GET_FORMAT_SIZE
      for ( auto k : vcf_line.format ){
	std::cerr << k << std::endl;
      }
      std::cerr << "\n";
#endif
      return vcf_line.format.size();
      break;
    }
    break;
  }
  return 0;
}

STR VCF_LINE::binary_call_set(){
  STR call_set;
  for ( int i = 1 ; i < format.size() ; i++ ){
    format[i]=="./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN"?call_set+="_0":call_set+="_1";
  }
  return call_set;
}

void VCF::add_binary_calls_to_id(STR file_name){
  int format_size = get_format_size();
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      vcf_line.id += vcf_line.binary_call_set();
      vcf_map[i.first][j.first]=vcf_line;
    }
  }
}

/**
 * @brief Extract ID SURVIVOR merge FORMAT
 * @details Replace FORMAT with ID
 */
void VCF_LINE::extract_id_from_SURVIVOR_merge_format(){
  STR update;
  for ( int i = 1 ; i < format.size() ; i++ ){
    if ( format[i]=="./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN"){
      update = "NONE";
    }
    else{
      STS ss(format[i]);
      for ( int i = 0 ; i < 8 ; i++ ){
	std::getline (ss, update, ':');
      }
    }
    format[i] = update;
  }
}

/**
 * @brief Extract ID SURVIVOR merge FORMAT
 * @details Replace FORMAT with ID
 */
void VCF::extract_id_from_SURVIVOR_merge_format(){
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      vcf_line.extract_id_from_SURVIVOR_merge_format();
      vcf_map[i.first][j.first]=vcf_line;
    }
  }
}


/**
 * @brief Update FORMAT name
 * 
 */
void VCF::update_format_name(){
  int ID_NUMBER=1;
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      vcf_line.update_format_name(ID_NUMBER++);
      vcf_map[i.first][j.first]=vcf_line;
    }
  }
}

void VCF_LINE::update_format_name(int ID_NUMBER){
  format[0] = "ID:SCORE";
    
  if ( format[1] == "NONE") format[1] += ":0";
  else format[1] += ":1";
    
  for ( int i = 2 ; i < format.size() ; i++ ){          
    STR format_backup = format[i];
    // find last "_" in format[i]
    int last_underscore = format_backup.find_last_of("_");
    // Before last_underscore is the ID
    STR id = format_backup.substr(0, last_underscore);
    // After last_underscore is the score
    STR score = format_backup.substr(last_underscore+1);
    score=="NONE"?score="-1":score=score;
    format[i] = id + ":" + score;
  }
}


/**
 * @brief Update FORMAT name for silver standard
 * 
 */
void VCF::update_format_name_for_silver_standard(){
  int ID_NUMBER=1;
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      vcf_line.update_format_name_for_silver_standard(ID_NUMBER++);
      vcf_map[i.first][j.first]=vcf_line;
    }
  }
}

void VCF_LINE::update_format_name_for_silver_standard(int ID_NUMBER){
  format[0] = "ID:SCORE";
    
  for ( int i = 1 ; i < format.size() ; i++ ){          
    STR format_backup = format[i];
    // find last "_" in format[i]
    int last_underscore = format_backup.find_last_of("_");
    // Before last_underscore is the ID
    STR id = format_backup.substr(0, last_underscore);
    // After last_underscore is the score
    STR score = format_backup.substr(last_underscore+1);
    score=="NONE"?score="-1":score=score;
    format[i] = id + ":" + score;
  }
}


void VCF::remove_low_qual_records(){
  VCF_MAP replace;
  for ( auto i : vcf_map ){
    for ( auto j : i.second ){
      VCF_LINE vcf_line = j.second;
      if ( stod(vcf_line.qual) >= cutoff ){
	replace[i.first][j.first] = vcf_line;
      }
    }
  }
  vcf_map = replace;
  // CERR << tool_name << "\t" << vcf_file << "\t" << size() << "\tafter_remove_low_quality_cutoff_" << cutoff << "\n";
}



#endif
