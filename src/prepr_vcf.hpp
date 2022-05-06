/**
 * @file prepr_vcf.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.12
 * @date 2022-05-06
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef PRERP_VCF
#define PRERP_VCF

#include <string>
#include <vector>
#include <map>
#include "args.hpp"

/**
 * @brief Position type
 * 
 * @details Position::first = (chromosome number)
 * @details Position::second = (position)
 */
using Position = std::pair < int , int > ;


/// Position type
typedef int64_t pos_type;

/// Breakpoint type
typedef std::tuple<pos_type,pos_type,pos_type> BP_type;

/// Breakend (BND) type
typedef std::pair<BP_type,BP_type> BND_type;


/**
 * @brief Class for a node in the BP graph
 * 
 * @details Each vcf line has at least 10 columns
 * @details 1. CHROM
 * @details 2. POS
 * @details 3. ID
 * @details 4. REF
 * @details 5. ALT
 * @details 6. QUAL
 * @details 7. FILTER
 * @details 8. INFO
 * @details 9. FORMAT
 * @details 10. SAMPLES in FORMAT
 */


class VCF_LINE{
private:
  /**
   * @brief Parsing tool to find a string with keyword
   * @note If input = "SCORE=13;CHR2=13;END=5253219;" and key = "CHR2",
   * @note it returns "13" and modify input to "END=5253219;" by cut
   * @param [in,out] &input Input string
   * @param [in] key Keyword
   * @return String related with the key
   * @return Cut input string after the keyword
   */
  std::string cut_str(std::string & input, std::string key);

  /**
   * @brief Find value of a key in a string
   * @param [in] input Input string
   * @param [in] key Keyword
   * @return String related with the key
   */

public:
  /// Constructor 
  VCF_LINE();
  /** 
      @brief Constructor with an input
      @param [in] line A vcf line in string
  */
  VCF_LINE(std::string line);
  /** 
      @brief Destructor 
  */
  ~VCF_LINE();


  /** 
      @defgroup VCF_LINE_ATTRIBUTES VCF attributes
      @details Each vcf line has at least 10 columns
      @tableofcontents
      | Number | Field | Variable | Description |
      | -------|--------|----------|-------------|
      | 1 | CHROM  | @ref chr1 | Chromosome of the mutation |
      | 2 | POS | @ref pos1 | Position of the mutation |
      | 3 | ID | @ref id | Mutation ID |
      | 4 | REF | @ref ref | Nucleotide (or sequence) in refeference genome |
      | 5 | ALT | @ref alt | Alteration of the mutation |
      | 6 | QUAL | @ref qual | Prediction quality |
      | 7 | FILTER | @ref filter | Quality tag(s) such as PASS, LOWQUAL, ... |
      | 8 | INFO | @ref info | Information of the mutation |
      | 9 | FORMAT | @ref FORMAT | Description of FORMAT |
      | 10- | SAMPLES | @ref SAMPLES | In FORMAT |
  */

  /**
     @brief Chromosome ID
     @ingroup VCF_LINE_ATTRIBUTES
  */
  std::string chr1;
  /**
     @brief Position of the mutation in chr1
     @ingroup VCF_LINE_ATTRIBUTES
  */
  pos_type pos1;
  /**
     @brief Direction of BP
     @ingroup VCF_LINE_ATTRIBUTES
  */
  pos_type dir1;
  
  /**
     @brief SV id
     @ingroup VCF_LINE_ATTRIBUTES
  */
  std::string id;

  /**
     @brief Nucleotide (or sequence) in refeference genome
     @ingroup VCF_LINE_ATTRIBUTES
  */
  std::string ref;
  /**
     @brief Alteration of the mutation
     @details
     @tableofcontents
     | ID | Description|
     |----|------------|
     |DEL | Deletion |
     |DUP | Tandem duplication|
     |INV | Inversion |
     |TRA | Inter-chromosomal translocation|
     |BND | Breakend (a pair of BPs), or undetermined SV |
     |SND | Single breakend (without mate-BP) |
     @see svtype
     @ingroup VCF_LINE_ATTRIBUTES
  */
  std::string alt;
 
  /**
     @brief Prediction quality
     @ingroup VCF_LINE_ATTRIBUTES
  */
  std::string qual;

  /**
     @brief Quality tag(s) 
     @details
     @tableofcontents
     | ID | Description|
     |----|------------|
     |PASS|High confidential variation greater than or equal to cut-off|
     |LOWQAUL|Low quality variation lower than cut-off|
     @ingroup VCF_LINE_ATTRIBUTES
  */
  std::string filter;

  /**
     @brief Information of the mutation
     @details
     @tableofcontents
     | Entry  | Variable  | Description |
     | -------|-----------|-------------|
     | CHR2   | @ref chr2 | Mate chromosome |
     | END    | @ref pos2 | Mate position |
     | MATEID | @ref mate_id | Mate ID |
     | (none) | @ref dir2 | Mate clip direction |
     | STRAND | @ref strand | Combined direction |
     | SVTYPE | @ref svtype | SV-types |
     | SVLEN  | @ref svlen  | SV length (0 for inter-chromosomal translocation) |
     @ingroup VCF_LINE_ATTRIBUTES
     @see @ref INFO_ENTRIES
  */
  std::string info;
  
  /**
     @brief String template vector of FORMAT
  */
  std::vector < std::string > format;

  
  /**
     @defgroup INFO_ENTRIES Entries for INFO
     @brief INFO (@ref info) field entries 
     @details
     @tableofcontents
     | Entry  | Variable | Description |
     | -------|----------|-------------|
     | CHR2   | @ref chr2 | Mate chromosome |
     | END    | @ref pos2 | Mate position |
     | MATEID | @ref mate_id | Mate ID |
     | (none) | @ref dir2 | Mate clip direction |
     | STRAND | @ref strand | Combined direction |
     | SVTYPE | @ref svtype | SV-types |
     | SVLEN  | @ref svlen  | SV length (0 for inter-chromosomal translocation) |
  */
  
  /// Chromosome of mate-BP
  std::string chr2;
  /// Position of mate-BP
  pos_type pos2;
  /// Direction of mate-BP
  pos_type dir2;

  /**
     @brief Combined strand (or directions of BPs)
     @details
     @tableofcontents
     | Strand | Description| Rearrangement path (REPATH) |
     | -------|------------|-----------------------------|
     | FF     | Tail-tail  | Chr10:45394261(+)-Chr5:38253193(-) |
     | FR     | Tail-head  | Chr12:32229443(+)-Chr12:32209976(+) |
     | RF     | Head-tail  | Chr12:32209976(-)-Chr12:32229443(-) |
     | RR     | Head-head  | Chr3:93511849(-)-Chr21:10746646(+) |
     | F      | Tail-(unknown)  | Chr1:70072208(+)-UNKNOWN |
     | R      | Head-(unknown)  | Chr8:124091696(-)-UNKNOWN |
  */
  std::string strand;  
  /// mate-BP id
  std::string mate_id;
  /**
     @brief SV-types such as DEL, DUP, INV, BND
     @details
     @tableofcontents
     | ID | Description|
     |----|------------|
     |DEL | Deletion |
     |DUP | Tandem duplication|
     |INV | Inversion |
     |TRA | Inter-chromosomal translocation|
     |BND | Breakend (a pair of BPs), or undetermined SV |
     |SND | Single breakend (without mate-BP) |
     @see alt
  */
  std::string svtype; 
  /// SV length (0 for inter-chromosomal translocation)
  pos_type svlen;

  

  /**
     @brief Reading input vcf line
     @param [in] line Input vcf lie
     @return Update @ref VCF_LINE_ATTRIBUTES
  */
  void read_vcf_line(std::string vcf_string);

  /**
     @brief Update INFO field with svtype
  */
  void update_svtype_in_info();

  /**
     @brief Update INFO field
  */
  void update_info();

  /**
     @brief Compile a vcf line
     @return Vcf line
  */
  std::string to_string();

  /**
     @brief Compile a vcf line
     @return Vcf line without FORMAT...
  */
  std::string to_string_short();

  /**
   * @brief Flip BPs
   * 
   */
  void flip();


  /**
   * @brief Check if second BP is empty
   * @return True if second BP is empty
   */
  bool is_empty_second();
  
  std::string binary_call_set();

  /**
   * @brief Extract ID SURVIVOR merge FORMAT
   * @details Replace FORMAT with ID
   */
  void extract_id_from_SURVIVOR_merge_format();

  /**
   * @brief Update FORMAT name
   * 
   */
  void update_format_name(int number);

  /**
   * @brief Update FORMAT name
   * 
   */
  void update_format_name_for_silver_standard(int number);
};

/**
   @brief Generate mate vcf line
   @param [in] input VCF_LINE
   @details Copy input to output
   @details swap ( chr1 , chr2 )
   @details swap ( pos1 , pos2 )
   @details strand "FR" <--> "FR" (Note: "FF" and "RR" are not altered)
   @return Mate VCF line
*/
VCF_LINE return_mate(VCF_LINE input);
std::string mate_alt ( VCF_LINE input );
std::string mate_info ( VCF_LINE input );

/**
   @class BP types.hpp "types.hpp"
   @brief A class of BP 
   @details BP is defined by a vector such that
   @details \f$ BP = ( c , x , s ) \f$
   @details where \f$c\f$, \f$x\f$, and \f$s\f$ are chromosome, position in a chromosome, and direction, respectively.
   @details This vecter is stored in a tuple @ref bp = std::tuple< pos_type, pos_type, pos_type > 
*/


class BP{
private:
public:
  /// Constructor
  BP(){init();}
  /**
     @brief Constructor
     @param [in,out] a Chromosome
     @param [in,out] b Position in a chromosome
     @param [in,out] c Direction
     @result Initialize tuple
  */
  BP(pos_type a, pos_type b, pos_type c){init();insert(a,b,c);}
  /**
     @brief Constructor
     @param [in,out] input Initialize bp
     @result Initialize @ref bp 
  */
  BP(BP_type input){init();insert(input);}

  /// Destructor
  ~BP(){};

  /**
     @brief A tuple for a BP
     @note bp = ( c , x , s )
     @note where
     @note c is chromosome
     @note x is position
     @note s is direction
  */
  BP_type bp;

  /**
     @brief Initialize bp with (-1,-1,0) for default
  */
  void init(){
    std::get<0>(this->bp)=-1;
    std::get<1>(this->bp)=-1;
    std::get<2>(this->bp)=0;
  }

  /// Remove data from bp
  void clear(){init();};

  /// Insert (a,b,c)
  void insert (pos_type a, pos_type b, pos_type c){this->bp=std::make_tuple (a,b,c);}

  /// Insert BP_type input to bp
  void insert (BP_type input){this->bp = input;}

  /// Return chr (c)
  pos_type chr() const {return std::get<0>(this->bp);}
  /// Return pos (x)
  pos_type pos() const {return std::get<1>(this->bp);}
  /// Return dir (s)
  pos_type dir() const {return std::get<2>(this->bp);}

  /// Update chr (c) with a
  void chr(pos_type a) {std::get<0>(this->bp)=a;}
  /// Update pos (x) with a
  void pos(pos_type b) {std::get<1>(this->bp)=b;}
  /// Update dir (s) with a
  void dir(pos_type c) {std::get<2>(this->bp)=c;}

  /// Return (c,x) without s
  Position Pos() const {return std::make_pair(std::get<0>(this->bp),std::get<1>(this->bp));};

  /// Filp direction (c,x,s) --> (c,x,-s)
  void flip(){std::get<2>(this->bp)*=-1;};

  /// Return 1 if bp is default (-1,-1,0 ), else 0
  bool is_empty(){return std::get<0>(this->bp)!=-1?0:std::get<1>(this->bp)!=-1?0:std::get<2>(this->bp)?0:1;}

  /// Return (c,x,s) in a format Chr:Pos:(+/-)
  std::string to_string() const {
    return std::to_string(chr()) + ":" + std::to_string(pos()) + ":" + std::to_string(dir());
  }

  BP operator + (pos_type x) const {BP a(*this); std::get<1>(a.bp) += x; return a;}
  BP operator - (pos_type x) const {BP a(*this); std::get<1>(a.bp) -= x; return a;}

  BP & operator += (pos_type x) {std::get<1>(this->bp) += x; return *this;}
  BP & operator -= (pos_type x) {std::get<1>(this->bp) -= x; return *this;}

  BP & operator =  (BP_type bp) {this->bp = bp;return *this;}
  
  bool operator == (BP input) const {return this->bp==input.bp;}
  bool operator != (BP input) const {return this->bp!=input.bp;}
  bool operator <  (BP input) const {return this->bp< input.bp;}
  bool operator >  (BP input) const {return this->bp> input.bp;}
  bool operator <= (BP input) const {return this->bp<=input.bp;}
  bool operator >= (BP input) const {return this->bp>=input.bp;}

};



/**
 * @brief A class of BP pair (BND) 
 * @details BND is defined by
 * @details \f$ BND = ( BP_1 , BP_2 ) \f$
 * @details where \f$BP_1\f$ and \f$BP_2\f$ are BP and mate-BP, respectively.
 * @note BND::first = BP_1
 * @note BND::second= BP_2
 * @details If BP_2 is empty, we call it single BND (SND).
 */

class BND{
private:
public:
  /// Constructor
  BND(){};
  /**
     @brief Constructor
     @param [in] first 
     @param [in] second
     @return BND::first = first
     @return BND::second= second
  */
  BND(BP first, BP second){
    this->first = first;
    this->second = second;
  };
  
  /**
     @brief Constructor
     @param [in] input BND
     @return BND::first = input.first
     @return BND::second= input.second
  */
  BND(BND_type input) {this->first=input.first;this->second=input.second;};

  /// Destructor
  ~BND(){};

  /// BP_1 (BP)
  BP first;
  /// BP_2 (mate-BP)
  BP second;
  

  /// 1 if BP_1 and BP_2 are empty simutaneously, else 0
  bool is_empty(){return this->first.is_empty() && this->second.is_empty() ; }

  /// Insert (bp1, bp2) 
  void insert (BP bp1, BP bp2){this->first = bp1; this->second = bp2;} 

  /// Insert bnd
  void insert (BND bnd){ *this = bnd;} 
  
  const bool operator == (BND input) const {return (first.bp==input.first.bp&&second.bp==input.second.bp);}
  const bool operator != (BND input) const {return ! (*this == input) ;}
  const bool operator <  (BND input) const {
    return this->first<input.first?1:this->first>input.first?0:this->second<input.second;
  }
  const bool operator >  (BND input) const {
    return this->first>input.first?1:this->first<input.first?0:this->second>input.second;
  }
  const bool operator <= (BND input) const {return ! (*this < input);}
  const bool operator >= (BND input) const {return ! (*this > input);}


  /// Return ( c , x ) for BP_1
  Position Pos1() const {
    return std::make_pair(this->first.chr(),this->first.pos());
  }
  
  /// Return ( c , x ) for BP_2
  Position Pos2() const {
    return std::make_pair(this->second.chr(),this->second.pos());
  }

  /// Return a string in a format chr1:pos1:dir1-chr2:pos2:dir2
  std::string to_string(){
    return first.to_string() + "-" + second.to_string() ;
  }
  
  pos_type chr1() const {return this->first.chr();}
  pos_type pos1() const {return this->first.pos();}
  pos_type dir1() const {return this->first.dir();}
  pos_type chr2() const {return this->second.chr();}
  pos_type pos2() const {return this->second.pos();}
  pos_type dir2() const {return this->second.dir();}
  
};



/**
   @brief Map for VCF_LINE
*/
using BP_MAP = std::map < BP, VCF_LINE >;

/**
   @brief Map for VCF_LINE
   @note VCF_MAP [ BP ] [ BP ] = VCF_LINE
*/
using VCF_MAP = std::map < BP , std::map < BP , VCF_LINE > >;



/**
 * @brief A class for vcf file
 */
class VCF{
private:
  /// Vector arrary of chromosome ID and size
  std::vector < std::pair < std::string , std::size_t > > genome_info;

  /**
   * @brief Build id-to-number and number-to-id maps
   * @param [in] infile VCF file with "contig" entry in header
   */
  void build_id_ref_map(const std::string infile);

  /**
   * @brief Update id-to-number and number-to-id maps
   * @param [in] chr Chromosome ID
   * @note It updates id_ref_map and ref_id_map if chr ID was not in id_ref_map.
   */
  void update_id_ref_map(std::string chr);  
  
  /**
   * @brief Return the nucleotide at the position
   * @param [in] chr Chromosome ID
   * @param [in] pos Position in the chromosome
   */
  char GetNucl(const std::string & chr , const std::size_t & pos );



  /**
     @brief Return first word from input string
     @param [in] input Fasta ID line
  */
  std::string get_first ( std::string input );

public:
  VCF();
  VCF(const VCF & vcf) ;
  VCF(const std::string input_file);
  VCF(const source_info source_info);
  VCF(const args CONF, const int i);
  ~VCF();

  void init();


  /**
     @brief Main container
     @details VCF_MAP = std::map < BP , std::map < BP , VCF_LINE > >
     @note vcf_map [ BP1 ] [ BP2 ] = vcf_line 
  */
  VCF_MAP vcf_map;
  
  /// @brief File name of vcf file
  std::string vcf_file;
  std::string tool_name;
  std::string outfile;

  std::string metainfo;
  std::string header;
  std::string short_header;

  std::map < std::string , int > id_ref_map;
  std::map < int , std::string > ref_id_map;
  std::map < std::string , std::string > genome;

  // int gen_TRA_mate;

  double cutoff;
  std::string cutoff_file;
  
  VCF_MAP::iterator begin();
  VCF_MAP::iterator end();
  BND find (BND bnd) ;
  BND find (BP bp1, BP bp2) ;

  std::size_t size();
  std::size_t above_cutoff_size();

  void get_header(const std::string infile);
  void read_vcf_file(const std::string infile);
  void read_vcf_file();

  void add_metainfo ( std::string key, std::string input, bool force);
  
  void clear();

  void get_genome(std::string);

  void write();
  void write_no_header();
  void fwrite(std::string);
  void fwrite();

  /**
   * @brief Print without FORMATs
   * 
   */
  void write_short();

  /**
   * @brief Write a file without FORMATs
   * @param [in] outfile Output file name
   */
  void fwrite_short(std::string);
  void fwrite_short();

  void rename_id(std::string );
  void insert ( VCF_LINE vcf_line );

  void add_missed_header_lines();

  BP_MAP & operator []( BP bp);

  /** remove imprecise SVs */
  void remove_imprecise();
  int noimprecise;

  /** remove germline SVs and add score (only for lumpy)*/
  void prepr_lumpy ();
  /** add score (only for delly) */
  void prepr_delly ();
  /** add score (only for manta) */
  void prepr_manta ();
  int manta_PR_plus_SR ( std::vector < std::string > format_vec );

  /** collect of "PASS" and "qual" */
  void prepr_gridss ();

  /**
   * @brief Add ID for novoBreak
   * @param [in] score_file Score file
   */
  void prepr_novoBreak();

  void prepr_tools ();

  void add_score_to_id();

  void add_binary_calls_to_id(std::string file_name);
  
  int get_format_size();
  
  void extract_id_from_SURVIVOR_merge_format();

  void update_format_name();

  void update_format_name_for_silver_standard();

  void remove_low_qual_records();

  void typing_SV ();
  void remove_mate ();
  void rearrange_bp ();
  void make_mate_for_TRA();
  /**
   * @brief Remove short SVs 
   * @param [in] min_sv_size Minimum SV size
   */
  void remove_short (int min_sv_size);
  void remove_large (int max_sv_size);
  void remove_inter_chromosome (int intra_SV);
  void remove_intra_chromosome (int inter_SV);

  /**
     @brief Return current data in local time
     @note You can replace it with current_date_time() in VCF::make_header_short() 
     @note to display vcf generating time in YYYY-MM-DD.HH:mm:ss format
     @return YYYY-MM-DD
  */
  std::string current_date();

  /** 
      @brief Return current data and time in local time
      @return YYYY-MM-DD.HH:mm:ss
  */
  std::string current_date_time();

};

#endif
