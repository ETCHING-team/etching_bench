/**
 * @file define.hpp
 * @author Jang-il Sohn (sohnjangil@gmail.com)
 * @brief 
 * @version 0.2.11
 * @date 2022-04-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#define PROGRAM_NAME "etching_bench"
#define PROGRAM_VERSION "0.2.11"

#define DEFAULT_PREFIX PROGRAM_NAME

#define ETCHING_TOOL "etching"
#define DELLY_TOOL "delly"
#define LUMPY_TOOL "lumpy"
#define MANTA_TOOL "manta"
#define SVABA_TOOL "svaba"
#define NOVOBREAK_TOOL "novobreak"
#define GRIDSS_TOOL "gridss"

#define TOOL _TOOL
#define CUTOFF _CUTOFF
#define SCORE_TYPE _SCORE_TYPE

#define ETCHING_CUTOFF 0.4
#define DELLY_CUTOFF 18
#define LUMPY_CUTOFF 12
#define MANTA_CUTOFF 40
#define SVABA_CUTOFF 0
#define NOVOBREAK_CUTOFF 40
#define GRIDSS_CUTOFF 400

#define CONFIG_FILE ""
#define OUTFILE_PREFIX ""
#define TRUTH_FILE ""
#define MERGE_LIST_FILE ""
#define MERGE_WINDOW 10
#define MIN_SV_SIZE 100
#define MAX_SV_SIZE 0
#define HAS_TRUTH false

#define CONSENSUS_CUTOFF 3

#define ETCHING_SCORE_TYPE "float"
#define DELLY_SCORE_TYPE "int"
#define LUMPY_SCORE_TYPE "int"
#define MANTA_SCORE_TYPE "int"
#define SVABA_SCORE_TYPE "int"
#define NOVOBREAK_SCORE_TYPE "float"
#define GRIDSS_SCORE_TYPE "int"

#define CHECK_HELP(A, B) {if ( argc == 1 || any_help(argc, argv)){USAGE();}}
#define CERR_ERROR(msg) {std::cerr << "ERROR: " << msg << std::endl;}
#define CLOSE(A) {if (A.is_open()){A.close();}}
#define IFS std::ifstream
#define OFS std::ofstream
#define STS std::stringstream
#define STR std::string
#define VEC std::vector
#define COUT std::cout
#define CERR std::cerr
#define GETLINE std::getline
#define CHECK_FIN(A) {IFS fin(A); if (!fin.is_open()){std::cerr << "ERROR: Cannot open the file: " << A << "\n"; printf("CHECK...%s(%d)::%s\n", __FILE__, __LINE__, __FUNCTION__);exit(1);} else CLOSE(fin);}
#define CHECK_FOUT(A) {OFS fout(A); if (!fout.is_open()){std::cerr << "ERROR: Cannot open the file: " << A << "\n"; printf("CHECK...%s(%d)::%s\n", __FILE__, __LINE__, __FUNCTION__);exit(1);} else CLOSE(fout);}

#define CHECK_CONFIG {std::cerr << "Please check the configuration file: " << config_file << "\n"; exit(1);}
