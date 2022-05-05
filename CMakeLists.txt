cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 11)

project(etching_bench
  VERSION 0.2.11
  DESCRIPTION "Benchmarking program for ETCHING project"
  LANGUAGES CXX)

set(EXPORT_NAME ${PROJECT_NAME}-config)

add_library (
benchlib
STATIC
${CMAKE_SOURCE_DIR}/src/args.cpp
${CMAKE_SOURCE_DIR}/src/bench.cpp
${CMAKE_SOURCE_DIR}/src/instance.cpp
${CMAKE_SOURCE_DIR}/src/main.cpp
${CMAKE_SOURCE_DIR}/src/merge.cpp
${CMAKE_SOURCE_DIR}/src/prepr.cpp
${CMAKE_SOURCE_DIR}/src/prepr_tools.cpp
${CMAKE_SOURCE_DIR}/src/prepr_typer.cpp
${CMAKE_SOURCE_DIR}/src/prepr_vcf.cpp
)


add_library (
survivorlib
STATIC
${CMAKE_SOURCE_DIR}/src/instance.cpp
${CMAKE_SOURCE_DIR}/src/SURVIVOR/src/merge_vcf/IntervallTree.cpp
${CMAKE_SOURCE_DIR}/src/SURVIVOR/src/merge_vcf/combine_svs.cpp
${CMAKE_SOURCE_DIR}/src/SURVIVOR/src/vcfs/Merge_VCF.cpp
${CMAKE_SOURCE_DIR}/src/SURVIVOR/src/simulator/Eval_vcf.cpp
)

add_executable (etching_bench ${CMAKE_SOURCE_DIR}/src/main.cpp)

target_link_libraries ( etching_bench benchlib survivorlib )
