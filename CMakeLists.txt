cmake_minimum_required(VERSION 3.8)
project(ReadSamFile)

set(CMAKE_CXX_STANDARD 14)

# add extra include directories
include_directories(/usr/local/lib/pybind11/include
                    /usr/local/lib/spdlog/include
                    /usr/include/python2.7
                    /usr/local/lib/seqan/include
                    /usr/include)

# add boost library directories
set(BOOST_ROOT /usr/local/include/boost)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
include_directories(${Boost_INCLUDE_DIRS})
find_package( Boost COMPONENTS system filesystem timer chrono REQUIRED )

#setup common source directories
set(SOURCE_FILES kgl_logging.h kgl_logging.cc
        kgl_mt_queue.h kgl_mt_data.h kgl_genome_types.h
        kgl_lock.h kgl_mt_numpy.h
        kgl_mt_contig.h kgl_exec_env.cc kgl_mt_insert.h)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pg -pthread")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall -pg")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall")


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)

#generate shared python library.
add_library(read_sam SHARED ${SOURCE_FILES} kgl_python.cc)

#generate kgl_genome executable
add_executable (kgl_genome ${SOURCE_FILES} kgl_main.cc  kgl_gff_fasta.h kgl_gff_fasta.cc
                kgl_genome_db.h kgl_genome_db.cc kgl_patterns.h kgl_genome_feature.h kgl_genome_feature.cc
                kgl_exec_env_boost.cc kgl_variant.h kgl_variant.cc
                kgl_filter.h kgl_filter.cc kgl_variant_set.cc kgl_minority_env.cc
                kgl_minority_env.h kgl_minority_app.cc kgl_table.h kgl_table.cc
                kgl_variant_evidence.h kgl_variant_evidence.cc kgl_sam_parse.h kgl_sam_parse.cc
                kgl_sam_read.h kgl_sam_consume.h kgl_sam_process.h kgl_table_ncbi.h kgl_table_ncbi.cc
                kgl_sequence_base.h kgl_sequence_base.cc kgl_sequence_amino.h kgl_sequence_amino.cc
                kgl_sequence_mutant.h kgl_sequence_mutant.cc kgl_alphabet_base.h kgl_alphabet_base.cc
                kgl_alphabet_amino.h kgl_variant_db.h kgl_variant_db.cc kgl_variant_snp.h kgl_variant_snp.cc
                kgl_variant_compound.h kgl_variant_compound.cc)

target_link_libraries(kgl_genome ${Boost_LIBRARIES})
