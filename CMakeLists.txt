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
set(SOURCE_FILES kgl_logging.h kgl_logging.cc kgl_read_sam.h
        kgl_mt_queue.h kgl_mt_data.h kgl_genome_types.h kgl_consume_sam.h
        kgl_process_sam.h kgl_lock.h kgl_mt_numpy.h kgl_nucleotide.h kgl_parse_sam.h
        kgl_parse_sam.cc kgl_mt_contig.h kgl_exec_env.cc kgl_mt_insert.h)

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
                kgl_minority_env.h kgl_minority_app.cc kgl_nucleotide.cc kgl_amino.h kgl_table.h kgl_base_sequence.h kgl_table.cc kgl_mutant_sequence.h kgl_mutant_sequence.cc kgl_base_sequence.cc kgl_ncbi_table.h kgl_ncbi_table.cc kgl_amino_sequence.h kgl_amino_sequence.cc kgl_variant_evidence.h kgl_variant_evidence.cc)
target_link_libraries(kgl_genome ${Boost_LIBRARIES})
