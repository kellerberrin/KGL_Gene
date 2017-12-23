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
#set(BOOST_ROOT /usr/lib/boost)
set(BOOST_INCLUDEDIR /usr/include/boost)
set(BOOST_LIBRARYDIR /usr/lib/boost)
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

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -s -DNDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall")


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)

#generate shared python library.
add_library(read_sam SHARED ${SOURCE_FILES} kgl_python.cc)

#generate kgl_genome executable
add_executable (kgl_genome ${SOURCE_FILES} kgl_main.cc  kgl_gff_fasta.h kgl_gff_fasta.cc
                kgl_genome_db.h kgl_patterns.h kgl_genome_feature.h kgl_genome_feature.cc
                kgl_exec_env_boost.cc kgl_variant.h kgl_variant.cc
                kgl_filter.h kgl_filter.cc kgl_variant_set.cc
                kgl_table.h kgl_table.cc kgl_sam_parse.h kgl_sam_parse.cc
                kgl_sam_read.h kgl_sam_consume.h kgl_sam_process.h kgl_table_ncbi.h kgl_table_ncbi.cc
                kgl_sequence_base.h kgl_sequence_base.cc kgl_sequence_amino.h kgl_sequence_amino.cc
                kgl_alphabet_amino.h kgl_variant_db.h
                kgl_variant_db.cc kgl_variant_snp.cc kgl_variant_compound.h
                kgl_variant_compound.cc kgl_variant_compound_snp.cc kgl_genome_prelim.h
                kgl_genome_prelim.cc kgl_attributes.h kgl_attributes.cc kgl_phylogenetic_app.cc
                kgl_phylogenetic_env.h kgl_phylogenetic_env.cc kgl_sequence_virtual.h kgl_genome_contig.cc
                kgl_genome_verify.cc kgl_genome_genome.cc kgl_variant_io.cc kgl_sequence_codon.h
                kgl_sequence_codon.cc kgl_alphabet_string.h kgl_alphabet_amino.cc kgl_alphabet_dna5.h
                kgl_alphabet_dna5.cc kgl_alphabet_extend.h kgl_alphabet_extend.cc kgl_phylogenetic_analysis.h
                kgl_phylogenetic_analysis.cc kgl_variant_factory.h
        kgl_variant_factory_compound.cc kgl_variant_factory_compound_delete.cc
                kgl_variant_factory_compound_snp.cc kgl_variant_factory_compound.h
                kgl_variant_factory_compound_insert.cc kgl_variant_compound_insert.cc
                kgl_variant_compound_delete.cc kgl_alphabet_coding_dna5.h kgl_alphabet_coding_dna5.cc
                kgl_variant_factory.cc kgl_statistics.h kgl_statistics.cc kgl_variant_evidence.h
                kgl_variant_evidence.cc kgl_statistics_upgma.h kgl_statistics_upgma.cc kgl_variant_mutation.cc
                kgl_sequence_virtual_compare.h
        kgl_variant_single.h kgl_variant_delete.cc kgl_variant_insert.cc kgl_variant_single.cc kgl_variant_factory_single.cc)

target_link_libraries(kgl_genome ${Boost_LIBRARIES})
