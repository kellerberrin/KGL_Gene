cmake_minimum_required(VERSION 3.8)
project(KGL)

set(CMAKE_CXX_STANDARD 14)


# add extra include directories
include_directories(/usr/local/lib/spdlog/include /usr/local/lib/seqan/include /usr/include ./kgl_library)

# add boost library directories
#set(BOOST_ROOT /usr/lib/boost)
set(BOOST_INCLUDEDIR /usr/include/boost)
set(BOOST_LIBRARYDIR /usr/lib/boost)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
include_directories(${Boost_INCLUDE_DIRS})
find_package( Boost COMPONENTS system filesystem timer chrono REQUIRED )

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pg -pthread")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall -pg")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -s -DNDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -W -Wstrict-aliasing -pedantic -Wno-variadic-macros")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall -pthread -W -Wstrict-aliasing -pedantic -Wno-variadic-macros")


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)

set(SOURCE_FILES
        kgl_library/kgl_logging.h
        kgl_library/kgl_logging.cc
        kgl_library/kgl_mt_queue.h
        kgl_library/kgl_mt_data.h
        kgl_library/kgl_genome_types.h
        kgl_library/kgl_lock.h
        kgl_library/kgl_mt_numpy.h
        kgl_library/kgl_mt_contig.h
        kgl_library/kgl_exec_env.cc
        kgl_library/kgl_mt_insert.h
        kgl_library/kgl_gff_fasta.h
        kgl_library/kgl_gff_fasta.cc
        kgl_library/kgl_genome_db.h
        kgl_library/kgl_patterns.h
        kgl_library/kgl_genome_feature.h
        kgl_library/kgl_genome_feature.cc
        kgl_library/kgl_variant.h
        kgl_library/kgl_variant.cc
        kgl_library/kgl_filter.h
        kgl_library/kgl_filter.cc
        kgl_library/kgl_variant_set.cc
        kgl_library/kgl_table.h
        kgl_library/kgl_table.cc
        kgl_library/kgl_sam_parse.h
        kgl_library/kgl_sam_parse.cc
        kgl_library/kgl_sam_read.h
        kgl_library/kgl_sam_consume.h
        kgl_library/kgl_sam_process.h
        kgl_library/kgl_table_ncbi.h
        kgl_library/kgl_sequence_base.h
        kgl_library/kgl_sequence_base.cc
        kgl_library/kgl_sequence_amino.h
        kgl_library/kgl_sequence_amino.cc
        kgl_library/kgl_alphabet_amino.h
        kgl_library/kgl_variant_db.h
        kgl_library/kgl_variant_db.cc
        kgl_library/kgl_variant_snp.cc
        kgl_library/kgl_variant_compound.h
        kgl_library/kgl_variant_compound.cc
        kgl_library/kgl_variant_compound_snp.cc
        kgl_library/kgl_genome_prelim.h
        kgl_library/kgl_genome_prelim.cc
        kgl_library/kgl_attributes.h
        kgl_library/kgl_attributes.cc
        kgl_library/kgl_sequence_virtual.h
        kgl_library/kgl_genome_contig.cc
        kgl_library/kgl_genome_verify.cc
        kgl_library/kgl_genome_genome.cc
        kgl_library/kgl_variant_io.cc
        kgl_library/kgl_sequence_codon.h
        kgl_library/kgl_sequence_codon.cc
        kgl_library/kgl_alphabet_string.h
        kgl_library/kgl_alphabet_amino.cc
        kgl_library/kgl_alphabet_dna5.h
        kgl_library/kgl_alphabet_dna5.cc
        kgl_library/kgl_variant_factory.h
        kgl_library/kgl_variant_factory_compound.cc
        kgl_library/kgl_variant_factory_compound_delete.cc
        kgl_library/kgl_variant_factory_compound_snp.cc
        kgl_library/kgl_variant_factory_compound.h
        kgl_library/kgl_variant_factory_compound_insert.cc
        kgl_library/kgl_variant_compound_insert.cc
        kgl_library/kgl_variant_compound_delete.cc
        kgl_library/kgl_alphabet_coding_dna5.h
        kgl_library/kgl_alphabet_coding_dna5.cc
        kgl_library/kgl_variant_factory.cc
        kgl_library/kgl_statistics.h
        kgl_library/kgl_statistics.cc
        kgl_library/kgl_variant_evidence.h
        kgl_library/kgl_statistics_upgma.h
        kgl_library/kgl_statistics_upgma.cc
        kgl_library/kgl_variant_mutation.cc
        kgl_library/kgl_sequence_virtual_compare.h
        kgl_library/kgl_variant_single.h
        kgl_library/kgl_variant_delete.cc
        kgl_library/kgl_variant_insert.cc
        kgl_library/kgl_variant_single.cc
        kgl_library/kgl_variant_factory_single.cc
        kgl_library/kgl_alphabet_readcount.cc
        kgl_library/kgl_alphabet_readcount.h
        kgl_library/kgl_alphabet_extend.h
        kgl_library/kgl_alphabet_extend.cc
        kgl_library/kgl_variant_factory_single.h
        kgl_library/kgl_variant_factory_vcf.h
        kgl_library/kgl_utility.h
        kgl_library/kgl_utility.cc
        kgl_library/kgl_variant_factory_vcf.cc
        kgl_library/kgl_variant_factory_bam.h
        kgl_library/kgl_variant_factory_bam.cc
        kgl_library/kgl_variant_factory_vcf_impl.h
        kgl_library/kgl_variant_factory_vcf_impl.cc
        kgl_library/kgl_variant_factory_vcf_utils.cc
        kgl_library/kgl_variant_mutation.h
        kgl_library/kgl_variant_db_mutation.cc
        kgl_library/kgl_variant_db_contig.h
        kgl_library/kgl_variant_db_genome.h
        kgl_library/kgl_sequence_offset.h
        kgl_library/kgl_sequence_offset.cc
        kgl_library/kgl_variant_mutation_offset.h
        kgl_library/kgl_variant_mutation_offset.cc
        kgl_library/kgl_sequence_manip.h
        kgl_library/kgl_sequence_manip.cc
        kgl_library/kgl_variant_factory_compound_noncoding.cc
        kgl_library/kgl_variant_db_population.h
        kgl_library/kgl_table_impl.h
        kgl_library/kgl_table_impl.cc
        kgl_library/kgl_table_organism.h)

#generate shared python library.
add_library(kgl_lib SHARED ${SOURCE_FILES})

target_link_libraries(kgl_lib ${Boost_LIBRARIES})

#generate kgl_genome executable
add_executable (kgl_genome kgl_main.cc kgl_phylogenetic_env.h kgl_phylogenetic_env.cc kgl_phylogenetic_app.cc
                kgl_phylogenetic_analysis.h kgl_phylogenetic_analysis.cc kgl_phylogenetic_gene.h kgl_phylogenetic_gene.cc)

target_link_libraries(kgl_genome kgl_lib)
