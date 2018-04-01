cmake_minimum_required(VERSION 3.8)
project(KGL)

set(CMAKE_CXX_STANDARD 14)


# add extra include directories
include_directories(/usr/include
                    kgl_library
                    kgl_app
                    /usr/local/lib/spdlog/include   # 3rd party logging library.
                    /usr/local/lib/seqan/include    # 3rd party sequence analysis library.
                    contrib/edlib) # 3rd party library for fast Leveshtein sequence matching.

# add boost library directories
#set(BOOST_ROOT /usr/lib/boost)
set(BOOST_INCLUDEDIR /usr/include/boost)
set(BOOST_LIBRARYDIR /usr/lib/boost)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
include_directories(${Boost_INCLUDE_DIRS})
find_package( Boost COMPONENTS system filesystem timer chrono REQUIRED )

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -s -DNDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -W -Wstrict-aliasing -pedantic -Wno-variadic-macros")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall -pthread -W -Wstrict-aliasing -pedantic -Wno-variadic-macros")


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)

set(LIBRARY_SOURCE_FILES
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
        kgl_library/kgl_variant_factory_compound.h
        kgl_library/kgl_variant_factory_compound_insert.cc
        kgl_library/kgl_variant_compound_insert.cc
        kgl_library/kgl_variant_compound_delete.cc
        kgl_library/kgl_alphabet_coding_dna5.h
        kgl_library/kgl_alphabet_coding_dna5.cc
        kgl_library/kgl_variant_factory.cc
        kgl_library/kgl_variant_evidence.h
        kgl_library/kgl_statistics_upgma.h
        kgl_library/kgl_statistics_upgma.cc
        kgl_library/kgl_variant_mutation.cc
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
        kgl_library/kgl_variant_mutation.h
        kgl_library/kgl_variant_db_mutation.cc
        kgl_library/kgl_variant_db_contig.h
        kgl_library/kgl_variant_db_genome.h
        kgl_library/kgl_sequence_offset.h
        kgl_library/kgl_sequence_offset.cc
        kgl_library/kgl_variant_mutation_offset.h
        kgl_library/kgl_variant_mutation_offset.cc
        kgl_library/kgl_variant_db_population.h
        kgl_library/kgl_table_impl.h
        kgl_library/kgl_table_impl.cc
        kgl_library/kgl_table_organism.h
        kgl_library/kgl_variant_factory_fbvcf_impl.h
        kgl_library/kgl_variant_factory_fbvcf_impl.cc
        kgl_library/kgl_variant_factory_gatkvcf_impl.h
        kgl_library/kgl_variant_factory_gatkvcf_impl.cc
        kgl_library/kgl_variant_factory_vcf_impl.h
        kgl_library/kgl_variant_factory_vcf_impl.cc
        kgl_library/kgl_gaf_parser.h
        kgl_library/kgl_gaf_parser.cc
        kgl_library/kgl_sequence_compare.h
        kgl_library/kgl_sequence_compare_impl.h
        kgl_library/kgl_sequence_compare_impl.cc
        kgl_library/kgl_sequence_distance.h
        kgl_library/kgl_genome_types.cc
        kgl_library/kgl_sequence_distance_impl.h
        kgl_library/kgl_sequence_distance_impl.cc
        kgl_library/kgl_variant_factory_pf3k_impl.h
        kgl_library/kgl_variant_factory_pf3k_impl.cc
        kgl_library/kgl_variant_factory_readvcf_impl.h
        kgl_library/kgl_variant_factory_genome_vcf_impl.h
        kgl_library/kgl_variant_factory_genome_vcf_impl.cc
        kgl_library/kgl_variant_factory_vcf_parse_impl.h
        kgl_library/kgl_variant_factory_vcf_parse_impl.cc
        kgl_library/kgl_variant_factory_record_vcf_impl.h
        kgl_library/kgl_variant_factory_record_vcf_impl.cc
        kgl_library/kgl_ploidy_vcf_Pf3k_impl.cc
        kgl_library/kgl_ploidy_vcf_Pf3k_impl.h
        kgl_library/kgl_variant_factory_vcf_cigar_impl.h
        kgl_library/kgl_variant_factory_vcf_cigar_impl.cc
        kgl_library/kgl_genome_aux_csv.h
        kgl_library/kgl_genome_aux_csv_impl.cc)

# 3rd party library for fast Leveshtein sequence matching.
set(EDLIB_SOURCE_FILES contrib/edlib/edlib.cpp kgl_library/kgl_variant_db_genotype.h)

set(APPLICATION_SOURCE_FILES
        kgl_app/kgl_main.cc
        kgl_app/kgl_phylogenetic_env.h
        kgl_app/kgl_phylogenetic_env.cc
        kgl_app/kgl_phylogenetic_app.cc
        kgl_app/kgl_phylogenetic_analysis.h
        kgl_app/kgl_phylogenetic_analysis.cc
        kgl_app/kgl_phylogenetic_gene.h
        kgl_app/kgl_phylogenetic_gene.cc
        kgl_app/kgl_upgma.h
        kgl_app/kgl_upgma.cc
        kgl_app/kgl_phylogenetic_app_analysis.cc
        kgl_app/kgl_phylogenetic_app.h
        kgl_app/kgl_rna_search.cc
        kgl_app/kgl_rna_search.h
        kgl_app/kgl_ploidy_analysis.h
        kgl_app/kgl_ploidy_analysis.cc)


## DOXYGEN
#set(BUILD_DOC OFF CACHE BOOL "Use Doxygen to create the HTML/Latex based API documentation")
#option(BUILD_DOC "Use Doxygen to create the HTML/Latex based API documentation" ON)
if(BUILD_DOC)

    find_package(Doxygen)

    if(DOXYGEN_FOUND)


        # set input and output files
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.out)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM )

#        add_custom_target( doc_doxygen
#                    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                    COMMENT "Generating API documentation with Doxygen"
#                    VERBATIM )

    endif(DOXYGEN_FOUND)

endif(BUILD_DOC)


#generate shared library.
add_library(kgl_lib SHARED ${LIBRARY_SOURCE_FILES} ${EDLIB_SOURCE_FILES})

target_link_libraries(kgl_lib  ${Boost_LIBRARIES})

#generate kgl_genome executable
add_executable (kgl_genome ${APPLICATION_SOURCE_FILES})

target_link_libraries(kgl_genome kgl_lib)
