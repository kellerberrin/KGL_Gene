cmake_minimum_required(VERSION 3.17)
project(KGL)

set(CMAKE_CXX_STANDARD 23)


# add extra include directories
include_directories(/usr/local/include/spdlog/include/
                    contrib/edlib # 3rd party library for fast Levenshtein sequence matching.
                    /usr/local/include/nlopt
                    /usr/local/include/libhmsbeagle-1  # the beagle mcmc library used by kpl
                    /usr/include/eigen3
                    kel_utility         # Basic Infrastructure library
                    kel_thread          # multi-threading library
                    kel_io              # file and network io library
                    kgl_genomics        # Basic Genetic library
                    kel_app             # Application and logging
                    kel_math            # Mathematical Code
                    kgl_mixture                  # Deconvolution and mixture library
                    kgl_app
                    kgl_genomics/kgl_parser
                    kgl_genomics/kgl_evidence
                    kgl_genomics/kgl_sequence
                    kgl_genomics/kgl_database
                    kgl_genomics/kgl_genome
                    kgl_genomics/kgl_genome_io
                    kgl_genomics/kgl_variant_db
                    kgl_genomics/kgl_variant_filter
                    kgl_genomics/kgl_variant_analysis
                    kgl_genomics/kgl_literature
                    kgl_genomics/kgl_mutation
                    kol_ontology                  # Experimental Ontology Library
                    kol_ontology/kgl_ontology     # Higher level ontology functionality
                    kga_analytic/kga_inbreed
                    kga_analytic/kga_sequence_analysis
                    kga_analytic/kga_info
                    kga_analytic/kga_PfEMP
                    kga_analytic/kga_legacy
                    kga_analytic/kga_literature
                    kga_analytic/kga_mutation
                    kga_analytic/kga_template_analysis
                    kga_analytic                  # Genetic analysis library (top level)
                    contrib/rapidjson/include # 3rd party Json parser
                    contrib/rapidxml) # 3rd party XML parser

link_directories(/usr/local/lib
                 /usr/local/lib/ncl
                 /usr/local/lib/boost
                 /usr/local/lib/libhmsbeagle-1/lib
                 /usr/local/lib/nlopt)


# Add boost library directories
set(BOOST_INCLUDEDIR /usr/local/include)
include_directories(${Boost_INCLUDE_DIRS})
set(BOOST_LIBRARYDIR /usr/local/lib/boost)
# Specify only static libraries for portability
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
# Find the Boost packages used.
find_package( Boost COMPONENTS filesystem system timer chrono program_options iostreams REQUIRED )

# zlib for compressed VCF files.
find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

# Debug g++
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++23 -fPIC -fconcepts -ggdb")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-ggdb")

# Release g++
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++23 -O3 -s -DNDEBUG")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE " -O3 -s -DNDEBUG")

# Profile g++
set(CMAKE_CXX_FLAGS_PROFILE "-std=c++23 -fPIC -fconcepts -pg -g -ggdb -O3 -Xlinker -Map=output.map ")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g -ggdb -O3 -Xlinker -Map=output.map ")

# All possible warnings
#set(CMAKE_ERROR_WARNING_FLAGS "-std=c++20 -Wall -pthread -Wextra -Werror -Wpedantic -pedantic-errors -Wstrict-aliasing -Wno-variadic-macros")
# Relax slightly to temporarily allow unused parameters.
set(CMAKE_ERROR_WARNING_FLAGS "-std=c++23 -Wall -pthread -Wextra -Wpedantic -Wstrict-aliasing -Wno-variadic-macros")

# Actually set the g++ flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)

#Basic application setup with logger.
set(RUNTIME_SOURCE_FILES
        kel_app/kel_exec_env.h
        kel_app/kel_exec_env.cpp
        kel_app/kel_exec_env_app.h
        kel_app/kel_logging.h
        kel_app/kel_logging.cpp
        kel_app/kel_property_tree.cpp
        kel_app/kel_property_tree.h
)

set(MATHEMATICAL_SOURCE_FILES
        kel_math/kel_distribution.h
        kel_math/kel_distribution.cpp
        kel_math/kel_percentile.h
        kel_math/kel_optimize.cpp
        kel_math/kel_optimize.h)

# Basic infrastructure library
set(UTILITY_SOURCE_FILES
        kel_utility/kel_utility.h
        kel_utility/kel_utility.cpp
        kel_utility/kel_mem_alloc.h
        kel_utility/kel_mem_alloc.cpp
        kel_utility/kel_date_time.cpp
        kel_utility/kel_date_time.h
        kel_utility/kel_interval_type.h
        kel_utility/kel_interval_unsigned.h
        kel_utility/kel_interval_map.h
        kel_utility/kel_interval_set.h
        kel_utility/kel_interval_set.cpp)

# Multi-threading library
set(THREAD_SOURCE_FILES
        kel_thread/kel_queue_mt_safe.h
        kel_thread/kel_workflow_threads.h
        kel_thread/kel_queue_tidal.h
        kel_thread/kel_workflow_async.h
        kel_thread/kel_queue_monitor.h
        kel_thread/kel_workflow_example.cpp
        kel_thread/kel_workflow_pipeline.h
        kel_thread/kel_movefunction.h)

# The basic I/O library.
set(IO_SOURCE_FILES
        kel_io/kel_basic_io.cpp
        kel_io/kel_basic_io.h
        kel_io/kel_rest_api.cpp
        kel_io/kel_rest_api.h
        kel_io/kel_bzip_workflow.cpp
        kel_io/kel_bzip_workflow.h
        kel_io/kel_bzip_workflow_verify.cpp
        kel_io/kel_mt_buffer.cpp
        kel_io/kel_mt_buffer.h
        kel_io/kel_file_io.cpp
        kel_io/kel_file_io.h)

# Basic Genetic infrastructure library
set(GENETIC_LIBRARY_SOURCE_FILES
        kgl_genomics/kgl_genome/kgl_genome_types.h
        kgl_genomics/kgl_genome_io/kgl_io_gff_fasta.h
        kgl_genomics/kgl_genome_io/kgl_io_gff_fasta.cpp
        kgl_genomics/kgl_genome/kgl_genome_feature.h
        kgl_genomics/kgl_genome/kgl_genome_feature.cpp
        kgl_genomics/kgl_variant_db/kgl_variant_db.h
        kgl_genomics/kgl_variant_db/kgl_variant_db.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_variant.h
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_variant.cpp
        kgl_genomics/kgl_sequence/kgl_table.h
        kgl_genomics/kgl_sequence/kgl_table.cpp
        kgl_genomics/kgl_sequence/kgl_table_ncbi.h
        kgl_genomics/kgl_sequence/kgl_sequence_base.h
        kgl_genomics/kgl_sequence/kgl_sequence_base.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_amino.h
        kgl_genomics/kgl_sequence/kgl_sequence_amino.cpp
        kgl_genomics/kgl_sequence/kgl_alphabet_amino.h
        kgl_genomics/kgl_genome/kgl_genome_prelim.h
        kgl_genomics/kgl_genome/kgl_genome_prelim.cpp
        kgl_genomics/kgl_genome/kgl_genome_attributes.h
        kgl_genomics/kgl_genome/kgl_genome_attributes.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_virtual.h
        kgl_genomics/kgl_genome/kgl_genome_contig.cpp
        kgl_genomics/kgl_genome/kgl_genome_verify.cpp
        kgl_genomics/kgl_genome/kgl_genome_genome.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_codon.h
        kgl_genomics/kgl_sequence/kgl_sequence_codon.cpp
        kgl_genomics/kgl_sequence/kgl_alphabet_string.h
        kgl_genomics/kgl_sequence/kgl_alphabet_amino.cpp
        kgl_genomics/kgl_sequence/kgl_alphabet_dna5.h
        kgl_genomics/kgl_sequence/kgl_alphabet_dna5.cpp
        kgl_genomics/kgl_sequence/kgl_alphabet_coding_dna5.h
        kgl_genomics/kgl_sequence/kgl_alphabet_coding_dna5.cpp
        kgl_genomics/kgl_evidence/kgl_variant_evidence.h
        kgl_genomics/kgl_evidence/kgl_variant_evidence.cpp
        kgl_genomics/kgl_sequence/kgl_statistics_upgma.h
        kgl_genomics/kgl_sequence/kgl_statistics_upgma.cpp
        kgl_genomics/kgl_sequence/kgl_table_impl.h
        kgl_genomics/kgl_sequence/kgl_table_impl.cpp
        kgl_genomics/kgl_sequence/kgl_table_organism.h
        kgl_genomics/kgl_parser/kgl_gaf_parser.h
        kgl_genomics/kgl_parser/kgl_gaf_parser.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_compare.h
        kgl_genomics/kgl_sequence/kgl_sequence_compare_impl.h
        kgl_genomics/kgl_sequence/kgl_sequence_compare_impl.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_distance.h
        kgl_genomics/kgl_sequence/kgl_sequence_distance_impl.h
        kgl_genomics/kgl_sequence/kgl_sequence_distance_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_pf_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_pf_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_readvcf_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_cigar.h
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_cigar.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_record_vcf_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_record_vcf_impl.cpp
        kgl_genomics/kgl_parser/kgl_pfgenome_aux.h
        kgl_genomics/kgl_parser/kgl_pfgenome_aux.cpp
        kgl_genomics/kgl_variant_db/kgl_variant_db_genome.h
        kgl_genomics/kgl_variant_db/kgl_variant_db_population.h
        kgl_genomics/kgl_genome/kgl_genome_contig_feature.h
        kgl_genomics/kgl_genome/kgl_genome_contig_feature.cpp
        kgl_genomics/kgl_genome/kgl_genome_contig_aux.cpp
        kgl_genomics/kgl_genome/kgl_genome_contig_aux.h
        kgl_genomics/kgl_sequence/kgl_phylogenetic_tree.cpp
        kgl_genomics/kgl_sequence/kgl_phylogenetic_tree.h
        kgl_genomics/kgl_parser/kgl_variant_vcf_record.h
        kgl_genomics/kgl_parser/kgl_variant_factory_readvcf_impl.cpp
        kgl_genomics/kgl_parser/kgl_data_file_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_grch_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_grch_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_info.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_info.h
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_header.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_header.h
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence.h
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_data.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_data.h
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_memory.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_memory.h
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_data_blk.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_data_blk.h
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_data_blk_read.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_analysis.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_analysis.h
        kgl_genomics/kgl_parser/kgl_variant_factory_1000_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_1000_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_population.h
        kgl_genomics/kgl_variant_db/kgl_variant_db_offset.h
        kgl_genomics/kgl_variant_db/kgl_variant_db_contig.h
        kgl_genomics/kgl_parser/kgl_hsgenealogy_parser.h
        kgl_genomics/kgl_parser/kgl_hsgenealogy_parser.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_parsers.h
        kgl_genomics/kgl_parser/kgl_variant_factory_gnomad_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_gnomad_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_parsers.cpp
        kgl_genomics/kgl_variant_db/kgl_variant_db_contig.cpp
        kgl_genomics/kgl_variant_db/kgl_variant_db_genome.cpp
        kgl_genomics/kgl_variant_db/kgl_variant_db_population.cpp
        kgl_genomics/kgl_parser/kgl_data_file_type.h
        kgl_genomics/kgl_genome/kgl_genome_contig.h
        kgl_genomics/kgl_parser/kgl_data_file_type.cpp
        kgl_genomics/kgl_parser/kgl_pf3k_coi.cpp
        kgl_genomics/kgl_parser/kgl_pf3k_coi.h
        kgl_genomics/kgl_genome/kgl_genome_genome.h
        kgl_genomics/kgl_parser/kgl_square_parser.cpp
        kgl_genomics/kgl_parser/kgl_square_parser.h
        kgl_genomics/kgl_variant_db/kgl_variant_db_offset.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_virtual.h
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_mem_alloc.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_mem_alloc.h
        kgl_genomics/kgl_parser/kgl_ensembl_id_parser.cpp
        kgl_genomics/kgl_parser/kgl_ensembl_id_parser.h
        kgl_genomics/kgl_genome/kgl_genome_collection.h
        kgl_genomics/kgl_genome/kgl_genome_collection.cpp
        kgl_genomics/kgl_variant_analysis/kgl_variant_sort.cpp
        kgl_genomics/kgl_variant_analysis/kgl_variant_sort.h
        kgl_genomics/kgl_parser/kgl_hsgenome_aux.cpp
        kgl_genomics/kgl_parser/kgl_hsgenome_aux.h
        kgl_genomics/kgl_parser/kgl_uniprot_parser.cpp
        kgl_genomics/kgl_parser/kgl_uniprot_parser.h
        kgl_genomics/kgl_variant_analysis/kgl_variant_sort_analysis.cpp
        kgl_genomics/kgl_variant_analysis/kgl_variant_sort_analysis.h
        kgl_genomics/kgl_parser/kgl_json_parser.cpp
        kgl_genomics/kgl_parser/kgl_json_parser.h
        kgl_genomics/kgl_parser/kgl_citation_parser.cpp
        kgl_genomics/kgl_parser/kgl_citation_parser.h
        kgl_genomics/kgl_variant_db/kgl_variant_db_freq.cpp
        kgl_genomics/kgl_variant_db/kgl_variant_db_freq.h
        kgl_genomics/kgl_parser/kgl_entrez_parser.cpp
        kgl_genomics/kgl_parser/kgl_entrez_parser.h
        kgl_genomics/kgl_parser/kgl_bio_pmid_parser.cpp
        kgl_genomics/kgl_parser/kgl_bio_pmid_parser.h
        kgl_genomics/kgl_parser/kgl_pubmed_resource.h
        kgl_genomics/kgl_literature/kgl_literature.cpp
        kgl_genomics/kgl_literature/kgl_literature.h
        kgl_genomics/kgl_parser/kgl_pubmed_xml_parser.cpp
        kgl_genomics/kgl_parser/kgl_pubmed_xml_parser.h
        kgl_genomics/kgl_parser/kgl_pubmed_cache.cpp
        kgl_genomics/kgl_parser/kgl_pubmed_cache.h
        kgl_genomics/kgl_parser/kgl_pubmed_api.cpp
        kgl_genomics/kgl_parser/kgl_pubmed_api.h
        kgl_genomics/kgl_literature/kgl_literature_filter.cpp
        kgl_genomics/kgl_literature/kgl_literature_analysis.cpp
        kgl_genomics/kgl_literature/kgl_literature_analysis.h
        kgl_genomics/kgl_literature/kgl_literature_filter.h
        kgl_genomics/kgl_literature/kgl_literature_filter_virtual.h
        kgl_genomics/kgl_literature/kgl_literature_data.h
        kgl_genomics/kgl_literature/kgl_literature_data.cpp
        kgl_genomics/kgl_genome_io/kgl_io_gff3.cpp
        kgl_genomics/kgl_genome_io/kgl_io_gff3.h
        kgl_genomics/kgl_genome_io/kgl_io_fasta.cpp
        kgl_genomics/kgl_genome_io/kgl_io_fasta.h
        kgl_genomics/kgl_parser/kgl_variant_vcf_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_vcf_impl.h
        kgl_genomics/kgl_parser/kgl_pf7_sample_parser.cpp
        kgl_genomics/kgl_parser/kgl_pf7_sample_parser.h
        kgl_genomics/kgl_parser/kgl_pf7_fws_parser.cpp
        kgl_genomics/kgl_parser/kgl_pf7_fws_parser.h
        kgl_genomics/kgl_parser/kgl_pf7_genetic_distance_parser.cpp
        kgl_genomics/kgl_parser/kgl_pf7_genetic_distance_parser.h
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_vep.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_vep.h
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_genome.h
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_genome.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_type.h
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_info.h
        kgl_genomics/kgl_parser/kgl_Pf7_physical_distance.cpp
        kgl_genomics/kgl_parser/kgl_Pf7_physical_distance.h
        kgl_genomics/kgl_variant_db/kgl_variant_db_population_filter.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_Pf7.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_Pf7.h
        kgl_genomics/kgl_variant_db/kgl_variant_db_variant.cpp
        kgl_genomics/kgl_variant_db/kgl_variant_db_variant.h
        kgl_genomics/kgl_mutation/kgl_mutation_db.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_db.h
        kgl_genomics/kgl_mutation/kgl_mutation_analysis.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_analysis.h
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_coding.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_coding.h
        kgl_genomics/kgl_mutation/kgl_mutation_variant_filter.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_variant_filter.h
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_contig.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_contig.h
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_offset.cpp
        kgl_genomics/kgl_variant_filter/kgl_variant_filter_db_offset.h
        kgl_genomics/kgl_mutation/kgl_mutation_interval.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_interval.h
        kgl_genomics/kgl_mutation/kgl_mutation_variant_map.h
        kgl_genomics/kgl_mutation/kgl_mutation_sequence.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_sequence.h
        kgl_genomics/kgl_mutation/kgl_mutation_transcript.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_transcript.h
        kgl_genomics/kgl_mutation/kgl_mutation_translate.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_translate.h
        kgl_genomics/kgl_mutation/kgl_mutation_sequence_update.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_coding.cpp
        kgl_genomics/kgl_mutation/kgl_mutation_coding.h
)

# 3rd party library for fast Levenshtein sequence matching.
set(EDLIB_SOURCE_FILES contrib/edlib/edlib.cpp)

# Gene ontology library
set(ONTOLOGY_SOURCE_FILES
        kol_ontology/contrib/kol_Accumulators.h
        kol_ontology/kol_ParserGoObo.h
        kol_ontology/kol_ParserGoXml.h
        kol_ontology/kol_PolicyEvidence.h
        kol_ontology/kol_PolicyRelationship.h
        kol_ontology/kol_SetSimilarityAllPairsAverage.h
        kol_ontology/kol_SetSimilarityAllPairsMax.h
        kol_ontology/kol_InformationAncestorMean.h
        kol_ontology/kol_TermAnnotation.h
        kol_ontology/kol_SetSimilarityBestMatchAverage.h
        kol_ontology/kol_InformationCoutoGraSMAdjusted.h
        kol_ontology/kol_InformationCoutoGraSM.h
        kol_ontology/kol_EnrichmentTools.h
        kol_ontology/kol_InformationExclusiveInherited.h
        kol_ontology/kol_InformationFrontier.h
        kol_ontology/kol_SetSimilarityGentlemanSimUI.h
        kol_ontology/kol_ParserAnnotationGaf.h
        kol_ontology/kol_library.h
        kol_ontology/kol_GoEnums.h
        kol_ontology/contrib/kol_GoGraphImpl.h
        kol_ontology/kol_ParserGoFactory.h
        kol_ontology/kol_ParserGoInterface.h
        kol_ontology/kol_SetSimilarityJaccard.h
        kol_ontology/kol_SetSimilarityMazanduSimDIC.h
        kol_ontology/kol_SetSimilarityMazanduSimUIC.h
        kol_ontology/kol_SimilarityJiangConrath.h
        kol_ontology/kol_SimilarityLin.h
        kol_ontology/kol_SimilarityResnik.h
        kol_ontology/kol_SimilarityPekarStaab.h
        kol_ontology/kol_SetSimilarityPesquitaSimGIC.h
        kol_ontology/kol_SimilarityMatrix.h
        kol_ontology/kol_ParserGoRapidXml.h
        kol_ontology/kol_SetUtilities.h
        kol_ontology/kol_InformationInterface.h
        kol_ontology/kol_InformationDepthMap.h
        kol_ontology/kol_InformationContentDAG.h
        kol_ontology/kol_SetSimilarityInterface.h
        kol_ontology/kol_SimilarityInterface.h
        kol_ontology/kol_SimilarityWriter.h
        kol_ontology/kol_SimilarityCache.h
        contrib/rapidxml/rapidxml.h
        contrib/rapidxml/rapidxml_iterators.h
        contrib/rapidxml/rapidxml_print.h
        contrib/rapidxml/rapidxml_utils.h
        kol_ontology/kgl_ontology/kgl_ontology_database.h
        kol_ontology/contrib/kol_GoGraphImpl.cpp
        kol_ontology/kol_InformationFrontier.cpp
        kol_ontology/kol_SimilarityCache.cpp
        kol_ontology/kol_SimilarityCacheAsymmetric.cpp
        kol_ontology/kol_SimilarityCacheAsymmetric.h
        kol_ontology/kol_SetSimilarityBestMatchAverage.cpp
        kol_ontology/kol_SetSimilarityAllPairsAverage.cpp
        kol_ontology/kol_SetSimilarityAllPairsMax.cpp
        kol_ontology/kol_SetSimilarityGentlemanSimUI.cpp
        kol_ontology/kol_SetSimilarityJaccard.cpp
        kol_ontology/kol_SetSimilarityMazanduSimDIC.cpp
        kol_ontology/kol_SetSimilarityMazanduSimUIC.cpp
        kol_ontology/kol_SetSimilarityPesquitaSimGIC.cpp
        kol_ontology/kol_InformationExclusiveInherited.cpp
        kol_ontology/kol_InformationCoutoGraSM.cpp
        kol_ontology/kol_InformationCoutoGraSMAdjusted.cpp
        kol_ontology/kol_ParserGoObo.cpp
        kol_ontology/kol_ParserGoXml.cpp
        kol_ontology/kol_InformationAncestorMean.cpp
        kol_ontology/kol_InformationDepthMap.cpp
        kol_ontology/kol_SimilarityPekarStaab.cpp
        kol_ontology/kol_SimilarityWriter.cpp
        kol_ontology/kol_ParserGoRapidXml.cpp
        kol_ontology/kol_SimilarityMatrix.cpp
        kol_ontology/kol_SimilarityImpl.cpp
        kol_ontology/kol_EnrichmentTools.cpp
        kol_ontology/kgl_ontology/kgl_ontology_database.cpp
        kol_ontology/kol_OntologyTypes.h
        kol_ontology/kgl_ontology/kgl_ontology_database.h
        kol_ontology/kgl_ontology/kgl_gene_cache_ontology.cpp
        kol_ontology/kgl_ontology/kgl_gene_cache_ontology.h
        kol_ontology/kgl_ontology/kgl_ontology_database_test.cpp
        kol_ontology/kgl_ontology/kgl_ontology_database_test.h
        kol_ontology/kol_SimilarityRelevance.h
        kol_ontology/kol_InformationContent.h
        kol_ontology/kol_InformationContentDAG.cpp
        kol_ontology/kol_InformationContentImpl.cpp
        kol_ontology/kol_InformationContentImpl.h
        kol_ontology/kol_InformationContent.cpp
        kol_ontology/kol_GoEnums.cpp
        kol_ontology/kol_ParserAnnotationGaf.cpp
        kol_ontology/kol_TermAnnotation.cpp
        kol_ontology/kol_ParserGafRecord.cpp
        kol_ontology/kol_ParserGafRecord.h
        kol_ontology/kol_SetSimilarityAverageBestMatch.cpp
        kol_ontology/kol_SetSimilarityAverageBestMatch.h
        kol_ontology/kol_OntologyFactory.cpp
        kol_ontology/kol_OntologyFactory.h
        kol_ontology/kol_GoGraph.cpp
        kol_ontology/kol_GoGraph.h)

# Gene ontology unit test suite
set(ONTOLOGY_UNIT_TEST_FILES
        kol_ontology/unit_test/kol_test.cpp
        kol_ontology/unit_test/kol_test_GoParsers.cpp
        kol_ontology/unit_test/kol_test_GoGraph.cpp
        kol_ontology/unit_test/kol_test.h
        kol_ontology/unit_test/kol_test_AnnotateParsers.cpp
        kol_ontology/unit_test/kol_test_Annotate.cpp
        kol_ontology/unit_test/kol_test_EnrichmentTools.cpp
        kol_ontology/unit_test/kol_test_JiangConrathSimilarity.cpp
        kol_ontology/unit_test/kol_test_LinSimilarity.cpp
        kol_ontology/unit_test/kol_test_similarity.h
        kol_ontology/unit_test/kol_test_ResnikSimilarity.cpp
        kol_ontology/unit_test/kol_test_PekarStaabSimilarity.cpp
        kol_ontology/unit_test/kol_test_RelevanceSimilarity.cpp
        kol_ontology/unit_test/kol_test_AllPairsMaxSet.cpp
        kol_ontology/unit_test/kol_test_AllPairAvgSet.cpp
        kol_ontology/unit_test/kol_test_setsimilarity.h
        kol_ontology/unit_test/kol_test_BestPairAvgSet.cpp
        kol_ontology/unit_test/kol_test_JaccardSet.cpp
        kol_ontology/unit_test/kol_test_SimDICSet.cpp
        kol_ontology/unit_test/kol_test_SimGICSet.cpp
        kol_ontology/unit_test/kol_test_SimUISet.cpp
        kol_ontology/unit_test/kol_test_SimUICSet.cpp
        kol_ontology/unit_test/kol_test_SharedInfoAdjCouto.cpp
        kol_ontology/unit_test/kol_test_sharedinfo.h
        kol_ontology/unit_test/kol_test_SharedAncestorMean.cpp
        kol_ontology/unit_test/kol_test_SharedCouto.cpp
        kol_ontology/unit_test/kol_test_SharedExInherit.cpp
        kol_ontology/unit_test/kol_test_SharedFrontier.cpp
        kol_ontology/unit_test/kol_test_DepthMap.cpp
        kol_ontology/unit_test/kol_test_TermInfoMap.cpp
        kol_ontology/unit_test/kol_test_TermSimMatrix_JC_BP.cpp
        kol_ontology/unit_test/kol_test_TermSimMatrix_LIN_MF.cpp
        kol_ontology/unit_test/kol_test_TermSimMatrix_RES_CC.cpp
        kol_ontology/unit_test/kol_test_TermSimCache_LIN_BP.cpp
        kol_ontology/unit_test/kol_test_data.h
        kol_ontology/unit_test/kol_test_SymbolicSet.cpp
        kol_ontology/unit_test/kol_test_symbolicset.h
        kol_ontology/unit_test/kol_test_OntologyDatabase.cpp)

# Genetic analysis library
set(ANALYTIC_SOURCE_FILES
        kga_analytic/kga_PfEMP/kga_phylogenetic_analysis.h
        kga_analytic/kga_PfEMP/kga_phylogenetic_analysis.cpp
        kga_analytic/kga_PfEMP/kga_analysis_gene_sequence.h
        kga_analytic/kga_PfEMP/kga_analysis_gene_sequence.cpp
        kga_analytic/kga_PfEMP/kga_upgma_node.h
        kga_analytic/kga_PfEMP/kga_upgma_node.cpp
        kga_analytic/kga_legacy/kga_rna_search.cpp
        kga_analytic/kga_legacy/kga_rna_search.h
        kga_analytic/kga_legacy/kga_ploidy_analysis.h
        kga_analytic/kga_PfEMP/kga_upgma_unphased.h
        kga_analytic/kga_PfEMP/kga_upgma_unphased.cpp
        kga_analytic/kga_legacy/kga_sequence_complexity.h
        kga_analytic/kga_PfEMP/kga_upgma.h
        kga_analytic/kga_analysis_all.h
        kga_analytic/kga_template_analysis/kga_analysis_null.cpp
        kga_analytic/kga_template_analysis/kga_analysis_null.h
        kga_analytic/kga_info/kga_analysis_interval.cpp
        kga_analytic/kga_info/kga_analysis_interval.h
        kga_analytic/kga_info/kga_analysis_info_filter.cpp
        kga_analytic/kga_info/kga_analysis_info_filter.h
        kga_analytic/kga_info/kga_analysis_age.cpp
        kga_analytic/kga_info/kga_analysis_age.h
        kga_analytic/kga_legacy/kga_legacy_Pf_analysis.cpp
        kga_analytic/kga_legacy/kga_legacy_Pf_analysis.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_diploid.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_synthetic.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_calc.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_diploid.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_locus.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_locus.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_calc.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_freq.cpp
        kga_analytic/kga_sequence_analysis/kga_analysis_sequence.cpp
        kga_analytic/kga_sequence_analysis/kga_analysis_sequence.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_freq.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_synthetic.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_output.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_output.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_syngen.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_syngen.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_execute.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_execute.h
        kga_analytic/kga_analysis_virtual.h
        kga_analytic/kga_inbreed/kga_analysis_inbreed_args.cpp
        kga_analytic/kga_inbreed/kga_analysis_inbreed_args.h
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP.cpp
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP.h
        kga_analytic/kga_mutation/kga_analysis_mutation.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_output.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_ethnic.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_ethnic.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_clinvar.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_clinvar.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_variant.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_variant.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_stats.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_stats.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_ontology.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_ontology.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_allele.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_allele.h
        kga_analytic/kga_mutation/kga_analysis_mutation_data.h
        kga_analytic/kga_mutation/kga_analysis_mutation_data.cpp
        kga_analytic/kga_info/kga_analysis_json.cpp
        kga_analytic/kga_info/kga_analysis_json.h
        kga_analytic/kga_literature/kga_analysis_literature.cpp
        kga_analytic/kga_literature/kga_analysis_literature.h
        kga_analytic/kga_literature/kga_analysis_literature_gene.cpp
        kga_analytic/kga_literature/kga_analysis_literature_gene.h
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_allele_pop.cpp
        kga_analytic/kga_mutation/kga_analysis_mutation_gene_allele_pop.h
        kga_analytic/kga_literature/kga_analysis_literature_publication.cpp
        kga_analytic/kga_literature/kga_analysis_literature_publication.h
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_distance.cpp
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_perform.cpp
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_variant.cpp
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_variant.h
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_heterozygous.cpp
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_heterozygous.h
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_filter.cpp
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_FWS.cpp
        kga_analytic/kga_PfEMP/kga_analysis_PfEMP_FWS.h
)

# Deconvolution and mixture files.
set(MIXTURE_SOURCE_FILES
        kgl_mixture/kgl_variant_classify.h
        kgl_mixture/kgl_variant_classify.cpp
        kgl_mixture/kgl_variant_classify_io.cpp
        kgl_mixture/kgl_read_phasing.h
        kgl_mixture/kgl_read_phasing.cpp)

# Application level objects
set(APPLICATION_SOURCE_FILES
        kgl_app/kgl_main.cpp
        kgl_app/kgl_gene_app.cpp
        kgl_app/kgl_gene_app.h
        kgl_app/kgl_gene_parse.cpp
        kgl_app/kgl_properties.h
        kgl_app/kgl_properties.cpp
        kgl_app/kgl_package.cpp
        kgl_app/kgl_package.h
        kgl_app/kgl_package_analysis.cpp
        kgl_app/kgl_package_analysis.h
        kgl_app/kgl_runtime.cpp
        kgl_app/kgl_runtime.h
        kgl_app/kgl_runtime_config.cpp
        kgl_app/kgl_runtime_config.h
        kgl_app/kgl_package_resource.cpp
        kgl_app/kgl_runtime_resource.h
        kgl_app/kgl_runtime_resource.cpp
        kgl_app/kgl_properties_resource.cpp
        kgl_app/kgl_properties_resource.h
        kgl_app/kgl_properties_resource_items.cpp
        kgl_app/kgl_properties_resource_pf.cpp
        kgl_app/kgl_package_resource_pf.cpp)



add_library(kel_utility STATIC ${UTILITY_SOURCE_FILES} ${MATHEMATICAL_SOURCE_FILES} ${EDLIB_SOURCE_FILES} )

add_library(kel_app STATIC ${RUNTIME_SOURCE_FILES})

add_library(kel_thread STATIC ${THREAD_SOURCE_FILES})

add_library(kel_io STATIC ${IO_SOURCE_FILES})

add_library(kgl_genomics STATIC ${GENETIC_LIBRARY_SOURCE_FILES})

target_link_libraries(kgl_genomics kel_utility kel_thread)

add_library(kgl_analysis STATIC ${ANALYTIC_SOURCE_FILES})

target_link_libraries(kgl_analysis kgl_genomics kel_utility kel_thread)

#generate kgl_genome executable
add_executable (kgl_genome ${APPLICATION_SOURCE_FILES} )

# Specify the static libraries.
target_link_libraries(kgl_genome kgl_analysis kgl_genomics kol_ontology kel_app kel_utility kel_thread kel_io ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} nlopt bz2 curl)

add_library(kol_ontology STATIC ${ONTOLOGY_SOURCE_FILES})

#generate libraries.
add_executable(kol_test ${ONTOLOGY_UNIT_TEST_FILES})

target_link_libraries (kol_test kol_ontology kel_app kel_utility kel_thread kel_io ${Boost_LIBRARIES})


############################################################################################################
# Generate the kpl executable.



set(KPL_LIBRARY_SOURCE_FILES
        kpl_phylogenetic/kpl_node.h
        kpl_phylogenetic/kpl_tree.h
        kpl_phylogenetic/kpl_treemanip.h
        kpl_phylogenetic/kpl_xstrom.h
        kpl_phylogenetic/kpl_splittree.h
        kpl_phylogenetic/kpl_treesummary.h
        kpl_phylogenetic/kpl_strom.h
        kpl_phylogenetic/kpl_geneticcode.h
        kpl_phylogenetic/kpl_genetictype.h
        kpl_phylogenetic/kpl_partition.h
        kpl_phylogenetic/kpl_geneticdata.h
        kpl_phylogenetic/kpl_geneticcode.cpp
        kpl_phylogenetic/kpl_geneticdata.cpp
        kpl_phylogenetic/kpl_genetictype.cpp
        kpl_phylogenetic/kpl_node.cpp
        kpl_phylogenetic/kpl_partition.cpp
        kpl_phylogenetic/kpl_splittree.cpp
        kpl_phylogenetic/kpl_strom.cpp
        kpl_phylogenetic/kpl_tree.cpp
        kpl_phylogenetic/kpl_treemanip.cpp
        kpl_phylogenetic/kpl_treesummary.cpp
        kpl_phylogenetic/kpl_xstrom.cpp
        kpl_phylogenetic/kpl_likelihood.cpp
        kpl_phylogenetic/kpl_likelihood.h
        kpl_phylogenetic/kpl_model.cpp
        kpl_phylogenetic/kpl_model.h
        kpl_phylogenetic/kpl_qmatrix.cpp
        kpl_phylogenetic/kpl_qmatrix.h
        kpl_phylogenetic/kpl_qmatrixnucleotide.cpp
        kpl_phylogenetic/kpl_qmatrixnucleotide.h
        kpl_phylogenetic/kpl_qmatrixcodon.cpp
        kpl_phylogenetic/kpl_qmatrixcodon.h
        kpl_phylogenetic/kpl_asrv.cpp
        kpl_phylogenetic/kpl_asrv.h
        kpl_phylogenetic/kpl_random.cpp
        kpl_phylogenetic/kpl_random.h
        kpl_phylogenetic/kpl_mcmc_updater.cpp
        kpl_phylogenetic/kpl_mcmc_updater.h
        kpl_phylogenetic/kpl_mcmc_gamma.cpp
        kpl_phylogenetic/kpl_mcmc_gamma.h
        kpl_phylogenetic/kpl_mcmc_chain.cpp
        kpl_phylogenetic/kpl_mcmc_chain.h
        kpl_phylogenetic/kpl_mcmc_output.cpp
        kpl_phylogenetic/kpl_mcmc_output.h
        kpl_phylogenetic/kpl_mcmc_dirichlet.cpp
        kpl_phylogenetic/kpl_mcmc_dirichlet.h
        kpl_phylogenetic/kpl_mcmc_statefrequpdater.cpp
        kpl_phylogenetic/kpl_mcmc_statefrequpdater.h
        kpl_phylogenetic/kpl_mcmc_exchangeupdater.cpp
        kpl_phylogenetic/kpl_mcmc_exchangeupdater.h
        kpl_phylogenetic/kpl_mcmc_subsetupdater.cpp
        kpl_phylogenetic/kpl_mcmc_subsetupdater.h
        kpl_phylogenetic/kpl_mcmc_pinvarupdater.cpp
        kpl_phylogenetic/kpl_mcmc_pinvarupdater.h
        kpl_phylogenetic/kpl_mcmc_omegaupdater.cpp
        kpl_phylogenetic/kpl_mcmc_omegaupdater.h
        kpl_phylogenetic/kpl_mcmc_treeupdater.cpp
        kpl_phylogenetic/kpl_mcmc_treeupdater.h
        kpl_phylogenetic/kpl_mcmc_treelength.cpp
        kpl_phylogenetic/kpl_mcmc_treelength.h
        kpl_phylogenetic/kpl_mcmc_polytomyupdater.cpp
        kpl_phylogenetic/kpl_mcmc_polytomyupdater.h
        kpl_phylogenetic/kpl_mcmc_polytomyprior.cpp
        kpl_phylogenetic/kpl_mcmc_polytomyprior.h
        kpl_phylogenetic/kpl_tree_io.cpp
        kpl_phylogenetic/kpl_tree_io.h
        kpl_phylogenetic/kpl_main_app.cpp
        kpl_phylogenetic/kpl_main_app.h)

add_library(kpl_library_lib STATIC ${KPL_LIBRARY_SOURCE_FILES})


#generate kpl executable
add_executable (kpl kpl_phylogenetic/kpl_main.cpp)

#specify link libraries.
target_link_libraries(kpl kpl_library_lib kgl_genomics kel_app kel_utility kel_thread libncl.a ${Boost_LIBRARIES} hmsbeagle)

## DOXYGEN
#set(BUILD_KGL_KGD_DOC ON CACHE BOOL "Use Doxygen to create the HTML/Latex based API documentation")
#option(BUILD_DOC "Use Doxygen to create the HTML/Latex based API documentation" ON)
if(BUILD_KGL_KGD_DOC)

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
        #                    FASTA_COMMENT_ "Generating API documentation with Doxygen"
        #                    VERBATIM )

    endif(DOXYGEN_FOUND)

endif(BUILD_KGL_KGD_DOC)

