cmake_minimum_required(VERSION 3.17)
project(KGL)

set(CMAKE_CXX_STANDARD 20)


# add extra include directories
include_directories( /usr/local/include/seqan3/submodules/range-v3/include   # the range library required by seqan3.
                    /usr/local/include/seqan3/submodules/sdsl-lite/include   # the sdsl library required by seqan3.
                    /usr/local/include/seqan3/submodules/cereal/include   # the cereal library required by seqan3.
                    /usr/local/include/seqan3/include
                    /usr/local/include/spdlog/include/
                    contrib/edlib # 3rd party library for fast Levenshtein seqence matching.
                    /usr/local/include/libhmsbeagle-1  # the beagle mcmc library used by kpl
                    /usr/include/eigen3
                    kel_utility         # Basic Infrastructure library
                    kgl_genomics                  # Basic Genetic library
                    kgl_genomics/kgl_parser
                    kgl_genomics/kgl_evidence
                    kgl_genomics/kgl_sequence
                    kgl_genomics/kgl_database
                    kgl_mixture                  # Deconvolution and mixture library
                    kgl_analytic                  # Genetic analysis library (top level)
                    kol_ontology                  # Experimental Ontology Library
                    kol_ontology/kgl_ontology     # Higher level ontology functionality
                    kgl_analytic/kgl_inbreed
                    kgl_analytic/kgl_verify
                    kgl_analytic/kgl_info
                    kgl_analytic/kgl_Pf
                    kgl_analytic/kgl_legacy
                    kgl_app)



# add boost library directories
set(BOOST_INCLUDEDIR /usr/local/include)
set(BOOST_LIBRARYDIR /usr/local/lib/boost)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
include_directories(${Boost_INCLUDE_DIRS})
find_package( Boost COMPONENTS system filesystem timer chrono program_options iostreams REQUIRED )

# zlib for compressed VCF files.
find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

# Debug g++
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++20 -fPIC -fconcepts -ggdb -DBOOST_MATH_DISABLE_FLOAT128")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-ggdb -DBOOST_MATH_DISABLE_FLOAT128")

# Release g++
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++20 -fPIC -fconcepts -O3 -s -DNDEBUG -DBOOST_MATH_DISABLE_FLOAT128")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE " -O3 -s -DNDEBUG -DBOOST_MATH_DISABLE_FLOAT128")

# Profile g++
set(CMAKE_CXX_FLAGS_PROFILE "-std=c++20 -fPIC -fconcepts -pg -g -ggdb -O3 -DBOOST_MATH_DISABLE_FLOAT128 -Xlinker -Map=output.map ")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g -ggdb -O3 -DBOOST_MATH_DISABLE_FLOAT128 -Xlinker -Map=output.map ")

# All possible warnings
#set(CMAKE_ERROR_WARNING_FLAGS "-std=c++20 -Wall -pthread -Wextra -Werror -Wpedantic -pedantic-errors -Wstrict-aliasing -Wno-variadic-macros")
# Relax slightly to temporarily allow unused parameters.
set(CMAKE_ERROR_WARNING_FLAGS "-std=c++20 -Wall -pthread -Wextra -Wpedantic -Wstrict-aliasing -Wno-variadic-macros")

# Actually set the g++ flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)

# Basic infrastructure library
set(LIBRARY_SOURCE_FILES
        kel_utility/kel_logging.h
        kel_utility/kel_logging.cpp
        kel_utility/kel_utility.h
        kel_utility/kel_utility.cpp
        kel_utility/kel_exec_env.h
        kel_utility/kel_exec_env.cpp
        kel_utility/kel_property_tree.cpp
        kel_utility/kel_property_tree.h
        kel_utility/kel_distribution.h
        kel_utility/kel_distribution.cpp
        kel_utility/kel_basic_io.cpp
        kel_utility/kel_basic_io.h
        kel_utility/kel_percentile.h
        kel_utility/kel_thread_pool.h
        kel_utility/kel_optimize.cpp
        kel_utility/kel_optimize.h
        kel_utility/kel_patterns.h
        kel_utility/kel_lock.h
        kel_utility/kel_mt_queue.h
        kel_utility/kel_exec_env_app.h
        kel_utility/kel_bzip.cpp
        kel_utility/kel_bzip.h
        kel_utility/kel_bound_queue.h
        kel_utility/kel_mem_alloc.h
        kel_utility/kel_mem_alloc.cpp)

# Basic Genetic infrastructure library
set(GENETIC_LIBRARY_SOURCE_FILES
        kgl_genomics/kgl_database/kgl_genome_types.h
        kgl_genomics/kgl_database/kgl_gff_fasta.h
        kgl_genomics/kgl_database/kgl_gff_fasta.cpp
        kgl_app/kgl_resource_db.h
        kgl_genomics/kgl_database/kgl_genome_feature.h
        kgl_genomics/kgl_database/kgl_genome_feature.cpp
        kgl_genomics/kgl_database/kgl_variant.h
        kgl_genomics/kgl_database/kgl_variant.cpp
        kgl_genomics/kgl_database/kgl_variant_filter.h
        kgl_genomics/kgl_database/kgl_variant_filter.cpp
        kgl_genomics/kgl_sequence/kgl_table.h
        kgl_genomics/kgl_sequence/kgl_table.cpp
        kgl_genomics/kgl_sequence/kgl_table_ncbi.h
        kgl_genomics/kgl_sequence/kgl_sequence_base.h
        kgl_genomics/kgl_sequence/kgl_sequence_base.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_amino.h
        kgl_genomics/kgl_sequence/kgl_sequence_amino.cpp
        kgl_genomics/kgl_sequence/kgl_alphabet_amino.h
        kgl_genomics/kgl_database/kgl_genome_prelim.h
        kgl_genomics/kgl_database/kgl_genome_prelim.cpp
        kgl_genomics/kgl_database/kgl_genome_attributes.h
        kgl_genomics/kgl_database/kgl_genome_attributes.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_virtual.h
        kgl_genomics/kgl_database/kgl_genome_contig.cc
        kgl_genomics/kgl_database/kgl_genome_verify.cpp
        kgl_genomics/kgl_database/kgl_genome_genome.cpp
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
        kgl_genomics/kgl_database/kgl_variant_mutation.cpp
        kgl_genomics/kgl_database/kgl_variant_mutation.h
        kgl_genomics/kgl_database/kgl_variant_db_mutation.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_offset.h
        kgl_genomics/kgl_sequence/kgl_sequence_offset.cpp
        kgl_genomics/kgl_database/kgl_variant_mutation_offset.h
        kgl_genomics/kgl_database/kgl_variant_mutation_offset.cpp
        kgl_genomics/kgl_sequence/kgl_table_impl.h
        kgl_genomics/kgl_sequence/kgl_table_impl.cpp
        kgl_genomics/kgl_sequence/kgl_table_organism.h
        kgl_genomics/kgl_parser/kgl_gaf_parser.h
        kgl_genomics/kgl_parser/kgl_gaf_parser.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_compare.h
        kgl_genomics/kgl_sequence/kgl_sequence_compare_impl.h
        kgl_genomics/kgl_sequence/kgl_sequence_compare_impl.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_distance.h
        kgl_genomics/kgl_database/kgl_genome_types.cpp
        kgl_genomics/kgl_sequence/kgl_sequence_distance_impl.h
        kgl_genomics/kgl_sequence/kgl_sequence_distance_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_pf_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_pf_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_readvcf_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_cigar.h
        kgl_genomics/kgl_parser/kgl_variant_factory_vcf_parse_cigar.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_record_vcf_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_record_vcf_impl.cpp
        kgl_genomics/kgl_parser/kgl_genome_aux_csv.h
        kgl_genomics/kgl_parser/kgl_genome_aux_csv_impl.cpp
        kgl_genomics/kgl_database/kgl_variant_db_genome.h
        kgl_genomics/kgl_database/kgl_variant_db_population.h
        kgl_genomics/kgl_database/kgl_genome_contig_feature.h
        kgl_genomics/kgl_database/kgl_genome_contig_feature.cpp
        kgl_genomics/kgl_database/kgl_genome_contig_aux.cc
        kgl_genomics/kgl_database/kgl_genome_contig_aux.h
        kgl_genomics/kgl_sequence/kgl_phylogenetic_tree.cpp
        kgl_genomics/kgl_sequence/kgl_phylogenetic_tree.h
        kgl_genomics/kgl_parser/kgl_variant_file_vcf_record.h
        kgl_genomics/kgl_parser/kgl_variant_factory_readvcf_impl.cpp
        kgl_genomics/kgl_parser/kgl_data_file_impl.cpp
        kgl_genomics/kgl_parser/kgl_data_file_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_grch_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_grch_impl.h
        kgl_genomics/kgl_parser/kgl_variant_file_vcf_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_file_vcf_impl.h
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
        kgl_genomics/kgl_database/kgl_variant_db_offset.h
        kgl_genomics/kgl_database/kgl_variant_db_contig.h
        kgl_genomics/kgl_database/kgl_variant_db_mutation.h
        kgl_genomics/kgl_parser/kgl_ped_parser.h
        kgl_genomics/kgl_parser/kgl_ped_parser.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_parsers.h
        kgl_genomics/kgl_parser/kgl_variant_factory_gnomad_impl.cpp
        kgl_genomics/kgl_parser/kgl_variant_factory_gnomad_impl.h
        kgl_genomics/kgl_parser/kgl_variant_factory_parsers.cpp
        kgl_genomics/kgl_database/kgl_variant_db_contig.cpp
        kgl_genomics/kgl_database/kgl_variant_db_genome.cpp
        kgl_genomics/kgl_database/kgl_variant_db_population.cpp
        kgl_genomics/kgl_database/kgl_variant_db_type.h
        kgl_genomics/kgl_database/kgl_genome_contig.h
        kgl_genomics/kgl_database/kgl_variant_db_type.cpp
        kgl_genomics/kgl_parser/kgl_Pf3k_COI.cpp
        kgl_genomics/kgl_parser/kgl_Pf3k_COI.h
        kgl_genomics/kgl_database/kgl_genome_genome.h
        kgl_genomics/kgl_parser/kgl_square_parser.cpp
        kgl_genomics/kgl_parser/kgl_square_parser.h
        kgl_genomics/kgl_database/kgl_variant_db_offset.cpp
        kgl_genomics/kgl_database/kgl_variant_filter_virtual.h
        kgl_genomics/kgl_database/kgl_variant_db_offset_diploid.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_mem_alloc.cpp
        kgl_genomics/kgl_evidence/kgl_variant_factory_vcf_evidence_mem_alloc.h
        kgl_genomics/kgl_parser/kgl_ensembl_id_parser.cpp
        kgl_genomics/kgl_parser/kgl_ensembl_id_parser.h)

# 3rd party library for fast Levenshtein sequence matching.
set(EDLIB_SOURCE_FILES contrib/edlib/edlib.cpp )

# Gene ontology library
set(ONTOLOGY_SOURCE_FILES
        kol_ontology/kol_Accumulators.h
        kol_ontology/kol_AllowedRelationshipOboGoParser.h
        kol_ontology/kol_AllowedRelationshipXmlGoParser.h
        kol_ontology/kol_AllowedSetEvidencePolicy.h
        kol_ontology/kol_AllowedSetRelationshipPolicy.h
        kol_ontology/kol_AllPairsAverageSetSimilarity.h
        kol_ontology/kol_AllPairsMaxSetSimilarity.h
        kol_ontology/kol_AncestorMeanSharedInformation.h
        kol_ontology/kol_AnnotationData.h
        kol_ontology/kol_AnnotationParserFactory.h
        kol_ontology/kol_AnnotationParserInterface.h
        kol_ontology/kol_BestMatchAverageSetSimilarity.h
        kol_ontology/kol_CoutoGraSMAdjustedSharedInformation.h
        kol_ontology/kol_CoutoGraSMSharedInformation.h
        kol_ontology/kol_DisallowedSetEvidencePolicy.h
        kol_ontology/kol_EnrichmentTools.h
        kol_ontology/kol_EntrezGene2GoAnnotationParser.h
        kol_ontology/kol_EvidencePolicyInterface.h
        kol_ontology/kol_ExclusivelyInheritedSharedInformation.h
        kol_ontology/kol_ExperimentalEvidencePolicy.h
        kol_ontology/kol_FrontierSharedInformation.h
        kol_ontology/kol_GafAnnotationParser.h
        kol_ontology/kol_GentlemanSimUISetSimilarity.h
        kol_ontology/kol_GoaAnnotationParser.h
        kol_ontology/kol_library.h
        kol_ontology/kol_GoEnums.h
        kol_ontology/kol_GoGraph.h
        kol_ontology/kol_GoParserFactory.h
        kol_ontology/kol_GoParserInterface.h
        kol_ontology/kol_JaccardSetSimilarity.h
        kol_ontology/kol_MazanduSimDICSetSimilarity.h
        kol_ontology/kol_MazanduSimUICSetSimilarity.h
        kol_ontology/kol_MgiAnnotationParser.h
        kol_ontology/kol_MICASharedInformation.h
        kol_ontology/kol_SimilarityJiangConrath.h
        kol_ontology/kol_SimilarityLin.h
        kol_ontology/kol_SimilarityResnik.h
        kol_ontology/kol_SimilarityPekarStaab.h
        kol_ontology/kol_PesquitaSimGICSetSimilarity.h
        kol_ontology/kol_PrecomputedMatrixTermSimilarity.h
        kol_ontology/kol_RapidXmlGoParser.h
        kol_ontology/kol_RelationshipPolicyInterface.h
        kol_ontology/kol_SetUtilities.h
        kol_ontology/kol_SharedInformationInterface.h
        kol_ontology/kol_StandardOboGoParser.h
        kol_ontology/kol_StandardRelationshipPolicy.h
        kol_ontology/kol_StandardXmlGoParser.h
        kol_ontology/kol_TermDepthMap.h
        kol_ontology/kol_TermInformationContentMap.h
        kol_ontology/kol_TermSetSimilarityInterface.h
        kol_ontology/kol_TermSimilarityInterface.h
        kol_ontology/kol_TermSimilarityWriter.h
        kol_ontology/kol_TermSimilarityCache.h
        kol_ontology/xml/rapidxml.h
        kol_ontology/xml/rapidxml_iterators.h
        kol_ontology/xml/rapidxml_print.h
        kol_ontology/xml/rapidxml_utils.h
        kol_ontology/kgl_ontology/kgl_ontology_database.h
        kol_ontology/kol_AnnotationData.cpp
        kol_ontology/kol_GoGraph.cpp
        kol_ontology/kol_FrontierSharedInformation.cpp
        kol_ontology/kol_TermSimilarityCache.cpp
        kol_ontology/kol_AsymmetricSimilarityCache.cpp
        kol_ontology/kol_AsymmetricSimilarityCache.h
        kol_ontology/kol_BestMatchAverageSetSimilarity.cpp
        kol_ontology/kol_AllPairsAverageSetSimilarity.cpp
        kol_ontology/kol_AllPairsMaxSetSimilarity.cpp
        kol_ontology/kol_GentlemanSimUISetSimilarity.cpp
        kol_ontology/kol_JaccardSetSimilarity.cpp
        kol_ontology/kol_MazanduSimDICSetSimilarity.cpp
        kol_ontology/kol_MazanduSimUICSetSimilarity.cpp
        kol_ontology/kol_PesquitaSimGICSetSimilarity.cpp
        kol_ontology/kol_ExclusivelyInheritedSharedInformation.cpp
        kol_ontology/kol_CoutoGraSMSharedInformation.cpp
        kol_ontology/kol_CoutoGraSMAdjustedSharedInformation.cpp
        kol_ontology/kol_AllowedRelationshipOboGoParser.cpp
        kol_ontology/kol_AllowedRelationshipXmlGoParser.cpp
        kol_ontology/kol_AncestorMeanSharedInformation.cpp
        kol_ontology/kol_MICASharedInformation.cpp
        kol_ontology/kol_TermDepthMap.cpp
        kol_ontology/kol_SimilarityPekarStaab.cpp
        kol_ontology/kol_TermSimilarityWriter.cpp
        kol_ontology/kol_RapidXmlGoParser.cpp
        kol_ontology/kol_PrecomputedMatrixTermSimilarity.cpp
        kol_ontology/kol_SimilarityImpl.cpp
        kol_ontology/kol_EnrichmentTools.cpp
        kol_ontology/kol_EntrezGene2GoAnnotationParser.cpp
        kol_ontology/kgl_ontology/kgl_ontology_database.cpp
        kol_ontology/kol_OntologyTypes.h
        kol_ontology/kgl_ontology/kgl_ontology_database.h
        kol_ontology/kol_AnnotationData.cpp
        kol_ontology/kgl_ontology/kgl_gene_cache_ontology.cpp
        kol_ontology/kgl_ontology/kgl_gene_cache_ontology.h
        kol_ontology/kgl_ontology/kgl_ontology_database_test.cpp kol_ontology/kgl_ontology/kgl_ontology_database_test.h kol_ontology/kol_TermInformationInterface.h kol_ontology/kol_SimilarityRelevance.h kol_ontology/kol_TermInformationContentUnique.h kol_ontology/kol_TermInformationContentMap.cpp kol_ontology/kol_TermInformationContentImpl.cpp kol_ontology/kol_TermInformationContentImpl.h kol_ontology/kol_TermInformationContentUnique.cpp)

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
        kol_ontology/unit_test/kol_test_SharedMICA.cpp
        kol_ontology/unit_test/kol_test_DepthMap.cpp
        kol_ontology/unit_test/kol_test_TermInfoMap.cpp
        kol_ontology/unit_test/kol_test_TermSimMatrix_JC_BP.cpp
        kol_ontology/unit_test/kol_test_TermSimMatrix_LIN_MF.cpp
        kol_ontology/unit_test/kol_test_TermSimMatrix_RES_CC.cpp
        kol_ontology/unit_test/kol_test_TermSimCache_LIN_BP.cpp
        kol_ontology/unit_test/kol_test_data.h)

# Genetic analysis library
set(ANALYTIC_SOURCE_FILES
        kgl_analytic/kgl_Pf/kgl_phylogenetic_analysis.h
        kgl_analytic/kgl_Pf/kgl_phylogenetic_analysis.cpp
        kgl_analytic/kgl_Pf/kgl_analysis_gene_sequence.h
        kgl_analytic/kgl_Pf/kgl_analysis_gene_sequence.cpp
        kgl_analytic/kgl_Pf/kgl_upgma_node.h
        kgl_analytic/kgl_Pf/kgl_upgma_node.cpp
        kgl_analytic/kgl_legacy/kgl_rna_search.cpp
        kgl_analytic/kgl_legacy/kgl_rna_search.h
        kgl_analytic/kgl_legacy/kgl_ploidy_analysis.h
        kgl_analytic/kgl_Pf/kgl_upgma_unphased.h
        kgl_analytic/kgl_Pf/kgl_upgma_unphased.cpp
        kgl_analytic/kgl_legacy/kgl_sequence_complexity.h
        kgl_analytic/kgl_legacy/kgl_epigenetic_motif.cpp
        kgl_analytic/kgl_legacy/kgl_epigenetic_motif.cpp
        kgl_analytic/kgl_Pf/kgl_upgma.h
        kgl_analytic/kgl_analysis_all.h
        kgl_analytic/kgl_verify/kgl_analysis_null.cpp
        kgl_analytic/kgl_verify/kgl_analysis_null.h
        kgl_analytic/kgl_info/kgl_analysis_interval.cpp
        kgl_analytic/kgl_info/kgl_analysis_interval.h
        kgl_analytic/kgl_info/kgl_analysis_info_filter.cpp
        kgl_analytic/kgl_info/kgl_analysis_info_filter.h
        kgl_analytic/kgl_info/kgl_analysis_age.cpp
        kgl_analytic/kgl_info/kgl_analysis_age.h
        kgl_analytic/kgl_legacy/kgl_legacy_Pf_analysis.cpp
        kgl_analytic/kgl_legacy/kgl_legacy_Pf_analysis.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_diploid.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_synthetic.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_calc.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_diploid.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_locus.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_locus.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_calc.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_freq.cpp
        kgl_analytic/kgl_verify/kgl_analysis_verify.cpp
        kgl_analytic/kgl_verify/kgl_analysis_verify.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_freq.h
        kgl_genomics/kgl_database/kgl_variant_db_freq.cpp
        kgl_genomics/kgl_database/kgl_variant_db_freq.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_synthetic.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_output.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_output.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_syngen.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_syngen.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_execute.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_execute.h
        kgl_analytic/kgl_analysis_virtual.h
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_args.cpp
        kgl_analytic/kgl_inbreed/kgl_analysis_inbreed_args.h
        kgl_analytic/kgl_Pf/kgl_analysis_PfEMP.cpp
        kgl_analytic/kgl_Pf/kgl_analysis_PfEMP.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation.h
        kgl_analytic/kgl_Pf/kgl_phase_Pf.cpp
        kgl_analytic/kgl_Pf/kgl_phase_Pf.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_output.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_clinvar.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_clinvar.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_ethnic.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_ethnic.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_clinvar.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_clinvar.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_variant.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_variant.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_stats.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_stats.h
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_ontology.cpp
        kgl_analytic/kgl_mutation/kgl_analysis_mutation_gene_ontology.h)

# Deconvolution and mixture files.

set(MIXTURE_SOURCE_FILES
        kgl_mixture/kgl_variant_classify.h
        kgl_mixture/kgl_variant_classify.cpp
        kgl_mixture/kgl_variant_classify_io.cpp
        kgl_mixture/kgl_read_phasing.h
        kgl_mixture/kgl_read_phasing.cpp)


# Application level objects
set(APPLICATION_SOURCE_FILES
        kgl_app/kgl_main.cc
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
        kgl_app/kgl_package_resources.cpp)





add_library(kel_utility STATIC ${LIBRARY_SOURCE_FILES} ${EDLIB_SOURCE_FILES})

add_library(kgl_genomics STATIC ${GENETIC_LIBRARY_SOURCE_FILES})

target_link_libraries(kgl_genomics kel_utility )

add_library(kgl_analysis STATIC ${ANALYTIC_SOURCE_FILES})

target_link_libraries(kgl_analysis kgl_genomics kel_utility )

#generate kgl_genome executable
add_executable (kgl_genome ${APPLICATION_SOURCE_FILES} )

# Specify the static libraries.
target_link_libraries(kgl_genome kgl_analysis kgl_genomics kol_ontology kel_utility ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} nlopt)

add_library(kol_ontology STATIC ${ONTOLOGY_SOURCE_FILES})

#generate libraries.
add_executable(kol_test ${ONTOLOGY_UNIT_TEST_FILES})

target_link_libraries (kol_test kol_ontology kel_utility ${Boost_LIBRARIES})



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
        kpl_phylogenetic/kpl_tree_io.h)

add_library(kpl_library_lib STATIC ${KPL_LIBRARY_SOURCE_FILES})

link_directories(/usr/local/lib)

link_directories(/usr/local/lib/ncl)

link_directories(/usr/local/lib/boost)

link_directories(/usr/local/lib/libhmsbeagle-1/lib)

#generate kpl executable
add_executable (kpl kpl_phylogenetic/kpl_main.cpp)

#specify link libraries.
target_link_libraries(kpl kpl_library_lib kgl_genomics kel_utility libncl.a ${Boost_LIBRARIES} hmsbeagle)

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
        #                    COMMENT_ "Generating API documentation with Doxygen"
        #                    VERBATIM )

    endif(DOXYGEN_FOUND)

endif(BUILD_KGL_KGD_DOC)

