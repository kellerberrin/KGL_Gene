cmake_minimum_required(VERSION 3.5)
project(KGL)

set(CMAKE_CXX_STANDARD 17)


# add extra include directories
include_directories( /usr/local/include/seqan3/submodules/range-v3/include   # the range library required by seqan3.
                    /usr/local/include/seqan3/submodules/sdsl-lite/include   # the sdsl library required by seqan3.
                    contrib/edlib # 3rd party library for fast Levenshtein seqence matching.
                     /usr/local/include/libhmsbeagle-1  # the beagle mcmc library used by kpl
                    kel_library                 # Basic Infrastructure library
                    kgl_library                  # Basic Genetic library
                    kgl_mixture                  # Deconvolution and mixture library
                    kgl_analytic                  # Genetic analysis library
                    kgl_app)



# add boost library directories
set(BOOST_INCLUDEDIR /usr/local/include)
set(BOOST_LIBRARYDIR /usr/local/lib/boost)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
include_directories(${Boost_INCLUDE_DIRS})
find_package( Boost COMPONENTS system filesystem timer chrono program_options iostreams REQUIRED )

# Debug g++
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fconcepts -ggdb -DBOOST_MATH_DISABLE_FLOAT128")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-ggdb -DBOOST_MATH_DISABLE_FLOAT128")
# Release g++
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fconcepts -O3 -s -DNDEBUG -DBOOST_MATH_DISABLE_FLOAT128")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE " -O3 -s -DNDEBUG -DBOOST_MATH_DISABLE_FLOAT128")

# Profile g++
set(CMAKE_CXX_FLAGS_PROFILE "-fconcepts -pg -g -ggdb -O3 -DBOOST_MATH_DISABLE_FLOAT128 -Xlinker -Map=output.map ")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg -g -ggdb -O3 -DBOOST_MATH_DISABLE_FLOAT128 -Xlinker -Map=output.map ")

# All possible warnings
#set(CMAKE_ERROR_WARNING_FLAGS "-std=c++17 -Wall -pthread -Wextra -Werror -Wpedantic -pedantic-errors -Wstrict-aliasing -Wno-variadic-macros")
# Relax slightly to temporarily allow unused parameters.
set(CMAKE_ERROR_WARNING_FLAGS "-std=c++17 -Wall -pthread -Wextra -Wpedantic -Wstrict-aliasing -Wno-variadic-macros")

# Actually set the g++ flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_LINKER_FLAGS_RELEASE}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_ERROR_WARNING_FLAGS} ${CMAKE_LINKER_FLAGS_RELEASE}")


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)

# Basic infrastructure library
set(LIBRARY_SOURCE_FILES
        kel_library/kel_logging.h
        kel_library/kel_logging.cc
        kel_library/kel_utility.h
        kel_library/kel_utility.cc
        kel_library/kel_exec_env.h
        kel_library/kel_exec_env.cc
        kel_library/kel_property_tree.cc
        kel_library/kel_property_tree.h
        kel_library/kel_distribution.h
        kel_library/kel_distribution.cc
        kel_library/kel_basic_io.cc
        kel_library/kel_basic_io.h)

# Basic Genetic infrastructure library
set(GENETIC_LIBRARY_SOURCE_FILES
        kel_library/kel_mt_queue.h
        kgl_library/kgl_genome_types.h
        kel_library/kel_lock.h
        kgl_library/kgl_gff_fasta.h
        kgl_library/kgl_gff_fasta.cc
        kgl_library/kgl_genome_db.h
        kel_library/kel_patterns.h
        kgl_library/kgl_genome_feature.h
        kgl_library/kgl_genome_feature.cc
        kgl_library/kgl_variant.h
        kgl_library/kgl_variant.cc
        kgl_library/kgl_filter.h
        kgl_library/kgl_filter.cc
        kgl_library/kgl_variant_set.cc
        kgl_library/kgl_table.h
        kgl_library/kgl_table.cc
        kgl_library/kgl_table_ncbi.h
        kgl_library/kgl_sequence_base.h
        kgl_library/kgl_sequence_base.cc
        kgl_library/kgl_sequence_amino.h
        kgl_library/kgl_sequence_amino.cc
        kgl_library/kgl_alphabet_amino.h
        kgl_library/kgl_variant_db.h
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
        kgl_library/kgl_alphabet_coding_dna5.h
        kgl_library/kgl_alphabet_coding_dna5.cc
        kgl_library/kgl_variant_evidence.h
        kgl_library/kgl_variant_evidence.cc
        kgl_library/kgl_statistics_upgma.h
        kgl_library/kgl_statistics_upgma.cc
        kgl_library/kgl_variant_mutation.cc
        kgl_library/kgl_variant_factory_vcf.h
        kgl_library/kgl_variant_factory_vcf.cc
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
        kgl_library/kgl_variant_factory_vcf_parse_cigar.h
        kgl_library/kgl_variant_factory_vcf_parse_cigar.cc
        kgl_library/kgl_variant_factory_record_vcf_impl.h
        kgl_library/kgl_variant_factory_record_vcf_impl.cc
        kgl_library/kgl_genome_aux_csv.h
        kgl_library/kgl_genome_aux_csv_impl.cc
        kgl_library/kgl_variant_db_homologous.h
        kgl_library/kgl_variant_db_homologous.cc
        kgl_library/kgl_variant_db_genome.cc
        kgl_library/kgl_variant_db_population.cc
        kgl_library/kgl_variant_db_contig.cc
        kgl_library/kgl_variant_db_unphased.h
        kgl_library/kgl_variant_db_unphased.cc
        kgl_library/kgl_variant_db_unphased_population.h
        kgl_library/kgl_variant_db_unphased_population.cc
        kgl_library/kgl_variant_vcf.h
        kgl_library/kgl_variant_vcf.cc
        kgl_library/kgl_variant_vcf_mutation.cc
        kgl_library/kgl_genome_contig_feature.h
        kgl_library/kgl_genome_contig_feature.cc
        kgl_library/kgl_genome_contig_aux.cc
        kgl_library/kgl_genome_contig_aux.h
        kgl_library/kgl_phylogenetic_tree.cc
        kgl_library/kgl_phylogenetic_tree.h
        kgl_library/kgl_variant_file_vcf_record.h
        kgl_library/kgl_variant_factory_readvcf_impl.cc
        kgl_library/kgl_variant_file_impl.cc
        kgl_library/kgl_variant_file_impl.h
        kgl_library/kgl_variant_factory_grch_impl.cc
        kgl_library/kgl_variant_factory_grch_impl.h
        kgl_library/kgl_variant_file_vcf_impl.cc
        kgl_library/kgl_variant_file_vcf_impl.h
        kgl_library/kgl_variant_factory_vcf_parse_info.cc
        kgl_library/kgl_variant_factory_vcf_parse_info.h
        kgl_library/kgl_variant_factory_vcf_parse_header.cc
        kgl_library/kgl_variant_factory_vcf_parse_header.h
        kgl_library/kgl_variant_factory_vcf_evidence.cc
        kgl_library/kgl_variant_factory_vcf_evidence.h
        kgl_library/kgl_variant_factory_vcf_evidence_data.cc
        kgl_library/kgl_variant_factory_vcf_evidence_data.h
        kgl_library/kgl_variant_factory_vcf_evidence_memory.cc
        kgl_library/kgl_variant_factory_vcf_evidence_memory.h
        kgl_library/kgl_variant_factory_vcf_evidence_data_blk.cc
        kgl_library/kgl_variant_factory_vcf_evidence_data_blk.h
        kgl_library/kgl_variant_factory_vcf_evidence_data_blk_read.cc)

# 3rd party library for fast Leveshtein sequence matching.
set(EDLIB_SOURCE_FILES contrib/edlib/edlib.cpp )

# Genetic analysis library
set(ANALYTIC_SOURCE_FILES
        kgl_analytic/kgl_phylogenetic_analysis.h
        kgl_analytic/kgl_phylogenetic_analysis.cc
        kgl_analytic/kgl_phylogenetic_gene.h
        kgl_analytic/kgl_phylogenetic_gene.cc
        kgl_analytic/kgl_upgma_node.h
        kgl_analytic/kgl_upgma_node.cc
        kgl_analytic/kgl_rna_search.cc
        kgl_analytic/kgl_rna_search.h
        kgl_analytic/kgl_ploidy_analysis.h
        kgl_analytic/kgl_distribution_analysis.h
        kgl_analytic/kgl_distribution_analysis.cc
        kgl_analytic/kgl_upgma_unphased.h
        kgl_analytic/kgl_upgma_unphased.cc
        kgl_analytic/kgl_finestructure_analysis.h
        kgl_analytic/kgl_finestructure_analysis.cc
        kgl_analytic/kgl_sequence_complexity.h
        kgl_analytic/kgl_sequence_complexity.cc
        kgl_analytic/kgl_epigenetic_motif.cc
        kgl_analytic/kgl_epigenetic_motif.cc
        kgl_analytic/kgl_upgma.h
        kgl_analytic/kgl_analysis_all.h
        kgl_analytic/kgl_analysis_null.cc
        kgl_analytic/kgl_analysis_null.h
        kgl_analytic/kgl_analysis_interval.cc
        kgl_analytic/kgl_analysis_interval.h kgl_analytic/kgl_analysis_info_filter.cc kgl_analytic/kgl_analysis_info_filter.h)

# Deconvolution and mixture files.

set(MIXTURE_SOURCE_FILES
        kgl_mixture/kgl_variant_phase.h
        kgl_mixture/kgl_variant_phase.cc
        kgl_mixture/kgl_variant_classify.h
        kgl_mixture/kgl_variant_classify.cc
        kgl_mixture/kgl_variant_classify_io.cc
        kgl_mixture/kgl_read_phasing.h
        kgl_mixture/kgl_read_phasing.cc kgl_library/kgl_variant_factory_vcf_evidence_analysis.cc kgl_library/kgl_variant_factory_vcf_evidence_analysis.h)


# Application level objects
set(APPLICATION_SOURCE_FILES
        kgl_app/kgl_main.cc
        kgl_app/kgl_gene_app.cc
        kgl_app/kgl_legacy_Pf_analysis.cc
        kgl_app/kgl_legacy_Pf_analysis.h
        kgl_app/kgl_gene_app.h
        kgl_app/kgl_gene_parse.cc
        kgl_app/kgl_properties.cc
        kgl_app/kgl_package.cc
        kgl_app/kgl_package.h
        kgl_app/kgl_package_analysis.cc
        kgl_app/kgl_package_analysis.h
        kgl_app/kgl_runtime.cc
        kgl_app/kgl_runtime.h
        )





############kgd_deconvolv########################################################

find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()


include_directories(kgd_deconvolv/kgd_deconvolv_app
        kgd_deconvolv/kgd_deconvolv_lib
        kgd_deconvolv/contrib/codeCogs
        kgd_deconvolv/contrib/gzstream
        kgd_deconvolv/kgd_random)


set(DECONVOLV_LIB_CONTRIB_FILES
        kgd_deconvolv/contrib/codeCogs/asympt_expn.h
        kgd_deconvolv/contrib/codeCogs/errorfn.h
        kgd_deconvolv/contrib/codeCogs/errorfnc.h
        kgd_deconvolv/contrib/codeCogs/errorfnc_exp.h
        kgd_deconvolv/contrib/codeCogs/expx2.h
        kgd_deconvolv/contrib/codeCogs/gamma.h
        kgd_deconvolv/contrib/codeCogs/ln_add1.h
        kgd_deconvolv/contrib/codeCogs/log_gamma.h
        kgd_deconvolv/contrib/codeCogs/logbeta.h
        kgd_deconvolv/contrib/codeCogs/loggammafrac.h
        kgd_deconvolv/contrib/codeCogs/loggammasum.h
        kgd_deconvolv/contrib/codeCogs/machine_epsilon.h
        kgd_deconvolv/contrib/codeCogs/mersenne.hpp
        kgd_deconvolv/contrib/codeCogs/poly_eval.h
        kgd_deconvolv/contrib/codeCogs/randomSample.hpp
        kgd_deconvolv/contrib/codeCogs/stirling.h
        kgd_deconvolv/contrib/codeCogs/xsub_ln_add1.h
        kgd_deconvolv/contrib/gzstream/gzstream.cpp
        kgd_deconvolv/contrib/gzstream/gzstream.h
        kel_library/kel_distribution.h)

set(DECONVOLV_LIB_SOURCE_FILES
        kgd_deconvolv/kgd_deconvolv_lib/kgd_vcf_reader_debug.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_deploid_io_export.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_write_mcmc_related.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_deploid_io.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_deploid_io.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_exceptions.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_global.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ibdconfig.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ibdconfig.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_panel.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_panel.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_txt_reader.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_txt_reader.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_update_haplotype.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_update_haplotype.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_utility.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_utility.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_variant_index.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_variant_index.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_vcf_reader.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_vcf_reader.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_update_single_haplotype.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_update_single_haplotype.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_update_pair_haplotype.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_update_pair_haplotype.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_hprior.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_hprior.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ibdpath.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ibdpath.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_combinatorial.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_combinatorial.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_sample.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ibdrecombprobs.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ibdrecombprobs.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_deconvolv_types.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_allelefreq.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_allelefreq.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_alleleread.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_alleleread.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_ibd.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_ibd.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_base.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_base.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_hap.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_hap.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_virtual.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_virtual.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_titre.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_mcmc_titre.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_data.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_data.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_parameter.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_parameter.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_function.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_function.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_output.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_output.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_input.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_ctl_input.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_deconvolv_app.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_deconvolv_app.cc
        kgd_deconvolv/kgd_deconvolv_lib/kgd_prob_cache.h
        kgd_deconvolv/kgd_deconvolv_lib/kgd_prob_cache.cc)

# Application level objects
set(DECONVOLV_APPLICATION_SOURCE_FILES
        kgd_deconvolv/kgd_deconvolv_app/kgd_deconvolv_main.cc
        kgd_deconvolv/kgd_deconvolv_app/kgd_deconvolv_parse.cc)

#generate libraries.

add_library(kel_library_lib STATIC ${LIBRARY_SOURCE_FILES} ${EDLIB_SOURCE_FILES})

add_library(kgl_genome_lib STATIC ${GENETIC_LIBRARY_SOURCE_FILES} ${ANALYTIC_SOURCE_FILES} ${MIXTURE_SOURCE_FILES})

target_link_libraries(kgl_genome_lib  kel_library_lib ${Boost_LIBRARIES})

add_library(kgd_contrib_lib STATIC ${DECONVOLV_LIB_CONTRIB_FILES})

add_library(kgd_deconvolv_lib STATIC ${DECONVOLV_LIB_SOURCE_FILES})

target_link_libraries(kgd_deconvolv_lib kgd_contrib_lib kel_library_lib ${ZLIB_LIBRARIES})

#generate kgd_deconvolv executable
add_executable (kgd_deconvolv ${DECONVOLV_APPLICATION_SOURCE_FILES})

target_link_libraries(kgd_deconvolv kgd_deconvolv_lib kgl_genome_lib kel_library_lib)

#generate kgl_genome executable
add_executable (kgl_genome ${APPLICATION_SOURCE_FILES} kgl_app/kgl_properties.h)

target_link_libraries(kgl_genome kgl_genome_lib kgd_deconvolv_lib kel_library_lib)

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
        kpl_phylogenetic/kpl_geneticcode.cc
        kpl_phylogenetic/kpl_geneticdata.cc
        kpl_phylogenetic/kpl_genetictype.cc
        kpl_phylogenetic/kpl_node.cc
        kpl_phylogenetic/kpl_partition.cc
        kpl_phylogenetic/kpl_splittree.cc
        kpl_phylogenetic/kpl_strom.cc
        kpl_phylogenetic/kpl_tree.cc
        kpl_phylogenetic/kpl_treemanip.cc
        kpl_phylogenetic/kpl_treesummary.cc
        kpl_phylogenetic/kpl_xstrom.cc
        kpl_phylogenetic/kpl_likelihood.cc
        kpl_phylogenetic/kpl_likelihood.h
        kpl_phylogenetic/kpl_model.cc
        kpl_phylogenetic/kpl_model.h
        kpl_phylogenetic/kpl_qmatrix.cc
        kpl_phylogenetic/kpl_qmatrix.h
        kpl_phylogenetic/kpl_qmatrixnucleotide.cc
        kpl_phylogenetic/kpl_qmatrixnucleotide.h
        kpl_phylogenetic/kpl_qmatrixcodon.cc
        kpl_phylogenetic/kpl_qmatrixcodon.h
        kpl_phylogenetic/kpl_asrv.cc
        kpl_phylogenetic/kpl_asrv.h
        kpl_phylogenetic/kpl_random.cc
        kpl_phylogenetic/kpl_random.h
        kpl_phylogenetic/kpl_mcmc_updater.cc
        kpl_phylogenetic/kpl_mcmc_updater.h
        kpl_phylogenetic/kpl_mcmc_gamma.cc
        kpl_phylogenetic/kpl_mcmc_gamma.h
        kpl_phylogenetic/kpl_mcmc_chain.cc
        kpl_phylogenetic/kpl_mcmc_chain.h
        kpl_phylogenetic/kpl_mcmc_output.cc
        kpl_phylogenetic/kpl_mcmc_output.h
        kpl_phylogenetic/kpl_mcmc_dirichlet.cc
        kpl_phylogenetic/kpl_mcmc_dirichlet.h
        kpl_phylogenetic/kpl_mcmc_statefrequpdater.cc
        kpl_phylogenetic/kpl_mcmc_statefrequpdater.h
        kpl_phylogenetic/kpl_mcmc_exchangeupdater.cc
        kpl_phylogenetic/kpl_mcmc_exchangeupdater.h
        kpl_phylogenetic/kpl_mcmc_subsetupdater.cc
        kpl_phylogenetic/kpl_mcmc_subsetupdater.h
        kpl_phylogenetic/kpl_mcmc_pinvarupdater.cc
        kpl_phylogenetic/kpl_mcmc_pinvarupdater.h
        kpl_phylogenetic/kpl_mcmc_omegaupdater.cc
        kpl_phylogenetic/kpl_mcmc_omegaupdater.h
        kpl_phylogenetic/kpl_mcmc_treeupdater.cc
        kpl_phylogenetic/kpl_mcmc_treeupdater.h
        kpl_phylogenetic/kpl_mcmc_treelength.cc
        kpl_phylogenetic/kpl_mcmc_treelength.h
        kpl_phylogenetic/kpl_mcmc_polytomyupdater.cc
        kpl_phylogenetic/kpl_mcmc_polytomyupdater.h
        kpl_phylogenetic/kpl_mcmc_polytomyprior.cc
        kpl_phylogenetic/kpl_mcmc_polytomyprior.h
        kpl_phylogenetic/kpl_tree_io.cc
        kpl_phylogenetic/kpl_tree_io.h)

add_library(kpl_library_lib STATIC ${KPL_LIBRARY_SOURCE_FILES})

link_directories(/usr/local/lib)

link_directories(/usr/local/lib/ncl)

link_directories(/usr/local/lib/boost)

link_directories(/usr/local/lib/beagle/lib)

#generate kpl executable
add_executable (kpl kpl_phylogenetic/kpl_main.cpp)

#specify link libraries.
target_link_libraries(kpl kpl_library_lib kgl_genome_lib kel_library_lib libncl.a ${Boost_LIBRARIES} hmsbeagle)

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
        #                    COMMENT "Generating API documentation with Doxygen"
        #                    VERBATIM )

    endif(DOXYGEN_FOUND)

endif(BUILD_KGL_KGD_DOC)

