//
// Created by kellerberrin on 5/1/21.
//

#include "kga_analysis_mutation_data.h"
#include "kga_analysis_mutation.h"
#include "kgl_analysis_gene_sequence.h"
#include "kgl_variant_filter_db_offset.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const ActiveParameterList& named_parameters,
                                               const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
                            ident_work_directory_);

  }

  // Get the analysis parameters.
  ExecEnv::log().info("Default Analysis Id: {} initialized with analysis results directory: {}", ident(), ident_work_directory_);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  if (not getParameters(named_parameters)) {

    ExecEnv::log().critical("Analysis Id: {}, problem parsing parameters, program ends.", ident());

  }

  // Get the analysis resources.
  ref_genome_ptr_ = resource_ptr->getSingleResource<const GenomeReference>(ResourceProperties::GENOME_RESOURCE_ID_);
  ontology_db_ptr_ = resource_ptr->getSingleResource<const kol::OntologyDatabase>(ResourceProperties::ONTOLOGY_RESOURCE_ID_);
  genome_aux_ptr_ = resource_ptr->getSingleResource<const HsGenomeAux>(ResourceProperties::GENOMEAUX_RESOURCE_ID_);
  uniprot_nomenclature_ptr_ = resource_ptr->getSingleResource<const UniprotResource>(ResourceProperties::GENE_NOMENCLATURE_RESOURCE_ID_,
                                                                                     ResourceProperties::NOMENCLATURE_UNIPROTID);
  ensembl_nomenclature_ptr_ = resource_ptr->getSingleResource<const EnsemblHGNCResource>(ResourceProperties::GENE_NOMENCLATURE_RESOURCE_ID_,
                                                                                         ResourceProperties::NOMENCLATURE_ENSEMBL);
  entrez_nomenclature_ptr_ = resource_ptr->getSingleResource<const EntrezResource>(ResourceProperties::ENTREZ_RESOURCE_ID_);
  allele_citation_ptr_ = resource_ptr->getSingleResource<const CitationResource>(ResourceProperties::CITATION_RESOURCE_ID_);
  pubmed_requestor_ptr_ = resource_ptr->getSingleResource<const PubmedRequester>(ResourceProperties::PUBMED_API_RESOURCE_ID_);

 // Update the template populations.
  gene_mutation_.genomeAnalysis( MutationAnalysisData::OMIMGeneSymbol(),
                                 ref_genome_ptr_,
                                 genome_aux_ptr_,
                                 ontology_db_ptr_,
                                 uniprot_nomenclature_ptr_,
                                 ensembl_nomenclature_ptr_,
                                 entrez_nomenclature_ptr_);


  // Initialize the gene allele analysis objects with resources.
  gene_alleles_.initialize( MutationAnalysisData::OMIMGeneSymbol(),
                            uniprot_nomenclature_ptr_,
                            entrez_nomenclature_ptr_,
                            allele_citation_ptr_,
                            pubmed_requestor_ptr_);

  all_pmid_alleles_.initialize( MutationAnalysisData::OMIMGeneSymbol(),
                                uniprot_nomenclature_ptr_,
                                entrez_nomenclature_ptr_,
                                allele_citation_ptr_,
                                pubmed_requestor_ptr_);

  pop_pub_alleles_.initialize( genome_aux_ptr_,
                               uniprot_nomenclature_ptr_,
                               entrez_nomenclature_ptr_,
                               allele_citation_ptr_,
                               pubmed_requestor_ptr_);

  return true;

}




bool kgl::MutationAnalysis::getParameters(const ActiveParameterList& named_parameters) {

  for (auto const& named_block : named_parameters.getMap()) {

    auto [block_name, block_vector] = named_block.second;

    if (block_vector.size() != 1) {

      ExecEnv::log().error("MutationAnalysis::getParameters; parameter block: {} vector size: {}, expected size = 1",
                           block_name, block_vector.size());
      return false;

    }

    ExecEnv::log().info("Analysis: {} parsing parameter block: {}", ident(), block_name);

    for (auto const& xml_vector : block_vector) {


      auto output_opt = xml_vector.getString(OUTPUT_FILE_);
      if (output_opt) {

        output_file_name_ = output_opt.value().front() + std::string(CSV_FILE_EXT_);
        output_file_name_ = Utility::filePath(output_file_name_, ident_work_directory_);
        ExecEnv::log().info("Analysis: {} outputfile: {}", ident(), output_file_name_);

      } else {

        ExecEnv::log().error("MutationAnalysis::getParameters; bad value for parameter: {}", OUTPUT_FILE_);
        return false;

      }

    }

  }

  // Recreate the Population Literature Subdirectory.
  population_lit_allele_directory_ = ident_work_directory_ + std::string("/") + POPULATION_LIT_SUBDIRECTORY_;
  if (not Utility::recreateDirectory(population_lit_allele_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::getParameters, unable to recreate gene analysis results directory: {}",
                            population_lit_allele_directory_);

  }

  // Recreate the Allele Literature Subdirectory.
  allele_literature_directory_ = ident_work_directory_ + std::string("/") + ALLELE_LIT_SUBDIRECTORY_;
  if (not Utility::recreateDirectory(allele_literature_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::getParameters, unable to recreate gene analysis results directory: {}",
                            allele_literature_directory_);

  }

  // Recreate the Allele VEP Subdirectory.
  allele_vep_directory_ = ident_work_directory_ + std::string("/") + ALLELE_VEP_SUBDIRECTORY_;
  if (not Utility::recreateDirectory(allele_vep_directory_)) {

    ExecEnv::log().critical("MutationAnalysis::getParameters, unable to recreate gene analysis results directory: {}",
                            allele_vep_directory_);

  }

  return true;

}


// Perform pre-processing (generally just type casting) for each file read into the analysis object.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called with file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();

  if (file_characteristic.data_implementation == DataImplEnum::PopulationVariant) {

    if (file_characteristic.data_source == DataSourceEnum::Genome1000
        or file_characteristic.data_source == DataSourceEnum::GnomadGenome3_1) {

      auto const_population = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

      if (not const_population) {

        ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast Genome1000 data file to population, severe error.");

      }

      // Must cast to a non-count pointer because we want to do some in-Situ filtering.
      // In-Situ filtering is used to minimise memory usage.
      auto population = std::const_pointer_cast<PopulationDB>(const_population);

      ExecEnv::log().info("Begin uniqueness filter for population: {} variant count: {}", population->populationId(), population->variantCount());
      auto pass_results = population->selfFilter(PassFilter());
      auto diploid_results = population->selfFilter(DiploidFilter());
      ExecEnv::log().info("Filtered Population: {} 'SNP and Pass' count: {}, 'Diploid' count: {}",
                          population->populationId(), pass_results.second, diploid_results.second);

      // Filtering is done so assign to the object pointer.
      population_ptr_ = population;
      // Generate another population that just contains unique variants.
      unphased_population_ptr_ = population_ptr_->uniqueUnphasedGenome();

    } else if ( file_characteristic.data_source == DataSourceEnum::Gnomad2_1
               or file_characteristic.data_source == DataSourceEnum::Gnomad3_1) {

      unphased_population_ptr_ = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

      if (not unphased_population_ptr_) {

        ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast Gnomad 2.1 data file to population, severe error.");

      }

    } else if (file_characteristic.data_source == DataSourceEnum::Clinvar) {

      clinvar_population_ptr_ = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

      if (not clinvar_population_ptr_) {

        ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast Clinvar data file to population, severe error.");

      }

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis unknown file type: {}", file_characteristic.source_text);

    }

  } else if (file_characteristic.data_source == DataSourceEnum::BioPMID) {

    auto bio_pmid_ptr = std::dynamic_pointer_cast<const BioPMIDFileData>(data_ptr);

    if (not bio_pmid_ptr) {

      ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast data file to Bio PMID object, severe error.");

    }

    ExecEnv::log().info("Bio PMID MeSH Disease Map Size: {}", bio_pmid_ptr->diseaseMeSHMap().size());
    ExecEnv::log().info("Bio PMID Entrez Gene Map Size: {}", bio_pmid_ptr->entrezMap().size());

    auto const disease_pmid_set = bio_pmid_ptr->selectDiseaseBioPMID(MutationAnalysisData::malariaMeSHList());
    auto disease_allele_map = allele_citation_ptr_->filteredAlleleIndexed(disease_pmid_set);
    gene_alleles_.addDiseaseAlleles(disease_allele_map);
    all_pmid_alleles_.addDiseaseAlleles(disease_allele_map);
    pop_pub_alleles_.addDiseaseAlleles(disease_allele_map);

    gene_mutation_.updatePMIDStatistics(disease_pmid_set, bio_pmid_ptr);

  } else {

    ExecEnv::log().error("MutationAnalysis::fileReadAnalysis unknown file type: {}", file_characteristic.source_text);

  }

  return true;

}


// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  // Check that everything is defined for the iteration analysis.
  if ( population_ptr_
       and unphased_population_ptr_
       and clinvar_population_ptr_
       and genome_aux_ptr_
       and ontology_db_ptr_) {

    // Sort variants by Gene Ensembl Code
    auto sorted_variants_ptr = std::make_shared<SortedVariantAnalysis>(unphased_population_ptr_);
    ExecEnv::log().info("Unphased filter Sorted by VEP Gene Identifier: {}, number sorted by non-ensembl identifier: {}",
                        sorted_variants_ptr->ensemblMap()->size(), VariantSort::nonEnsemblIdentifiers(*(sorted_variants_ptr->ensemblMap())));
    // Perform the analysis
    gene_mutation_.variantAnalysis( population_ptr_,
                                    unphased_population_ptr_,
                                    clinvar_population_ptr_,
                                    genome_aux_ptr_,
                                    allele_citation_ptr_,
                                    sorted_variants_ptr->ensemblMap());
    // Add the sorted variants to the gene allele analysis.
    gene_alleles_.addGeneCitedVariants(sorted_variants_ptr);
    all_pmid_alleles_.addDiseaseCitedVariants(sorted_variants_ptr);
    pop_pub_alleles_.processPopulation(population_ptr_, sorted_variants_ptr);


    // Write chromosome output.
    std::string population_lit_allele_file = population_ptr_->populationId() + std::string(POP_LIT_ALLELE_FILE_) + std::string(TEXT_FILE_EXT_);
    population_lit_allele_file = Utility::filePath(population_lit_allele_file, population_lit_allele_directory_);
    pop_pub_alleles_.writePopLiterature(population_lit_allele_file);

    // Write allele literature output.
    std::string literature_allele_file = population_ptr_->populationId() + std::string(LIT_ALLELE_FILE_) + std::string(TEXT_FILE_EXT_);
    literature_allele_file = Utility::filePath(literature_allele_file, allele_literature_directory_);
    all_pmid_alleles_.writeLiteratureAlleleSummary(literature_allele_file);

    std::string allele_literature_file = population_ptr_->populationId() + std::string(ALLELE_LIT_FILE_) + std::string(TEXT_FILE_EXT_);
    allele_literature_file = Utility::filePath(allele_literature_file, allele_literature_directory_);
    all_pmid_alleles_.writeAlleleLiteratureSummary(allele_literature_file);

    // Gene Allele Specific Analysis.
    std::string gene_allele_file = population_ptr_->populationId() + std::string(GENE_ALLELE_OUTPUT_FILE_) + std::string(CSV_FILE_EXT_);
    gene_allele_file = Utility::filePath(gene_allele_file, allele_vep_directory_);
    gene_alleles_.writeOutput(gene_allele_file, OUTPUT_DELIMITER_);

    std::string all_allele_file = population_ptr_->populationId() + std::string(ALL_ALLELE_OUTPUT_FILE_) + std::string(CSV_FILE_EXT_);
    all_allele_file = Utility::filePath(all_allele_file, allele_vep_directory_);
    all_pmid_alleles_.writeOutput(all_allele_file, OUTPUT_DELIMITER_);

    // Reset for the next chromosome.
    pop_pub_alleles_.clear();
    gene_alleles_.clear();
    all_pmid_alleles_.clear();

  }

  std::pair<size_t, size_t> mem_pair = Utility::process_mem_usage2(); // pair.first is process vm_usage, pair.second is resident memory set.
  ExecEnv::log().info("Before Clear(), Variant Objects:{}, Data Blocks:{}, VM Usage: {}, Resident Memory: {}",
                      Variant::objectCount(), DataMemoryBlock::objectCount(), mem_pair.first, mem_pair.second);
  // Explicitly clean up the populations to recover memory.
  population_ptr_ = nullptr;
  unphased_population_ptr_ = nullptr;
  // Some memory auditing code to look for memory leaks.
  AuditMemory::trimFreeStore();
  mem_pair = Utility::process_mem_usage2();
  ExecEnv::log().info("After Clear(), Variant Objects:{}, Data Blocks:{}, VM Usage: {}, Resident Memory: {}",
                      Variant::objectCount(), DataMemoryBlock::objectCount(), mem_pair.first, mem_pair.second);


  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  gene_mutation_.writeOutput(genome_aux_ptr_,  output_file_name_,  OUTPUT_DELIMITER_);

  return true;

}


