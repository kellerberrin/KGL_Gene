//
// Created by kellerberrin on 4/5/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_variant_db_freq.h"
#include "kgl_analysis_null.h"
#include "kgl_variant_sort.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::NullAnalysis::initializeAnalysis(const std::string& work_directory,
                                           const ActiveParameterList& named_parameters,
                                           const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  for (auto const& genome_resource_ptr : resource_ptr->getResources(RuntimeResourceType::GENOME_DATABASE)) {

    auto genome_ptr = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_ptr);
    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome_ptr->genomeId());

  }

  work_directory_ = work_directory;

/*
  std::string pubmed_efetch_url{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"};
  std::string pubmed_pmid_args{"db=pubmed&id=19281305"};

  std::string pubmed_elink_url{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"};
  std::string pubmed_citations_args{"dbfrom=pubmed&linkname=pubmed_pubmed_citedin&id=19281305&id=20401335&id=21029472&id=21790707&id=21867552&id=21929748"};

  std::string api_key {"&api_key=8cd3dde4cbf1eeb71b5ae469ae8a99247609"};

  auto [summary_result, summary_text] = test_rest_api_.synchronousRequest(pubmed_efetch_url, pubmed_pmid_args + api_key);
  auto [citation_result, citation_text] = test_rest_api_.synchronousRequest(pubmed_elink_url, pubmed_citations_args + api_key);

  if (not citation_result) {

    ExecEnv::log().error("NullAnalysis::initializeAnalysis; problem with Pubmed API: {}", citation_text);

  } else {

    ExecEnv::log().info("NullAnalysis::initializeAnalysis; executed request, result: {}", citation_text);

  }

  // Sleep.
  const size_t seconds = 10;
  ExecEnv::log().info("*********** Sleep for seconds: {} ******************", seconds);
  std::chrono::seconds timespan(seconds); // or whatever
  std::this_thread::sleep_for(timespan);
*/

  std::vector<std::string> pmidids{ "19281305", "20401335", "21029472", "21790707", "21867552", "21929748"};

  auto citation_map = test_api_.getCitations(pmidids);

  for (auto const& [pmid, citations] : citation_map) {

    std::string citations_list;
    for(auto const& cite : citations) {

      citations_list += cite;
      if (cite != citations.back()) {

        citations_list += "&";

      }

    }

    ExecEnv::log().info("pmid_: {}, citations: {}, cite list: {}", pmid, citations.size(), citations_list);

  }

  auto reference_map = test_api_.getReferences(pmidids);

  for (auto const& [pmid, citations] : reference_map) {

    std::string citations_list;
    for(auto const& cite : citations) {

      citations_list += cite;
      if (cite != citations.back()) {

        citations_list += "&";

      }

    }

    ExecEnv::log().info("pmid_: {}, reference: {}, reference list: {}", pmid, citations.size(), citations_list);

  }

  auto publication_map = test_api_.getPublicationDetails(pmidids);

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called with file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();

  // If the file is a population then list genomes and all available info fields.
  if (file_characteristic.data_implementation == DataImplEnum::PopulationVariant
      and file_characteristic.data_organism == DataOrganism::HomoSapien) {

    std::shared_ptr<const PopulationDB> population = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

  if (not population) {

    ExecEnv::log().critical("NullAnalysis::fileReadAnalysis; File: {} is not a Homo Sapien Population", data_ptr->fileId());

  }

  ExecEnv::log().info("NullAnalysis::fileReadAnalysis; Diploid Population: ,{}, Total Genomes: {}",
                      population->populationId(), population->getMap().size());

    // Iterate any available info fields
    std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = population->getVCFInfoEvidenceHeader();
    if (not info_header_opt) {

      ExecEnv::log().warn("InfoFilterAnalysis::fileReadAnalysis could not get Info Field Header from VCF population: {}",
                          population->populationId());
      return true; // Not necessarily an error.

    }

    std::shared_ptr<const InfoEvidenceHeader> evidence_header_ptr = info_header_opt.value();

    ExecEnv::log().info("Analysis Id: {}, VCF File vcf_population: {}, Obtained header for: {} Info Fields (listed below)",
                        ident(), evidence_header_ptr->getConstMap().size(), population->populationId());

    // List the available INFO fields
    for (auto const&[ident, field_item] : evidence_header_ptr->getConstMap()) {

      ExecEnv::log().info("Field Id: {}, Type: {}, Number: {}, Description: {}",
                          ident, field_item.infoVCF().type, field_item.infoVCF().number,
                          field_item.infoVCF().description);

    }

    ensemblIdIndex(population);

  }

  return true;

}

void kgl::NullAnalysis::ensemblIdIndex(const std::shared_ptr<const PopulationDB>& population) {

  // Used in adhoc polymorphism analysis (chromosome 19).
  const std::vector<std::string> LILRB1_ensembl_symbol = {

      "ENSG00000167613",   // Alt  LAIR1
      "ENSG00000104972", // Alt LILRB1
      "ENSG00000204577", // Alt LILRB3
      "ENSG00000244482", // Alt LILRA6
      "ENSG00000105609", // Alt LILRB5	leukocyte immunoglobulin like receptor B5
      "ENSG00000131042", // Alt LILRB2	leukocyte immunoglobulin like receptor B2
      "ENSG00000187116", //	Alt LILRA5	leukocyte immunoglobulin like receptor A5
      "ENSG00000239961", // Alt LILRA4	leukocyte immunoglobulin like receptor A4
      "ENSG00000167618", //	Alt LAIR2	leukocyte associated immunoglobulin like receptor 2
      "ENSG00000239998", //	Alt LILRA2	leukocyte immunoglobulin like receptor A2
      "ENSG00000104974", //	Alt LILRA1	leukocyte immunoglobulin like receptor A1
      "ENSG00000186818",  // 	Alt LILRB4	leukocyte immunoglobulin like receptor B4
      "ENSG00000242019", //	Alt KIR3DL3	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail 3
      "ENSG00000243772", //	Alt KIR2DL3	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail 3
      "ENSG00000125498", //	Alt KIR2DL1	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail 1
      "ENSG00000189013", //	Alt KIR2DL4	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail 4
      "ENSG00000167633", //	Alt KIR3DL1	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail 1
      "ENSG00000221957", //	Alt KIR2DS4	killer cell immunoglobulin like receptor%2C two Ig domains and short cytoplasmic tail 4
      "ENSG00000240403",	// Alt KIR3DL2	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail 2


      "ENSG00000276452",   // LILRB1
      "ENSG00000276163",   //  LAIR1
      "ENSG00000277816",  //   LILRB3
      "ENSG00000275584",	// LILRA6	leukocyte immunoglobulin like receptor A6
      "ENSG00000277414", 	// LILRB5	leukocyte immunoglobulin like receptor B5
      "ENSG00000274513",  //	LILRB2	leukocyte immunoglobulin like receptor B2
      "ENSG00000278355", //	LILRA5	leukocyte immunoglobulin like receptor A5
      "ENSG00000276798", //	LILRA4	leukocyte immunoglobulin like receptor A4
      "ENSG00000274084", //	LAIR2	leukocyte associated immunoglobulin like receptor 2
      "ENSG00000274000", //	LILRA2	leukocyte immunoglobulin like receptor A2
      "ENSG00000274935", //	LILRA1	leukocyte immunoglobulin like receptor A1
      "ENSG00000278555",  // 	LILRB4	leukocyte immunoglobulin like receptor B4
      "ENSG00000276433", //	KIR3DL3	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail 3
      "ENSG00000273947", //	KIR2DL3	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail 3
      "ENSG00000276820", //	KIR2DL1	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail 1
      "ENSG00000276779", //	KIR2DL4	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail 4
      "ENSG00000273775", //	KIR3DL1	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail 1
      "ENSG00000274324", //	KIR2DS4	killer cell immunoglobulin like receptor%2C two Ig domains and short cytoplasmic tail 4
      "ENSG00000273735"	// KIR3DL2	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail 2

  };

  // Used in adhoc FCGR gene  polymorphism analysis (chromosome 1).
  const std::vector<std::string> FCGR_ensembl_symbol = {

      "ENSG00000150337",   // FCGR1A
      "ENSG00000143226",   // FCGR2A
      "ENSG00000072694",   // FCGR2B
      "ENSG00000162747",  // FCGR3B
      "ENSG00000203747"   // FCGR3A

  };


  ExecEnv::log().info("Starting Genome Ensembl Variant sort ...");

  auto all_ensembl_map_ptr = VariantSort::ensemblIndex(population);
  ExecEnv::log().info("All Ensembl Variant sort, variants found: {}", all_ensembl_map_ptr->size());

  std::shared_ptr<EnsemblIndexMap> ensembl_index_ptr(std::make_shared<EnsemblIndexMap>());
  VariantSort::ensemblAddIndex( population,LILRB1_ensembl_symbol, ensembl_index_ptr);

  ExecEnv::log().info("Completed Genome Ensembl Variant sort, variants found: {}", ensembl_index_ptr->size());


  for (auto const& [ensembl_id, variant_ptr] : *ensembl_index_ptr) {

    auto AN_opt = FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_ALL_);
    auto AC_opt = FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_ALL_);
    auto AN_afr_opt = FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_AFR_);
    auto AC_afr_opt = FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_AFR_);
    std::string variant_text = variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false);
    ExecEnv::log().info("Ensembl: {}, AN: {}, AC:{}, AN_Afr: {}, AC_Afr: {}, variant: {}",
                        ensembl_id, AN_opt.value(), AC_opt.value(), AN_afr_opt.value(), AC_afr_opt.value(), variant_text);

  }

}


void kgl::NullAnalysis::genomeIdIndex(const std::shared_ptr<const PopulationDB>& population) {

  // See how many genomes have the variant rs62418762.
  ExecEnv::log().info("Starting Genome Variant sort ...");

  auto genome_variant_index_ptr = VariantSort::variantGenomeIndexMT(population);

  ExecEnv::log().info("Completed Genome Variant sort");

  //    const std::vector<std::string> variant_id_vector{"rs62418762"};
  const std::vector<std::string> variant_id_vector { "rs10944597", "rs114240342", "rs114373344", "rs116065218", "rs12524302", "rs12527921",
                                                     "rs12528802", "rs12529529", "rs12661085", "rs12661202", "rs12664281", "rs12665626", "rs1328491", "rs141787228", "rs145866672",
                                                     "rs147631439", "rs149821037", "rs1999063", "rs2183001", "rs34821188", "rs34967714", "rs35673902", "rs4707739", "rs4707740",
                                                     "rs58551026", "rs58832941", "rs62418762", "rs62418800", "rs62418818", "rs62418819", "rs62420860", "rs6904002", "rs6925094",
                                                     "rs6939249", "rs6939736", "rs71897753", "rs72926662", "rs72928719", "rs75488031", "rs75869001", "rs7741533", "rs78233953",
                                                     "rs78702660", "rs78935233", "rs9353916" };

  size_t genome_count{0};
  for (auto const& [genome_id, sorted_variant_ptr] : *genome_variant_index_ptr) {

    for (auto const& variant_id : variant_id_vector) {

      auto result = sorted_variant_ptr->find(variant_id);
      if (result != sorted_variant_ptr->end()) {

        auto const& [variant_id, variant_ptr] = *result;

        ++genome_count;
        auto AN_opt = FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_ALL_);
        auto AC_opt = FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_ALL_);
        auto AN_afr_opt = FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_AFR_);
        auto AC_afr_opt = FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_AFR_);
        std::string variant_text = variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false);
        ExecEnv::log().info("Genome: {}, AN: {}, AC:{}, AN_Afr: {}, AC_Afr: {}, variant: {}",
                            genome_id, AN_opt.value(), AC_opt.value(), AN_afr_opt.value(), AC_afr_opt.value(), variant_text);

      }

    }

  }

}



void kgl::NullAnalysis::investigateVepFields(const std::shared_ptr<const PopulationDB>& population) {


  // Investigate vep field values.

    InfoEvidenceAnalysis::vepSubFieldCount("Gene", population);
    InfoEvidenceAnalysis::vepSubFieldValues("Consequence", population);
    InfoEvidenceAnalysis::vepSubFieldValues("IMPACT", population);
    InfoEvidenceAnalysis::vepSubFieldValues("Feature_type", population);
    InfoEvidenceAnalysis::vepSubFieldValues("BIOTYPE", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("EXON", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("INTRON", population);
    InfoEvidenceAnalysis::vepSubFieldValues("LoF", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("LoF_filter", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("LoF_flags", population);
    InfoEvidenceAnalysis::vepSubFieldValues("CLIN_SIG", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("LoF_info", population);
//    InfoEvidenceAnalysis::vepSubFieldValues("Gene", population);
    InfoEvidenceAnalysis::vepSubFieldValues("SOMATIC", population);


}


// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::NullAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}

