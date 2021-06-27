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
    /*
    for (auto const& [genome_id, genome_ptr] : population->getMap()) {

      ExecEnv::log().info("NullAnalysis::fileReadAnalysis; Diploid Population: {}, GenomeID: {}, Contigs: {}",
                          population->populationId(), genome_id, genome_ptr->getMap().size());

      for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

        ExecEnv::log().info("NullAnalysis::fileReadAnalysis; ContigID: {}, Contig Offsets: {}, Contig Variants: {}",
                            contig_id, contig_ptr->getMap().size(), contig_ptr->variantCount());

      }

    }

    ExecEnv::log().info("NullAnalysis::fileReadAnalysis; Diploid Population: ,{}, Total Genomes: {}",
                        population->populationId(), genome_count_);

*/

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
//    ExecEnv::log().info("Genomes: {} containing variants: {}", genome_count, variant_id);

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


  return true;

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

