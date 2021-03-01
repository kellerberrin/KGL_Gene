//
// Created by kellerberrin on 4/5/20.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_analysis_null.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::NullAnalysis::initializeAnalysis(const std::string& work_directory,
                                           const ActiveParameterList& named_parameters,
                                           std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called with file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();

  // If the file is a population then list genomes and all available info fields.
  if (file_characteristic.data_implementation == DataImplEnum::PopulationVariant) {

    std::shared_ptr<const PopulationDB> population = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);
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
                        population->populationId(), genome_count);

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

      // Investigate vep field values.

//      InfoEvidenceAnalysis::vepSubFieldValues("Consequence", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("IMPACT", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("Feature_type", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("BIOTYPE", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("EXON", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("INTRON", population);
    InfoEvidenceAnalysis::vepSubFieldValues("LoF", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("LoF_filter", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("LoF_flags", population);
//     InfoEvidenceAnalysis::vepSubFieldValues("CLIN_SIG", population);
//      InfoEvidenceAnalysis::vepSubFieldValues("LoF_info", population);



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

