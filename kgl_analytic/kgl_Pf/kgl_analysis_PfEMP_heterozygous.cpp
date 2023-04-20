//
// Created by kellerberrin on 18/04/23.
//

#include "kgl_analysis_PfEMP_heterozygous.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include <fstream>


namespace kgl = kellerberrin::genome;



void kgl::HeteroHomoZygous::analyzeVariantPopulation(const std::shared_ptr<const PopulationDB> &gene_population_ptr) {

  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    auto& genome_count_map = variant_analysis_map_[genome_id];

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      auto& contig_count = genome_count_map[contig_id];

      if (contig_ptr->variantCount() == 0) {

        continue;

      }

      contig_count.total_variants_ = contig_ptr->variantCount();

      for (auto const& [offset, offset_array_ptr] : contig_ptr->getMap()) {

        for (auto const& variant_ptr : offset_array_ptr->getVariantArray()) {

          if (variant_ptr->isSNP()) {

            ++contig_count.snp_count_;

          } else {

            ++contig_count.indel_count_;

          }

        }

        if (offset_array_ptr->getVariantArray().size() == 1) {

          ++contig_count.single_variant_;

        } else if (offset_array_ptr->getVariantArray().size() == 2) {

          if (offset_array_ptr->getVariantArray().front()->HGVS() == offset_array_ptr->getVariantArray().back()->HGVS()) {

            ++contig_count.homozygous_count_;

          } else {

            ++contig_count.heterozygous_count_;

          }

        } else {

          ExecEnv::log().warn("HeteroHomoZygous::analyzeVariantPopulation; genome: {}, contig: {}, offset: {}, unexpected number of variants: {}",
                              genome_id, contig_id, offset, offset_array_ptr->getVariantArray().size());
          for (auto const& variant_ptr : offset_array_ptr->getVariantArray()) {

            ExecEnv::log().warn("HeteroHomoZygous::analyzeVariantPopulation; genome: {}, contig: {}, offset: {}, Variant: {}, Format: {}, VQSLOD: {}",
                                genome_id, contig_id, offset, variant_ptr->HGVS(),
                                variant_ptr->evidence().output(',', VariantOutputIndex::START_0_BASED),
                                getVQSLOD(variant_ptr));

          }

        }

      }

    }

  }

}


void kgl::HeteroHomoZygous::write_results(const std::string& file_name) const {

  std::ofstream analysis_file(file_name);

  if (not analysis_file.good()) {

    ExecEnv::log().error("GenomeGeneVariantAnalysis::setGeneVector; Unable to open gene variant results file: {}", file_name);
    return;

  }

  analysis_file << "Genome"
               << CSV_DELIMITER_
               << "Contig"
               << CSV_DELIMITER_
               << "Variant Count"
               << CSV_DELIMITER_
               << "Single Variant"
               << CSV_DELIMITER_
               << "Homozygous"
               << CSV_DELIMITER_
               << "Heterozygous"
               << CSV_DELIMITER_
               << "SNP"
               << CSV_DELIMITER_
               << "Indel"
               << '\n';

  for (auto const& [genome_id, contig_map] : variant_analysis_map_) {

    analysis_file << genome_id;

    for (auto const& [contig_id, variant_counts] : contig_map) {

      analysis_file << CSV_DELIMITER_
                    << contig_id
                    << CSV_DELIMITER_
                    << variant_counts.total_variants_
                    << CSV_DELIMITER_
                    << variant_counts.single_variant_
                    << CSV_DELIMITER_
                    << variant_counts.homozygous_count_
                    << CSV_DELIMITER_
                    << variant_counts.heterozygous_count_
                    << CSV_DELIMITER_
                    << variant_counts.snp_count_
                    << CSV_DELIMITER_
                    << variant_counts.indel_count_;

    }

    analysis_file << '\n';

  }

}

double kgl::HeteroHomoZygous::getVQSLOD(const std::shared_ptr<const Variant>& variant_ptr) {

  auto field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, VQSLOD_FILTER_NAME_);
  if (field_opt) {

    std::vector<double> float_vector = InfoEvidenceAnalysis::varianttoFloats(field_opt.value());
    if (float_vector.size() == 1) {

      return float_vector.front();
//      return std::pow(10.0, float_vector.front());

    }

  }

  return std::nan("n/a");

}

