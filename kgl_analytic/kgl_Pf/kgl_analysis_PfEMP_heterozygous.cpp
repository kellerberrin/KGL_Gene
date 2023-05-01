//
// Created by kellerberrin on 18/04/23.
//

#include "kgl_analysis_PfEMP_heterozygous.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include <fstream>


namespace kgl = kellerberrin::genome;



void kgl::HeteroHomoZygous::analyzeVariantPopulation(const std::shared_ptr<const PopulationDB> &gene_population_ptr,
                                                     const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr,
                                                     const std::shared_ptr<const Pf7SampleResource>& Pf7_sample_ptr,
                                                     const std::shared_ptr<const Pf7SampleLocation>& Pf7_physical_distance_ptr) {

  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    double FWS_statistic = Pf7_fws_ptr->getFWS(genome_id);
    std::string city;
    std::string country;
    std::string study;
    std::string year;

    if (Pf7_sample_ptr->getMap().contains(genome_id)) {

      auto iter = Pf7_sample_ptr->getMap().find(genome_id);
      auto const& [sample_id, sample_record] = *iter;
      city = sample_record.location1_;
      country = sample_record.country_;
      study = sample_record.study_;
      year = sample_record.year_;

    }

    auto [genome_iter, result] = variant_analysis_map_.try_emplace(genome_id, genome_id, FWS_statistic, city, country, study, year);
    auto& [_genome_id, analysis_obj] = *genome_iter;

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      auto [contig_iter, contig_result] = analysis_obj.getMap().try_emplace(contig_id);
      auto& [_contig_id, contig_count] = *contig_iter;

      if (contig_ptr->variantCount() == 0) {

        continue;

      }

      contig_count.total_variants_ = contig_ptr->variantCount();

      for (auto const& [offset, offset_array_ptr] : contig_ptr->getMap()) {

        if (offset_array_ptr->getVariantArray().empty()) {

          continue;

        }

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

            ExecEnv::log().warn("HeteroHomoZygous::analyzeVariantPopulation; genome: {}, contig: {}, offset: {}, Variant: {}, Cigar: {}, Format: {}",
                                genome_id,
                                contig_id,
                                offset,
                                variant_ptr->HGVS(),
                                variant_ptr->alternateCigar(),
                                variant_ptr->evidence().output(',', VariantOutputIndex::START_0_BASED));


          }

        }

      }

    }

  }

}


void kgl::HeteroHomoZygous::write_results(const std::string& file_name) {

  std::ofstream analysis_file(file_name);

  if (not analysis_file.good()) {

    ExecEnv::log().error("GenomeGeneVariantAnalysis::setGeneVector; Unable to open gene variant results file: {}", file_name);
    return;

  }

  // Write the header.
  // Get the number of contigs.
  auto& [genome, contig_data] = *variant_analysis_map_.begin();
  size_t contig_count = contig_data.getMap().size();

  analysis_file << "Genome" << CSV_DELIMITER_
                << "FWS" << CSV_DELIMITER_
                << "City" << CSV_DELIMITER_
                << "Country" << CSV_DELIMITER_
                << "Study" << CSV_DELIMITER_
                << "Year";

  for (size_t i = 0; i < contig_count; ++i) {

    analysis_file << CSV_DELIMITER_
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
                  << "Indel";

  }

  analysis_file << '\n';

  for (auto& [genome_id, contig_map] : variant_analysis_map_) {

    analysis_file << genome_id << CSV_DELIMITER_
                  << contig_map.getFWS() << CSV_DELIMITER_
                  << contig_map.getCity() << CSV_DELIMITER_
                  << contig_map.getCountry() << CSV_DELIMITER_
                  << contig_map.getStudy() << CSV_DELIMITER_
                  << contig_map.getYear();

    for (auto& [contig_id, variant_counts] : contig_map.getMap()) {

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

