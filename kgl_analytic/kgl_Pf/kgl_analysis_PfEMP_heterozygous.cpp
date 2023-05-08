//
// Created by kellerberrin on 18/04/23.
//

#include "kgl_analysis_PfEMP_heterozygous.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_PfEMP.h"
#include <fstream>


namespace kgl = kellerberrin::genome;



void kgl::HeteroHomoZygous::analyzeVariantPopulation(const std::shared_ptr<const PopulationDB> &gene_population_ptr,
                                                     const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr,
                                                     const std::shared_ptr<const Pf7SampleResource>& Pf7_sample_ptr) {

  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {


    auto record_iter = Pf7_sample_ptr->getMap().find(genome_id);
    if (record_iter == Pf7_sample_ptr->getMap().end()) {

      ExecEnv::log().error("HeteroHomoZygous::analyzeVariantPopulation; Unexpected, could not find sample record for genome:{}", genome_id);
      continue;

    }
    auto const& [sample_id, sample_record] = *record_iter;

    double FWS_statistic = Pf7_fws_ptr->getFWS(genome_id);

    auto [genome_iter, result] = variant_analysis_map_.try_emplace(genome_id, genome_id, sample_record, FWS_statistic);
    auto& [_genome_id, analysis_obj] = *genome_iter;

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      auto [contig_iter, contig_result] = analysis_obj.getMap().try_emplace(contig_id);
      auto& [_contig_id, contig_count] = *contig_iter;

      if (contig_ptr->variantCount() == 0) {

        continue;

      }

      for (auto const& [offset, offset_array_ptr] : contig_ptr->getMap()) {

        if (offset_array_ptr->getVariantArray().empty()) {

          continue;

        }

        for (auto const& variant_ptr : offset_array_ptr->getVariantArray()) {

          ++contig_count.total_variants_;

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

        } // het or hom.

      } // offset

    } // contig

  } // genome (sample)

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
                << "Year" << CSV_DELIMITER_
                << "Hom/Het";

  for (size_t i = 0; i <= contig_count; ++i) {

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


    std::vector<GenomeId_t> single_vector{ genome_id };
    auto aggregated = aggregateResults(single_vector);

    double hom_het_ratio{0.0};
    size_t total_heterozygous = aggregated.single_variant_ + aggregated.heterozygous_count_;
    if (total_heterozygous > 0) {

      hom_het_ratio = static_cast<double>(aggregated.homozygous_count_) / static_cast<double>(total_heterozygous);

    }

    analysis_file << genome_id << CSV_DELIMITER_
                  << contig_map.getFWS() << CSV_DELIMITER_
                  << contig_map.getCity() << CSV_DELIMITER_
                  << contig_map.getCountry() << CSV_DELIMITER_
                  << contig_map.getStudy() << CSV_DELIMITER_
                  << contig_map.getYear() << CSV_DELIMITER_
                  << hom_het_ratio;

    analysis_file << CSV_DELIMITER_
                  << "Combined"
                  << CSV_DELIMITER_
                  << aggregated.total_variants_
                  << CSV_DELIMITER_
                  << aggregated.single_variant_
                  << CSV_DELIMITER_
                  << aggregated.homozygous_count_
                  << CSV_DELIMITER_
                  << aggregated.heterozygous_count_
                  << CSV_DELIMITER_
                  << aggregated.snp_count_
                  << CSV_DELIMITER_
                  << aggregated.indel_count_;

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


kgl::VariantAnalysisType kgl::HeteroHomoZygous::aggregateResults(const std::vector<GenomeId_t>& sample_vector) const {

  VariantAnalysisType analysis_summary;
  std::set<GenomeId_t> sample_set(sample_vector.begin(), sample_vector.end());

  for (auto const& genome_id : sample_set) {

    if (variant_analysis_map_.contains(genome_id)) {

      auto const& [id, analysis_obj] = *variant_analysis_map_.find(genome_id);

      for (auto const& [contig_id, het_hom_record] : analysis_obj.getConstMap()) {

        analysis_summary.total_variants_ += het_hom_record.total_variants_;
        analysis_summary.single_variant_ += het_hom_record.single_variant_;
        analysis_summary.homozygous_count_ += het_hom_record.homozygous_count_;
        analysis_summary.heterozygous_count_ += het_hom_record.heterozygous_count_;
        analysis_summary.snp_count_ += het_hom_record.snp_count_;
        analysis_summary.indel_count_ += het_hom_record.indel_count_;

      }

    }

  }

  return analysis_summary;

}


void kgl::HeteroHomoZygous::write_location_results(const std::string& file_name,
                                                   const std::shared_ptr<const Pf7SampleResource>& Pf7_sample_ptr,
                                                   const std::shared_ptr<const Pf7SampleLocation>& Pf7_physical_distance_ptr,
                                                   double radius_km,
                                                   const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr) {

  std::ofstream analysis_file(file_name);

  if (not analysis_file.good()) {

    ExecEnv::log().error("HeteroHomoZygous::write_location_results; Unable to open results file: {}", file_name);
    return;

  }

  // Generate a set of Pass genomes.
  std::set<GenomeId_t> pass_genomes;
  for (auto const& [genome_id, sample_record] : Pf7_sample_ptr->getMap()) {

    if (sample_record.pass()) {

      pass_genomes.insert(genome_id);

    }

  }

  analysis_file << "Location"
                << CSV_DELIMITER_
                << "Type"
                << CSV_DELIMITER_
                << "City"
                << CSV_DELIMITER_
                << "Country"
                << CSV_DELIMITER_
                << "Region"
                << CSV_DELIMITER_
                << "Radius KM"
                << CSV_DELIMITER_
                << "Genomes (samples)"
                << CSV_DELIMITER_
                << "Passed QC"
                << CSV_DELIMITER_
                << "Studies"
                << CSV_DELIMITER_
                << "QC Monoclonal"
                << CSV_DELIMITER_
                << "Hom/Het"
                << CSV_DELIMITER_
                << "Variant Count"
                << CSV_DELIMITER_
                << "Variant Rate"
                << CSV_DELIMITER_
                << "Single Variant"
                << CSV_DELIMITER_
                << "Homozygous"
                << CSV_DELIMITER_
                << "Heterozygous"
                << CSV_DELIMITER_
                << "SNP"
                << CSV_DELIMITER_
                << "Indel" << '\n';


  for (auto const& [location, location_record] : Pf7_physical_distance_ptr->locationMap()) {

    // All samples for location.
    auto radii_samples = Pf7_physical_distance_ptr->sampleRadius(location, radius_km);

    // Get a vector of samples that have FWS statistics (have passed QC).
    std::vector<GenomeId_t> radii_passed;
    for (auto const& sample : radii_samples) {

      if (pass_genomes.contains(sample)) {

        radii_passed.push_back(sample);

      }

    }

    auto aggregated = aggregateResults(radii_samples);

    auto const& [loc, location_type] = location_record.location();

    std::string type = location_type == LocationType::City ? "City" : "Country";

    double hom_het_ratio{0.0};
    size_t total_heterozygous = aggregated.single_variant_ + aggregated.heterozygous_count_;
    if (total_heterozygous > 0) {

      hom_het_ratio = static_cast<double>(aggregated.homozygous_count_) / static_cast<double>(total_heterozygous);

    }

    double variant_rate{0.0};
    if (not radii_samples.empty()) {

      variant_rate = static_cast<double>(aggregated.total_variants_) / static_cast<double>(radii_samples.size());

    }

    double monoclonal{0.0};
    if (not radii_passed.empty()){

      auto mono_samples = Pf7_fws_ptr->filterFWS(FwsFilterType::GREATER_EQUAL, Pf7FwsResource::MONOCLONAL_FWS_THRESHOLD, radii_passed);
      monoclonal = static_cast<double>(mono_samples.size()) / static_cast<double>(radii_passed.size());

    }

    analysis_file << location
                  << CSV_DELIMITER_
                  << type
                  << CSV_DELIMITER_
                  << location_record.city()
                  << CSV_DELIMITER_
                  << location_record.country()
                  << CSV_DELIMITER_
                  << location_record.region()
                  << CSV_DELIMITER_
                  << radius_km
                  << CSV_DELIMITER_
                  << radii_samples.size()
                  << CSV_DELIMITER_
                  << radii_passed.size()
                  << CSV_DELIMITER_
                  << location_record.locationStudies().size()
                  << CSV_DELIMITER_
                  << monoclonal
                  << CSV_DELIMITER_
                  << hom_het_ratio
                  << CSV_DELIMITER_
                  << aggregated.total_variants_
                  << CSV_DELIMITER_
                  << variant_rate
                  << CSV_DELIMITER_
                  << aggregated.single_variant_
                  << CSV_DELIMITER_
                  << aggregated.homozygous_count_
                  << CSV_DELIMITER_
                  << aggregated.heterozygous_count_
                  << CSV_DELIMITER_
                  << aggregated.snp_count_
                  << CSV_DELIMITER_
                  << aggregated.indel_count_ << '\n';

  }


}