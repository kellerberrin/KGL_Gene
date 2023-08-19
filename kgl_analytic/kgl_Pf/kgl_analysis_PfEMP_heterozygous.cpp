//
// Created by kellerberrin on 18/04/23.
//

#include "kgl_analysis_PfEMP_heterozygous.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_variant_filter_db_offset.h"
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

        updateVariantAnalysisType(offset_array_ptr, contig_count);

      } // offset

    } // contig

  } // genome (sample)

}


void kgl::HeteroHomoZygous::updateVariantAnalysisType( const std::shared_ptr<const OffsetDB>& offset_ptr,
                                                       VariantAnalysisType& analysis_record) {

  if (offset_ptr->getVariantArray().empty()) {

    return;

  }

  for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

    ++analysis_record.total_variants_;

    if (variant_ptr->isSNP()) {

      ++analysis_record.snp_count_;

    } else {

      ++analysis_record.indel_count_;

    }

  }

  if (offset_ptr->getVariantArray().size() == 1) {

    ++analysis_record.heterozygous_reference_minor_alleles_;

  } else {

    auto homozgygous_offset = offset_ptr->viewFilter(HomozygousFilter());
    if (not offset_ptr->getVariantArray().empty()) {

      auto unique_homozygous = offset_ptr->viewFilter(UniqueUnphasedFilter());
      analysis_record.homozygous_minor_alleles_ += unique_homozygous->getVariantArray().size();

    }

    auto heterozgygous_offset = offset_ptr->viewFilter(HeterozygousFilter());
    analysis_record.heterozygous_minor_alleles_ += heterozgygous_offset->getVariantArray().size();

  }

}


void kgl::HeteroHomoZygous::write_variant_results(const std::string& file_name, const LocationSummaryMap& location_summary) {

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
                << "FIS (inbreed)" << CSV_DELIMITER_
                << "City" << CSV_DELIMITER_
                << "Country" << CSV_DELIMITER_
                << "Region" << CSV_DELIMITER_
                << "Study" << CSV_DELIMITER_
                << "Year" << CSV_DELIMITER_
                << "Hom/Het";

  for (size_t i = 0; i <= contig_count; ++i) {

    analysis_file << CSV_DELIMITER_
                  << "Contig"
                  << CSV_DELIMITER_
                  << "Variant Count"
                  << CSV_DELIMITER_
                  << "Hom Ref (A;A)"
                  << CSV_DELIMITER_
                  << "Het Ref Minor (A;a)"
                  << CSV_DELIMITER_
                  << "Hom Minor (a;a)"
                  << CSV_DELIMITER_
                  << "Het Diff Minor (a;b)"
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
    size_t total_heterozygous = aggregated.heterozygous_reference_minor_alleles_ + aggregated.heterozygous_minor_alleles_;
    if (total_heterozygous > 0) {

      hom_het_ratio = static_cast<double>(aggregated.homozygous_minor_alleles_) / static_cast<double>(total_heterozygous);

    }

    std::string region;
    if (location_summary.contains(contig_map.getCity())) {

      auto iter = location_summary.find(contig_map.getCity());
      auto const& [location, location_record] = *iter;
      region = location_record.region_;

    }

    analysis_file << genome_id << CSV_DELIMITER_
                  << contig_map.getFWS() << CSV_DELIMITER_
                  << contig_map.getFIS() << CSV_DELIMITER_
                  << contig_map.getCity() << CSV_DELIMITER_
                  << contig_map.getCountry() << CSV_DELIMITER_
                  << region << CSV_DELIMITER_
                  << contig_map.getStudy() << CSV_DELIMITER_
                  << contig_map.getYear() << CSV_DELIMITER_
                  << hom_het_ratio;

    analysis_file << CSV_DELIMITER_
                  << "Combined"
                  << CSV_DELIMITER_
                  << aggregated.total_variants_
                  << CSV_DELIMITER_
                  << aggregated.homozygous_reference_alleles_
                  << CSV_DELIMITER_
                  << aggregated.heterozygous_reference_minor_alleles_
                  << CSV_DELIMITER_
                  << aggregated.homozygous_minor_alleles_
                  << CSV_DELIMITER_
                  << aggregated.heterozygous_minor_alleles_
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
                    << variant_counts.homozygous_reference_alleles_
                    << CSV_DELIMITER_
                    << variant_counts.heterozygous_reference_minor_alleles_
                    << CSV_DELIMITER_
                    << variant_counts.homozygous_minor_alleles_
                    << CSV_DELIMITER_
                    << variant_counts.heterozygous_minor_alleles_
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
        analysis_summary.heterozygous_reference_minor_alleles_ += het_hom_record.heterozygous_reference_minor_alleles_;
        analysis_summary.homozygous_minor_alleles_ += het_hom_record.homozygous_minor_alleles_;
        analysis_summary.heterozygous_minor_alleles_ += het_hom_record.heterozygous_minor_alleles_;
        analysis_summary.snp_count_ += het_hom_record.snp_count_;
        analysis_summary.indel_count_ += het_hom_record.indel_count_;
        analysis_summary.homozygous_reference_alleles_ += het_hom_record.homozygous_reference_alleles_;

      }

    }

  }

  return analysis_summary;

}


kgl:: LocationSummaryMap kgl::HeteroHomoZygous::location_summary( const std::shared_ptr<const Pf7SampleResource>& Pf7_sample_ptr,
                                                                  const std::shared_ptr<const Pf7SampleLocation>& Pf7_physical_distance_ptr,
                                                                  double radius_km,
                                                                  const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr) const {

  // Generate a set of Pass genomes.
  std::set<GenomeId_t> pass_genomes;
  for (auto const& [genome_id, sample_record] : Pf7_sample_ptr->getMap()) {

    if (sample_record.pass()) {

      pass_genomes.insert(genome_id);

    }

  }

  LocationSummaryMap summary_map;
  // Summarize each location.
  for (auto const& [location, location_record] : Pf7_physical_distance_ptr->locationMap()) {

    // All samples for location.
    auto radii_samples = Pf7_physical_distance_ptr->sampleRadius(location, radius_km);

    // Get a vector of samples that have FWS statistics (have passed QC).
    std::vector<GenomeId_t> radii_passed;
    for (auto const &sample: radii_samples) {

      if (pass_genomes.contains(sample)) {

        radii_passed.push_back(sample);

      }

    }

    auto aggregated = aggregateResults(radii_samples);

    auto const &[loc, location_type] = location_record.location();

    std::string type = location_type == LocationType::City ? "City" : "Country";

    double hom_het_ratio{0.0};
    size_t total_heterozygous = aggregated.heterozygous_reference_minor_alleles_ + aggregated.heterozygous_minor_alleles_;
    if (total_heterozygous > 0) {

      hom_het_ratio = static_cast<double>(aggregated.homozygous_minor_alleles_) / static_cast<double>(total_heterozygous);

    }

    double variant_rate{0.0};
    if (not radii_samples.empty()) {

      variant_rate = static_cast<double>(aggregated.total_variants_) / static_cast<double>(radii_samples.size());

    }

    double monoclonal{0.0};
    if (not radii_passed.empty()) {

      auto mono_samples = Pf7_fws_ptr->filterFWS(FwsFilterType::GREATER_EQUAL, Pf7FwsResource::MONOCLONAL_FWS_THRESHOLD, radii_passed);
      monoclonal = static_cast<double>(mono_samples.size()) / static_cast<double>(radii_passed.size());

    }

    LocationSummary location_summary;

    location_summary.location_ = location;
    location_summary.location_type_ = location_type;
    location_summary.city_ = location_record.city();
    location_summary.country_ = location_record.city();
    location_summary.region_ = location_record.region();
    location_summary.radius_km_ = radius_km;
    location_summary.radii_samples_ = radii_samples.size();
    location_summary.radii_samples_OK_ = radii_passed.size();
    location_summary.studies_ = location_record.locationStudies();
    location_summary.monoclonal_Fst_ = monoclonal;
    location_summary.hom_het_ratio_ = hom_het_ratio;
    location_summary.total_variants_ = aggregated.total_variants_;
    location_summary.variant_rate_ = variant_rate;
    location_summary.homozygous_reference_alleles_ = aggregated.homozygous_reference_alleles_;
    location_summary.heterozygous_reference_minor_alleles_ = aggregated.heterozygous_reference_minor_alleles_;
    location_summary.homozygous_minor_alleles_ = aggregated.homozygous_minor_alleles_;
    location_summary.heterozygous_minor_alleles_ = aggregated.heterozygous_minor_alleles_;
    location_summary.snp_count_ = aggregated.snp_count_;
    location_summary.indel_count_ = aggregated.indel_count_;

    summary_map[location] = location_summary;

  }

  return summary_map;

}


void kgl::HeteroHomoZygous::UpdateSampleLocation(const LocationSummaryMap& location_summary_map) {

  for (auto& [genome_id, contig_map] : variant_analysis_map_) {

    // Lookup the location summary

    auto record_iter = location_summary_map.find(contig_map.getCity());
    if (record_iter != location_summary_map.end()) {

      auto const& [city_location, city_location_summary] = *record_iter;

      // If the city has below minimum samples then check the country
      if (city_location_summary.radii_samples_OK_ < MINIMUM_LOCATION_SAMPLES_) {

        if (location_summary_map.contains(contig_map.getCountry())) {

          record_iter = location_summary_map.find(contig_map.getCountry());

        } else {

          ExecEnv::log().error("HeteroHomoZygous::UpdateSampleLocation; Unable to find the location record for sample/genome country: {}", contig_map.getCountry());
          continue;

        }

      }

    } else {

      ExecEnv::log().error("HeteroHomoZygous::UpdateSampleLocation; Unable to find the location record for sample/genome city: {}", contig_map.getCity());
      continue;

    }

    auto const& [location, location_summary] = *record_iter;
    std::vector<GenomeId_t> single_vector{ genome_id };
    auto aggregated = aggregateResults(single_vector);

    double wrights_inbreeding{0.0};
    if (location_summary.total_variants_ > 0 and aggregated.total_variants_ > 0) {

      double expected_heterozygosity = static_cast<double>(location_summary.heterozygous_minor_alleles_ + location_summary.heterozygous_reference_minor_alleles_) / static_cast<double>(location_summary.total_variants_);
      double observed_heterozygosity = static_cast<double>(aggregated.heterozygous_minor_alleles_ + aggregated.heterozygous_reference_minor_alleles_) / static_cast<double>(aggregated.total_variants_);

      wrights_inbreeding = (expected_heterozygosity - observed_heterozygosity) / expected_heterozygosity;

    }

    contig_map.setFIS(wrights_inbreeding);

  }

}



void kgl::HeteroHomoZygous::write_location_results(const std::string& file_name, const LocationSummaryMap& summary_map) const {

  std::ofstream analysis_file(file_name);

  if (not analysis_file.good()) {

    ExecEnv::log().error("HeteroHomoZygous::write_location_results; Unable to open results file: {}", file_name);
    return;

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
                << "Hom Ref (A;A)"
                << CSV_DELIMITER_
                << "Het Ref Minor (A;a)"
                << CSV_DELIMITER_
                << "Hom Minor (a;a)"
                << CSV_DELIMITER_
                << "Het Diff Minor (a;b)"
                << CSV_DELIMITER_
                << "SNP"
                << CSV_DELIMITER_
                << "Indel" << '\n';


  for (auto const& [location, summary] : summary_map) {


    analysis_file << location
                  << CSV_DELIMITER_
                  << (summary.location_type_ == LocationType::City ? "City" : "Country")
                  << CSV_DELIMITER_
                  << summary.city_
                  << CSV_DELIMITER_
                  << summary.country_
                  << CSV_DELIMITER_
                  << summary.region_
                  << CSV_DELIMITER_
                  << summary.radius_km_
                  << CSV_DELIMITER_
                  << summary.radii_samples_
                  << CSV_DELIMITER_
                  << summary.radii_samples_OK_
                  << CSV_DELIMITER_
                  << summary.studies_.size()
                  << CSV_DELIMITER_
                  << summary.monoclonal_Fst_
                  << CSV_DELIMITER_
                  << summary.hom_het_ratio_
                  << CSV_DELIMITER_
                  << summary.total_variants_
                  << CSV_DELIMITER_
                  << summary.variant_rate_
                  << CSV_DELIMITER_
                  << summary.homozygous_reference_alleles_
                  << CSV_DELIMITER_
                  << summary.heterozygous_reference_minor_alleles_
                  << CSV_DELIMITER_
                  << summary.homozygous_minor_alleles_
                  << CSV_DELIMITER_
                  << summary.heterozygous_minor_alleles_
                  << CSV_DELIMITER_
                  << summary.snp_count_
                  << CSV_DELIMITER_
                  << summary.indel_count_ << '\n';

  }


}