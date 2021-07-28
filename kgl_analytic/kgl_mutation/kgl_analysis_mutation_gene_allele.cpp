//
// Created by kellerberrin on 20/6/21.
//

#include "kgl_analysis_mutation_gene_allele.h"
#include "kgl_variant_db_freq.h"
#include "kgl_variant_sort.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kel_distribution.h"

#include <fstream>

namespace kgl = kellerberrin::genome;


void kgl::GenerateGeneAllele::addSortedVariants(const std::shared_ptr<const SortedVariantAnalysis>& sorted_variants) {

  auto filtered_variants = sorted_variants->filterEnsembl(ensembl_gene_list_);

  for (auto const& [ensembl_id, variant_ptr] : filtered_variants) {

    sorted_allele_map_.emplace(ensembl_id, variant_ptr);

  }

}

void kgl::GenerateGeneAllele::filterAlleleMap(const double ALL_frequency, const double AFR_frequency) {

  EnsemblIndexMap freq_sorted_allele_map;

  for (auto const& [ensembl_id, variant_ptr] : sorted_allele_map_) {

    auto frequency_AFR_opt = FrequencyDatabaseRead::superPopFrequency(*variant_ptr, AFR_SUPER_POP_);
    auto frequency_ALL_opt =  FrequencyDatabaseRead::superPopFrequency(*variant_ptr, ALL_SUPER_POP_);

    if (frequency_AFR_opt and frequency_ALL_opt) {

      double afr_frequency = frequency_AFR_opt.value();
      double all_frequency = frequency_ALL_opt.value();

      if (afr_frequency >= AFR_frequency and all_frequency >= ALL_frequency) {

        freq_sorted_allele_map.emplace(ensembl_id, variant_ptr);

      }

    }

  }

  ExecEnv::log().info("Allele map filtered for ALL Freq >= {}, AFR Freq >= {}, unfiltered size: {}, filtered size: {}",
                      ALL_frequency, AFR_frequency, sorted_allele_map_.size(), freq_sorted_allele_map.size());


  sorted_allele_map_ = std::move(freq_sorted_allele_map);

}


void kgl::GenerateGeneAllele::writeHeader(std::ofstream& outfile, const char delimiter) {

  outfile << "EnsemblGeneId" << delimiter
          << "VariantId" << delimiter
          << "ContigId" << delimiter
          << "Offset" << delimiter
          << "Reference" << delimiter
          << "Alternate" << delimiter
          << "CiteCount" << delimiter
          << "GlobalFreq" << delimiter
          << "AFRFreq"  << delimiter;

  for (auto const& field_id : VEP_FIELD_LIST_) {

    outfile << field_id << delimiter;

  }

  outfile << "TotalCount" << delimiter
          << "AltCount" << delimiter
          << "AFRCount" << delimiter
          << "AFRAltCount" << delimiter
          << "AFR_Freq%" << delimiter
          << "NonAFR_Freq%" << delimiter
          << "Citations" << '\n';

}


void kgl::GenerateGeneAllele::writeOutput(const std::string& output_file, const char delimiter) const {

  std::ofstream out_file(output_file);

  if (not out_file.good()) {

    ExecEnv::log().error("GenerateGeneAllele::writeOutput; cannot open output file: {}", output_file);
    return;

  }

  writeHeader(out_file, delimiter);

  for (auto const& [ensembl_id, variant_ptr] : sorted_allele_map_) {

    out_file << ensembl_id << delimiter
             << variant_ptr->identifier() << delimiter
             << variant_ptr->contigId() << delimiter
             << variant_ptr->offset() << delimiter
             << variant_ptr->reference().getSequenceAsString() << delimiter
             << variant_ptr->alternate().getSequenceAsString() << delimiter;

    if (not variant_ptr->identifier().empty()) {

      auto find_result = allele_citation_ptr_->citationMap().find(variant_ptr->identifier());
      if (find_result != allele_citation_ptr_->citationMap().end()) {

        auto const& [rsid, citations] = *find_result;
        out_file << citations.size() << delimiter;

      } else {

        out_file << 0 << delimiter;

      }

    } else {

      out_file << 0 << delimiter;

    }

    double global_freq{0.0};
    auto frequency_opt = FrequencyDatabaseRead::superPopFrequency(*variant_ptr, ALL_SUPER_POP_);
    if (frequency_opt) {

      global_freq = frequency_opt.value();

    }

    out_file << global_freq << delimiter;

    double afr_freq{0.0};
    auto afr_frequency_opt = FrequencyDatabaseRead::superPopFrequency(*variant_ptr, AFR_SUPER_POP_);
    if (afr_frequency_opt) {

      afr_freq = afr_frequency_opt.value();

    }

    out_file << afr_freq << delimiter;

    auto vep_fields =  retrieveVepFields( variant_ptr, VEP_FIELD_LIST_);
    for (auto const& [field_id, field_value] : vep_fields) {

      out_file << field_value << delimiter;

    }


    auto total_allele_count_opt =  FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, ALL_SUPER_POP_);
    auto alt_allele_count_opt =  FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, ALL_SUPER_POP_);
    auto afr_total_allele_count_opt =  FrequencyDatabaseRead::superPopTotalAlleles(*variant_ptr, AFR_SUPER_POP_);
    auto afr_alt_allele_count_opt =  FrequencyDatabaseRead::superPopAltAlleles(*variant_ptr, AFR_SUPER_POP_);

    size_t total_allele_count{0};
    size_t alt_allele_count{0};
    size_t afr_total_allele_count{0};
    size_t afr_alt_allele_count{0};

    if (total_allele_count_opt and alt_allele_count_opt and afr_total_allele_count_opt and afr_alt_allele_count_opt) {

      total_allele_count = total_allele_count_opt.value();
      alt_allele_count = alt_allele_count_opt.value();
      afr_total_allele_count = afr_total_allele_count_opt.value();
      afr_alt_allele_count = afr_alt_allele_count_opt.value();

    }

    out_file << total_allele_count << delimiter
             << alt_allele_count << delimiter
             << afr_total_allele_count << delimiter
             << afr_alt_allele_count << delimiter;

    if (afr_total_allele_count > 0) {

      double afr_percent = (static_cast<double>(afr_alt_allele_count) / static_cast<double>(afr_total_allele_count)) * 100.0;
      out_file << afr_percent << delimiter;

    } else {

      out_file << 0 << delimiter;

    }

    int64_t non_afr_count = total_allele_count - afr_total_allele_count;
    if (non_afr_count  > 0) {

      double non_afr_percent = ((static_cast<double>(alt_allele_count)-static_cast<double>(afr_alt_allele_count)) / static_cast<double>(non_afr_count)) * 100.0;
      out_file << non_afr_percent  << delimiter;

    } else {

      out_file << 0 << delimiter;

    }

    if (not variant_ptr->identifier().empty()) {

      auto find_result = allele_citation_ptr_->citationMap().find(variant_ptr->identifier());
      if (find_result != allele_citation_ptr_->citationMap().end()) {

        auto const& [rsid, citations] = *find_result;
        std::string cite_string;
        for (auto const& cite : citations) {

          cite_string += cite;
          if (cite != citations.back()) {

            cite_string += CONCATENATE_VEP_FIELDS_;

          }

        }

        out_file << cite_string;

      } else {

        out_file << "";

      }

    } else {

      out_file << "";

    }

    out_file << '\n';

  }

}


std::map<std::string, std::string> kgl::GenerateGeneAllele::retrieveVepFields( const std::shared_ptr<const Variant>& variant_ptr,
                                                                               const std::vector<std::string>& field_list) {

  static std::vector<std::pair<std::string, size_t>> field_indicies;
  static bool initialized{false};
  static bool error_flag{false};

  if (not initialized) {

    field_indicies = InfoEvidenceAnalysis::getVepIndexes(*variant_ptr, field_list);
    initialized = true;

  }

  auto field_vector = InfoEvidenceAnalysis::getVepData(*variant_ptr, field_indicies);

  std::map<std::string, std::set<std::string>> aggregated_fields;
  // There may be multiple VEP fields.
  for (auto& field : field_vector) {

    if (field.size() != field_list.size()) {

      if (not error_flag) {

        ExecEnv::log().warn("GenerateGeneAllele::retrieveVepFields; requested fields: {}, returned VEP fields: {}", field_list.size(), field.size());
        error_flag = true;

      }

      continue;

    }

    for (auto const& [field_id, field_value] : field) {

      auto result = aggregated_fields.find(field_id);
      if (result == aggregated_fields.end()) {

        aggregated_fields.emplace(field_id, std::set<std::string>{field_value});

      } else {

        if (not field_value.empty()) {

          auto& [agg_field_id, agg_set_value] = *result;
          agg_set_value.insert(field_value);

        }

      }

    }

  }

  // Concatenate multiple field values.
  std::map<std::string, std::string> vep_fields;
  for (auto const&[field_id, field_value_set] : aggregated_fields) {

    std::string concat_fields;
    for (auto const& field_value : field_value_set) {

      concat_fields += field_value + CONCATENATE_VEP_FIELDS_;

    }

    vep_fields[field_id] = concat_fields;

  }

  return vep_fields;

}
