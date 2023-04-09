//
// Created by kellerberrin on 9/04/23.
//


#include "kgl_analysis_PfEMP_variant.h"
#include "kgl_variant_filter.h"

#include <fstream>



void kgl::GenomeGeneVariantAnalysis::createGeneMap(const GeneVector& gene_vector) {

  variant_gene_map_.clear();

  for (auto const& gene_ptr : gene_vector) {

    std::shared_ptr<ContigDB> gene_contig_ptr(std::make_shared<ContigDB>(gene_ptr->id()));

    std::pair<std::shared_ptr<const GeneFeature>, std::shared_ptr<ContigDB>> map_value{ gene_ptr, gene_contig_ptr};

    auto [iter, result] = variant_gene_map_.insert({gene_ptr->id(), map_value});

    if(not result) {

      ExecEnv::log().warn("PGenomeGeneVariantAnalysis::createGeneMap; Duplicate Gene Id: {}", gene_ptr->id());

    }

  }

}


// Get variants only occurring within codingFeatures for all mRNA sequences.
void kgl::GenomeGeneVariantAnalysis::getGeneVariants(const std::shared_ptr<const PopulationDB>& population_ptr) {

  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      if (contig_ptr->getMap().empty()) {

        continue;

      }

      for (auto const &[gene_id, value_pair]: variant_gene_map_) {

        auto const &[gene_ptr, gene_contig_ptr] = value_pair;

        if (gene_ptr->contig()->contigId() == contig_id) {

          auto coding_sequence_array = GeneFeature::getTranscriptionSequences(gene_ptr);

          // Only analyze the first sequence.
          if (not coding_sequence_array->empty()) {

            auto sequence_ptr = coding_sequence_array->getFirst();

            for (const auto &[feature_id, feature_ptr]: sequence_ptr->getFeatureMap()) {

              gene_contig_ptr->merge(contig_ptr->subset(feature_ptr->sequence().begin(), feature_ptr->sequence().end()));

            } // For all cds

          } // if not empty

        } // if same contig.

      } // genemap

    } // contig

  } // genome

}


void kgl::GenomeGeneVariantAnalysis::writeGeneResults(const std::string& variant_file_name) {

  std::ofstream variant_file(variant_file_name);

  if (not variant_file.good()) {

    ExecEnv::log().error("GenomeGeneVariantAnalysis::createGeneMap; Unable to open gene variant results file: {}", variant_file_name);
    return;

  }

  for (auto const& [gene_id, value_pair] : variant_gene_map_) {

    auto const& [gene_ptr, contig_ptr] = value_pair;

    auto sequence_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);

    size_t coding_length{0};
    if (not sequence_array_ptr->empty()) {

      coding_length = sequence_array_ptr->getFirst()->codingNucleotides();

    }

    // Actual variants
    size_t variant_count = contig_ptr->variantCount();
    // Determine the number of unique variants.
    auto unique_variant_contig = contig_ptr->filterVariants(UniqueUnphasedFilter());
    size_t unique_variant_count = unique_variant_contig->variantCount();

    double variant_rate{0.0};
    double unique_variant_rate{0.0};
    if (coding_length > 0) {

      variant_rate = static_cast<double>(variant_count) / static_cast<double>(coding_length);
      unique_variant_rate = static_cast<double>(unique_variant_count) / static_cast<double>(coding_length);

    }

    variant_file << gene_id
                 << CSV_DELIMITER_
                 << "Coding Length:"
                 << CSV_DELIMITER_
                 << coding_length
                 << CSV_DELIMITER_
                 << "Variant Count:"
                 << CSV_DELIMITER_
                 << variant_count
                 << CSV_DELIMITER_
                 << "Variant Rate:"
                 << CSV_DELIMITER_
                 << variant_rate
                 << CSV_DELIMITER_
                 << "Unique Variants:"
                 << CSV_DELIMITER_
                 << unique_variant_count
                 << CSV_DELIMITER_
                 << "Unique Variant Rate:"
                 << CSV_DELIMITER_
                 << unique_variant_rate
                 << CSV_DELIMITER_
                 << gene_ptr->featureText(CSV_DELIMITER_)
                 << '\n';
  }

}
