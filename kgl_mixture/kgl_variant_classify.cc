//
// Created by kellerberrin on 15/06/18.
//


#include "kgl_variant_classify.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant_single.h"
#include <fstream>
#include <algorithm>


namespace kgl = kellerberrin::genome;

namespace kgd = kellerberrin::deconvolv;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects accept unphased variants and produces ref/alt read statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Classify the ref/alt counts by variant (row) and then genome (column)
kgl::VariantClassifier::VariantClassifier(std::shared_ptr<const UnphasedPopulation> vcf_population_ptr) {

  // Get a list of all genomes.
  genomes_ = vcf_population_ptr->genomeList();

  // for all variants index by contig/offset
  for (auto genome : vcf_population_ptr->getMap()) {

    for(auto contig : genome.second->getMap()) {

      for(auto offset : contig.second->getMap()) {

        for(auto variant : offset.second) {

          ContigOffsetVariants ordered_variant(variant);
          auto result = variant_map_.find(ordered_variant);

          if (result != variant_map_.end()) {


            // We ignore duplicate inserts, since this is just the variant caller calling
            // 2 identical variants under the ploidy = 2 setting. Both variants have identical ref/alt ratios.
            // with a small (generally zero) ref count and a much larger alt count.
            result->second.insert(GenomeOffsetMap::value_type(variant->genomeId(), variant));

          } else {

            GenomeOffsetMap genome_map;

            auto insert_result = variant_map_.insert(VariantMap::value_type(ordered_variant, genome_map));

            if (not insert_result.second) {

              ExecEnv::log().error("Attempt to insert duplicate variant: {}",
                                   variant->output(' ',VariantOutputIndex::START_0_BASED, false));

            } else {

              // We ignore duplicate inserts, since this is just the variant caller calling
              // 2 identical variants under the ploidy = 2 setting. Both variants have identical ref/alt ratios.
              // with a small (generally zero) ref count and a much larger alt count.
              insert_result.first->second.insert(GenomeOffsetMap::value_type(variant->genomeId(), variant));

            }

          } // if found

        } // for variant

      } // for offset

    } // for contig

  } // for genome

}

kgd::MixtureDataObj kgl::VariantClassifier::convertToMixture(const GenomeId_t& genome_count, size_t min_count) const {

  std::vector<std::string> contig_vector;
  std::vector<size_t> indexOfContigStarts;
  std::vector<double> plaf_vector;
  std::vector<std::vector<size_t>> position;
  std::vector<double> refCount;
  std::vector<double> altCount;

  kgd::MixtureDataObj mixture_data;
  std::vector<size_t> current_position;
  ContigId_t current_contig;
  bool first_contig = true;

  // for all variant offsets.
  for (auto variant_offset : getMap()) {

    size_t allele_count = 0;

    // Initialize the current contig.
    if (first_contig) {

      current_contig = variant_offset.first.variant()->contigId();
      contig_vector.push_back(current_contig);
      first_contig = false;

    }

    for (auto genome : getGenomes()) {

      auto find_result = variant_offset.second.find(genome);

      if (find_result != variant_offset.second.end()) {

        std::shared_ptr<const CountEvidence> count_evidence_ptr = std::dynamic_pointer_cast<const CountEvidence>(
        find_result->second->evidence());

        if (count_evidence_ptr) {

          size_t total_count = count_evidence_ptr->refCount() + count_evidence_ptr->altCount();

          if (total_count >= min_count) {

            ++allele_count;

          }

        } else {

          ExecEnv::log().error("Variant without count evidence: {}",
                               find_result->second->output(' ', VariantOutputIndex::START_0_BASED, false));

        }

      }

    } // all genomes;

    if (allele_count > 0) {

      if (variant_offset.first.variant()->contigId() == current_contig) {

        current_position.push_back(variant_offset.first.variant()->offset());

      } else {   // new contig.

        indexOfContigStarts.push_back(current_position.size());
        position.push_back(current_position);
        current_position.clear();
        current_position.push_back(variant_offset.first.variant()->offset());

        // Set to new contig.
        current_contig = variant_offset.first.variant()->contigId();
        contig_vector.push_back(current_contig);

      }

      double allele_population_freq = static_cast<double>(allele_count) / static_cast<double>(getGenomes().size());
      plaf_vector.push_back(allele_population_freq);

    } // allele count > 0

  }  // all variants.

  // insert final position count.
  if (not current_position.empty()) {

    indexOfContigStarts.push_back(current_position.size());
    position.push_back(current_position);

  }

  return mixture_data;

}
