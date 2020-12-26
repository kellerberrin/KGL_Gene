//
// Created by kellerberrin on 15/06/18.
//


#include "kgl_variant_classify.h"
#include "kgl_variant_evidence.h"
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
kgl::VariantClassifier::VariantClassifier(std::shared_ptr<const PopulationVariant> vcf_population_ptr) {

  // for all variants index by contig/offset
  for (auto const& [genome, genome_ptr] : vcf_population_ptr->getMap()) {

    // Get a list of all genomes.
    genomes_.push_back(genome);

    for(auto const& contig : genome_ptr->getMap()) {

      for(auto const& offset : contig.second->getMap()) {

        for(auto const& variant : offset.second->getVariantArray()) {

          ContigOffsetVariants ordered_variant(variant);
          auto result = variant_map_.find(ordered_variant);

          if (result != variant_map_.end()) {

            auto map_result = result->second.try_emplace(genome, variant);

            if (not map_result.second) {

              ExecEnv::log().error("VariantClassifier; Attempt to insert duplicate variant: {}",
                                   variant->output(' ',VariantOutputIndex::START_0_BASED, false));

            }


          } else {

            GenomeOffsetMap genome_map;

            auto insert_result = variant_map_.insert(VariantMap::value_type(ordered_variant, genome_map));

            if (not insert_result.second) {

              ExecEnv::log().error("VariantClassifier; Attempt to insert duplicate variant: {}",
                                   variant->output(' ',VariantOutputIndex::START_0_BASED, false));

            } else {

              auto map_result= insert_result.first->second.try_emplace(genome, variant);

              if (not map_result.second) {

                ExecEnv::log().error("VariantClassifier; Attempt to insert duplicate variant: {}",
                                     variant->output(' ',VariantOutputIndex::START_0_BASED, false));

              }

            }

          } // if found

        } // for variant

      } // for offset

    } // for contig

  } // for genome

}

kgd::MixtureDataObj kgl::VariantClassifier::convertToMixture(const GenomeId_t& analysis_genome, size_t min_count, size_t max_count) const {

  std::vector<std::string> contig_vector;
  std::vector<size_t> indexOfContigStarts;
  std::vector<double> plaf_vector;
  std::vector<std::vector<size_t>> position;
  std::vector<double> refCount;
  std::vector<double> altCount;

  std::vector<size_t> current_position;
  size_t chrom_start_offset = 0;
  ContigId_t current_contig;
  bool first_contig = true;

  ExecEnv::log().info("Converting Genome: {} to Mixture Format, Minimum (Alt+Ref) Count: {}, Maximum Count: {}",
                      analysis_genome, min_count, max_count);

  // for all variant offsets.
  for (auto variant_offset : getMap()) {

    // Initialize the current contig.
    if (first_contig) {

      current_contig = variant_offset.first.variant()->contigId();
      contig_vector.push_back(current_contig);
      indexOfContigStarts.push_back(chrom_start_offset);  // setup the contig start vector.
      first_contig = false;

    }

    // Check if this is an SNP.
    if (not variant_offset.first.variant()->isSNP()) continue; // If not then skip to next variant.

    size_t allele_count = 0;
    for (auto genome : getGenomes()) {

      // Does this genome have this variant?
      auto find_result = variant_offset.second.find(genome);

      // If yes then increment the allele counter.
      if (find_result != variant_offset.second.end()) {

        std::optional<std::shared_ptr<const FormatData>> count_evidence_opt = find_result->second->evidence().formatData();

        if (count_evidence_opt) {

          size_t total_count = count_evidence_opt.value()->refCount() + count_evidence_opt.value()->altCount();

          // check that minimum ref+alt count requirement is achieved.
          if (total_count >= min_count and total_count <= max_count) {


            ++allele_count; // increment the allele count.

          }

        } else {

          ExecEnv::log().error("Variant without count evidence: {}",
                               find_result->second->output(' ', VariantOutputIndex::START_0_BASED, false));

        }

      }

    } // all genomes;

    // If more than one genome had this allele.
    if (allele_count > 0) {

      // Check if there has been a change of contig.
      if (variant_offset.first.variant()->contigId() == current_contig) {

        // store the variant position.
        current_position.push_back(variant_offset.first.variant()->offset());

      } else {   // new contig.

        chrom_start_offset += current_position.size();
        indexOfContigStarts.push_back(chrom_start_offset);  // setup the contig start vector.
        position.push_back(current_position);  // store the positions for the the current contig.
        current_position.clear(); // clear the position vector fot the next contig
        current_position.push_back(variant_offset.first.variant()->offset()); // and add this variant to it.

        current_contig = variant_offset.first.variant()->contigId(); // Set to new contig.
        contig_vector.push_back(current_contig); // push the new contig onto the contig vector.

      }

      double allele_population_freq = static_cast<double>(allele_count) / static_cast<double>(getGenomes().size());
      plaf_vector.push_back(allele_population_freq);  // push the plaf for this variant.

    } // allele count > 0


    // If yes then setup the ref and alt fields.

    if (allele_count > 0) {

      // Does the analysis genome have this variant.
      auto find_result = variant_offset.second.find(analysis_genome);

      if (find_result != variant_offset.second.end()) {

        std::optional<std::shared_ptr<const FormatData>> count_evidence_opt = find_result->second->evidence().formatData();

        if (count_evidence_opt) {

          size_t total_count = count_evidence_opt.value()->refCount() + count_evidence_opt.value()->altCount();

          // check that minimum ref+alt count requirement is achieved.
          if (total_count >= min_count and total_count <= max_count) {

            refCount.push_back(count_evidence_opt.value()->refCount()); // Record counts
            altCount.push_back(count_evidence_opt.value()->altCount());


          } else {

            refCount.push_back(0.0); // No counts.
            altCount.push_back(0.0);

          }

        } else {

          ExecEnv::log().error("Analysis Genome Variant without count evidence: {}",
                               find_result->second->output(' ', VariantOutputIndex::START_0_BASED, false));

          refCount.push_back(0.0); // No counts
          altCount.push_back(0.0);

        }

      } else {

        refCount.push_back(0.0); // No counts.
        altCount.push_back(0.0);

      }


    }


  }  // all variants.

  // insert final position count.
  if (not current_position.empty()) {

    position.push_back(current_position); // store the positions for the the current contig.

  }


  // Return the kgd data object.
  kgd::MixtureDataObj mixture_data( contig_vector,
                                    indexOfContigStarts,
                                    position,
                                    plaf_vector,
                                    refCount,
                                    altCount);

  if (not mixture_data.verifyPrint(true /* true - print object details */)) {

    ExecEnv::log().critical("VariantClassifier::convertToMixture(); Mixture object did not pass integrity checks.");

  }

  return mixture_data;

}
