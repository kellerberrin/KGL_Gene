//
// Created by kellerberrin on 12/11/17.
//


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_filter.h"
#include "kgl_gff_fasta.h"

namespace kgl = kellerberrin::genome;


bool kgl::GenomeVariant::outputCSV(const std::string& file_name, VariantOutputIndex output_index, bool detail) const {

  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << output(',', output_index, detail);

  return out_file.good();

}

bool kgl::GenomeVariant::mutantProteins( const ContigId_t& contig_id,
                                         const FeatureIdent_t& gene_id,
                                         const FeatureIdent_t& sequence_id,
                                         const std::shared_ptr<const GenomeDatabase>& genome_db,
                                         std::vector<std::shared_ptr<AminoSequence>>& amino_sequence_vector) const {
  // Get the contig.
  std::shared_ptr<ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }


  // Get the DNA sequence
  std::shared_ptr<DNA5SequenceCoding> dna_sequence_ptr;
  if (not contig_ptr->getDNA5SequenceCoding(coding_sequence_ptr, dna_sequence_ptr))  {

    ExecEnv::log().warn("No valid DNA sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  // Extract the variants for processing.
  OffsetVariantMap coding_variant_map;
  getCodingSortedVariants(contig_id, coding_sequence_ptr->start(), coding_sequence_ptr->end(), coding_variant_map);

  // There maybe more than one variant specified per offset.
  // If this is the case then we create alternative mutation paths.
  // In other words, there can be more than one mutant protein alternative.
  // This function is exponential. Alternatives = 2 ^ (#equal offset variants) - assuming 2 alternatives.
  // A message is issued if there are more than 32 alternatives (5 equal offset variants).
  // For performance reasons, a hard limit of 128 alternatives is imposed (can be varied).
  std::vector<OffsetVariantMap> variant_map_vector;
  size_t alternative_count = 0;
  getMutationAlternatives(coding_variant_map, variant_map_vector, alternative_count);

  for (auto variant_map : variant_map_vector) {

    // Make a copy of the protein dna.
    std::shared_ptr<DNA5SequenceCoding> copy_dna_sequence_ptr(std::make_shared<DNA5SequenceCoding>(*dna_sequence_ptr));

    // And mutate it.
    if (not GenomeVariant::mutateDNA(variant_map, sequence_id, copy_dna_sequence_ptr)) {

      ExecEnv::log().warn("Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                          contig_id, gene_id, sequence_id);
      return false;

    }

    amino_sequence_vector.push_back(contig_ptr->getAminoSequence(copy_dna_sequence_ptr));

  }

  return true;

}

// Recursively generates mutation alternatives.
// Warning - not a long function, but the coding logic is convoluted.
// Study carefully before modification.
void kgl::GenomeVariant::getMutationAlternatives(const OffsetVariantMap coding_variant_map,
                                                 std::vector<OffsetVariantMap>& variant_map_vector,
                                                 size_t& alternative_count,
                                                 size_t soft_limit,
                                                 size_t hard_limit) {

  alternative_count += 1;

  if (alternative_count == soft_limit) {

    ExecEnv::log().info("Soft limit of {} coding sequence mutations reached", soft_limit);

  }

  if (alternative_count > hard_limit) {

    ExecEnv::log().info("Hard limit of {} coding sequence mutations reached, no additional mutations generated", hard_limit);
    return;

  }

  OffsetVariantMap alternative_map;
  for (auto it = coding_variant_map.begin(); it != coding_variant_map.end(); ++it) {

    if (not alternative_map.empty()) {

      if (alternative_map.rbegin()->second->contigOffset() == it->second->contigOffset()) {

        OffsetVariantMap copy_alternative_map = alternative_map;
        copy_alternative_map.erase(std::prev(copy_alternative_map.end()));  // Pop the last element

        for (auto copy_it = it; copy_it != coding_variant_map.end(); ++copy_it) { // Copy the rest.

          copy_alternative_map.insert(*copy_it);

        }
        // Recursive call.
        getMutationAlternatives(copy_alternative_map, variant_map_vector, alternative_count, soft_limit, hard_limit);

      } else { // no matching offset variant.

        alternative_map.insert(*it);

      }

    } else { // empty

      alternative_map.insert(*it);

    }

  }

  variant_map_vector.push_back(alternative_map);

}



// Perform the actual mutation.
bool kgl::GenomeVariant::mutateDNA(const OffsetVariantMap& variant_map,
                                   const FeatureIdent_t& sequence_id,
                                   std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr) {

  // Mutate the base sequence.
  for (const auto& variant : variant_map) {

    if (not variant.second->mutateCodingSequence(sequence_id, dna_sequence_ptr)) {

      ExecEnv::log().warn("mutateDNA(), problem with variant contig: {}, offset: {}",
                          variant.second->contigId(), variant.second->offset());
      return false;

    }

  }

  return true;

}


std::ostream& operator<<(std::ostream &os, const kgl::GenomeVariant& genome_variant) {

  os << genome_variant.output(' ', kgl::VariantOutputIndex::START_1_BASED, true);
  os.flush();

  return os;

}

