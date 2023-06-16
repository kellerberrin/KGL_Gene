//
// Created by kellerberrin on 29/05/23.
//

#include "kgl_genome_interval.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>



namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::GeneIntervalStructure::codingInterval(const std::shared_ptr<const GeneFeature>& gene_feature) {

  gene_feature_ = gene_feature;

  // Get the gene transciptions.
  auto coding_sequence_array = GeneFeature::getTranscriptionSequences(gene_feature);

  // Process all the gene transcriptions.
  for (auto const& [sequence_name, coding_sequence] : coding_sequence_array->getMap()) {

    // Process all the CDS intervals within a transcript.
    IntervalSet sequence_intervals;
    for (auto const& [feature_offset, coding_feature] : coding_sequence->getFeatureMap()) {

      auto const& sequence = coding_feature->sequence();
      OpenRightInterval sequence_interval(sequence.begin(), sequence.end());
      // Add to the interval set.
      auto [insert_iter, result] = sequence_intervals.insert(sequence_interval);
      if (not result) {

        ExecEnv::log().warn("GeneIntervalStructure::codingInterval; Gene: {}, Transcript: {} has duplicate coding regions: [{}, {})",
                            gene_feature->id(), sequence_name, sequence_interval.lower(), sequence_interval.upper());

      }

    } // For each coding CDS feature within a transcript.

    // Insert the transcript interval set and identifier into the transcripts map.
    auto [insert_iter, result] = gene_coding_transcripts_.try_emplace(sequence_name, sequence_intervals);
    if (not result) {

      ExecEnv::log().warn("GeneIntervalStructure::codingInterval; Gene: {}, unable to insert Transcript: {} (duplicate)",
                          gene_feature->id(), sequence_name);

    }

  } // For each coding sequence.

}


bool kgl::GeneIntervalStructure::isMemberCoding(ContigOffset_t offset) const {

  for (auto const& [transcript_id, interval_set] : gene_coding_transcripts_) {

    if (interval_set.containsOffset(offset)) {

      return true;

    }

  }

  return false;

}


bool kgl::GeneIntervalStructure::codingModifier(const Variant& variant) const {

  // Check the contig id.
  if (not isSameContig(variant.contigId())) {

    return false;

  }

  // Check if the variant modifies the gene interval.
 if (not variant.sequenceModifier(gene_interval_.lower(), gene_interval_.size())) {

   return false;

 }

 // For all transcripts check if the variant modifies any of the coding sequences.
 for (auto const& [transcript_id, interval_set] : gene_coding_transcripts_) {

   auto [variant_offset, variant_size] = variant.extentOffset();
   auto interval_iter = interval_set.findUpperEqualIter(OpenRightInterval(variant_offset, variant_offset + variant_size));
   if (interval_iter == interval_set.end()) {

     return false;

   }

   auto const& coding_interval = *interval_iter;
   if (variant.sequenceModifier(coding_interval.lower(), coding_interval.size())) {

     return true;

   }

 }

  return false;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// GeneIntervalStructure objects are stored in an IntervalMap by their gene intervals.,
// The IntervalMaps are further indexed by contig in the ContigIntervalMap container.
// This enables the entire gene set in a genome to be stored as GeneIntervalStructure objects.
// Thus, we can test if a variant is within a gene coding region, and if so, return the gene feature that it belongs to.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::IntervalCodingVariants::IntervalCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr) {


  std::vector<std::shared_ptr<const GeneFeature>> gene_vector;
  for (auto const& [contig_id, contig_ptr] : reference_ptr->getMap()) {

    for (auto const& [gene_offset, gene_ptr] : contig_ptr->getGeneMap()) {

      gene_vector.push_back(gene_ptr);

    }

  }

  InitializeGeneVector(gene_vector);

}

void kgl::IntervalCodingVariants::InitializeGeneVector(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector) {

  for (auto const& gene_feature : gene_vector) {

    // Find the gene contig, insert if not present.
    auto contig_iter = contig_interval_map_.find(gene_feature->contig()->contigId());
    if (contig_iter == contig_interval_map_.end()) {

      auto [inserted_iter, result] = contig_interval_map_.try_emplace(gene_feature->contig()->contigId(), IntervalMap<GeneIntervalStructure>());
      if (not result) {

        ExecEnv::log().error("IntervalCodingVariants::InitializeGeneVector; cannot insert contig: {} (duplicate)", gene_feature->contig()->contigId());
        continue;

      }

      contig_iter = inserted_iter;

    }

    // We have a valid contig IntervalMap so insert the GeneIntervalStructure object by the gene interval.
    auto& [contig_id, interval_map] = *contig_iter;

    GeneIntervalStructure gene_coding_intervals(gene_feature);
    auto gene_interval = gene_coding_intervals.geneInterval();
    auto [insert_iter, result] = interval_map.try_emplace(gene_interval, gene_coding_intervals);
    if (not result) {

      ExecEnv::log().warn("IntervalCodingVariants::IntervalCodingVariants; cannot insert Gene: {}, Contig: {}, Interval: [{}, {}) (duplicate interval)"
                           , gene_feature->id(), contig_id, gene_interval.lower(), gene_interval.upper());
      auto iter = interval_map.find(gene_interval);
      if (iter != interval_map.end()) {

        const auto& [duplicate_interval, gene_coding] = *iter;
        const auto& sequence = gene_coding.getGene()->sequence();
        ExecEnv::log().warn("IntervalCodingVariants::IntervalCodingVariants; Previously inserted gene: {}, Contig: {}, Interval: [{}, {})"
                             , gene_coding.getGene()->id(), gene_coding.getGene()->contig()->contigId(), sequence.begin(), sequence.end());

      }

    }

  } // For each gene.

}

// Returns true if the variant is within a gene coding region.
bool kgl::IntervalCodingVariants::codingRegionVariant(const Variant &variant) const {

  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return false;

  }

  auto const& [contig_id, interval_map] = *contig_iter;

  // Lookup the IntervalMap to see if there is a candidate gene for this variant.
  auto candidate_gene_iter = interval_map.findUpperOffsetIter(variant.offset());
  if (candidate_gene_iter != interval_map.end()) {

    auto const& [interval_key, gene_interval_struct] = *candidate_gene_iter;
    if (gene_interval_struct.codingModifier(variant)) {

      return true;

    }

  }

  return false;

}


std::optional<std::shared_ptr<const kgl::GeneFeature>> kgl::IntervalCodingVariants::getGeneCoding(const Variant &variant) const {



  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return std::nullopt;

  }

  auto const& [contig_id, interval_map] = *contig_iter;

  // Lookup the IntervalMap to see if there is a candidate gene for this variant.
  auto candidate_gene_iter = interval_map.findUpperOffsetIter(variant.offset());
  if (candidate_gene_iter != interval_map.end()) {

    auto const& [interval_key, gene_interval_struct] = *candidate_gene_iter;
    if (gene_interval_struct.codingModifier(variant)) {

      return gene_interval_struct.getGene();

    }

  }

  return std::nullopt;

}

