//
// Created by kellerberrin on 29/05/23.
//

#include "kgl_seq_interval.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::GeneIntervalStructure::codingInterval(const std::shared_ptr<const GeneFeature>& gene_feature) {

  gene_feature_ = gene_feature;

  // Get the gene transciptions.
  auto coding_sequence_array = GeneFeature::getTranscriptionSequences(gene_feature);
  if (not coding_sequence_array->getMap().empty()) {

    strand_ = coding_sequence_array->getFirst()->strand();

  }

  // Process all the gene transcriptions.
  for (auto const& [transcript_id, transcript_ptr] : coding_sequence_array->getMap()) {

    // Process all the CDS intervals within a transcript.
    auto sequence_intervals = transcriptIntervals(transcript_ptr);

    // Insert the transcript interval set and identifier into the transcripts map.
    auto [insert_iter, result] = gene_coding_transcripts_.try_emplace(transcript_id, sequence_intervals);
    if (not result) {

      ExecEnv::log().warn("GeneIntervalStructure::codingInterval; Gene: {}, unable to insert Transcript: {} (duplicate)",
                          gene_feature->id(), transcript_id);

    }

  } // For each coding sequence.

  // Create a union of all the coding transcripts
  for (auto const& [transcript_name, transcript_set] : gene_coding_transcripts_) {

    transcript_union_ = transcript_union_.intervalSetUnion(transcript_set);

  }

}

// Using the transcript interval map, generate a map of introns.
kgl::GeneCodingTranscriptMap kgl::GeneIntervalStructure::createIntronMap() const {

  GeneCodingTranscriptMap intron_map;

  for (auto const& [transcript_id, exon_interval_set] : gene_coding_transcripts_) {

    IntervalSetLower intron_set;

    auto iter_interval = exon_interval_set.begin();
    while (iter_interval != exon_interval_set.end()) {

      auto iter_next_interval = std::ranges::next(iter_interval, 1, exon_interval_set.end());
      if (iter_next_interval == exon_interval_set.end()) {

        break;

      }

      if (iter_interval->disjoint(*iter_next_interval) and not iter_interval->adjacent(*iter_next_interval)) {

        OpenRightUnsigned intron_interval{iter_interval->upper(), iter_next_interval->lower()};
        auto [insert_iter, result] = intron_set.insert(intron_interval);
        if (not result) {

          ExecEnv::log().warn("Unable to enter insert intron interval: {} (duplicate)", intron_interval.toString());

        }

      }

      iter_interval = iter_next_interval;

    } // Next exon pair.

    auto [insert_iter, result] = intron_map.try_emplace(transcript_id, intron_set);
    if (not result) {

      ExecEnv::log().warn("Unable to enter insert intron set for transcript: {} (duplicate)", transcript_id);

    }

  }

  return intron_map;

}

kel::IntervalSetLower
kgl::GeneIntervalStructure::transcriptIntervals(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) {

  // Process all the CDS intervals within a transcript.
  IntervalSetLower sequence_intervals;
  for (auto const& [feature_offset, coding_feature] : transcript_ptr->getFeatureMap()) {

    auto const& sequence = coding_feature->sequence();
    OpenRightUnsigned sequence_interval(sequence.begin(), sequence.end());
    // Add to the interval set.
    auto [insert_iter, result] = sequence_intervals.insert(sequence_interval);
    if (not result) {

      ExecEnv::log().warn("IntervalCodingVariants::transcriptIntervals; Gene: {}, Transcript: {} has duplicate coding region: {}",
                          transcript_ptr->getGene()->id(), transcript_ptr->getGene()->id(), sequence_interval.toString());

    }

  } // For each coding CDS feature within a transcript.

  return sequence_intervals;

}

kel::IntervalSetLower
kgl::GeneIntervalStructure::transcriptIntronIntervals(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) {

  IntervalSetLower intron_set;

  auto exon_interval_set = transcriptIntervals(transcript_ptr);

  auto iter_interval = exon_interval_set.begin();
  while (iter_interval != exon_interval_set.end()) {

    auto iter_next_interval = std::ranges::next(iter_interval, 1, exon_interval_set.end());
    if (iter_next_interval == exon_interval_set.end()) {

      break;

    }

    if (iter_interval->disjoint(*iter_next_interval) and not iter_interval->adjacent(*iter_next_interval)) {

      OpenRightUnsigned intron_interval{iter_interval->upper(), iter_next_interval->lower()};
      auto [insert_iter, result] = intron_set.insert(intron_interval);
      if (not result) {

        ExecEnv::log().warn("Unable to enter insert intron interval: {} (duplicate)", intron_interval.toString());

      }

    }

    iter_interval = iter_next_interval;

  } // Next exon pair.

  return intron_set;;

}



bool kgl::GeneIntervalStructure::codingModifier(const Variant& variant) const {

  // Check the contig id.
  if (not isSameContig(variant.contigId())) {

    return false;

  }

  // Check if the variant modifies the gene interval.
  auto const [variant_type, variant_interval] = variant.memberInterval();

  if (not gene_interval_.intersects(variant_interval)) {

   return false;

 }

  return transcript_union_.intersectsInterval(variant_interval);

}


bool kgl::GeneIntervalStructure::transcriptModifier(const Variant& variant, const FeatureIdent_t& transcript) const {

  // Check the contig id.
  if (not isSameContig(variant.contigId())) {

    return false;

  }

  // Check if the variant modifies the gene interval.
  auto const [variant_type, variant_interval] = variant.memberInterval();
  if (not variant_interval.intersects(gene_interval_)) {

    return false;

  }

  // Finally check if the variant modifies the transcript.
  auto find_iter = gene_coding_transcripts_.find(transcript);
  if (find_iter == gene_coding_transcripts_.end()) {

    return false;

  }

  auto const& [found_transcipt, transcript_intervals] = *find_iter;
  return transcript_intervals.intersectsInterval(variant_interval);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// std::multimap implementation
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

      auto [inserted_iter, result] = contig_interval_map_.try_emplace(gene_feature->contig()->contigId(), IntervalLowerMultiMap<std::shared_ptr<const GeneIntervalStructure>>());
      if (not result) {

        ExecEnv::log().error("IntervalCodingVariants::InitializeGeneVector; cannot insert contig: {} (duplicate)", gene_feature->contig()->contigId());
        continue;

      }

      contig_iter = inserted_iter;

    }

    // We have a valid contig IntervalMap so insert the GeneIntervalStructure object by the gene interval.
    auto& [contig_id, interval_map] = *contig_iter;

    auto gene_struct_ptr = std::make_shared<const GeneIntervalStructure>(gene_feature);
    auto gene_interval = gene_struct_ptr->geneInterval();
    interval_map.emplace(gene_interval, gene_struct_ptr);

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
  auto const [variant_type, variant_interval] = variant.memberInterval();

  auto gene_ptr_vector = interval_map.findIntersectsIntervals(variant_interval);
  for (auto const& gene_struct_ptr : gene_ptr_vector) {

    if (gene_struct_ptr->codingModifier(variant)) {

      return true;

    }

  }

  return false;

}

