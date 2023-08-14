//
// Created by kellerberrin on 29/05/23.
//

#include "kgl_genome_interval.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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

  // Create a union of all the coding transcripts
  for (auto const& [transcript_name, transcript_set] : gene_coding_transcripts_) {

    transcript_union_ = transcript_union_.intervalSetUnion(transcript_set);

  }

}


bool kgl::GeneIntervalStructure::codingModifier(const Variant& variant) const {

  // Check the contig id.
  if (not isSameContig(variant.contigId())) {

    return false;

  }

  // Check if the variant modifies the gene interval.
  auto [variant_offset, extent] = variant.extentOffset();
  OpenRightInterval variant_interval(variant_offset, variant_offset + extent);

 if (not variant_interval.intersects(gene_interval_)) {

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
  auto [variant_offset, extent] = variant.extentOffset();
  OpenRightInterval variant_interval(variant_offset, variant_offset + extent);

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

      auto [inserted_iter, result] = contig_interval_map_.try_emplace(gene_feature->contig()->contigId(), IntervalMultiMap<std::shared_ptr<const GeneIntervalStructure>>());
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
  auto const [offset, extent] = variant.extentOffset();
  auto gene_ptr_vector = interval_map.findIntersectsIntervals(OpenRightInterval(offset, offset + extent));
  for (auto const& gene_struct_ptr : gene_ptr_vector) {

    if (gene_struct_ptr->codingModifier(variant)) {

      return true;

    }

  }

  return false;

}


std::vector<std::shared_ptr<const kgl::GeneFeature>> kgl::IntervalCodingVariants::getGeneCoding(const Variant &variant) const {


  std::vector<std::shared_ptr<const kgl::GeneFeature>> gene_vector;
  // Check if the contig is defined.
  auto contig_iter = contig_interval_map_.find(variant.contigId());
  if (contig_iter == contig_interval_map_.end()) {

    return gene_vector; // Return the empty vector.

  }

  auto const& [contig_id, interval_map] = *contig_iter;

  // Lookup the IntervalMap to see if there are candidate genes for this variant.
  // Variants may (and do) map to more than one gene.
  auto const [offset, extent] = variant.extentOffset();
  auto gene_ptr_vector = interval_map.findIntersectsIntervals(OpenRightInterval(offset, offset + extent));
  for (auto const& gene_struct_ptr : gene_ptr_vector) {

    if (gene_struct_ptr->codingModifier(variant)) {

      gene_vector.push_back(gene_struct_ptr->getGene());

    }

  }

  return gene_vector;

}


