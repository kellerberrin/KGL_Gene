//
// Created by kellerberrin on 27/09/23.
//

#include "kgl_seq_coding.h"
#include "kgl_seq_interval.h"


namespace kgl = kellerberrin::genome;



std::optional<kgl::DNA5SequenceCoding>
kgl::CodingTranscript::codingSequence( const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                                       const std::shared_ptr<const ContigReference>& contig_ptr) {

  auto cds_interval_set = GeneIntervalStructure::transcriptIntervals(transcript_ptr);
  std::vector<OpenRightUnsigned> interval_vector(cds_interval_set.begin(), cds_interval_set.end());

  auto concat_sequence_opt = contig_ptr->sequence().concatSequences(interval_vector);
  if (not concat_sequence_opt) {

    ExecEnv::log().warn("Unable to concat sequence intervals for Gene: {}, Transcript: {}",
                        transcript_ptr->getGene()->id(), transcript_ptr->getParent()->id());

    return std::nullopt; // Return an empty coding sequence.

  }
  auto& concat_sequence = concat_sequence_opt.value();

  return concat_sequence.codingSequence(transcript_ptr->strand());

}


std::optional<kgl::IntronMap>
kgl::CodingTranscript::intronSequence( const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                                       const std::shared_ptr<const ContigReference>& contig_ptr) {

  auto intron_interval_set = GeneIntervalStructure::transcriptIntronIntervals(transcript_ptr);

  IntronMap intron_map;
  for (auto const& intron_interval : intron_interval_set) {

    auto intron_sequence_opt = contig_ptr->sequence().subOptSequence(intron_interval);
    if (not intron_sequence_opt) {

      ExecEnv::log().warn("Unable to generate intron sequence, interval {}, Gene: {}, Transcript: {}",
                          intron_interval.toString(),
                          transcript_ptr->getGene()->id(),
                          transcript_ptr->getParent()->id());

      return std::nullopt;

    }
    auto& intron_sequence = intron_sequence_opt.value();
    auto intron_coding = intron_sequence.codingSequence(transcript_ptr->strand());

    auto intron_seq_ptr = std::make_shared<const DNA5SequenceCoding>(std::move(intron_coding));
    auto key_interval = intron_seq_ptr->interval();
    intron_map.insert({key_interval, intron_seq_ptr});

  }

  return intron_map;

}
