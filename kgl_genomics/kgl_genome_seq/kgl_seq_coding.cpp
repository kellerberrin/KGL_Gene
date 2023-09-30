//
// Created by kellerberrin on 27/09/23.
//

#include "kgl_seq_coding.h"
#include "kgl_seq_interval.h"


namespace kgl = kellerberrin::genome;


std::optional<kgl::IntronMap>
kgl::CodingTranscript::intronSequence( const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                                       const std::shared_ptr<const ContigReference>& contig_ref_ptr) {

  IntronMap intron_map;

  auto intron_interval_set = GeneIntervalStructure::transcriptIntronIntervals(transcript_ptr);
  for (auto const& intron_interval : intron_interval_set) {

    auto intron_sequence_opt = contig_ref_ptr->sequence().subSequence(intron_interval);
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
