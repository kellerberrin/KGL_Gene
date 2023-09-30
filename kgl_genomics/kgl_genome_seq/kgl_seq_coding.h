//
// Created by kellerberrin on 27/09/23.
//

#ifndef KGL_SEQ_CODING_H
#define KGL_SEQ_CODING_H


#include "kgl_genome_contig.h"
#include "kgl_seq_coding.h"
#include "kgl_sequence_base.h"
#include "kel_interval_map.h"


namespace kellerberrin::genome {   //  organization::project level namespace


using IntronMap = IntervalLowerMultiMapType<std::shared_ptr<const DNA5SequenceCoding>>;

class CodingTranscript {

public:

  CodingTranscript() = default;
  ~CodingTranscript() = default;

  // The entire sequence defined by the TranscriptionSequence is returned.

  [[nodiscard]] static std::optional<IntronMap>
    intronSequence( const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                    const std::shared_ptr<const ContigReference>& contig_ref_ptr);

private:

};


} // Namespace.


#endif //KGL_SEQ_CODING_H
