//
// Created by kellerberrin on 4/7/20.
//

#ifndef KGL_SEQ_VARIANT_DB_H
#define KGL_SEQ_VARIANT_DB_H

#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_genome_attributes.h"
#include "kgl_variant_db.h"
#include "kgl_variant_db_population.h"
#include "kgl_mutation_db.h"


namespace kellerberrin::genome {   //  organization level namespace


class GenomeMutation {

public:

  // A namespace / object for mutation functionality.
  GenomeMutation() = delete;
  ~GenomeMutation() = delete;

  // All contig_id variants use the zero-based half-open convention [start, end).
  // End points past the last variant; end = (last + 1).

  // Returns reference and protein mutations.
  [[nodiscard]] static bool mutantProteins( const ContigId_t &contig_id,
                                            const FeatureIdent_t &gene_id,
                                            const FeatureIdent_t &sequence_id,
                                            const std::shared_ptr<const GenomeReference> &genome_ref_ptr,
                                            const OffsetVariantMap &variant_map,
                                            AminoSequence &reference_sequence,
                                            AminoSequence &sequence_vector);

  // Returns reference and mutant stranded DNA.
  [[nodiscard]] static bool mutantCodingDNA(const ContigId_t &contig_id,
                                            const FeatureIdent_t &gene_id,
                                            const FeatureIdent_t &transcript_id,
                                            const std::shared_ptr<const GenomeReference> &genome_ref_ptr,
                                            const OffsetVariantMap &variant_map,
                                            DNA5SequenceCoding &reference_sequence,
                                            DNA5SequenceCoding &mutant_sequence);

  // Returns reference and mutant unstranded DNA region
  [[nodiscard]] static bool mutantRegion(const ContigId_t &contig_id,
                                         ContigOffset_t region_offset,
                                         ContigSize_t region_size,
                                         const std::shared_ptr<const GenomeReference> &genome_ref_ptr,
                                         const OffsetVariantMap &variant_map,
                                         DNA5SequenceLinear &reference_sequence,
                                         DNA5SequenceLinear &mutant_sequence);


};


} // namespace


#endif //KGL_SEQ_VARIANT_DB_H
