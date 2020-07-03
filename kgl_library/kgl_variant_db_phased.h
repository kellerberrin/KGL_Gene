//
// Created by kellerberrin on 8/01/18.
//

#ifndef KGL_VARIANT_DB_POPULATION_H
#define KGL_VARIANT_DB_POPULATION_H




#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ContigVariant - All the variant features that map onto that region/sequence.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using OffsetVariantMap = std::map<ContigOffset_t, std::shared_ptr<const Variant>>;

class ContigVariant : public ContigOffsetVariant<UnphasedContigListOffset> {

public:

  explicit ContigVariant(const ContigId_t &contig_id) : ContigOffsetVariant<UnphasedContigListOffset>(contig_id) {}
  ContigVariant(const ContigVariant &) = delete; // Use deep copy.
  ~ContigVariant() override = default;

  [[nodiscard]] bool getSortedVariants( PhaseId_t phase,
                                        ContigOffset_t start,
                                        ContigOffset_t end,
                                        OffsetVariantMap& variant_map) const;

  static constexpr PhaseId_t HAPLOID_HOMOLOGOUS_INDEX = 0;  // The first (and only) haploid homologous contig.

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The phased variant database Genome class
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GenomeVariant : public GenomeVariantArray<ContigVariant> {

public:

  explicit GenomeVariant(const GenomeId_t &genome_id) : GenomeVariantArray<ContigVariant>(genome_id) {}

  GenomeVariant(const GenomeVariant &) = delete; // Use deep copy.
  ~GenomeVariant() override = default;

  GenomeVariant &operator=(const GenomeVariant &) = delete; // Use deep copy.

  // All contig_id variants use the zero-based half-open convention [start, end).
  // End points past the last variant; end = (last + 1).
  [[nodiscard]] bool getSortedVariants(ContigId_t contig_id,
                                       PhaseId_t phase,
                                       ContigOffset_t start,
                                       ContigOffset_t end,
                                       OffsetVariantMap &variant_map) const;

  // Returns reference and protein mutations.
  [[nodiscard]] bool mutantProteins(const ContigId_t &contig_id,
                                    PhaseId_t phase,
                                    const FeatureIdent_t &gene_id,
                                    const FeatureIdent_t &sequence_id,
                                    const std::shared_ptr<const GenomeReference> &genome_db,
                                    OffsetVariantMap &variant_map,
                                    AminoSequence &reference_sequence,
                                    AminoSequence &sequence_vector) const;

  // Returns reference and mutant stranded DNA.
  [[nodiscard]] bool mutantCodingDNA(const ContigId_t &contig_id,
                                     PhaseId_t phase,
                                     const FeatureIdent_t &gene_id,
                                     const FeatureIdent_t &sequence_id,
                                     const std::shared_ptr<const GenomeReference> &genome_db,
                                     OffsetVariantMap &variant_map,
                                     DNA5SequenceCoding &reference_sequence,
                                     DNA5SequenceCoding &mutant_sequence) const;

  // Returns reference and mutant unstranded DNA region
  [[nodiscard]] bool mutantRegion(const ContigId_t &contig_id,
                                  PhaseId_t phase,
                                  ContigOffset_t region_offset,
                                  ContigSize_t region_size,
                                  const std::shared_ptr<const GenomeReference> &genome_db,
                                  OffsetVariantMap &variant_map,
                                  DNA5SequenceLinear &reference_sequence,
                                  DNA5SequenceLinear &mutant_sequence) const;

  // Returns pointer reference to the contig and mutant unstranded contig.
  [[nodiscard]] bool mutantContig(const ContigId_t &contig_id,
                                  PhaseId_t phase,
                                  const std::shared_ptr<const GenomeReference> &genome_db,
                                  std::shared_ptr<const DNA5SequenceContig> &reference_contig_ptr,
                                  DNA5SequenceContig &mutant_contig_ptr) const;

// Ploidy constants.
  constexpr static const PhaseId_t HAPLOID_GENOME = 1;
  constexpr static const PhaseId_t DIPLOID_GENOME = 2;

};




// General phased population.
using PhasedPopulation = PopulationVariant<GenomeVariant>;

// Haploid.
using HaploidContig = ContigOffsetVariant<HaploidOffset>;
using HaploidGenome = ContigOffsetVariant<HaploidContig>;
using HaploidPopulation = PopulationVariant<HaploidGenome>;

// Diploid.
using DiploidContig = ContigOffsetVariant<DiploidOffset>;
using DiploidGenome = GenomeVariantArray<DiploidContig>;
using DiploidPopulation = PopulationVariant<DiploidGenome>;


}   // end namespace



#endif //KGL_VARIANT_DB_POPULATION_H
