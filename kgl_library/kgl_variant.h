///
// Created by kellerberrin on 13/10/17.
//

#ifndef KGL_VARIANT_H
#define KGL_VARIANT_H

#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_genome_types.h"
#include "kgl_alphabet_amino.h"
#include "kgl_genome_db.h"
#include "kgl_sequence_base.h"
#include "kgl_variant_evidence.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin::genome {   //  organization level namespace


class Variant;   // Forward decl.

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  [[nodiscard]] virtual bool applyFilter(const Variant& variant) const = 0;

  [[nodiscard]] virtual std::shared_ptr<VariantFilter> clone() const = 0;

  [[nodiscard]] std::string filterName() const { return filter_name_; }

  void filterName(std::string filter_name) { filter_name_ = std::move(filter_name); }


private:

  std::string filter_name_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant offset output convention.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class VariantSequenceType { CDS_CODING, INTRON, NON_CODING };
class VariantSequence {

public:

  VariantSequence(const GenomeId_t& genome_id,
                  const ContigId_t& contig_id,
                  PhaseId_t phase_id,
                  ContigOffset_t contig_offset,
                  const VariantEvidence& evidence) : genome_id_(genome_id),
                                                     contig_id_(contig_id),
                                                     phase_id_(phase_id),
                                                     contig_offset_(contig_offset),
                                                     evidence_(evidence) {}
  virtual ~VariantSequence() = default;

  [[nodiscard]] const GenomeId_t& genomeId() const { return genome_id_; }
  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  [[nodiscard]] ContigOffset_t offset() const { return contig_offset_; }
  [[nodiscard]] PhaseId_t phaseId() const { return phase_id_; }
  virtual void updatePhaseId(PhaseId_t phase_id) { protectedPhaseId(phase_id); }

  [[nodiscard]] const VariantEvidence& evidence() const { return evidence_; }

  [[nodiscard]] std::string genomeOutput(char delimiter, VariantOutputIndex output_index) const;  // Genome information text.

  inline static constexpr PhaseId_t UNPHASED = 255;

protected:

  void protectedPhaseId(PhaseId_t phase_id) { phase_id_ = phase_id; }

private:

  GenomeId_t genome_id_;                          // The source of this variant
  ContigId_t contig_id_;                          // The contig of this variant
  PhaseId_t phase_id_;                            // The phase of this variant (which homologous contig)
  ContigOffset_t contig_offset_;                  // Location on the contig.
  VariantEvidence evidence_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Base class for a genome variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Defined Variant Types.
enum class VariantType { VCF_VARIANT };


class Variant : public VariantSequence {

public:


  Variant(const GenomeId_t& genome_id,
          const ContigId_t& contig_id,
          PhaseId_t phase_id,
          ContigOffset_t contig_offset,
          const VariantEvidence& evidence) : VariantSequence(genome_id,
                                             contig_id,
                                             phase_id,
                                             contig_offset,
                                             evidence) {}
  ~Variant() override = default;
  [[nodiscard]] bool filterVariant(const VariantFilter& filter) const { return applyFilter(filter); }

  [[nodiscard]] std::string name() const;

  [[nodiscard]] virtual VariantType variantType() const = 0;
  [[nodiscard]] virtual size_t size() const = 0;
  [[nodiscard]] virtual size_t referenceSize() const { return 1; }
  [[nodiscard]] virtual std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const = 0;
  [[nodiscard]] virtual std::string mutation(char delimiter, VariantOutputIndex output_index) const = 0;
  [[nodiscard]] virtual bool mutateSequence(SignedOffset_t offset_adjust,
                                            DNA5SequenceLinear& dna_sequence,
                                            SignedOffset_t& sequence_size_modify) const = 0;

  [[nodiscard]] virtual bool isSNP() const { return false; }

  [[nodiscard]] virtual bool equivalent(const Variant& cmp_var) const = 0;
  [[nodiscard]] bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };

  [[nodiscard]] virtual bool lessThan(const Variant& cmp_var) const = 0;
  [[nodiscard]] bool operator<(const Variant& cmp_var) const { return lessThan(cmp_var); };

private:

  [[nodiscard]] virtual bool applyFilter(const VariantFilter& filter) const = 0;

};


}   // end namespace



#endif //KGL_VARIANT_H
