//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_VARIANT_DB_GENOME_H
#define KGL_VARIANT_DB_GENOME_H



#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"
#include "kgl_variant_db_contig.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The variant contig map
using GenomeVariantMap = std::map<ContigId_t, std::shared_ptr<ContigVariant>>;

class GenomeVariant {

public:

  explicit GenomeVariant(const GenomeId_t& genome_id, PhaseId_t ploidy) : genome_id_(genome_id), ploidy_(ploidy) {}
  GenomeVariant(const GenomeVariant&) = delete; // Use deep copy.
  ~GenomeVariant() = default;

  GenomeVariant& operator=(const GenomeVariant&) = delete; // Use deep copy.

  // Always use deep copy when modifying this object (filter and set operations).
  std::shared_ptr<GenomeVariant> deepCopy() const;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& contig_id) { genome_id_ = contig_id; }

  PhaseId_t ploidy() const{ return ploidy_; }

  size_t contigCount() const { return genome_variant_map_.size(); }

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  bool getContigVariant(const ContigId_t& contig_id, std::shared_ptr<ContigVariant>& contig_variant) const;

  bool addVariant(std::shared_ptr<const Variant> variant);

  size_t variantCount() const;

  std::shared_ptr<GenomeVariant> filterVariants(const VariantFilter& filter) const;

  const GenomeVariantMap& getMap() const { return genome_variant_map_; }

  bool isElement(const Variant& variant) const;
  std::shared_ptr<GenomeVariant> Union(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Intersection(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Difference(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;

  static std::shared_ptr<GenomeVariant> emptyGenomeVariant(const GenomeId_t& genome_id,
                                                           PhaseId_t ploidy,
                                                           const std::shared_ptr<const GenomeDatabase>& genome_db);

  std::string output(char field_delimiter, VariantOutputIndex output_index, bool detail) const;
  bool outputCSV(const std::string& file_name, VariantOutputIndex output_index, bool detail) const;

  void getVariants(std::vector<std::shared_ptr<const Variant>>& variant_vector) const;

  // All contig_id variants use the zero-based half-open convention [start, end).
  // End points past the last variant; end = (last + 1).
  bool getSortedVariants(ContigId_t contig_id,
                         PhaseId_t phase,
                         ContigOffset_t start,
                         ContigOffset_t end,
                         OffsetVariantMap& variant_map) const;


  // Returns reference and protein mutations.
  bool mutantProteins( const ContigId_t& contig_id,
                       PhaseId_t phase,
                       const FeatureIdent_t& gene_id,
                       const FeatureIdent_t& sequence_id,
                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                       OffsetVariantMap& variant_map,
                       std::shared_ptr<AminoSequence>& reference_sequence,
                       std::shared_ptr<AminoSequence>& sequence_vector) const;

  // Returns reference and mutant stranded DNA.
  bool mutantCodingDNA( const ContigId_t& contig_id,
                        PhaseId_t phase,
                        const FeatureIdent_t& gene_id,
                        const FeatureIdent_t& sequence_id,
                        const std::shared_ptr<const GenomeDatabase>& genome_db,
                        OffsetVariantMap& variant_map,
                        std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                        std::shared_ptr<DNA5SequenceCoding>& mutant_sequence) const;

  // Returns reference and mutant unstranded DNA region
  bool mutantRegion( const ContigId_t& contig_id,
                     PhaseId_t phase,
                     ContigOffset_t region_offset,
                     ContigSize_t region_size,
                     const std::shared_ptr<const GenomeDatabase>& genome_db,
                     OffsetVariantMap& variant_map,
                     std::shared_ptr<DNA5SequenceLinear>& reference_sequence,
                     std::shared_ptr<DNA5SequenceLinear>& mutant_sequence) const;

  // Returns reference and mutant unstranded contig.
  bool mutantContig( const ContigId_t& contig_id,
                     PhaseId_t phase,
                     std::shared_ptr<const GenomeDatabase> genome_db,
                     std::shared_ptr<const DNA5SequenceContig>& reference_contig_ptr,
                     std::shared_ptr<const DNA5SequenceContig>& mutant_contig_ptr) const;


  static constexpr PhaseId_t HAPLOID_GENOME = 1;
  static constexpr PhaseId_t DIPLOID_GENOME = 2;

private:

  GenomeId_t genome_id_;
  PhaseId_t ploidy_;
  GenomeVariantMap genome_variant_map_;

};


}   // namespace genome
}   // namespace kellerberrin


// Not in kgl:: namespace.
std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::GenomeVariant& genome_variant);
std::ostream & operator<<(std::ostream &os, std::shared_ptr<const kellerberrin::genome::GenomeVariant> genome_variant_ptr);
std::ostream & operator<<(std::ostream &os, std::shared_ptr<kellerberrin::genome::GenomeVariant> genome_variant_ptr);


#endif //KGL_VARIANT_DB_GENOME_H
