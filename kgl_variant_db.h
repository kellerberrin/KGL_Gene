//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_VARIANT_DB_H
#define KGL_VARIANT_DB_H


#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<const Variant>>;
class ContigVariant {

public:

  explicit ContigVariant(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  ContigVariant(const ContigVariant&) = delete; // Use deep copy.
  ~ContigVariant() = default;

  ContigVariant& operator=(const ContigVariant&) = delete; // Use deep copy.

  // Always use deep copy when modifying this object (filter and set operations).
  std::shared_ptr<ContigVariant> deepCopy() const;

  bool addVariant(std::shared_ptr<const Variant>& variant_ptr);

  const ContigId_t& contigId() const { return contig_id_; }
  size_t variantCount() const { return offset_variant_map_.size(); }

  // Set functions.
  bool isElement(const Variant& variant) const;
  std::shared_ptr<ContigVariant> Union(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Intersection(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Difference(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;

  std::shared_ptr<ContigVariant> filterVariants(const VariantFilter& filter) const;

  const OffsetVariantMap& getMap() const { return offset_variant_map_; }

  size_t size() const { return offset_variant_map_.size(); }

  // All variants in [start, end) - note that end points past the last variant; end = (last + 1).
  bool getSortedVariants(ContigOffset_t start,
                         ContigOffset_t end,
                         OffsetVariantMap& variant_map) const;


private:

  ContigId_t contig_id_;
  OffsetVariantMap offset_variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The delete offset accounting map records where deletions occur so that inserts and deletes can be properly aligned.
using IndelAccountingMap = std::map<ContigOffset_t, SignedOffset_t>;

// The variant contig map
using GenomeVariantMap = std::map<ContigId_t, std::shared_ptr<ContigVariant>>;

class GenomeVariant {

public:

  explicit GenomeVariant(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  GenomeVariant(const GenomeVariant&) = delete; // Use deep copy.
  ~GenomeVariant() = default;

  GenomeVariant& operator=(const GenomeVariant&) = delete; // Use deep copy.

  // Always use deep copy when modifying this object (filter and set operations).
  std::shared_ptr<GenomeVariant> deepCopy() const;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& contig_id) { genome_id_ = contig_id; }

  size_t contigCount() const { return genome_variant_map_.size(); }

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  bool getContigVariant(const ContigId_t& contig_id, std::shared_ptr<ContigVariant>& contig_variant) const;

  bool addVariant(std::shared_ptr<const Variant> variant);

  size_t size() const;

  std::shared_ptr<GenomeVariant> filterVariants(const VariantFilter& filter) const;

  const GenomeVariantMap& getMap() const { return genome_variant_map_; }

  bool isElement(const Variant& variant) const;
  std::shared_ptr<GenomeVariant> Union(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Intersection(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Difference(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;

  static std::shared_ptr<GenomeVariant> emptyGenomeVariant(const GenomeId_t& genome_id,
                                                           const std::shared_ptr<const GenomeDatabase>& genome_db);

  std::string output(char field_delimiter, VariantOutputIndex output_index, bool detail) const;
  bool outputCSV(const std::string& file_name, VariantOutputIndex output_index, bool detail) const;

  void getVariants(std::vector<std::shared_ptr<const Variant>>& variant_vector) const;

  // All contig_id variants use the zero-based half-open convention [start, end).
  // End points past the last variant; end = (last + 1).
  bool getSortedVariants(ContigId_t contig_id,
                         ContigOffset_t start,
                         ContigOffset_t end,
                         OffsetVariantMap& variant_map) const;

  // Convenience routine only returns coding variants.
  bool getCodingSortedVariants(ContigId_t contig_id,
                               ContigOffset_t start,
                               ContigOffset_t end,
                               OffsetVariantMap& variant_map) const;

  const Attributes& attributes() const { return attributes_; }
  Attributes& attributes() { return attributes_; }
  void attributes(const Attributes& attributes) { attributes_ = attributes; }

  // Returns a maximum of MUTATION_HARD_LIMIT_ alternative protein mutations.
  bool mutantProteins( const ContigId_t& contig_id,
                       const FeatureIdent_t& gene_id,
                       const FeatureIdent_t& sequence_id,
                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                       std::shared_ptr<AminoSequence>& reference_sequence,
                       std::vector<std::shared_ptr<AminoSequence>>& mutant_sequence_vector) const;

  // Returns a maximum of MUTATION_HARD_LIMIT_ alternative protein mutations.
  bool mutantCodingDNA( const ContigId_t& contig_id,
                        const FeatureIdent_t& gene_id,
                        const FeatureIdent_t& sequence_id,
                        const std::shared_ptr<const GenomeDatabase>& genome_db,
                        std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                        std::vector<std::shared_ptr<DNA5SequenceCoding>>& mutant_sequence_vector) const;

  // Returns a maximum of MUTATION_HARD_LIMIT_ alternative mutations.
  bool mutantRegion( const ContigId_t& contig_id,
                     const ContigOffset_t & region_offset,
                     const ContigSize_t region_size,
                     const std::shared_ptr<const GenomeDatabase>& genome_db,
                     std::shared_ptr<DNA5SequenceLinear>& reference_sequence,
                     std::vector<std::shared_ptr<DNA5SequenceLinear>>& mutant_sequence_vector) const;



private:

  GenomeId_t genome_id_;
  GenomeVariantMap genome_variant_map_;
  Attributes attributes_;

  constexpr static size_t MUTATION_SOFT_LIMIT_ = 32;
  constexpr static size_t MUTATION_HARD_LIMIT_ = 128;

  // The coding variants in the variant_map are used to mutate the dna_sequence.
  static bool mutateDNA(const OffsetVariantMap& variant_map,
                        std::shared_ptr<const ContigFeatures> contig_ptr,
                        std::shared_ptr<const CodingSequence> coding_sequence_ptr,
                        std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr);

  // Returns a vector of alternative mutation paths based on the number of equal offset mutations in the coding variants.
  // There maybe more than one variant specified per offset.
  // If there are equal offset variants then we create alternative mutation paths.
  // This function is exponential. Alternative Mutations = 2 ^ (#equal offset variants).
  // A warning is issued if the soft limit is reached; default 32 alternatives (5 equal offset variants).
  // The number of variant paths is capped by the hard limit; default 128 alternatives (9 equal offset variants)
  static void getMutationAlternatives(std::shared_ptr<const OffsetVariantMap> variant_map_ptr,
                                      std::vector<OffsetVariantMap>& variant_map_vector,
                                      size_t& alternative_count,
                                      size_t soft_limit,
                                      size_t hard_limit);

  // Split the variant map into SNP, Delete and Insert Variants.
  static void SplitVariantMap(const OffsetVariantMap& variant_map,
                              OffsetVariantMap& snp_variant_map,
                              OffsetVariantMap& delete_variant_map,
                              OffsetVariantMap& insert_variant_map);

  static bool mutateDNA(const OffsetVariantMap& insert_variant_map,
                        ContigOffset_t contig_offset,
                        std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr,
                        IndelAccountingMap& indel_accounting_map);

// Mutate the DNA sequence using SNP variants
  static bool mutateSNPs(const OffsetVariantMap& snp_variant_map,
                         ContigOffset_t contig_offset,
                         std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr);

// Mutate the DNA sequence using a single SNP.
  static bool mutateSingleSNP(std::shared_ptr<const Variant> variant_ptr,
                              ContigOffset_t contig_offset,
                              std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr);

// Mutate the DNA sequence using Delete variants
  static bool mutateDeletes(const OffsetVariantMap& snp_variant_map,
                            ContigOffset_t contig_offset,
                            std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr,
                            IndelAccountingMap& indel_accounting_map);

// Mutate the DNA sequence using Insert variants
  static bool mutateInserts(const OffsetVariantMap& snp_variant_map,
                            ContigOffset_t contig_offset,
                            std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr,
                            IndelAccountingMap& indel_accounting_map);


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container to hold genome variants for populations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using PopulationVariantMap = std::map<GenomeId_t, std::shared_ptr<const GenomeVariant>>;
class PopulationVariant {

public:

  explicit PopulationVariant(const std::string& population_id) : population_id_(population_id) {}
  PopulationVariant(const PopulationVariant&) = default;
  ~PopulationVariant() = default;

  bool addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant);

  bool getGenomeVariant(const GenomeId_t& genome_id, std::shared_ptr<const GenomeVariant>& genome_variant);

  const PopulationVariantMap& getMap() const { return population_variant_map_; }

private:

  PopulationVariantMap population_variant_map_;
  std::string population_id_;

};


}   // namespace genome
}   // namespace kellerberrin


// Not in kgl:: namespace.
std::ostream & operator<<(std::ostream &os, const kellerberrin::genome::GenomeVariant& genome_variant);


#endif //KGL_VARIANT_DB_H
