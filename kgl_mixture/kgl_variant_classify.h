//
// Created by kellerberrin on 15/06/18.
//

#ifndef KGL_VARIANT_CLASSIFY_H
#define KGL_VARIANT_CLASSIFY_H


#include <kgd_ctl_data.h>
#include "kel_utility.h"
#include "kgl_variant_db.h"
#include "kgd_ctl_data.h"

namespace kellerberrin::genome {   //  organization::project level namespace

namespace kgd = kellerberrin::deconvolv;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects accept unphased variants and produce out ref/alt read statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigOffsetVariants {

public:

  explicit ContigOffsetVariants(std::shared_ptr<const Variant> variant) : variant_(variant) {}
  ContigOffsetVariants(const ContigOffsetVariants& copy) = default;
  ~ContigOffsetVariants() = default;

  std::shared_ptr<const Variant> variant() const { return variant_; }

  bool operator<(const ContigOffsetVariants& cmp) const { return variant_->lessThan(*cmp.variant()); }

private:

  std::shared_ptr<const Variant> variant_;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Orders by variants and then by genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Structures used to order the variant info into variant (row) by genome (column) ref/alt count info.
using GenomeOffsetMap = std::map<GenomeId_t, std::shared_ptr<const Variant>>;
using VariantMap = std::map<ContigOffsetVariants, GenomeOffsetMap>;

// Alt count ordered structure to hold the variant count info.
// This is used in testing the J. O'brien 'R' MCMC mixture code.
// pair.first is ref count, pair.second is alt count.
using OrderedAltCount = std::multimap<size_t, std::pair<std::vector<size_t>, std::vector<size_t>>>;

class VariantClassifier {

public:

  explicit VariantClassifier(std::shared_ptr<const PopulationDB> vcf_population_ptr);
  ~VariantClassifier() = default;

  const VariantMap& getMap() const { return variant_map_; }
  const std::vector<GenomeId_t>& getGenomes() const { return genomes_; }

  // Produce a kgd data object.
  kgd::MixtureDataObj convertToMixture(const GenomeId_t& genome_count, size_t min_count, size_t max_count) const;

  // IO routines - to produce files for testing.
  // Write ref/alt counts into separate files by variant order.
  bool writeVariants(char delimiter,
                     const std::string& ref_file_name,
                     const std::string& alt_file_name,
                     size_t min_count) const;

// Write population allele frequency.
  bool writePlaf(char delimiter,
                 const std::string& plaf_file_name,
                 const std::string& ref_file_prefix,
                 const std::string& alt_file_prefix,
                 size_t min_count) const;

// Write ref/alt counts simultaneously into separate files by total genome alt count order.
// The maximum number of genomes with non-zero alt counts.
// This is used in initial 'R' testing of the the J. O'brien MCMC mixture algorithm
  bool writeOrderedVariants(char delimiter,
                            const std::string& ref_file_name,
                            const std::string& alt_file_name,
                            size_t min_count,
                            size_t max_lines) const;

private:

  VariantMap variant_map_;
  std::vector<GenomeId_t> genomes_;

  bool orderVariantCount(OrderedAltCount& ordered_count, size_t min_count) const;
  bool writeOrderedCount(OrderedAltCount& ordered_count,
                         char delimiter,
                         const std::string& ref_file_name,
                         const std::string& alt_file_name,
                         size_t max_lines) const;

};


}   // end::namespace


#endif //KGL_VARIANT_CLASSIFY_H
