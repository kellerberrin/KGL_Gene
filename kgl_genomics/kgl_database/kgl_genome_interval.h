//
// Created by kellerberrin on 29/05/23.
//

#ifndef KGL_GENOME_INTERVAL_H
#define KGL_GENOME_INTERVAL_H


#include "kgl_genome_genome.h"
#include "kgl_variant.h"

#include <set>
#include <map>


namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Defines a simple right open interval [lower_, upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class OpenRightInterval {

public:

  OpenRightInterval(ContigOffset_t lower, ContigOffset_t upper);
  ~OpenRightInterval() = default;
  OpenRightInterval(const OpenRightInterval& copy) { lower_ = copy.lower_; upper_ = copy.upper_; }

  OpenRightInterval& operator=(const OpenRightInterval& copy) { lower_ = copy.lower_; upper_ = copy.upper_; return *this; }

  [[nodiscard]] ContigOffset_t lower() const { return lower_; }
  [[nodiscard]] ContigOffset_t upper() const { return upper_; }

  [[nodiscard]] size_t size() const { return upper_ - lower_; }
  [[nodiscard]] bool containsOffset(ContigOffset_t offset) const { return offset >= lower_ and offset < upper_; }
  [[nodiscard]] bool containsInterval(const OpenRightInterval& interval) const { return interval.lower_ >= lower_ and (interval.lower_ + interval.size()) <= upper_; }

private:

  ContigOffset_t upper_{0};
  ContigOffset_t lower_{0};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Adapters for std::set and std::map to use OpenRightInterval.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Comparison operator used to order intervals within indexed containers std::set or std::map.
struct CompareInterval {

  bool operator()(const OpenRightInterval& lhs, const OpenRightInterval& rhs) const
  {
    return lhs.lower() < rhs.lower();
  }

};

// Interval adapter for std::set.
class IntervalSet : public std::set<OpenRightInterval, CompareInterval> {

public:

  IntervalSet() = default;
  ~IntervalSet() = default;

  // Find the interval that contains the argument OR the interval immediately greater (lower > arg.lower) than the argument interval.
  [[nodiscard]] auto findUpperEqualIter(const OpenRightInterval& interval) const {

    auto iter = this->lower_bound(interval);
    // Return the previous map interval if it contains the argument interval.
    auto prev_iter = std::prev(iter, 1);
    if (prev_iter != this->end()) {

      auto const& interval_key = *prev_iter;
      if (interval_key.containsInterval(interval)) {

        return prev_iter;

      }

    }

    // Else just return the upper interval
    return iter;

  }

  [[nodiscard]] bool containsInterval(const OpenRightInterval& interval) const;
  [[nodiscard]] bool containsOffset(ContigOffset_t offset) const { return contains(OpenRightInterval(offset, offset+1)); }

};

// Interval adapter for std::map.
template<typename ValueType>
using IntervalMapType = std::map<OpenRightInterval, ValueType, CompareInterval>;
template<typename ValueType>
class IntervalMap : public IntervalMapType<ValueType> {

public:

  IntervalMap() = default;
  ~IntervalMap() = default;

  // Returns any the map interval that contains the argument or end().
  [[nodiscard]] auto findIntervalIter(const OpenRightInterval& interval) const {

    auto iter = this->lower_bound(interval);
    if (iter != this->end()) {

      // Do the lower bound of the intervals match.
      auto const& [interval_key, value] = *iter;
      if (interval_key.lower() == interval.lower()) {

        if (interval_key.containsInterval(interval)) {

          return iter;

        } else {

          return this->end();

        }

      }

    }

    // Look at the previous interval.
    iter = std::prev(iter, 1);
    if (iter != this->end()) {

      auto const& [interval_key, value] = *iter;
      if (interval_key.containsInterval(interval)) {

        return iter;

      } else {

        return this->end();

      }

    }

    return this->end();

  }

  // Find the interval that contains the argument OR the interval immediately greater (lower > arg.lower) than the argument interval.
  [[nodiscard]] auto findUpperEqualIter(const OpenRightInterval& interval) const {

    auto iter = this->lower_bound(interval);
    // Return the previous map interval if it contains the argument interval.
    auto prev_iter = std::prev(iter, 1);
    if (prev_iter != this->end()) {

      auto const& [interval_key, value] = *prev_iter;
      if (interval_key.containsInterval(interval)) {

        return prev_iter;

      }

    }

    // Else just return the upper interval
    return iter;

  }

  [[nodiscard]] std::optional<ValueType> findUpperInterval(const OpenRightInterval& interval) const {

    auto iter = findUpperEqualIter(interval);
    if (iter != this->end()) {

      auto const& [interval_key, value] = *iter;
      return value;

    }

    return std::nullopt;

  }

  [[nodiscard]] auto findUpperOffsetIter(ContigOffset_t offset) const { return findUpperEqualIter(OpenRightInterval(offset, offset + 1)); }
  [[nodiscard]] std::optional<ValueType> findUpperOffset(ContigOffset_t offset) const { return findUpperInterval(OpenRightInterval(offset, offset + 1)); }

  [[nodiscard]] bool containsInterval(const OpenRightInterval& interval) const { return findIntervalIter(interval) != this->end(); }
  [[nodiscard]] bool containsOffset(ContigOffset_t offset) const { return containsInterval(OpenRightInterval(offset, offset + 1)); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A gene interval node.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// A map of named gene coding transcripts.
using GeneCodingTranscriptMap = std::map<std::string, IntervalSet>;
class GeneIntervalStructure {

public:

  explicit GeneIntervalStructure(const std::shared_ptr<const GeneFeature>& gene_feature);
  ~GeneIntervalStructure() = default;

  GeneIntervalStructure(const GeneIntervalStructure& copy) {

    gene_feature_ = copy.getGene();
    gene_coding_transcripts_ = copy.codingTranscripts();
    gene_interval_ = copy.geneInterval();

  }

  // Object access
  [[nodiscard]] const std::shared_ptr<const GeneFeature>& getGene() const { return gene_feature_; }
  [[nodiscard]] const GeneCodingTranscriptMap& codingTranscripts() const { return gene_coding_transcripts_; }
  [[nodiscard]] const OpenRightInterval& geneInterval() const { return gene_interval_; }

  // Given an offset, does the offset fall within a defined gene interval (can include 5 prime, coding intervals, introns, 3 prime).
  [[nodiscard]] bool isMemberGene(ContigOffset_t offset) const { return gene_interval_.containsOffset(offset); }
  // Given an offset, does the offset fall within a gene interval coding region.
  [[nodiscard]] bool isMemberCoding(ContigOffset_t offset) const;
  // Used with the functions above to determine if a contig + offset resides within a gene interval or the coding intervals of a gene.
  [[nodiscard]] bool isSameContig(const ContigId_t& contig) const { return contig == (gene_feature_->contig()->contigId()); }
  // Test if a variant modifies any of the coding transcripts of this gene structure.
  [[nodiscard]] bool codingModifier(const Variant& variant) const;

private:

  std::shared_ptr<const GeneFeature> gene_feature_;
  GeneCodingTranscriptMap gene_coding_transcripts_;
  OpenRightInterval gene_interval_{0, 0};

  void codingInterval(const std::shared_ptr<const GeneFeature>& gene_vector);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// GeneIntervalStructure objects are stored in an IntervalMap by their gene intervals.,
// The IntervalMaps are further indexed by contig in the ContigIntervalMap container.
// This enables the entire gene set in a genome to be stored as GeneIntervalStructure objects.
// Thus, we can test if a variant is within a gene coding region, and if so, return the gene feature that it belongs to.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ContigIntervalMap = std::map<ContigId_t, IntervalMap<GeneIntervalStructure>>;
class IntervalCodingVariants {

public:

  explicit IntervalCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector) { InitializeGeneVector(gene_vector); }
  explicit IntervalCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr);
  ~IntervalCodingVariants() = default;

  // Returns true if the variant is within a gene coding region.
  [[nodiscard]] bool codingRegionVariant(const Variant& variant) const;
  // Returns std::nullopt if the variant is not within a gene coding region.
  [[nodiscard]] std::optional<std::shared_ptr<const GeneFeature>> getGeneCoding(const Variant &variant) const;

private:

  ContigIntervalMap contig_interval_map_;

  void InitializeGeneVector(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector);

};



} // Namespace



#endif //KGL_GENOME_INTERVAL_H
