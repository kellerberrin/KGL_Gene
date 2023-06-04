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

struct CompareInterval {

  bool operator()(const OpenRightInterval& lhs, const OpenRightInterval& rhs) const
  {
    return lhs.lower() < rhs.lower();
  }

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Adapters for std::set and std::map to use OpenRightInterval.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class IntervalSet : public std::set<OpenRightInterval, CompareInterval> {

public:

  IntervalSet() = default;
  ~IntervalSet() = default;

  [[nodiscard]] bool containsInterval(const OpenRightInterval& interval) const;
  [[nodiscard]] bool containsOffset(ContigOffset_t offset) const { return contains(OpenRightInterval(offset, offset+1)); }

};


template<typename ValueType>
class IntervalMap : public std::map<OpenRightInterval, ValueType, CompareInterval> {

public:

  IntervalMap() = default;
  ~IntervalMap() = default;

  [[nodiscard]] auto findInterval(const OpenRightInterval& interval) const {

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
  [[nodiscard]] bool containsInterval(const OpenRightInterval& interval) const { return findInterval(interval) != this->end(); }
  [[nodiscard]] bool containsOffset(ContigOffset_t offset) const { return containsInterval(OpenRightInterval(offset, offset+1)); }

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

  [[nodiscard]] const std::shared_ptr<const GeneFeature>& getGene() const { return gene_feature_; }
  [[nodiscard]] const GeneCodingTranscriptMap& codingTranscripts() const { return gene_coding_transcripts_; }
  [[nodiscard]] const OpenRightInterval& geneInterval() const { return gene_interval_; }

  [[nodiscard]] bool isMemberGene(ContigOffset_t offset) const { return gene_interval_.containsOffset(offset); }
  [[nodiscard]] bool isMemberCoding(ContigOffset_t offset) const;

private:

  std::shared_ptr<const GeneFeature> gene_feature_;
  GeneCodingTranscriptMap gene_coding_transcripts_;
  OpenRightInterval gene_interval_{0, 0};

  void codingInterval(const std::shared_ptr<const GeneFeature>& gene_vector);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ContigIntervalMap = std::map<ContigId_t, IntervalMap<GeneIntervalStructure>>;
class IntervalCodingVariants {

public:

  explicit IntervalCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector);
  ~IntervalCodingVariants() = default;

  // Returns true if the variant is within a gene coding region.
  [[nodiscard]] bool codingRegionVariant(const Variant &variant) const;
  // Returns std::nullopt if the variant is not within a gene coding region.
  [[nodiscard]] std::optional<std::shared_ptr<const GeneFeature>> getGeneCoding(const Variant &variant) const;
  // Returns std::nullopt if the variant is not within a gene interval (5 prime, intron, coding, 3 prime).
  [[nodiscard]] std::optional<std::shared_ptr<const GeneFeature>> getGeneInterval(const Variant &variant) const;

private:

  ContigIntervalMap contig_interval_map_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FeatureInterval {

public:

  explicit FeatureInterval(const std::vector<std::shared_ptr<const Feature>>& feature_vector);
  ~FeatureInterval() = default;

  [[nodiscard]] bool variantContained(const Variant &variant) const;

private:

  ContigIntervalMap contig_interval_map_;

};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A gene interval node.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class IntervalAllCodingVariants {

public:

  explicit IntervalAllCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr);
  ~IntervalAllCodingVariants();

  IntervalAllCodingVariants(const IntervalAllCodingVariants&);

  [[nodiscard]] bool contains(const Variant& variant) const;

private:


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter variants that are (possibly partially if an indel) within the coding sequence of gene.
// The vector specifies which genes are filtered.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
class IntervalCodingVariants {

public:

  explicit IntervalCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector);
  ~IntervalCodingVariants();

  IntervalCodingVariants(const IntervalCodingVariants&);

  [[nodiscard]] bool contains(const Variant& variant) const;

private:


};
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter a vector of arbitrary features.
// If, for example, a feature is a gene, this can (and mostly does) include 5 prime, 3 prime and intron non-coding regions.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ImplementFeatureInterval; // Pimpl implementation for arbitrary feature intervals only.

class IntervalFeatures {

public:

  explicit IntervalFeatures(const std::vector<std::shared_ptr<const Feature>>& feature_vector);
  ~IntervalFeatures();

  IntervalFeatures(const IntervalFeatures&);

  [[nodiscard]] bool contains(const Variant& variant) const;

private:

  std::shared_ptr<const ImplementFeatureInterval> pimpl_feature_interval_ptr_;

};



} // Namespace



#endif //KGL_GENOME_INTERVAL_H
