//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_GENOME_DB_H
#define KGL_GENOME_DB_H


#include <memory>
#include <string>
#include <vector>
#include <map>
#include "kgl_genome_types.h"
#include "kgl_genome_feature.h"
#include "kgl_mt_data.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures - A contiguous region, the associated sequence,  and all features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetFeatureMap = std::multimap<ContigOffset_t, std::shared_ptr<Feature>>;
using IdFeatureMap = std::multimap<FeatureIdent_t, std::shared_ptr<Feature>>;
using CDSTable = std::vector<std::shared_ptr<CDSFeature>>;

class ContigFeatures {

public:

  ContigFeatures(ContigId_t contig_id, Sequence_t sequence) : contig_id_(std::move(contig_id)),
                                                            sequence_(std::move(sequence)) {}
  ContigFeatures(const ContigFeatures&) = default;
  ~ContigFeatures() = default;

  ContigFeatures& operator=(const ContigFeatures&) = default;

  bool addFeature(std::shared_ptr<Feature>& feature_ptr);
  // false if not found.
  bool findFeatureId(FeatureIdent_t& feature_id, std::vector<std::shared_ptr<Feature>>& feature_ptr_vec);
  // false if offset is not in an cds else returns a vector of cds (these will be in different genes).
  bool findOffsetCDS(ContigOffset_t offset, std::vector<std::shared_ptr<CDSFeature>>& cds_ptr_vec);

  const ContigId_t& contigId() const { return contig_id_; }
  const Sequence_t& sequence() const { return sequence_; }
  ContigSize_t contigSize() { return sequence_.length(); }

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();
  void verifyCDSPhasePeptide();

private:

  ContigId_t contig_id_;
  Sequence_t sequence_;
  OffsetFeatureMap offset_feature_map_;
  IdFeatureMap id_feature_map_;
  CDSTable cds_table_;

  void verifyContigOverlap();
  void verifySubFeatureSuperFeatureDimensions();
  void verifySubFeatureDuplicates();
  void verifySuperFeatureDuplicates();
  void removeSubFeatureDuplicates();
  void removeSuperFeatureDuplicates();
  void createCDSTable();

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase - A map of contigs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenomeSequenceMap = std::map<ContigId_t, std::shared_ptr<ContigFeatures>>;
class GenomeDatabase {

public:

  explicit GenomeDatabase() {}
  GenomeDatabase(const GenomeDatabase&) = default;
  ~GenomeDatabase() = default;

  GenomeDatabase& operator=(const GenomeDatabase&) = default;

  // Return false if contig already exists.
  bool addContigSequence(const ContigId_t& contig, Sequence_t sequence);
  // Returns false if key not found.
  bool getContigSequence(const ContigId_t& contig, std::shared_ptr<ContigFeatures>& contig_ptr) const;
  GenomeSequenceMap& getGenomeSequenceMap() { return genome_sequence_map_; }

  void createVerifyGenomeDatabase();

  void registerContigData(std::shared_ptr<ContigCountData>& contig_data_ptr) const;

private:

  GenomeSequenceMap genome_sequence_map_;

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();
  void verifyCDSPhasePeptide();

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_DB_H
