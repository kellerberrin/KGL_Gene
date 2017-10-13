// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
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
// ContigRecord - A contiguous region, the associated sequence,  and all features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetFeatureMap = std::multimap<ContigOffset_t, std::shared_ptr<FeatureRecord>>;
using IdFeatureMap = std::multimap<FeatureIdent_t, std::shared_ptr<FeatureRecord>>;
class ContigRecord {

public:

  ContigRecord(ContigId_t contig_id, Sequence_t sequence) : contig_id_(std::move(contig_id)),
                                                            sequence_(std::move(sequence)) {}
  ContigRecord(const ContigRecord&) = default;
  ~ContigRecord() = default;

  ContigRecord& operator=(const ContigRecord&) = default;

  bool addFeature(std::shared_ptr<FeatureRecord>& feature_ptr);
  // false if not found.
  bool findFeatureId(FeatureIdent_t& feature_id, std::vector<std::shared_ptr<FeatureRecord>>& feature_ptr_vec);

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

  void verifyContigOverlap();
  void verifySubFeatureSuperFeatureDimensions();
  void verifySubFeatureDuplicates();
  void verifySuperFeatureDuplicates();
  void removeSubFeatureDuplicates();
  void removeSuperFeatureDuplicates();

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase - A map of contigs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenomeSequenceMap = std::map<ContigId_t, std::shared_ptr<ContigRecord>>;
class GenomeDatabase {

public:

  explicit GenomeDatabase() {}
  GenomeDatabase(const GenomeDatabase&) = default;
  ~GenomeDatabase() = default;

  GenomeDatabase& operator=(const GenomeDatabase&) = default;

  // Return false if contig already exists.
  bool addContigSequence(const ContigId_t& contig, Sequence_t sequence);
  // Returns false if key not found.
  bool getContigSequence(const ContigId_t& contig, std::shared_ptr<ContigRecord>& contig_ptr) const;
  GenomeSequenceMap& getGenomeSequenceMap() { return genome_sequence_map_; }

  void createVerifyGenomeDatabase();

  void registerContigData(std::shared_ptr<ContigCountData>& contig_data_ptr);

private:

  GenomeSequenceMap genome_sequence_map_;

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();
  void verifyCDSPhasePeptide();

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_DB_H
