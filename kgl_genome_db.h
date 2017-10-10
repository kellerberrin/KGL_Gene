

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


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FeatureAttributes Object to hold { key=value } pairs.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using AttributeMap = std::multimap<const std::string, std::string>;
class FeatureAttributes {

public:

  explicit FeatureAttributes() = default;
  FeatureAttributes(const FeatureAttributes&) = default;
  ~FeatureAttributes() = default;

  FeatureAttributes& operator=(const FeatureAttributes&) = default;

  // General access routines.
  bool getAttributes(const std::string& key, std::vector<std::string>& value_vec) const; // false if no key.
  void insertAttribute(const std::string& key, const std::string& value); // Always succeeds; keys are uppercase.
  void getAllAttributes(std::vector<std::pair<std::string, std::string>>& all_key_value_pairs) const;

  // Attribute keys.
  constexpr static const char* ID_KEY = "ID";
  constexpr static const char* SUPER_FEATURE_KEY = "PARENT";

  // Convenience access routines.
  bool getIds(std::vector<std::string> &value_vec) const { return getAttributes(ID_KEY, value_vec); }
  bool getSuperFeatureIds(std::vector<std::string> &value_vec) const { return getAttributes(SUPER_FEATURE_KEY, value_vec); }


private:

  AttributeMap attributes_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FeatureSequence - Feature location and sense (if applicable)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class StrandSense { FORWARD = '+', REVERSE = '-', UNKNOWN = '.'};
class FeatureSequence {

public:

  FeatureSequence(const ContigOffset_t& begin_offset,
                  const ContigOffset_t& end_offset,
                  const StrandSense& strand_sense)
      : begin_offset_(begin_offset),
        end_offset_(end_offset),
        strand_sense_(strand_sense) {}
  ~FeatureSequence() = default;
  FeatureSequence(const FeatureSequence&) = default;

  FeatureSequence& operator=(const FeatureSequence&) = default;

  ContigOffset_t begin() const { return begin_offset_; }
  ContigOffset_t end() const { return end_offset_; }
  StrandSense sense() const { return strand_sense_; }
  void begin(ContigOffset_t begin) { begin_offset_ = begin; }
  void end(ContigOffset_t end) { end_offset_ = end; }
  void sense(StrandSense strand) { strand_sense_ = strand; }


private:

  ContigOffset_t begin_offset_;    // [O, size-1] based contig sequence offset.
  ContigOffset_t end_offset_;      // [O, size-1] based contig sequence offset.
  StrandSense strand_sense_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FeatureRecord - Annotated contig features
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ContigRecord; // Forward decl;
class FeatureRecord; // Forward decl.
using FeatureIdent_t = std::string;
using FeatureType_t = std::string;
using SubFeatureMap = std::multimap<const FeatureIdent_t, std::shared_ptr<FeatureRecord>>;
using SuperFeatureMap = std::multimap<const FeatureIdent_t, std::shared_ptr<FeatureRecord>>;
class FeatureRecord {

public:

  FeatureRecord(const FeatureIdent_t& id,
                const FeatureType_t& type,
                const std::shared_ptr<ContigRecord>& contig_ptr,
                const FeatureSequence& sequence): id_(id), type_(type), contig_ptr_(contig_ptr), sequence_(sequence) {}
  FeatureRecord(const FeatureRecord&) = default;
  virtual ~FeatureRecord() = default;

  FeatureRecord& operator=(const FeatureRecord&) = default;

  const FeatureIdent_t& id() const { return id_; }
  const FeatureSequence& sequence() const { return sequence_; }
  void sequence(const FeatureSequence& new_sequence) { sequence_ = new_sequence; }
  void setAttributes(const FeatureAttributes& attributes) { attributes_ = attributes; }
  const FeatureAttributes& getAttributes() const { return attributes_; }

  // Hierarchy routines.
  void clearHierachy() { sub_features_.clear(); super_features_.clear(); }
  void addSuperFeature(const FeatureIdent_t &super_feature_id, const std::shared_ptr<FeatureRecord> &super_feature_ptr);
  void addSubFeature(const FeatureIdent_t& sub_feature_id, const std::shared_ptr<FeatureRecord>& sub_feature_ptr);

  SuperFeatureMap& superFeatures() { return super_features_; }
  SubFeatureMap& subFeatures() { return sub_features_; }

private:

  FeatureIdent_t id_;
  FeatureType_t type_;
  std::shared_ptr<ContigRecord> contig_ptr_;
  FeatureSequence sequence_;
  SubFeatureMap sub_features_;
  SuperFeatureMap super_features_;
  FeatureAttributes attributes_;

};


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
// GenomeSequences - A map of contigs and associated features.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenomeSequenceMap = std::map<ContigId_t, std::shared_ptr<ContigRecord>>;
class GenomeSequences {

public:

  explicit GenomeSequences() {}
  GenomeSequences(const GenomeSequences&) = default;
  ~GenomeSequences() = default;

  GenomeSequences& operator=(const GenomeSequences&) = default;

  // Return false if contig already exists.
  bool addContigSequence(const ContigId_t& contig, Sequence_t sequence);
  // Returns false if key not found.
  bool getContigSequence(const ContigId_t& contig, std::shared_ptr<ContigRecord>& contig_ptr) const;
  GenomeSequenceMap& getGenomeSequenceMap() { return genome_sequence_map_; }

  void createVerifyGenomeDatabase();


private:

  GenomeSequenceMap genome_sequence_map_;

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_DB_H
