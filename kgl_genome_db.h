

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
#include <map>
#include "kgl_genome_types.h"
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Object to hold { index : item } pairs.
using AttributeMap = std::map<const std::string, const std::string>;
class GenomeAttributes {

public:

  explicit GenomeAttributes() = default;
  GenomeAttributes(const GenomeAttributes&) = default;
  ~GenomeAttributes() = default;

  GenomeAttributes& operator=(const GenomeAttributes&) = default;

  bool getAttributes(const std::string& key, std::string& value);   // Returns bool false if the key not found.
  void getAllAttributes(std::vector<std::pair<std::string, std::string>>& key_value_pairs);
  void insertAttribute(const std::string& key, const std::string& value); // Always succeeds - overwrites any existing.

private:

  AttributeMap attributes_;

};



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

private:

  ContigOffset_t begin_offset_;    // O based contig sequence offset.
  ContigOffset_t end_offset_;      // 0 based contig sequence offset.
  StrandSense strand_sense_;

};

class ContigRecord; // Forward decl;
class SubFeatureRecord; // Forward decl.
using SubFeaturePtrMap = std::map<ContigOffset_t, std::unique_ptr<SubFeatureRecord>>;
enum class FeatureType { GENE, EXON, CDS, INTRON };
class FeatureRecord {

public:

  FeatureRecord(ContigRecord& contig,
                const FeatureSequence& sequence): contig_(contig), sequence_(sequence) {}
  FeatureRecord(const FeatureRecord&) = default;
  virtual ~FeatureRecord() = default;

  FeatureRecord& operator=(const FeatureRecord&) = default;

  bool isSubFeature() { return (bool) getParent(); }
  virtual std::shared_ptr<FeatureRecord> getParent() { std::shared_ptr<FeatureRecord> null; return null; }
  virtual FeatureType featureType() { return FeatureType::GENE; }

private:

  ContigRecord& contig_;
  FeatureSequence sequence_;
  SubFeaturePtrMap sub_features_;
  GenomeAttributes attributes_;

};


class SubFeatureRecord: public FeatureRecord {

public:

  SubFeatureRecord(ContigRecord& contig,
                   const FeatureSequence& sequence,
                   std::shared_ptr<FeatureRecord> parent,
                   FeatureType feature_type) : FeatureRecord(contig, sequence),
                                               parent_(std::move(parent)),
                                               feature_type_(feature_type){}
  ~SubFeatureRecord() override = default;
  SubFeatureRecord(const SubFeatureRecord&) = default;

  std::shared_ptr<FeatureRecord> getParent() final { return parent_; }
  FeatureType featureType() final { return feature_type_; }

private:

  std::shared_ptr<FeatureRecord> parent_;
  FeatureType feature_type_;

};


using FeatureMap = std::map<ContigOffset_t, std::unique_ptr<FeatureRecord>>;
class ContigRecord {

public:

  ContigRecord(ContigId_t contig_id, Sequence_t sequence) : contig_id_(std::move(contig_id)),
                                                            sequence_(std::move(sequence)) {}
  ContigRecord(const ContigRecord&) = default;
  ~ContigRecord() = default;

  ContigRecord& operator=(const ContigRecord&) = default;

  const ContigId_t& contigID() const { return contig_id_; }
  const Sequence_t& sequence() const { return sequence_; }

private:

  ContigId_t contig_id_;
  Sequence_t sequence_;
  FeatureMap feature_map_;
  GenomeAttributes attributes_;

};


using GenomeSequenceMap = std::map<ContigId_t, std::unique_ptr<ContigRecord>>;
class GenomeSequences {

public:

  explicit GenomeSequences(Logger& logger) : log(logger) {}
  GenomeSequences(const GenomeSequences&) = default;
  ~GenomeSequences() = default;

  GenomeSequences& operator=(const GenomeSequences&) = default;

  void addContigSequence(ContigId_t& contig, Sequence_t sequence);
  GenomeSequenceMap& getGenomeSequenceMap() { return genome_sequence_map_; }

private:

  Logger& log;

  GenomeSequenceMap genome_sequence_map_;
  GenomeAttributes attributes_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_DB_H
