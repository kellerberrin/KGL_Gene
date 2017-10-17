///
// Created by kellerberrin on 10/10/17.
//

#ifndef KGL_GENOME_FEATURE_H
#define KGL_GENOME_FEATURE_H


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

class CDSRecord : public FeatureRecord {

public:

  CDSRecord(const FeatureIdent_t &id,
            const CDSPhaseType_t phase,
            const std::shared_ptr<ContigRecord> &contig_ptr,
            const FeatureSequence &sequence) : FeatureRecord(id, CDS_TYPE, contig_ptr, sequence), phase_(phase) {}

  CDSRecord(const CDSRecord &) = default;
  ~CDSRecord() override = default;

  CDSRecord &operator=(const CDSRecord &) = default;

  CDSPhaseType_t phase() const { return phase_; }

  // CDS Type.
  constexpr static const char* CDS_TYPE = "CDS";

private:

  CDSPhaseType_t phase_;

};


class GeneRecord : public FeatureRecord {

public:

  GeneRecord(const FeatureIdent_t &id,
             const std::shared_ptr<ContigRecord> &contig_ptr,
             const FeatureSequence &sequence) : FeatureRecord(id, GENE_TYPE, contig_ptr, sequence) {}

  GeneRecord(const GeneRecord &) = default;
  ~GeneRecord() override = default;

  GeneRecord &operator=(const GeneRecord &) = default;

  // CDS Type.
  constexpr static const char* GENE_TYPE = "GENE";

private:


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_FEATURE_H
