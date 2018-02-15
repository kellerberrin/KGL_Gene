//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_DISTANCE_H
#define KGL_SEQUENCE_DISTANCE_H

#include "kgl_genome_types.h"
#include "kgl_sequence_virtual.h"
#include "kgl_sequence_compare_impl.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




class SequenceDistance {

public:

  SequenceDistance() = default;
  virtual ~SequenceDistance() = default;

  virtual std::string distanceType() const = 0;

  virtual CompareDistance_t distance(std::shared_ptr<const VirtualSequence> sequenceA,
                                     std::shared_ptr<const VirtualSequence> sequenceB) const = 0;

private:


};


class SequenceLocalDistance : public SequenceDistance {

public:

  SequenceLocalDistance() = default;
  virtual ~SequenceLocalDistance() = default;


private:


};


class SequenceGlobalDistance : public SequenceDistance {

public:

  SequenceGlobalDistance() = default;
  virtual ~SequenceGlobalDistance() = default;


private:


};


class DNASequenceLocalDistance : public SequenceLocalDistance {

public:

  DNASequenceLocalDistance() = default;
  ~DNASequenceLocalDistance() override = default;


private:


};


class AminoSequenceLocalDistance  : public SequenceLocalDistance {

public:

  AminoSequenceLocalDistance() = default;
  ~AminoSequenceLocalDistance() override = default;


private:


};


class DNASequenceGlobalDistance : public SequenceGlobalDistance {

public:

  DNASequenceGlobalDistance() = default;
  ~DNASequenceGlobalDistance() override = default;


private:


};


class AminoSequenceGlobalDistance  : public SequenceGlobalDistance {

public:

  AminoSequenceGlobalDistance() = default;
  ~AminoSequenceGlobalDistance() override = default;


private:


};






class LevenshteinGlobal: public SequenceDistance {

public:

  LevenshteinGlobal() = default;
  virtual ~LevenshteinGlobal() = default;

  std::string distanceType() const override { return "Levenshtein Global"; }

  CompareDistance_t distance(std::shared_ptr<const VirtualSequence> sequenceA,
                             std::shared_ptr<const VirtualSequence> sequenceB) const override {

    return SequenceManipulation().LevenshteinGlobal(sequenceA->getSequenceAsString(), sequenceB->getSequenceAsString());

  }


private:


};




class LevenshteinLocal: public SequenceDistance {

public:

  LevenshteinLocal() = default;
  virtual ~LevenshteinLocal() = default;

  std::string distanceType() const override { return "Levenshtein Local"; }

  CompareDistance_t distance(std::shared_ptr<const VirtualSequence> sequenceA,
                             std::shared_ptr<const VirtualSequence> sequenceB) const override {

    return SequenceManipulation().LevenshteinLocal(sequenceA->getSequenceAsString(), sequenceB->getSequenceAsString());

  }


private:


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_KGL_SEQUENCE_DISTANCE_H
