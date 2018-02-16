//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_DISTANCE_H
#define KGL_SEQUENCE_DISTANCE_H

#include "kgl_genome_types.h"
#include "kgl_sequence_virtual.h"
#include "kgl_sequence_compare_impl.h"
#include "kgl_sequence_amino.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance is conceptually different from comparison.
// 1. Distances are always normalized for sequence length.
// For local distances, this is the size of the match.
// For global distances, this is the length of the sequences.
// 2. Distances are always positive or zero.
// 3. Distances are symmetric, d(x,y) = d(y,x)
// 4. Distances observe the triangle inequality, d(x,y) + d(y,z) >= d(x,z)
// 5. Distances are returned as a CompareDistance_t which is a double.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This virtual class hierarchy enforces correct typing for different combinations of
// local or global, dna or amino distances algorithms.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SequenceDistance {

public:

  SequenceDistance() = default;
  virtual ~SequenceDistance() = default;

  virtual std::string distanceType() const = 0;

protected:

  virtual CompareDistance_t distanceImpl(std::shared_ptr<const VirtualSequence> sequenceA,
                                         std::shared_ptr<const VirtualSequence> sequenceB) const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local distance algorithms for both DNA and Amino acids should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LocalSequenceDistance : public virtual SequenceDistance {

public:

  LocalSequenceDistance() = default;
  ~LocalSequenceDistance() override = default;


private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global distance algorithms for both DNA and Amino acids should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GlobalSequenceDistance : public virtual SequenceDistance {

public:

  GlobalSequenceDistance() = default;
  ~GlobalSequenceDistance() override = default;


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA specific distance algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNASequenceDistance : public virtual SequenceDistance {

public:

  DNASequenceDistance() = default;
  ~DNASequenceDistance() override = default;

  CompareDistance_t distance(std::shared_ptr<const DNA5SequenceLinear> sequenceA,
                             std::shared_ptr<const DNA5SequenceLinear> sequenceB) const {

    return distanceImpl(sequenceA, sequenceB);

  }

  CompareDistance_t distance(std::shared_ptr<const DNA5SequenceCoding> sequenceA,
                             std::shared_ptr<const DNA5SequenceCoding> sequenceB) const {

    return distanceImpl(sequenceA, sequenceB);

  }


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino specific distance algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AminoSequenceDistance  : public virtual SequenceDistance {

public:

  AminoSequenceDistance() = default;
  ~AminoSequenceDistance() override = default;


  CompareDistance_t distance(std::shared_ptr<const AminoSequence> sequenceA,
                             std::shared_ptr<const AminoSequence> sequenceB) const {

    return distanceImpl(sequenceA, sequenceB);

  }


private:


};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino specific distance algorithms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LocalAminoSequenceDistance  : public virtual AminoSequenceDistance, public virtual LocalSequenceDistance {

public:

  LocalAminoSequenceDistance() = default;
  ~LocalAminoSequenceDistance() override = default;


private:


};



class GlobalAminoSequenceDistance  : public virtual AminoSequenceDistance, public virtual GlobalSequenceDistance {

public:

  GlobalAminoSequenceDistance() = default;
  ~GlobalAminoSequenceDistance() override = default;


private:


};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA specific distance algorithms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LocalDNASequenceDistance  : public virtual DNASequenceDistance, public virtual LocalSequenceDistance {

public:

  LocalDNASequenceDistance() = default;
  ~LocalDNASequenceDistance() override = default;


private:


};



class GlobalDNASequenceDistance  : public virtual DNASequenceDistance, public virtual GlobalSequenceDistance {

public:

  GlobalDNASequenceDistance() = default;
  ~GlobalDNASequenceDistance() override = default;


private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of Distance algorithms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global Levenshtein (Edit) distance - any size difference is counted as an edit operation (DNA and Amino)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LevenshteinGlobal: public GlobalAminoSequenceDistance, public GlobalDNASequenceDistance {

public:

  LevenshteinGlobal() = default;
  virtual ~LevenshteinGlobal() = default;

  std::string distanceType() const override { return "Levenshtein Global"; }



private:

  CompareDistance_t distanceImpl(std::shared_ptr<const VirtualSequence> sequenceA,
                                 std::shared_ptr<const VirtualSequence> sequenceB) const override {

    return SequenceManipulation().LevenshteinGlobal(sequenceA->getSequenceAsString(), sequenceB->getSequenceAsString());

  }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local Levenshtein (Edit) distance - only the local match generates edit distance (DNA and Amino)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LevenshteinLocal: public LocalAminoSequenceDistance, public LocalDNASequenceDistance {

public:

  LevenshteinLocal() = default;
  virtual ~LevenshteinLocal() = default;

  std::string distanceType() const override { return "Levenshtein Local"; }

private:

  CompareDistance_t distanceImpl(std::shared_ptr<const VirtualSequence> sequenceA,
                                 std::shared_ptr<const VirtualSequence> sequenceB) const override {

    return SequenceManipulation().LevenshteinLocal(sequenceA->getSequenceAsString(), sequenceB->getSequenceAsString());

  }

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_KGL_SEQUENCE_DISTANCE_H
