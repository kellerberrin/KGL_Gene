//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_DISTANCE_H
#define KGL_SEQUENCE_DISTANCE_H

#include "kgl_genome_types.h"
#include "kgl_sequence_virtual.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_sequence_amino.h"


namespace kellerberrin::genome {   //  organization::project level namespace


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

  [[nodiscard]] virtual std::string distanceType() const = 0;

protected:

  [[nodiscard]] virtual CompareDistance_t distanceImpl( const VirtualSequence& sequenceA,
                                                        const VirtualSequence& sequenceB) const = 0;

};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Linear DNA specific calculateDistance algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LinearDNASequenceDistance : public SequenceDistance {

public:

  LinearDNASequenceDistance() = default;
  ~LinearDNASequenceDistance() override = default;

  [[nodiscard]] CompareDistance_t linear_distance( const DNA5SequenceLinear& sequenceA,
                                            const DNA5SequenceLinear& sequenceB) const {

    return distanceImpl(sequenceA, sequenceB);

  }

protected:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Linear DNA specific calculateDistance algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CodingDNASequenceDistance : public SequenceDistance {

public:

  CodingDNASequenceDistance() = default;
  ~CodingDNASequenceDistance() override = default;

  [[nodiscard]] CompareDistance_t coding_distance( const DNA5SequenceCoding& sequenceA,
                                            const DNA5SequenceCoding& sequenceB) const {

    return distanceImpl(sequenceA, sequenceB);

  }


private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA specific calculateDistance algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNASequenceDistance : public CodingDNASequenceDistance, public LinearDNASequenceDistance  {

public:

  DNASequenceDistance() = default;
  ~DNASequenceDistance() override = default;



private:


};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino specific calculateDistance algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AminoSequenceDistance  : public SequenceDistance {

public:

  AminoSequenceDistance() = default;
  ~AminoSequenceDistance() override = default;


  [[nodiscard]] CompareDistance_t amino_distance( const AminoSequence& sequenceA,
                                            const AminoSequence& sequenceB) const {

    return distanceImpl(sequenceA, sequenceB);

  }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA specific calculateDistance algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AnySequenceDistance : public DNASequenceDistance, public AminoSequenceDistance  {

public:

  AnySequenceDistance() = default;
  ~AnySequenceDistance() override = default;



private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global Levenshtein (Edit) calculateDistance - any size difference is counted as an edit operation (DNA and Amino)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LevenshteinGlobal: public AnySequenceDistance {

public:

  LevenshteinGlobal() = default;
  ~LevenshteinGlobal() override = default;

  [[nodiscard]] std::string distanceType() const override { return "Levenshtein Global"; }



protected:

  [[nodiscard]] CompareDistance_t distanceImpl( const VirtualSequence& sequenceA,
                                                const VirtualSequence& sequenceB) const override {

    return SequenceDistanceImpl().LevenshteinGlobal(sequenceA.getSequenceAsString(), sequenceB.getSequenceAsString());

  }


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local Levenshtein (Edit) distance - only the local match generates edit calculateDistance (DNA and Amino)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LevenshteinLocal: public AnySequenceDistance {

public:

  LevenshteinLocal() = default;
  ~LevenshteinLocal() override = default;

  [[nodiscard]] std::string distanceType() const override { return "Levenshtein Local"; }

protected:

  [[nodiscard]] CompareDistance_t distanceImpl( const VirtualSequence& sequenceA,
                                                const VirtualSequence& sequenceB) const override {

    return SequenceDistanceImpl().LevenshteinLocal(sequenceA.getSequenceAsString(), sequenceB.getSequenceAsString());

  }


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global Blosum80
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Blosum80Global: public AminoSequenceDistance {

public:

  Blosum80Global() = default;
  ~Blosum80Global() override = default;

  [[nodiscard]] std::string distanceType() const override { return "Blosum80 Global"; }



protected:

  [[nodiscard]] CompareDistance_t distanceImpl( const VirtualSequence& sequenceA,
                                                const VirtualSequence& sequenceB) const override {

    return SequenceDistanceImpl().globalblosum80Distance(sequenceA.getSequenceAsString(), sequenceB.getSequenceAsString());

  }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local Blosum80
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Blosum80Local: public AminoSequenceDistance {

public:

  Blosum80Local() = default;
  ~Blosum80Local() override = default;

  [[nodiscard]] std::string distanceType() const override { return "Blosum80 Local"; }



protected:

  [[nodiscard]] CompareDistance_t distanceImpl( const VirtualSequence& sequenceA,
                                                const VirtualSequence& sequenceB) const override {

    return SequenceDistanceImpl().localblosum80Distance(sequenceA.getSequenceAsString(), sequenceB.getSequenceAsString());

  }

};



}   // end namespace


#endif //KGL_KGL_SEQUENCE_DISTANCE_H
