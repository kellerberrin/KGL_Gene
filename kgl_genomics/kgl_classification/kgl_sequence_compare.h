//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_COMPARE_H
#define KGL_SEQUENCE_COMPARE_H

#include "kgl_genome_types.h"
#include "kgl_sequence_virtual.h"
#include "kgl_sequence_compare_impl.h"
#include "kgl_sequence_amino.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This virtual class hierarchy enforces correct typing for different combinations of
// local or global, dna or amino distances algorithms.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SequenceCompare {

public:

  SequenceCompare() = default;
  virtual ~SequenceCompare() = default;

  [[nodiscard]] virtual std::string compareType() const = 0;

protected:

  [[nodiscard]]virtual CompareScore_t compareImpl( const VirtualSequence& sequenceA,
                                                   const VirtualSequence& sequenceB,
                                                   std::string& compare_str) const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local calculateDistance algorithms for both DNA and Amino acids should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LocalSequenceCompare : public virtual SequenceCompare {

public:

  LocalSequenceCompare() = default;
  ~LocalSequenceCompare() override = default;


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global calculateDistance algorithms for both DNA and Amino acids should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GlobalSequenceCompare : public virtual SequenceCompare {

public:

  GlobalSequenceCompare() = default;
  ~GlobalSequenceCompare() override = default;


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA specific comparison algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DNASequenceCompare : public virtual SequenceCompare {

public:

  DNASequenceCompare() = default;
  ~DNASequenceCompare() override = default;

  [[nodiscard]] CompareScore_t compare( const DNA5SequenceLinear& sequenceA,
                                        const DNA5SequenceLinear& sequenceB,
                                        std::string& compare_str) const {

    return compareImpl(sequenceA, sequenceB, compare_str);

  }

  [[nodiscard]] CompareScore_t compare( const DNA5SequenceCoding& sequenceA,
                                        const DNA5SequenceCoding& sequenceB,
                                        std::string& compare_str) const {

    return compareImpl(sequenceA, sequenceB, compare_str);

  }


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino specific comparison algorithms should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AminoSequenceCompare  : public virtual SequenceCompare {

public:

  AminoSequenceCompare() = default;
  ~AminoSequenceCompare() override = default;


  [[nodiscard]] CompareDistance_t distance( const AminoSequence& sequenceA,
                                            const AminoSequence& sequenceB,
                                            std::string& compare_str) const {

    return compareImpl(sequenceA, sequenceB, compare_str);

  }


private:


};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino specific comparison algorithms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LocalAminoSequenceCompare  : public virtual AminoSequenceCompare, public virtual LocalSequenceCompare {

public:

  LocalAminoSequenceCompare() = default;
  ~LocalAminoSequenceCompare() override = default;


private:


};



class GlobalAminoSequenceCompare  : public virtual AminoSequenceCompare, public virtual GlobalSequenceCompare {

public:

  GlobalAminoSequenceCompare() = default;
  ~GlobalAminoSequenceCompare() override = default;


private:


};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA specific comparison algorithms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LocalDNASequenceCompare  : public virtual DNASequenceCompare, public virtual LocalSequenceCompare {

public:

  LocalDNASequenceCompare() = default;
  ~LocalDNASequenceCompare() override = default;

private:


};



class GlobalDNASequenceCompare  : public virtual DNASequenceCompare, public virtual GlobalSequenceCompare {

public:

  GlobalDNASequenceCompare() = default;
  ~GlobalDNASequenceCompare() override = default;


private:


};




class DNALocalAffineGap: public LocalDNASequenceCompare {

public:

  DNALocalAffineGap() = default;
  virtual ~DNALocalAffineGap() = default;

  [[nodiscard]] std::string compareType() const override { return "DNALocalAffineGap"; }

private:

  [[nodiscard]] CompareScore_t compareImpl( const VirtualSequence& sequenceA,
                                            const VirtualSequence& sequenceB,
                                            std::string& compare_str) const override {

    return SequenceComparison().DNALocalAffineGap(sequenceA.getSequenceAsString(),
                                                    sequenceB.getSequenceAsString(),
                                                    compare_str);

  }

};



}   // end namespace


#endif //KGL_SEQUENCE_COMPARE_H
