//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_COMPARE_H
#define KGL_SEQUENCE_COMPARE_H

#include "kgl_genome_types.h"
#include "kgl_sequence_virtual.h"
#include "kgl_sequence_compare_impl.h"
#include "kgl_sequence_amino.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This virtual class hierarchy enforces correct typing for different combinations of
// local or global, dna or amino distances algorithms.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SequenceCompare {

public:

  SequenceCompare() = default;
  virtual ~SequenceCompare() = default;

  virtual std::string compareType() const = 0;

protected:

  virtual CompareScore_t compareImpl(std::shared_ptr<const VirtualSequence> sequenceA,
                                     std::shared_ptr<const VirtualSequence> sequenceB,
                                     std::string& compare_str) const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Local distance algorithms for both DNA and Amino acids should inherit from here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LocalSequenceCompare : public virtual SequenceCompare {

public:

  LocalSequenceCompare() = default;
  ~LocalSequenceCompare() override = default;


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global distance algorithms for both DNA and Amino acids should inherit from here
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

  CompareScore_t compare(std::shared_ptr<const DNA5SequenceLinear> sequenceA,
                         std::shared_ptr<const DNA5SequenceLinear> sequenceB,
                         std::string& compare_str) const {

    return compareImpl(sequenceA, sequenceB, compare_str);

  }

  CompareScore_t compare(std::shared_ptr<const DNA5SequenceCoding> sequenceA,
                         std::shared_ptr<const DNA5SequenceCoding> sequenceB,
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


  CompareDistance_t distance(std::shared_ptr<const AminoSequence> sequenceA,
                             std::shared_ptr<const AminoSequence> sequenceB,
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



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of Comparison algorithms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class MyerHirschbergGlobal: public GlobalAminoSequenceCompare, public GlobalDNASequenceCompare {

public:

  MyerHirschbergGlobal() = default;
  virtual ~MyerHirschbergGlobal() = default;

  std::string compareType() const override { return "MyerHirschberg Global"; }

private:

  CompareScore_t compareImpl(std::shared_ptr<const VirtualSequence> sequenceA,
                             std::shared_ptr<const VirtualSequence> sequenceB,
                             std::string& compare_str) const override {

    return SequenceComparison().MyerHirschbergGlobal(sequenceA->getSequenceAsString(),
                                                       sequenceB->getSequenceAsString(),
                                                       compare_str);

  }

};




class MyerHirschbergLocal: public LocalAminoSequenceCompare, public LocalDNASequenceCompare {

public:

  MyerHirschbergLocal() = default;
  virtual ~MyerHirschbergLocal() = default;

  std::string compareType() const override { return "MyerHirschberg Local"; }

private:

  CompareScore_t compareImpl(std::shared_ptr<const VirtualSequence> sequenceA,
                             std::shared_ptr<const VirtualSequence> sequenceB,
                             std::string& compare_str) const override {

    return SequenceComparison().MyerHirschbergLocal(sequenceA->getSequenceAsString(),
                                                      sequenceB->getSequenceAsString(),
                                                      compare_str);

  }

};



class DNALocalAffineGap: public LocalDNASequenceCompare {

public:

  DNALocalAffineGap() = default;
  virtual ~DNALocalAffineGap() = default;

  std::string compareType() const override { return "DNALocalAffineGap"; }

private:

  CompareScore_t compareImpl(std::shared_ptr<const VirtualSequence> sequenceA,
                             std::shared_ptr<const VirtualSequence> sequenceB,
                             std::string& compare_str) const override {

    return SequenceComparison().DNALocalAffineGap(sequenceA->getSequenceAsString(),
                                                    sequenceB->getSequenceAsString(),
                                                    compare_str);

  }

};





}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_COMPARE_H
