//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_COMPARE_IMPL_H
#define KGL_SEQUENCE_COMPARE_IMPL_H


#include <memory>
#include <string>
#include <map>
#include <sstream>
#include "kel_logging.h"
#include "kgl_genome_types.h"
#include "kgl_sequence_base.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class EditItem{

public:

  EditItem() = default;
  EditItem(const EditItem&) = default;
  ~EditItem() = default;

  EditItem& operator=(const EditItem&) = default;

  char reference_char;
  ContigOffset_t reference_offset;
  char mutant_char;

  constexpr static const char INSERTION = '+';
  constexpr static const char DELETION = '-';

  bool isIndel() const { return mutant_char == INSERTION or mutant_char == DELETION; }

};


using EditVector = std::vector<EditItem>;

class SequenceComparison {

public:

  explicit SequenceComparison();
  ~SequenceComparison();

  // Comparison

  CompareScore_t MyerHirschbergGlobal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  CompareScore_t MyerHirschbergLocal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  CompareScore_t DNALocalAffineGap(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  // Compare sequences in mutation format.
  void editDNAItems(const std::string& reference,
                    const std::string& mutant,
                    EditVector& edit_vector) const;

private:

  class SequenceManipImpl;       // Forward declaration of the Sequence Manipulation implementation class
  std::unique_ptr<SequenceManipImpl> sequence_manip_impl_ptr_;    // Sequence Manipulation PIMPL


};



}   // end namespace



#endif //KGL_KGL_SEQUENCE_COMPARE_IMPL_H
