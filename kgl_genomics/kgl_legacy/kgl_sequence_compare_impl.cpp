//
// Created by kellerberrin on 15/02/18.
//


#include "kgl_sequence_compare_impl.h"
#include "kel_exec_env.h"
#include "kgl_genome_types.h"

#include "edlib.h"


namespace kgl = kellerberrin::genome;


class kgl::SequenceComparison::SequenceManipImpl {

public:

  SequenceManipImpl() = default;
  ~SequenceManipImpl() = default;



  void editDNAItems(const std::string& reference,
                    const std::string& mutant,
                    EditVector& edit_vector) const;
private:

  EditVector createEditItems(const std::string& reference_str, const std::string& mutant_str, EditVector& edit_vector) const;
  EditVector createEditItemsEdlib(const std::string& reference_str, const std::string& mutant_str, EditVector& edit_vector) const;

};


void kgl::SequenceComparison::SequenceManipImpl::editDNAItems(const std::string& reference,
                                                              const std::string& mutant,
                                                              EditVector& edit_vector) const {
  EditVector check_edit_vector;
  createEditItems(reference, mutant, edit_vector);
  createEditItemsEdlib(reference, mutant, check_edit_vector);

  if (edit_vector.size() != check_edit_vector.size()) {

    ExecEnv::log().error("Edit vector size mis-match, edit_vector size: {}, check edit vector size: {}", edit_vector.size(), check_edit_vector.size());

  } else {

    for (size_t index = 0; index < edit_vector.size(); ++index) {

      bool compare = edit_vector[index].mutant_char == check_edit_vector[index].mutant_char
                     and edit_vector[index].reference_offset == check_edit_vector[index].reference_offset
                     and edit_vector[index].reference_char == check_edit_vector[index].reference_char;


      if (not compare) {

        ExecEnv::log().error("Edit item mis-match, edit_vector: {}{}{}, check edit vector: {}{}{}",
                             edit_vector[index].reference_char, edit_vector[index].reference_offset, edit_vector[index].mutant_char,
                             check_edit_vector[index].reference_char, check_edit_vector[index].reference_offset, check_edit_vector[index].mutant_char);
      }

    }

  }



}



kgl::EditVector kgl::SequenceComparison::SequenceManipImpl::createEditItemsEdlib(const std::string& reference,
                                                                            const std::string& mutant,
                                                                            EditVector& edit_vector) const {

  edit_vector.clear();

  EdlibAlignResult result = edlibAlign(mutant.c_str(), mutant.size(),reference.c_str(), reference.size(),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
  if (result.status == EDLIB_STATUS_OK) {

    EditItem edit_item;
    int mut_index = 0;
    for (int ref_index = 0; ref_index < result.alignmentLength; ++ref_index) {

      switch(result.alignment[ref_index]) {

        case 0:
          break;

        case 1:
          edit_item.reference_char = reference[ref_index];
          edit_item.reference_offset = ref_index;
          edit_item.mutant_char = EditItem::INSERTION;
          edit_vector.push_back(edit_item);
          ++mut_index;
          break;

        case 2:
          edit_item.reference_char = reference[ref_index];
          edit_item.reference_offset = ref_index;
          edit_item.mutant_char = EditItem::DELETION;
          edit_vector.push_back(edit_item);
          --mut_index;
          break;

        case 3:
          edit_item.reference_char = reference[ref_index];
          edit_item.reference_offset = ref_index;
          edit_item.mutant_char = mutant[mut_index];
          edit_vector.push_back(edit_item);
          break;

      }

      ++mut_index;

    }


  } else {

    ExecEnv::log().error("Edlib - problem generating cigar reference:{}, alternate: {}", reference, mutant);

  }

  edlibFreeAlignResult(result);

  return edit_vector;

}


kgl::EditVector kgl::SequenceComparison::SequenceManipImpl::createEditItems(const std::string&,
                                                                            const std::string&,
                                                                            EditVector&) const {

  return std::vector<EditItem>();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::FreeBayesVCFImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::SequenceComparison::SequenceComparison() : sequence_manip_impl_ptr_(std::make_unique<kgl::SequenceComparison::SequenceManipImpl>()) {}
kgl::SequenceComparison::~SequenceComparison() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


void kgl::SequenceComparison::editDNAItems(const std::string& reference,
                                           const std::string& mutant,
                                           EditVector& edit_vector) const {

  sequence_manip_impl_ptr_->editDNAItems(reference, mutant, edit_vector);

}



