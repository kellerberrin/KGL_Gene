//
// Created by kellerberrin on 23/12/17.
//


#include "kgl_variant_single.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// InsertVariant - generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::DeleteVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << quality() << delimiter;
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail) {

    ss << evidence()->output(delimiter, output_index);

  }

  ss << '\n';

  return ss.str();

}


// This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
bool kgl::DeleteVariant::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                              SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                                              ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                                              SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                                              std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {


  CodingSequenceArray coding_sequence_array = codingSequences();

  // Check that we have a variant in a coding sequence.
  if (coding_sequence_array.empty()) {

    ExecEnv::log().warn("mutateCodingSequence(), variant: {} not in a coding sequence",
                        output(' ', VariantOutputIndex::START_0_BASED, true));
    return true; // just ignored.

  }

  std::shared_ptr<const CodingSequence> coding_sequence = coding_sequence_array.getFirst();

  // Check the sequence id.
  if (coding_sequence->getCDSParent()->id() != sequence_id) {

    ExecEnv::log().warn("mutateCodingSequence(), variant: {} does not mutate sequence id: {}",
                        output(' ', VariantOutputIndex::START_0_BASED, true), sequence_id);
    return true; // just ignored.

  }

  // Get the variant sequence offset
  ContigOffset_t sequence_offset;
  ContigSize_t sequence_length;
  if(not contig()->sequence().offsetWithinSequence(coding_sequence, contigOffset(), sequence_offset, sequence_length)) {

    ExecEnv::log().error("mutateCodingSequence(), unable to retrieve coding sequence offset for variant: {}",
                         output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;

  }

  // Check that the sequence lengths match.
  if (mutated_sequence->length() != sequence_size) {

    ExecEnv::log().error("delete::mutateCodingSequence(), sequence length: {} not equal to expected length: {}",
                         mutated_sequence->length(), sequence_length);
    return false;

  }

  SignedOffset_t check_offset = static_cast<SignedOffset_t>(sequence_offset) + offset_adjust;

  std::cout << " adjusted offset:" << check_offset << " seq offset:" << sequence_offset << " adjust:" <<  offset_adjust << std::endl;

  // Check the adjusted offset
  if (check_offset < 0 or check_offset >= static_cast<SignedOffset_t>(mutated_sequence->length())) {

    ExecEnv::log().error("mutateCodingSequence(), adjusted offset: {} out of range for sequence length: {}",
                         check_offset, mutated_sequence->length());
    return false;

  }

  ContigOffset_t adjusted_offset = static_cast<ContigOffset_t>(check_offset);

  // check the adjusted sequence offset + the delete size does not exceed the length of the sequence.
  if (adjusted_offset + size() > mutated_sequence->length()) {

    ExecEnv::log().error("mutateCodingSequence(), adjusted offset: {}, with delete size: {} out of range for sequence length: {}",
                         check_offset, size(), mutated_sequence->length());
    return false;

  }

  // Check that the sequence base code matches the original strand adjusted base code recorded in the variant.
  if (mutated_sequence->at(adjusted_offset) != strandReference()) {

    ExecEnv::log().warn("mutateCodingSequence(), unexpected; base: {} at seq. offset: {} (strand) reference: {}, probable duplicate variant",
                        CodingDNA5::convertToChar(mutated_sequence->at(sequence_offset)), sequence_offset,
                        CodingDNA5::convertToChar(strandReference()));

  }

  sequence_size_adjust = -1 * size(); // reduce the sequence size.

  // All is good, so mutate the sequence.
  return mutated_sequence->deleteSubSequence(sequence_offset, size());

}


bool kgl::DeleteVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const DeleteVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and type() == cmp_snp->type()
         and variantType() == cmp_snp->variantType()
         and codingSequenceId() == cmp_snp->codingSequenceId()
         and reference() == cmp_snp->reference();

}

std::string kgl::DeleteVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  if (type() == VariantSequenceType::CDS_CODING) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;


    ContigSize_t base_in_codon;
    ContigOffset_t codon_offset;

    codonOffset(codon_offset, base_in_codon);

    ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
    ss << offsetOutput(base_in_codon, output_index) << delimiter;
    ss << "-(" << size() << ")";
    ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
    ss << offsetOutput(base_in_codon, output_index) << delimiter;
    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << mutantChar() << delimiter;

  } else if (type() == VariantSequenceType::INTRON) {

    std::shared_ptr<const GeneFeature> gene_ptr = geneMembership().front();
    ss << gene_ptr->id() << delimiter;
    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << mutantChar() << delimiter;

  } else { // else non coding (non-gene) variant or unknown

    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << mutantChar() << delimiter;

  }

  return ss.str();

}
