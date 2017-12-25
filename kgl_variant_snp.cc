//
// Created by kellerberrin on 31/10/17.
//

#include "kgl_variant_single.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SNPVariant - SNPs generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::string kgl::SNPVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
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


bool kgl::SNPVariant::mutateCodingSequence(const FeatureIdent_t& sequence_id,
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
  if (sequence_length != sequence_size) {

    ExecEnv::log().error("mutateCodingSequence(), unexpected; variant sequence length: {} not equal mutate length: {}",
                         sequence_length, mutated_sequence->length());
    return false;

  }

  SignedOffset_t check_offset = static_cast<SignedOffset_t>(sequence_offset) + offset_adjust;

  // Check the adjusted offset
  if (check_offset < 0 or check_offset >= static_cast<SignedOffset_t>(mutated_sequence->length())) {

    ExecEnv::log().error("mutateCodingSequence(), adjusted offset: {} out of range for sequence length: {}",
                        check_offset, mutated_sequence->length());
    return false;

  }

  ContigOffset_t adjusted_offset = static_cast<ContigOffset_t>(check_offset);

  // Check that the sequence base code matches the original strand adjusted base code recorded in the variant.
  if (mutated_sequence->at(adjusted_offset) != strandReference()) {

    ExecEnv::log().warn("mutateCodingSequence(), unexpected; base: {} at seq. offset: {} (strand) reference: {}, probable duplicate variant",
    CodingDNA5::convertToChar(mutated_sequence->at(sequence_offset)), sequence_offset,
    CodingDNA5::convertToChar(strandReference()));

  }

  sequence_size_adjust = 0;

  // All is good, so mutate the sequence.
  return mutated_sequence->modifyBase(sequence_offset, strandMutant());


}


bool kgl::SNPVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const SNPVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and type() == cmp_snp->type()
         and variantType() == cmp_snp->variantType()
         and codingSequenceId() == cmp_snp->codingSequenceId()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


bool kgl::SNPVariant::codonMutation(ContigOffset_t &codon_offset,
                                    ContigSize_t &base_in_codon,
                                    AminoAcid::Alphabet &reference_amino,
                                    AminoAcid::Alphabet &mutant_amino) const {


  if (not codonOffset(codon_offset, base_in_codon)) {

    ExecEnv::log().error("codonMutation() called for non coding variant");
    reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
    mutant_amino = AminoAcid::AMINO_UNKNOWN;
    return false;

  }

  auto sequence_offset = static_cast<ContigOffset_t>(codon_offset * Codon::CODON_SIZE);

  const std::shared_ptr<const CodingSequence> coding_seq_ptr = codingSequences().getFirst();
  std::shared_ptr<DNA5SequenceCoding> codon_sequence = contig()->sequence().codingSubSequence(coding_seq_ptr,
                                                                                              sequence_offset,
                                                                                              Codon::CODON_SIZE);
  if (codon_sequence->length() != Codon::CODON_SIZE) {

    ExecEnv::log().error("codonMutation(), expected codon sequence size 3: got size: {}", codon_sequence->length());
    reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
    mutant_amino = AminoAcid::AMINO_UNKNOWN;
    return false;

  }

  Codon codon(codon_sequence, 0);  // Create the codon.

  reference_amino = contig()->getAminoAcid(codon);
  codon.modifyBase(base_in_codon, strandMutant());
  mutant_amino = contig()->getAminoAcid(codon);

  return true;

}


std::string kgl::SNPVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  if (type() == VariantSequenceType::CDS_CODING) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    ContigOffset_t codon_offset;
    ContigSize_t base_in_codon;
    AminoAcid::Alphabet reference_amino;
    AminoAcid::Alphabet mutant_amino;

    if (codonMutation(codon_offset, base_in_codon, reference_amino, mutant_amino)) {

      ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
      ss << offsetOutput(base_in_codon, output_index) << delimiter;
      ss << AminoAcid::convertToChar(reference_amino) << offsetOutput(codon_offset, output_index);
      ss << AminoAcid::convertToChar(mutant_amino) << delimiter;
      ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
      ss << mutantChar() << delimiter;
//        sequence->getGene()->recusivelyPrintsubfeatures();

    }


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
