//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include <memory>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the minimum read count SNP generation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ReadCountFilter::implementFilter(const SNPVariant& variant) const {

  if (variant.evidence()->isReadCount())  {

    auto read_count_ptr = std::dynamic_pointer_cast<const ReadCountEvidence>(variant.evidence());

    if (not read_count_ptr) {

      return true;

    }

    return read_count_ptr->readCount() >= read_count_;

  }

  return true;

}


std::string kgl::ReadCountFilter::filterName() const {

  std::ostringstream oss;
  oss << "Minimum Read Count >= " << read_count_;
  return oss.str();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the minimum mutant read proportion in a candidate SNP.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::MutantProportionFilter::implementFilter(const SNPVariant& variant) const {

  if (variant.evidence()->isReadCount())  {

    auto read_count_ptr = std::dynamic_pointer_cast<const ReadCountEvidence>(variant.evidence());

    if (not read_count_ptr) {

      return true;

    }

    return read_count_ptr->proportion() >= mutant_proportion_;

  }

  return true;

}


std::string kgl::MutantProportionFilter::filterName() const {

  std::ostringstream oss;
  oss << "Minimum Read Proportion >= " << mutant_proportion_;
  return oss.str();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variantss to coding sequences only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::CodingFilter::implementFilter(const Variant& variant) const {

  return variant.type() == VariantSequenceType::CDS_CODING;

}


std::string kgl::CodingFilter::filterName() const {

  return "Variant in Coding Sequence";

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ContigFilter::implementFilter(const Variant& variant) const {

    return variant.contigId() == contig_ident_;

}

std::string kgl::ContigFilter::filterName() const {

  return "Variant in Contig: " + contig_ident_;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular gene.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GeneFilter::implementFilter(const Variant& variant) const {

  GeneVector gene_ptr_vec;
  if (not variant.geneMembership().empty()) {

    std::shared_ptr<const GeneFeature> gene_ptr = variant.geneMembership().front();
    if (gene_ptr->id() == gene_ident_) {

      return true;

    }

  }

  return false;

}


std::string kgl::GeneFilter::filterName() const {

  return "Variant in Gene: " + gene_ident_;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::SequenceFilter::implementFilter(const Variant& variant) const {

  if (not variant.codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence_ptr = variant.codingSequences().getFirst();
    if (sequence_ptr->getCDSParent()->id() == sequence_ident_) {

      return true;

    }

  }

  return false;

}


std::string kgl::SequenceFilter::filterName() const {

  return "Variant in Sequence: " + sequence_ident_;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remove synonymous coding SNPs and compound SNPs.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::SynonymousSNPFilter::implementFilter(const Variant& variant) const {

  if (variant.variantType() == VariantType::SNP) {

    if (variant.type() == VariantSequenceType::CDS_CODING) {

      auto SNP = dynamic_cast<const SNPVariant*>(&variant);

      if (not SNP) {

        ExecEnv::log().error("SynonymousFilter(); Bad variant type for variant: {}",
                             variant.output(' ', VariantOutputIndex::START_0_BASED, true));
        return true;

      }

      ContigOffset_t codon_offset;
      ContigSize_t base_in_codon;
      AminoAcid::Alphabet reference_amino;
      AminoAcid::Alphabet mutant_amino;
      SNP->codonMutation(codon_offset, base_in_codon, reference_amino, mutant_amino);

      return not (reference_amino == mutant_amino);

    }


  } else if (variant.variantType() == VariantType::COMPOUND_SNP) {

    if (variant.type() == VariantSequenceType::CDS_CODING) {

      auto cmp_SNP = dynamic_cast<const CompoundSNP*>(&variant);

      if (not cmp_SNP) {

        ExecEnv::log().error("SynonymousFilter(); Bad variant type for variant: {}",
                             variant.output(' ', VariantOutputIndex::START_0_BASED, true));
        return true;

      }

      ContigOffset_t codon_offset;
      AminoAcid::Alphabet reference_amino;
      AminoAcid::Alphabet mutant_amino;
      cmp_SNP->codonMutation(codon_offset, reference_amino, mutant_amino);

      return not (reference_amino == mutant_amino);

    }

  }

  return true;

}



std::string kgl::SynonymousSNPFilter::filterName() const {

  return "Remove Synonymous Coding (single and compound) SNPs";

}
