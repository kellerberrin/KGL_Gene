//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kgl_patterns.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::VariantSequence::typestr() const {

  switch(genomeType()) {

    case VariantSequenceType::UNKNOWN: return "UNKNOWN";

    case VariantSequenceType::CDS_CODING : return "CDS Coding";

    case VariantSequenceType::INTRON : return "INTRON";

    case VariantSequenceType::NON_CODING : return "NON Coding";

  }

  return "ERROR"; // Should not happen.

}


std::string kgl::VariantSequence::genomeOutput() const {

  std:: stringstream ss;
// Contig.

  ss << contig()->contigId();
  ss << " " << typestr() << " ";

  if (type() != VariantSequenceType::NON_CODING) {

    ss << "Gene(s):";

    for (auto gene_ptr : geneMembership()) {

      ss << gene_ptr->id() << " " << (offset() - gene_ptr->sequence().begin()) << " ";

    }
  }

  return ss.str();

}

kgl::VariantSequenceType kgl::VariantSequence::genomeType() {

  if (variant_genome_type_ == VariantSequenceType::UNKNOWN) {

    if (contig_ptr_->findGenes(contig_offset_, gene_membership_)) {

      CDSArray cds_array;
      if(contig_ptr_->findOffsetCDS(contig_offset_, cds_array)) {

        variant_genome_type_ = VariantSequenceType::CDS_CODING;

      }
      else {

        variant_genome_type_ = VariantSequenceType::INTRON;

      }

    } else {

      variant_genome_type_ = VariantSequenceType::NON_CODING;

    }

  }

  return variant_genome_type_;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A compound variant. A collection of feature aligned and contiguous variants. Insertions and Deletions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::CompoundVariant::equivalent(const Variant& cmp_var) const {

  auto compound_var = dynamic_cast<const CompoundVariant*>(&cmp_var);

  if (compound_var == nullptr) return false;

  for (auto variant : variant_map_) {

    if (not ((*variant.second) == (*this))) return false;

  }

  return contigId() == compound_var->contigId() and contigOffset() == compound_var->contigOffset();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SNPVariant - SNPs generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::SNPVariantDNA5::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const SNPVariantDNA5*>(&cmp_var);

  if (cmp_snp == nullptr) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


std::string kgl::SNPVariantDNA5::output() const
{
  std::stringstream ss;
  ss << genomeOutput();
  ss << reference() << (contigOffset() + 1) << mutant() << " ";
  ss << mutantCount() << "/" << readCount() << " [";
  for (size_t idx = 0; idx < countArray().size(); ++idx) {
    ss << NucleotideColumn_DNA5::offsetToNucleotide(idx) << ":" << countArray()[idx] << " ";
  }
  ss << "]";
  return ss.str();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::ContigVariant> filtered_contig_ptr(std::make_shared<kgl::ContigVariant>(*this));
  // Complements the bool returned by filterVariant(filter) because the delete pattern expects bool true for deletion.
  auto predicate = [&](const OffsetVariantMap::const_iterator& it) { return not it->second->filterVariant(filter); };
  predicateIterableDelete(filtered_contig_ptr->offset_variant_map_,  predicate);

  return filtered_contig_ptr;

}


void kgl::ContigVariant::addVariant(ContigOffset_t contig_offset, std::shared_ptr<const kgl::Variant>& variant_ptr) {

  offset_variant_map_.insert(std::make_pair(contig_offset, variant_ptr));

}



std::ostream& kgl::operator<<(std::ostream &os, const kgl::ContigVariant& contig_variant) {

  for (auto& variant : contig_variant.offset_variant_map_) {

    os << *(variant.second) << '\n';

  }

  return os;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GenomeVariant::addContigVariant(std::shared_ptr<kgl::ContigVariant>& contig_variant) {

  auto result = genome_variant_map_.insert(std::make_pair(contig_variant->contigId(), contig_variant));

  return result.second;

}


std::shared_ptr<kgl::GenomeVariant> kgl::GenomeVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::GenomeVariant> filtered_genome_ptr(std::make_shared<kgl::GenomeVariant>(filter.filterName(),
                                                                                               genomeId()));

  ExecEnv::log().info("Applying filter: {}", filter.filterName());
  for (const auto& contig_variant : genome_variant_map_) {

    std::shared_ptr<kgl::ContigVariant> filtered_contig = contig_variant.second->filterVariants(filter);
    filtered_genome_ptr->addContigVariant(filtered_contig);
    ExecEnv::log().vinfo("Contig: {} has: {} filtered variants", contig_variant.first, filtered_contig->variantCount());

  }

  return filtered_genome_ptr;

}

std::ostream& kgl::operator<<(std::ostream &os, const kgl::GenomeVariant& genome_variant) {

  for (auto& contig_variant : genome_variant.genome_variant_map_) {

    os << contig_variant.second->contigId() << " variants: ";
    os << contig_variant.second->variantCount() << "\n";
    os << *(contig_variant.second);

  }

  os.flush();

  return os;

}
