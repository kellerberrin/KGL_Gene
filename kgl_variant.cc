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


std::string kgl::VariantGenome::typestr() const {

  switch(genomeType()) {

    case VariantGenomeType::UNKNOWN: return "UNKNOWN";

    case VariantGenomeType::CDS_CODING : return "CDS Coding";

    case VariantGenomeType::INTRON : return "INTRON";

    case VariantGenomeType::NON_CODING : return "NON Coding";

  }

  return "ERROR"; // Should not happen.

}

kgl::VariantGenomeType kgl::VariantGenome::genomeType() {

  if (variant_genome_type_ == VariantGenomeType::UNKNOWN) {

    std::shared_ptr<kgl::GeneFeature> gene_ptr;
    if (contig_ptr_->findGene(contig_offset_, gene_ptr)) {

      CDSArray cds_array;
      if(contig_ptr_->findOffsetCDS(contig_offset_, cds_array)) {

        variant_genome_type_ = VariantGenomeType::CDS_CODING;

      }
      else {

        variant_genome_type_ = VariantGenomeType::INTRON;

      }

    } else {

      variant_genome_type_ = VariantGenomeType::NON_CODING;

    }

  }

  return variant_genome_type_;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::ContigVariant> filtered_contig_ptr(std::make_shared<kgl::ContigVariant>(*this));
  // Inverts the bool returned by filterVariant(filter) because the delete pattern expects bool true for deletion.
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
