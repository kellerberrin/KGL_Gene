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

    case VariantGenomeType::NON_CODING : return "NON Coding";
  }

  return "ERROR"; // Should not happen.

}

kgl::VariantGenomeType kgl::VariantGenome::genomeType() {

  if (variant_genome_type_ == VariantGenomeType::UNKNOWN) {

    if (contig_ptr_->findOffsetCDS(contig_offset_, cds_ptr_vec_)) {

      variant_genome_type_ = VariantGenomeType::CDS_CODING;

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

std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::codingVariants(const std::shared_ptr<const GenomeDatabase> genome_db_ptr) const {

  // empty contig variants
  std::shared_ptr<ContigVariant> coding_contig_ptr(std::make_shared<ContigVariant>(contigId()));
  // Get the genome contig.
  std::shared_ptr<ContigFeatures> contig_ptr;
  if (genome_db_ptr->getContigSequence(contigId(), contig_ptr)) {

    for (auto cds_variant : offset_variant_map_) {

      std::vector<std::shared_ptr<CDSFeature>> cds_ptr_vec;
      if (contig_ptr->findOffsetCDS(cds_variant.second->contigOffset(), cds_ptr_vec)) {

        std::shared_ptr<Feature> gene_ptr = cds_ptr_vec.front()->getGene();
        if (gene_ptr) {

//          std::shared_ptr<Variant> coding_variant(std::make_shared<CodingSNPVariant>(cds_variant.second));

        } else {

          ExecEnv::log().error("Variant GeneFilter; Gene not found for CDS: {}", cds_ptr_vec.front()->id());

        }

      } else {

        ExecEnv::log().error("Contig: {} Offset: {} Variant CDS not found in genome database",
                             cds_variant.second->contigId(),
                             cds_variant.second->contigOffset());

      }

    }

  } else {

    ExecEnv::log().error("Variant contig: {} not found in genome database", contigId());

  }

  return coding_contig_ptr;

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


std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::codingVariants(const std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr) const {

  std::shared_ptr<GenomeVariant> variant_genome_ptr(std::make_shared<GenomeVariant>("Coding Variants", genomeId()));
  for (const auto& contig_variant : genome_variant_map_) {

    std::shared_ptr<ContigVariant> coding_contig_ptr = contig_variant.second->codingVariants(genome_db_ptr);
    variant_genome_ptr->addContigVariant(coding_contig_ptr);
    ExecEnv::log().vinfo("Contig: {} has: {} Coding variants", contig_variant.first, coding_contig_ptr->variantCount());

  }

  return variant_genome_ptr;

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
