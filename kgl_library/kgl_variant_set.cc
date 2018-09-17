// MIT License
//
// Created by kellerberrin on 17/10/17.
//

#include <algorithm>
#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_variant_db.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// HomologousVariant - All the variant features that map onto a homologous contig region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::HomologousVariant::isElement(const Variant& variant) const {

  auto result = offset_variant_map_.equal_range(variant.offset());

  for (auto it = result.first; it != result.second; ++it) {

    if (*(it->second) == variant) return true;

  }

  return false;

}


std::shared_ptr<kgl::HomologousVariant>
kgl::HomologousVariant::Union(std::shared_ptr<const kgl::HomologousVariant> homo_variant_ptr) const {

  std::shared_ptr<kgl::HomologousVariant> union_variant_ptr(std::make_shared<kgl::HomologousVariant>(phaseId()));

  std::set_union(offset_variant_map_.begin(), offset_variant_map_.end(),
                 homo_variant_ptr->offset_variant_map_.begin(), homo_variant_ptr->offset_variant_map_.end(),
                 std::inserter(union_variant_ptr->offset_variant_map_, union_variant_ptr->offset_variant_map_.begin()));

  return union_variant_ptr;

}


std::shared_ptr<kgl::HomologousVariant>
kgl::HomologousVariant::Intersection(std::shared_ptr<const kgl::HomologousVariant> homo_variant_ptr) const {

  std::shared_ptr<kgl::HomologousVariant> intersect_variant_ptr = deepCopy();

  auto pred = [](const OffsetVariantMap::const_iterator& mod_it,
                 const OffsetVariantMap::const_iterator& ref_it) { return *(mod_it->second) == *(ref_it->second); };

  intersectIterable(intersect_variant_ptr->offset_variant_map_, homo_variant_ptr->offset_variant_map_, pred);

  return intersect_variant_ptr;

}


std::shared_ptr<kgl::HomologousVariant>
kgl::HomologousVariant::Difference(std::shared_ptr<const kgl::HomologousVariant> homo_variant_ptr) const {

  std::shared_ptr<kgl::HomologousVariant> diff_variant_ptr = deepCopy();

  for (auto variant : homo_variant_ptr->getMap()) {

    diff_variant_ptr->eraseVariant(variant.second);

  }

  return diff_variant_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto a contig region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::ContigVariant::isElement(const Variant& variant) const {

  if (variant.phaseId() >= ploidy()) {

    ExecEnv::log().error("ContigVariant::isElement(); Contig ploidy: {}, Variant has incompatible phasing: {}",
                         ploidy(), variant.output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  return ploidy_vector_[variant.phaseId()]->isElement(variant);

}


std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::Union(std::shared_ptr<const kgl::ContigVariant> contig_variant_ptr) const {

  std::shared_ptr<kgl::ContigVariant> union_variant_ptr(std::make_shared<ContigVariant>(contigId(),ploidy()));

  if (ploidy() != contig_variant_ptr->ploidy()) {

    ExecEnv::log().error("ContigVariant::Union(); Mismatching contig ploidy: {}, union contig ploidy: {}",
                         ploidy(), contig_variant_ptr->ploidy());

    return union_variant_ptr; // empty contig.

  }

  size_t index = 0;
  for (auto homologous : getVector()) {

    union_variant_ptr->ploidy_vector_[index] = homologous->Union(contig_variant_ptr->ploidy_vector_[index]);
    index++;

  }

  return union_variant_ptr;

}


std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::Intersection(std::shared_ptr<const kgl::ContigVariant> contig_variant_ptr) const {

  std::shared_ptr<kgl::ContigVariant> intersect_variant_ptr(std::make_shared<ContigVariant>(contigId(),ploidy()));

  if (ploidy() != contig_variant_ptr->ploidy()) {

    ExecEnv::log().error("ContigVariant::Intersection(); Mismatching contig ploidy: {}, intersection contig ploidy: {}",
                         ploidy(), contig_variant_ptr->ploidy());

    return intersect_variant_ptr; // empty contig.

  }

  size_t index = 0;
  for (auto homologous : getVector()) {

    intersect_variant_ptr->ploidy_vector_[index] = homologous->Intersection(contig_variant_ptr->ploidy_vector_[index]);
    index++;

  }

  return intersect_variant_ptr;

}


std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::Difference(std::shared_ptr<const kgl::ContigVariant> contig_variant_ptr) const {

  std::shared_ptr<kgl::ContigVariant> difference_variant_ptr(std::make_shared<ContigVariant>(contigId(),ploidy()));

  if (ploidy() != contig_variant_ptr->ploidy()) {

    ExecEnv::log().error("ContigVariant::Difference(); Mismatching contig ploidy: {}, intersection contig ploidy: {}",
                         ploidy(), contig_variant_ptr->ploidy());

    return difference_variant_ptr; // empty contig.

  }

  size_t index = 0;
  for (auto homologous : getVector()) {

    difference_variant_ptr->ploidy_vector_[index] = homologous->Difference(contig_variant_ptr->ploidy_vector_[index]);
    index++;

  }

  return difference_variant_ptr;

}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::GenomeVariant::isElement(const Variant& variant) const {

  auto result = genome_variant_map_.find(variant.contigId());

  if (result == genome_variant_map_.end()) {

    return false;

  }

  return result->second->isElement(variant);

}

// Returns a GenomeVariant that is the union and returns (*this U *genome_variant_ptr)
std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::Union(std::shared_ptr<const kgl::GenomeVariant> genome_variant_ptr) const {

  std::shared_ptr<kgl::GenomeVariant> genome_union(std::make_shared<kgl::GenomeVariant>(genomeId(), ploidy()));

  for (auto contig_variant : genome_variant_map_) {

    auto result = genome_variant_ptr->genome_variant_map_.find(contig_variant.first);

    std::shared_ptr<kgl::ContigVariant> union_contig_ptr;

    if (result != genome_variant_map_.end()) {

      union_contig_ptr = contig_variant.second->Union(result->second);
      genome_union->addContigVariant(union_contig_ptr);

    } else {

      genome_union->addContigVariant(contig_variant.second);

    }

  }

  for (auto contig_variant: genome_variant_ptr->genome_variant_map_) {

    auto result = genome_union->genome_variant_map_.find(contig_variant.first);

    if (result == genome_union->genome_variant_map_.end()) {

      genome_union->addContigVariant(result->second);

    }

  }

  for (auto contig_variant : genome_union->genome_variant_map_) {

    kgl::ExecEnv::log().vinfo("Union Contig: {} contains: {} variants",
                              contig_variant.first,
                              contig_variant.second->variantCount());

  }

  return genome_union;

}


// Returns a GenomeVariant that is the intersection return (*this & *genome_variant)
std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::Intersection(std::shared_ptr<const kgl::GenomeVariant> genome_variant) const {


  std::shared_ptr<kgl::GenomeVariant> genome_intersection(std::make_shared<kgl::GenomeVariant>(genomeId(), ploidy()));

  for (auto contig_variant : genome_variant_map_) {

    auto result = genome_variant->genome_variant_map_.find(contig_variant.first);

    std::shared_ptr<kgl::ContigVariant> diff_contig_ptr;

    if (result != genome_variant_map_.end()) {

      diff_contig_ptr = contig_variant.second->Intersection(result->second);
      genome_intersection->addContigVariant(diff_contig_ptr);

    } else {

      diff_contig_ptr = std::make_shared<kgl::ContigVariant>(contig_variant.first, contig_variant.second->ploidy());
      genome_intersection->addContigVariant(diff_contig_ptr);

    }

    kgl::ExecEnv::log().vinfo("Intersection Contig: {} contains: {} variants",
                              contig_variant.first,
                              diff_contig_ptr->variantCount());

  }

  return genome_intersection;

}


// Returns a GenomeVariant that is the difference return (*this - *genome_variant)
std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::Difference(std::shared_ptr<const kgl::GenomeVariant> genome_variant) const {

  std::shared_ptr<kgl::GenomeVariant> genome_difference(std::make_shared<kgl::GenomeVariant>(genomeId(), ploidy()));

  for (auto contig_variant : genome_variant_map_) {

    auto result = genome_variant->genome_variant_map_.find(contig_variant.first);

    std::shared_ptr<kgl::ContigVariant> diff_contig_ptr;

    if (result != genome_variant_map_.end()) {

      diff_contig_ptr = contig_variant.second->Difference(result->second);
      genome_difference->addContigVariant(diff_contig_ptr);

    } else {

      genome_difference->addContigVariant(contig_variant.second);

    }

    kgl::ExecEnv::log().vinfo("Difference Contig: {} contains: {} variants",
                              contig_variant.first,
                              diff_contig_ptr->variantCount());

  }

  return genome_difference;

}

