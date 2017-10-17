// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// Created by kellerberrin on 13/10/17.
//

#include "kgl_variant.h"
#include "kgl_exec_env.h"
#include "kgl_patterns.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


size_t kgl::ContigVariant::filterVariants(const kgl::VariantFilter& filter) {

  // Inverts the bool returned by filterVariant(filter) because the delete pattern expects bool true for deletion.
  auto predicate = [&](const OffsetVariantMap::const_iterator& it) { return not it->second->filterVariant(filter); };
  predicateIterableDelete(offset_variant_map_,  predicate);

  return offset_variant_map_.size();

}


void kgl::ContigVariant::addVariant(ContigOffset_t contig_offset, std::shared_ptr<const kgl::Variant>& variant_ptr) {

  offset_variant_map_.insert(std::make_pair(contig_offset, variant_ptr));

}


bool kgl::ContigVariant::isElement(const Variant& variant) const {

  auto result = offset_variant_map_.equal_range(variant.contigOffset());

  for (auto it = result.first; it != result.second; ++it) {

    if (*(it->second) == variant) return true;

  }

  return false;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::GenomeVariant::addContigVariant(std::shared_ptr<kgl::ContigVariant>& contig_variant) {

  auto result = genome_variant_map_.insert(std::make_pair(contig_variant->contigId(), contig_variant));

  return result.second;

}

void kgl::GenomeVariant::filterVariants(const kgl::VariantFilter& filter) {

  ExecEnv::log().info("Applying filter: {}", filter.filterName());
  for (const auto& contig_variant : genome_variant_map_) {

    size_t v_size = contig_variant.second->filterVariants(filter);
    ExecEnv::log().info("Contig: {} has: {} filtered variants", contig_variant.first, v_size);

  }

}

bool kgl::GenomeVariant::isElement(const Variant& variant) const {

  auto result = genome_variant_map_.find(variant.contigId());

  if (result == genome_variant_map_.end()) {

    return false;

  }

  return result->second->isElement(variant);

}


std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::Union(std::shared_ptr<const kgl::GenomeVariant> genome_variant_ptr) const {

  std::shared_ptr<kgl::GenomeVariant> genome_union(std::make_shared<kgl::GenomeVariant>());

  return genome_union;

}

std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::Intersection(std::shared_ptr<const kgl::GenomeVariant> genome_variant) const {

  std::shared_ptr<kgl::GenomeVariant> genome_intersection(std::make_shared<kgl::GenomeVariant>());

  return genome_intersection;

}

std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::Difference(std::shared_ptr<const kgl::GenomeVariant> genome_variant) const {

  std::shared_ptr<kgl::GenomeVariant> genome_difference(std::make_shared<kgl::GenomeVariant>());

  return genome_difference;

}
