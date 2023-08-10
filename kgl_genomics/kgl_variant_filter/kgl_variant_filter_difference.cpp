//
// Created by kellerberrin on 10/08/23.
//

#include "kgl_variant_filter_difference.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the difference between two genomes.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unique_ptr<kgl::GenomeDB> kgl::GenomeDifferenceFilter::applyFilter(const GenomeDB& genome) const {

  auto diff_genome_ptr = std::make_unique<GenomeDB>(genome.genomeId());

  for (auto const& [contig_id, contig_ptr] : reference_ptr_->getMap()) {




  }

  return diff_genome_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the difference between two contigs.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::unique_ptr<kgl::ContigDB> kgl::ContigDifferenceFilter::applyFilter(const ContigDB& contig) const {

  auto diff_contig_ptr = std::make_unique<ContigDB>(contig.contigId());

  for (auto const& [offset, offset_ptr] : contig.getMap()) {

    if (reference_ptr_->getMap().contains(offset)) {

      auto const& [ref_offset, ref_offset_ptr] = *reference_ptr_->getMap().find(offset);
//      OffsetDifferenceFilter offset_diff_filter(ref_offset_ptr);
//      auto diff_offset = offset_diff_filter.applyFilter(*offset_ptr);

    } else {



    }


  }

  return diff_contig_ptr;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the difference between two offsets.
// Only variants that do not occur in the reference are returned.
// There is no homozygous or heterozygous difference information returned.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unique_ptr<kgl::OffsetDB> kgl::OffsetDifferenceFilter::applyFilter(const OffsetDB& offset) const {

  auto diff_offset_ptr = std::make_unique<OffsetDB>();

  std::set<std::string> reference_set;
  for (auto const& ref_variant_ptr : reference_ptr_->getVariantArray()) {

    reference_set.insert(ref_variant_ptr->HGVS());

  }

  std::map<std::string, std::shared_ptr<const Variant>> difference_map;
  for (auto const& variant_ptr : offset.getVariantArray()) {

    std::string variant_key = variant_ptr->HGVS();
    if (not reference_set.contains(variant_key) and not difference_map.contains(variant_key)) {

      auto [insert_iter, result] = difference_map.try_emplace(variant_key, variant_ptr);
      if (not result) {

        ExecEnv::log().error("OffsetDifferenceFilter::applyFilter; Unexpected, cannot insert variant: {} (duplicate)", variant_key);

      }

    }

  }

  for (auto const& [variant_key, variant_ptr] : difference_map) {

    diff_offset_ptr->addVariant(variant_ptr);

  }

  return diff_offset_ptr;

}



