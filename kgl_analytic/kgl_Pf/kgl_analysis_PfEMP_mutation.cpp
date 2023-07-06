//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_analysis_PfEMP_mutation.h"

namespace kgl = kellerberrin::genome;

void kgl::MutateGenes::initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [gene_offset, gene_ptr] : contig_ptr->getGeneMap()) {

      auto [iter, result] = gene_contig_map_.try_emplace(gene_ptr, contig_id);
      if (not result) {

        ExecEnv::log().warn("MutateGenes::initializeGeneContigMap; Unable to insert duplicate gene: {}", gene_ptr->id());

      }

    }

  }

}


std::vector<std::shared_ptr<const kgl::GeneFeature>> kgl::MutateGenes::contigGenes(const ContigId_t& contig_id) const {

  std::vector<std::shared_ptr<const GeneFeature>> gene_vector;

  for (auto const& [gene_ptr, gene_contig_id] : gene_contig_map_) {

    if (gene_contig_id == contig_id) {

      gene_vector.push_back(gene_ptr);

    }

  }

  return gene_vector;

}
