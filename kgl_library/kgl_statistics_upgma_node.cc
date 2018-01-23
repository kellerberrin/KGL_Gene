//
// Created by kellerberrin on 23/01/18.
//

#include "kgl_statistics_upgma_node.h"


namespace kgl = kellerberrin::genome;


kgl::DistanceType_t kgl::UPGMAContigDistance::distance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAContigDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAContigDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not up-cast node pointer");
    return 1.0;

  }

  DistanceType_t total_distance = 0;
  for (auto contig : mutated_contigs_) {

    auto result = node_ptr->mutated_contigs_.find(contig.first);

    if (result != node_ptr->mutated_contigs_.end()) {

      ExecEnv::log().info("distance(), Comparing Genome: {}, Contig: {} with Genome: {}",
                          genome_variant_ptr_->genomeId(), contig.first, node_ptr->genome_variant_ptr_->genomeId());
      CompareScore_t contig_score = contig.second->compareMyerHirschberg(result->second);
      total_distance += static_cast<DistanceType_t>(contig_score);
      ExecEnv::log().info("distance(), Calculated distance: {}", contig_score);

    } else {

      ExecEnv::log().error("distance(), Unexpected error, could not find contig: {}", contig.first);

    }

  }

  return total_distance;

}


std::shared_ptr<kgl::NodeVector<const kgl::UPGMADistanceNode>>
kgl::UPGMAContigDistance::upgma_matrix(std::shared_ptr<const PopulationVariant> pop_variant_ptr,
                                       std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  std::shared_ptr<NodeVector<const UPGMADistanceNode>> node_vector_ptr(std::make_shared<NodeVector<const UPGMADistanceNode>>());

  for (auto genome : pop_variant_ptr->getMap()) {

    std::shared_ptr<const UPGMAContigDistance> distance_ptr(std::make_shared<const UPGMAContigDistance>(genome.second, genome_db_ptr));
    std::shared_ptr<PhyloNode<const UPGMADistanceNode>> phylo_node_ptr(std::make_shared<PhyloNode<const UPGMADistanceNode>>(distance_ptr));
    node_vector_ptr->push_back(phylo_node_ptr);

  }

  return node_vector_ptr;

}


void kgl::UPGMAContigDistance::mutateContigs() {

  std::shared_ptr<const DNA5SequenceContig> reference_contig_ptr;
  std::shared_ptr<DNA5SequenceContig> mutant_contig_ptr;
  for (auto contig : genome_db_ptr_->getMap()) {

    if (not genome_variant_ptr_->mutantContig(contig.first, genome_db_ptr_, reference_contig_ptr, mutant_contig_ptr)) {

      ExecEnv::log().error("Unexpected error mutating contig; genome: {} , contig: {}", genome_variant_ptr_->genomeId(), contig.first);
      // Fail gracefully, just insert the reference contig.
      mutated_contigs_.insert(std::pair<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>(contig.first, contig.second->sequence_ptr()));

    } else {

      mutated_contigs_.insert(std::pair<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>(contig.first, mutant_contig_ptr));

    }

  }

}
