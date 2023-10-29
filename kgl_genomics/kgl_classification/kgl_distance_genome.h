//
// Created by kellerberrin on 14/08/18.
//

#ifndef KGL_UPGMA_UNPHASED_H
#define KGL_UPGMA_UNPHASED_H


#include "kgl_runtime_resource.h"
#include "kgl_variant_db_population.h"
#include "kgl_distance_tree_upgma.h"
#include "kgl_sequence_distance.h"


namespace kellerberrin::genome {   //  organization level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple distance measures the number of variant not co-occurring between different genomes.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GenomeDistance : public VirtualDistanceNode {

public:

  GenomeDistance(std::shared_ptr<const GenomeDB> genome_variant_ptr) : genome_variant_ptr_(std::move(genome_variant_ptr)) {}
  GenomeDistance(const GenomeDistance&) = default;
  ~GenomeDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  [[nodiscard]] DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

private:

  std::shared_ptr<const GenomeDB> genome_variant_ptr_;

};


}   // end namespace

#endif //KGL_UPGMA_UNPHASED_H
