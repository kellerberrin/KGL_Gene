//
// Created by kellerberrin on 14/08/18.
//

#ifndef KGL_UPGMA_UNPHASED_H
#define KGL_UPGMA_UNPHASED_H


#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_statistics_upgma.h"
#include "kgl_sequence_distance.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class UPGMAUnphasedDistance : public VirtualDistanceNode {

public:

  UPGMAUnphasedDistance(std::shared_ptr<const UnphasedGenome> genome_variant_ptr) : genome_variant_ptr_(genome_variant_ptr) {}
  UPGMAUnphasedDistance(const UPGMAUnphasedDistance&) = default;
  ~UPGMAUnphasedDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

private:



  std::shared_ptr<const UnphasedGenome> genome_variant_ptr_;

};


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_UPGMA_UNPHASED_H
