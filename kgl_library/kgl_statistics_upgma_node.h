//
// Created by kellerberrin on 23/01/18.
//




#ifndef KGL_STATISTICS_UPGMA_NODE_H
#define KGL_STATISTICS_UPGMA_NODE_H


#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_statistics_upgma.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual Distance class for the UPGMA graph.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UPGMADistanceNode {

public:

  UPGMADistanceNode() = default;
  UPGMADistanceNode(const UPGMADistanceNode&) = default;
  virtual ~UPGMADistanceNode() = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  virtual void write_node(std::ofstream& outfile) const = 0;
  // Pure Virtual calculates the distance between nodes.
  virtual DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const = 0;


private:

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Default Distance class for the UPGMA graph.
// Mutates the contigs and then compares the mutated contigs using the Myer Hirschberg sequence comparison
// algorthim (this is linear in space and can match large sequences such as contigs).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MutatedContigMap = std::map<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>;

class UPGMAContigDistance : public UPGMADistanceNode {

public:

  UPGMAContigDistance(std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                      std::shared_ptr<const GenomeDatabase> genome_db_ptr) : genome_variant_ptr_(genome_variant_ptr),
                                                                             genome_db_ptr_(genome_db_ptr) {
    mutateContigs();

  }
  UPGMAContigDistance(const UPGMAContigDistance&) = default;
  ~UPGMAContigDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void write_node(std::ofstream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;

  // Generate distance nodes from a population.
  static std::shared_ptr<NodeVector<const UPGMADistanceNode>> upgma_matrix(std::shared_ptr<const PopulationVariant> pop_variant_ptr,
                                                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr);


private:

  std::shared_ptr<const GenomeVariant> genome_variant_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  MutatedContigMap mutated_contigs_;

  void mutateContigs();

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_STATISTICS_UPGMA_NODE_H
