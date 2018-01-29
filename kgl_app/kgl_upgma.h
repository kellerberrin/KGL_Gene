//
// Created by kellerberrin on 30/01/18.
//

#ifndef KGL_UPGMA_H
#define KGL_UPGMA_H



#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_statistics_upgma.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates the contigs and then compares the mutated contigs using the Myer Hirschberg sequence comparison
// algorthim. This is linear in space - but quadratic in time. Need to find a faster comparison algorithm.
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates genome proteins and then compares the mutated proteins using the Myer Hirschberg sequence comparison
// algorthim. Myer Hirschberg is linear in space - but quadratic in time. Need to find a faster comparison algorithm.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MutatedProteinMap = std::map<FeatureIdent_t , std::shared_ptr<const DNA5SequenceCoding>>;

class UPGMAProteinDistance : public UPGMADistanceNode {

public:

  UPGMAProteinDistance(std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                       std::shared_ptr<const GenomeDatabase> genome_db_ptr) : genome_variant_ptr_(genome_variant_ptr),
                                                                              genome_db_ptr_(genome_db_ptr) {
    mutateProteins();

  }
  UPGMAProteinDistance(const UPGMAProteinDistance&) = default;
  ~UPGMAProteinDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void write_node(std::ofstream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;

  // Generate distance nodes from a population.
  static std::shared_ptr<NodeVector<const UPGMADistanceNode>> upgma_matrix(std::shared_ptr<const PopulationVariant> pop_variant_ptr,
                                                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr);

  const MutatedProteinMap& getMap() const { return  mutated_proteins_; }
  const GenomeId_t& genomeId() const { return genome_variant_ptr_->genomeId(); }
  const GeneOntology& ontology() const { return genome_db_ptr_->geneOntology(); }

private:

  std::shared_ptr<const GenomeVariant> genome_variant_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  MutatedProteinMap mutated_proteins_;

  void mutateProteins();

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance using protein families
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UPGMAFamilyDistance : public UPGMAProteinDistance{

public:

  UPGMAFamilyDistance(std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                      std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                      const std::string& family_code) : UPGMAProteinDistance(genome_variant_ptr, genome_db_ptr),
                                                        family_code_(family_code) {}

  UPGMAFamilyDistance(const UPGMAFamilyDistance&) = default;
  ~UPGMAFamilyDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;

  // Generate distance nodes from a population.
  static std::shared_ptr<NodeVector<const UPGMADistanceNode>> upgma_matrix(std::shared_ptr<const PopulationVariant> pop_variant_ptr,
                                                                           std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                           const std::string& family_code);

  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";

private:

  std::string family_code_;

};





}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_UPGMA_H
