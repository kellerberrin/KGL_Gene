//
// Created by kellerberrin on 30/01/18.
//

#ifndef KGL_UPGMA_H
#define KGL_UPGMA_H



#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_statistics_upgma.h"
#include "kgl_sequence_compare.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates the contigs and then compares the mutated contigs using the Myer Hirschberg sequence comparison
// algorthim. This is linear in space - but quadratic in time. Need to find a faster comparison algorithm.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MutatedContigMap = std::map<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>;

class UPGMAContigDistance : public UPGMADistanceNode {

public:

  UPGMAContigDistance(std::shared_ptr<const SequenceDistance> sequence_distance,
                      std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                      std::shared_ptr<const GenomeDatabase> genome_db_ptr) : sequence_distance_(sequence_distance),
                                                                             genome_variant_ptr_(genome_variant_ptr),
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

private:

  std::shared_ptr<const SequenceDistance> sequence_distance_;
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

  UPGMAProteinDistance(std::shared_ptr<const SequenceDistance> sequence_distance,
                       std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                       std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                       const std::string& protein_family) : sequence_distance_(sequence_distance),
                                                            protein_family_(protein_family),
                                                            genome_variant_ptr_(genome_variant_ptr),
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


  constexpr static const char* PROTEIN_FAMILY_WILDCARD = "*";
  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";
  constexpr static const char* SYMBOLIC_Na_H_FAMILY = "NHE";  // A 1 member metabolic family for comparison.
  constexpr static const char* SYMBOLIC_ATP4_FAMILY = "ATP4";  // PfATP4

private:

  std::shared_ptr<const SequenceDistance> sequence_distance_;
  std::string protein_family_;
  std::shared_ptr<const GenomeVariant> genome_variant_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  MutatedProteinMap mutated_proteins_;

  void mutateProteins();
  void getProtein(std::shared_ptr<const GeneFeature> gene_ptr);
  const MutatedProteinMap& getMap() const { return  mutated_proteins_; }
  const GenomeId_t& genomeId() const { return genome_variant_ptr_->genomeId(); }
  const GeneOntology& ontology() const { return genome_db_ptr_->geneOntology(); }

};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates a single gene and compares to other selected genes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UPGMAGeneDistance : public UPGMADistanceNode {

public:

  UPGMAGeneDistance(std::shared_ptr<const SequenceDistance> sequence_distance,
                    std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                    std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                    std::shared_ptr<const GeneFeature> gene_ptr,
                    const std::string& protein_family) : sequence_distance_(sequence_distance),
                                                         protein_family_(protein_family),
                                                         genome_variant_ptr_(genome_variant_ptr),
                                                         gene_ptr_(gene_ptr),
                                                         genome_db_ptr_(genome_db_ptr)  {
    mutateProtein();

  }
  UPGMAGeneDistance(const UPGMAGeneDistance&) = default;
  ~UPGMAGeneDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void write_node(std::ofstream& outfile) const override;
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;

  static bool geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         const std::string& protein_family);

  constexpr static const char* PROTEIN_FAMILY_WILDCARD = "*";
  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";
  constexpr static const char* SYMBOLIC_Na_H_FAMILY = "NHE";  // A 1 member metabolic family for comparison.
  constexpr static const char* SYMBOLIC_ATP4_FAMILY = "ATP4";  // PfATP4


protected:

  std::shared_ptr<const SequenceDistance> sequence_distance_;
  std::shared_ptr<const DNA5SequenceCoding> mutated_protein_;
  std::string protein_family_;
  std::shared_ptr<const GenomeVariant> genome_variant_ptr_;
  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

  void mutateProtein();

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares a single gene between isolates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UPGMAGenePhyloDistance : public UPGMAGeneDistance {

public:

  UPGMAGenePhyloDistance(std::shared_ptr<const SequenceDistance> sequence_distance,
                         std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         std::shared_ptr<const GeneFeature> gene_ptr,
                         const std::string& protein_family) :  UPGMAGeneDistance(sequence_distance,
                                                                                 genome_variant_ptr,
                                                                                 genome_db_ptr,
                                                                                 gene_ptr,
                                                                                 protein_family) {}

  UPGMAGenePhyloDistance(const UPGMAGenePhyloDistance&) = default;
  ~UPGMAGenePhyloDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void write_node(std::ofstream& outfile) const override;
  // Pure Virtual calculates the distance between nodes.

private:

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_UPGMA_H
