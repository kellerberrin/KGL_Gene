//
// Created by kellerberrin on 30/01/18.
//

#ifndef KGL_UPGMA_H
#define KGL_UPGMA_H


#include <kgl_sequence_distance.h>
#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_statistics_upgma.h"
#include "kgl_sequence_distance.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates the contigs and then compares the mutated contigs between genomes. Global DNA distance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MutatedContigMap = std::map<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>;

class UPGMAContigDistance : public UPGMADistanceNode {

public:

  UPGMAContigDistance(std::shared_ptr<const GlobalDNASequenceDistance> sequence_distance,
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
  void writeNode(std::ofstream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;

private:

  std::shared_ptr<const GlobalDNASequenceDistance> sequence_distance_;
  std::shared_ptr<const GenomeVariant> genome_variant_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  MutatedContigMap mutated_contigs_;

  void mutateContigs();

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares the mutated proteins across a single protein family (var, rifin, etc). Generates a comparison for each genome.
// Amino Local distance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MutatedProteinMap = std::map<FeatureIdent_t , std::shared_ptr<const AminoSequence>>;

class UPGMAProteinDistance : public UPGMADistanceNode {

public:

  UPGMAProteinDistance(std::shared_ptr<const LocalAminoSequenceDistance> sequence_distance,
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
  void writeNode(std::ofstream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;


  constexpr static const char* PROTEIN_FAMILY_WILDCARD = "*";
  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";
  constexpr static const char* SYMBOLIC_Na_H_FAMILY = "NHE";  // A 1 member metabolic family for comparison.
  constexpr static const char* SYMBOLIC_ATP4_FAMILY = "ATP4";  // PfATP4

private:

  std::shared_ptr<const LocalAminoSequenceDistance> sequence_distance_;
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
// For all genes that belong to a family, compares the same gene from different genomes.
// Generates a comparison for each gene. Local or Global Amino
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UPGMAGeneDistance : public UPGMADistanceNode {

public:

  UPGMAGeneDistance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
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
  void writeNode(std::ofstream& outfile) const override;
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

  std::shared_ptr<const AminoSequenceDistance> sequence_distance_;
  std::shared_ptr<const AminoSequence> mutated_protein_;
  std::string protein_family_;
  std::shared_ptr<const GenomeVariant> genome_variant_ptr_;
  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

  void mutateProtein();

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares a single gene between isolate genomes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UPGMAATP4Distance : public UPGMAGeneDistance {

public:

  UPGMAATP4Distance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                         std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         std::shared_ptr<const GeneFeature> gene_ptr,
                         const std::string& protein_family) :  UPGMAGeneDistance(sequence_distance,
                                                                                 genome_variant_ptr,
                                                                                 genome_db_ptr,
                                                                                 gene_ptr,
                                                                                 protein_family) {

    zero_metric_ptr_ = std::make_shared<const LevenshteinGlobal>();

  }

  UPGMAATP4Distance(const UPGMAATP4Distance&) = default;
  ~UPGMAATP4Distance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ofstream& outfile) const override;
  // virtual function ensures that identical nodes have zero distance.
  DistanceType_t zeroDistance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;

private:

  std::shared_ptr<const AminoSequenceDistance> zero_metric_ptr_;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares refernce (unmutated) genes across a single protein family (var, rifin, etc).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ReferenceGeneDistance : public UPGMADistanceNode {

public:

  ReferenceGeneDistance(std::shared_ptr<const DNASequenceDistance> sequence_distance,
                        std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                        std::shared_ptr<const GeneFeature> gene_ptr,
                        const std::string& protein_family) : sequence_distance_(sequence_distance),
                                                             genome_db_ptr_(genome_db_ptr),
                                                             gene_ptr_(gene_ptr),
                                                             protein_family_(protein_family) {
    getSequence();

  }

  explicit ReferenceGeneDistance(const ReferenceGeneDistance&) = default;
  ~ReferenceGeneDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ofstream& outfile) const override;
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const override;


  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";

  static bool geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         const std::string& protein_family);

private:

  std::shared_ptr<const DNASequenceDistance> sequence_distance_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::string protein_family_;

  std::shared_ptr<const DNA5SequenceLinear> sequence_ptr_;

  void getSequence();

};







}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_UPGMA_H
