//
// Created by kellerberrin on 30/01/18.
//

#ifndef KGL_UPGMA_NODE_H
#define KGL_UPGMA_NODE_H


#include "kgl_sequence_distance.h"
#include "kgl_runtime_resource.h"
#include "kgl_variant_db_population.h"
#include "kgl_distance_tree_upgma.h"
#include "kgl_sequence_distance.h"


namespace kellerberrin::genome {   //  organization level namespace




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares the mutated proteins across a single protein family (var, rifin, etc). Generates a comparison for each genome.
// Amino Local distance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MutatedProteinMap = std::map<FeatureIdent_t , std::shared_ptr<const AminoSequence>>;

class ProteinDistance : public VirtualDistanceNode {

public:

  ProteinDistance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                  std::shared_ptr<const GenomeDB> genome_variant_ptr,
                  std::shared_ptr<const GenomeReference> genome_db_ptr,
                  const std::string& protein_family) : sequence_distance_(sequence_distance),
                                                            protein_family_(protein_family),
                                                            genome_variant_ptr_(genome_variant_ptr),
                                                            genome_db_ptr_(genome_db_ptr) {
    mutateProteins();

  }
  ProteinDistance(const ProteinDistance&) = default;
  ~ProteinDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;


  constexpr static const char* PROTEIN_FAMILY_WILDCARD = "*";
  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";
  constexpr static const char* SYMBOLIC_Na_H_FAMILY = "NHE";  // A 1 member metabolic family for comparison.
  constexpr static const char* SYMBOLIC_ATP4_FAMILY = "ATP4";  // PfATP4

private:

  std::shared_ptr<const AminoSequenceDistance> sequence_distance_;
  std::string protein_family_;
  std::shared_ptr<const GenomeDB> genome_variant_ptr_;
  std::shared_ptr<const GenomeReference> genome_db_ptr_;
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

class GeneDistance : public VirtualDistanceNode {

public:

  GeneDistance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
               std::shared_ptr<const GenomeDB> genome_variant_ptr,
               std::shared_ptr<const GenomeReference> genome_db_ptr,
               std::shared_ptr<const GeneFeature> gene_ptr,
               const std::string& protein_family) : sequence_distance_(sequence_distance),
                                                         protein_family_(protein_family),
                                                         genome_variant_ptr_(genome_variant_ptr),
                                                         gene_ptr_(gene_ptr),
                                                         genome_db_ptr_(genome_db_ptr)  {
    mutateProtein();

  }
  ~GeneDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override;
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

  static bool geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                         std::shared_ptr<const GenomeReference> genome_db_ptr,
                         const std::string& protein_family);

  constexpr static const char* PROTEIN_FAMILY_WILDCARD = "*";
  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";
  constexpr static const char* SYMBOLIC_Na_H_FAMILY = "NHE";  // A 1 member metabolic family for comparison.
  constexpr static const char* SYMBOLIC_ATP4_FAMILY = "ATP4";  // PfATP4


protected:

  std::shared_ptr<const AminoSequenceDistance> sequence_distance_;
  AminoSequence mutated_protein_;
  std::string protein_family_;
  std::shared_ptr<const GenomeDB> genome_variant_ptr_;
  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const GenomeReference> genome_db_ptr_;

  void mutateProtein();

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares a single gene between isolate genomes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UPGMAATP4Distance : public GeneDistance {

public:

  UPGMAATP4Distance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                         std::shared_ptr<const GenomeDB> genome_variant_ptr,
                         std::shared_ptr<const GenomeReference> genome_db_ptr,
                         std::shared_ptr<const GeneFeature> gene_ptr,
                         const std::string& protein_family) : GeneDistance(sequence_distance,
                                                                           genome_variant_ptr,
                                                                           genome_db_ptr,
                                                                           gene_ptr,
                                                                           protein_family) {

  }

  ~UPGMAATP4Distance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override;

private:

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares reference (unmutated) genes across a single protein family (var, rifin, etc).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ReferenceGeneDistance : public VirtualDistanceNode {

public:

  ReferenceGeneDistance(std::shared_ptr<const GenomeReference> genome_db_ptr,
                        std::shared_ptr<const GeneFeature> gene_ptr,
                        std::string protein_family) : genome_db_ptr_(std::move(genome_db_ptr)),
                                                      gene_ptr_(std::move(gene_ptr)),
                                                      protein_family_(std::move(protein_family)) {}

  ReferenceGeneDistance(const ReferenceGeneDistance&) = default;
  ~ReferenceGeneDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override;


protected:

  std::shared_ptr<const GenomeReference> genome_db_ptr_;
  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::string protein_family_;

  std::shared_ptr<const TranscriptionSequence>  getCodingSequence();

};



class DNAGeneDistance : public ReferenceGeneDistance {

public:

  DNAGeneDistance(std::shared_ptr<const LinearDNASequenceDistance> sequence_distance,
                  std::shared_ptr<const GenomeReference> genome_db_ptr,
                  std::shared_ptr<const GeneFeature> gene_ptr,
                  const std::string& protein_family)
                  : ReferenceGeneDistance(genome_db_ptr, gene_ptr, protein_family),
                  sequence_distance_(sequence_distance) {

    getExonSequence(); // Distance of the

  }

  explicit DNAGeneDistance(const DNAGeneDistance&) = default;
  ~DNAGeneDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;


private:

  std::shared_ptr<const LinearDNASequenceDistance> sequence_distance_;
  DNA5SequenceLinear linear_sequence_;

  void  getExonSequence();

};


class AminoGeneDistance : public ReferenceGeneDistance {

public:

  AminoGeneDistance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                  std::shared_ptr<const GenomeReference> genome_db_ptr,
                  std::shared_ptr<const GeneFeature> gene_ptr,
                  const std::string& protein_family)
  : ReferenceGeneDistance(genome_db_ptr, gene_ptr, protein_family) {

    sequence_distance_ = std::dynamic_pointer_cast<const AminoSequenceDistance>(sequence_distance);

    if (not sequence_distance_) {

      ExecEnv::log().critical("AminoGeneDistance::AminoGeneDistance; distance metric: {} is not a superclass of 'LocalAminoSequenceDistance'",
                              sequence_distance->distanceType());

    }

    getAminoSequence();

  }

  ~AminoGeneDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

private:

  std::shared_ptr<const AminoSequenceDistance> sequence_distance_;
  AminoSequence amino_sequence_;

  void  getAminoSequence();

};



}   // end namespace genome


#endif //KGL_UPGMA_NODE_H
