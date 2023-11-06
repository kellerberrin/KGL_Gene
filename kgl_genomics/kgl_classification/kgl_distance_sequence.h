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
                  std::shared_ptr<const GenomeDB> genome_db_ptr,
                  std::shared_ptr<const GenomeReference> genome_ref_ptr,
                  std::string protein_family) : sequence_distance_(std::move(sequence_distance)),
                                                genome_db_ptr_(std::move(genome_db_ptr)),
                                                genome_ref_ptr_(std::move(genome_ref_ptr)),
                                                protein_family_(std::move(protein_family)) {
    mutateProteins();

  }
  ProteinDistance(const ProteinDistance&) = default;
  ~ProteinDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override { outfile << genome_db_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  [[nodiscard]] DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;


private:

  std::shared_ptr<const AminoSequenceDistance> sequence_distance_;
  std::shared_ptr<const GenomeDB> genome_db_ptr_;
  std::shared_ptr<const GenomeReference> genome_ref_ptr_;
  std::string protein_family_;

  MutatedProteinMap mutated_proteins_;

  void mutateProteins();
  void getProtein(std::shared_ptr<const GeneFeature> gene_ptr);
  [[nodiscard]] const MutatedProteinMap& getMap() const { return mutated_proteins_; }
  [[nodiscard]] const GenomeId_t& genomeId() const { return genome_db_ptr_->genomeId(); }
  [[nodiscard]] const GeneOntology& ontology() const { return genome_ref_ptr_->geneOntology(); }

};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// For all genes that belong to a family, compares the same gene from different genomes.
// Generates a comparison for each gene. Local or Global Amino
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneDistance : public VirtualDistanceNode {

public:

  GeneDistance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
               std::shared_ptr<const GenomeDB> genome_db_ptr,
               std::shared_ptr<const GenomeReference> genome_ref_ptr,
               std::shared_ptr<const GeneFeature> gene_ptr,
               std::string protein_family) : sequence_distance_(std::move(sequence_distance)),
                                             genome_db_ptr_(std::move(genome_db_ptr)),
                                             genome_ref_ptr_(std::move(genome_ref_ptr)),
                                             gene_ptr_(std::move(gene_ptr)),
                                             protein_family_(std::move(protein_family)) {

    mutateProtein();

  }
  ~GeneDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override;
  // Pure Virtual calculates the distance between nodes.
  [[nodiscard]] DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

  [[nodiscard]] static bool geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                                       std::shared_ptr<const GenomeReference> genome_db_ptr,
                                       const std::string& protein_family);


protected:

  std::shared_ptr<const AminoSequenceDistance> sequence_distance_;
  std::shared_ptr<const GenomeDB> genome_db_ptr_;
  std::shared_ptr<const GenomeReference> genome_ref_ptr_;
  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::string protein_family_;

  AminoSequence mutated_protein_;

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

  [[nodiscard]] std::shared_ptr<const TranscriptionSequence>  getCodingSequence();

};



class DNAGeneDistance : public ReferenceGeneDistance {

public:

  DNAGeneDistance(std::shared_ptr<const LinearDNASequenceDistance> sequence_distance,
                  std::shared_ptr<const GenomeReference> genome_db_ptr,
                  std::shared_ptr<const GeneFeature> gene_ptr,
                  std::string protein_family)
                  : ReferenceGeneDistance(std::move(genome_db_ptr), std::move(gene_ptr), std::move(protein_family)),
                  sequence_distance_(std::move(sequence_distance)) {

    getExonSequence(); // Distance of the

  }
  ~DNAGeneDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  [[nodiscard]] DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;


private:

  std::shared_ptr<const LinearDNASequenceDistance> sequence_distance_;
  DNA5SequenceLinear linear_sequence_;

  void getExonSequence();

};


class AminoGeneDistance : public ReferenceGeneDistance {

public:

  AminoGeneDistance(std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                    std::shared_ptr<const GenomeReference> genome_db_ptr,
                    std::shared_ptr<const GeneFeature> gene_ptr,
                    std::string protein_family)
  : ReferenceGeneDistance(std::move(genome_db_ptr), std::move(gene_ptr), std::move(protein_family)) {

    sequence_distance_ = std::dynamic_pointer_cast<const AminoSequenceDistance>(sequence_distance);
    if (not sequence_distance_) {

      ExecEnv::log().critical("AminoGeneDistance::AminoGeneDistance; distance metric: {} is not a superclass of 'LocalAminoSequenceDistance'",
                              sequence_distance->distanceType());

    }

    getAminoSequence();

  }

  ~AminoGeneDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  [[nodiscard]] DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

private:

  std::shared_ptr<const AminoSequenceDistance> sequence_distance_;
  AminoSequence amino_sequence_;

  void  getAminoSequence();

};





}   // end namespace genome


#endif //KGL_UPGMA_NODE_H
