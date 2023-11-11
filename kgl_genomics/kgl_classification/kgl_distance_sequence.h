//
// Created by kellerberrin on 11/11/23.
//

#ifndef KGL_KGL_DISTANCE_SEQUENCE_H
#define KGL_KGL_DISTANCE_SEQUENCE_H


#include "kgl_sequence_distance_impl.h"
#include "kgl_runtime_resource.h"
#include "kgl_variant_db_population.h"
#include "kgl_distance_tree_upgma.h"


namespace kellerberrin::genome {   //  organization level namespace




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares the mutated proteins across a single protein family (var, rifin, etc). Generates a comparison for each genome.
// Amino Local distance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using MutatedProteinMap = std::map<FeatureIdent_t , std::shared_ptr<const AminoSequence>>;

class ProteinDistance : public VirtualDistanceNode {

public:

  ProteinDistance(AminoDistanceMetric sequence_distance,
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

  AminoDistanceMetric sequence_distance_;
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

  GeneDistance(AminoDistanceMetric sequence_distance,
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

  AminoDistanceMetric sequence_distance_;
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

  UPGMAATP4Distance(AminoDistanceMetric sequence_distance,
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

  DNAGeneDistance(LinearDistanceMetric sequence_distance,
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

  LinearDistanceMetric sequence_distance_;
  DNA5SequenceLinear linear_sequence_;

  void getExonSequence();

};


class AminoGeneDistance : public ReferenceGeneDistance {

public:

  AminoGeneDistance(AminoDistanceMetric sequence_distance,
                    std::shared_ptr<const GenomeReference> genome_db_ptr,
                    std::shared_ptr<const GeneFeature> gene_ptr,
                    std::string protein_family)
  : ReferenceGeneDistance(std::move(genome_db_ptr), std::move(gene_ptr), std::move(protein_family)),
    sequence_distance_(sequence_distance) {

    getAminoSequence();

  }

  ~AminoGeneDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  [[nodiscard]] DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

private:

  AminoDistanceMetric sequence_distance_;
  AminoSequence amino_sequence_;

  void  getAminoSequence();

};



}   // end namespace genome



#endif //KGL_KGL_DISTANCE_SEQUENCE_H
