//
// Created by kellerberrin on 30/01/18.
//

#ifndef KGL_UPGMA_NODE_H
#define KGL_UPGMA_NODE_H


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

class UPGMAContigDistance : public VirtualDistanceNode {

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
  void writeNode(std::ostream& outfile) const override { outfile << genome_variant_ptr_->genomeId(); }
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

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

class UPGMAProteinDistance : public VirtualDistanceNode {

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

class UPGMAGeneDistance : public VirtualDistanceNode {

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
  void writeNode(std::ostream& outfile) const override;
  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

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

  }

  UPGMAATP4Distance(const UPGMAATP4Distance&) = default;
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

  ReferenceGeneDistance(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                        std::shared_ptr<const GeneFeature> gene_ptr,
                        const std::string& protein_family) : genome_db_ptr_(genome_db_ptr),
                                                             gene_ptr_(gene_ptr),
                                                             protein_family_(protein_family) {


  }

  explicit ReferenceGeneDistance(const ReferenceGeneDistance&) = default;
  ~ReferenceGeneDistance() override = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  void writeNode(std::ostream& outfile) const override;


  constexpr static const char* SYMBOLIC_VAR_FAMILY = "VAR";
  constexpr static const char* SYMBOLIC_RIFIN_FAMILY = "RIF";
  constexpr static const char* SYMBOLIC_MAURER_FAMILY = "MC-2TM";

  static bool geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         const std::string& protein_family);


protected:

  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::string protein_family_;

  std::shared_ptr<const CodingSequence>  getCodingSequence();

};



class DNAGeneDistance : public ReferenceGeneDistance {

public:

  DNAGeneDistance(std::shared_ptr<const LocalDNASequenceDistance> sequence_distance,
                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                  std::shared_ptr<const GeneFeature> gene_ptr,
                  const std::string& protein_family)
                  : ReferenceGeneDistance(genome_db_ptr, gene_ptr, protein_family),
                  sequence_distance_(sequence_distance) {

    getIntronSequence(); // Distance of the intron sequence.
//    getExonSequence(); // Distance of the

  }

  explicit DNAGeneDistance(const DNAGeneDistance&) = default;
  ~DNAGeneDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;


private:

  std::shared_ptr<const LocalDNASequenceDistance> sequence_distance_;
  std::shared_ptr<const DNA5SequenceLinear> sequence_ptr_;

  void  getExonSequence();
  void  getIntronSequence();

};


class AminoGeneDistance : public ReferenceGeneDistance {

public:

  AminoGeneDistance(std::shared_ptr<const LocalSequenceDistance> sequence_distance,
                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                  std::shared_ptr<const GeneFeature> gene_ptr,
                  const std::string& protein_family)
  : ReferenceGeneDistance(genome_db_ptr, gene_ptr, protein_family) {

    sequence_distance_ = std::dynamic_pointer_cast<const LocalAminoSequenceDistance>(sequence_distance);

    if (not sequence_distance_) {

      ExecEnv::log().critical("AminoGeneDistance::AminoGeneDistance; distance metric: {} is not a superclass of 'LocalAminoSequenceDistance'",
                              sequence_distance->distanceType());

    }

    getAminoSequence();

  }

  explicit AminoGeneDistance(const AminoGeneDistance&) = default;
  ~AminoGeneDistance() override = default;

  // Pure Virtual calculates the distance between nodes.
  DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const override;

private:

  std::shared_ptr<const LocalAminoSequenceDistance> sequence_distance_;
  std::shared_ptr<const AminoSequence> sequence_ptr_;

  void  getAminoSequence();

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_UPGMA_NODE_H
