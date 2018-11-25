//
// Created by kellerberrin on 16/12/17.
//

#ifndef KGL_STATISTICS_UPGMA_H
#define KGL_STATISTICS_UPGMA_H

#include <memory>
#include <map>
#include <vector>
#include <fstream>

#include "kgl_exec_env.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_population.h"
#include "kgl_utility.h"
#include "kgl_sequence_distance.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual Distance class implemented elsewhere that actually calculates the UPGMA distance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using DistanceType_t = double;
class PhyloNode;  // fwd.
using PhyloNodeVector = std::vector<std::shared_ptr<PhyloNode>>;

class UPGMADistanceNode {

public:

  UPGMADistanceNode() = default;
  UPGMADistanceNode(const UPGMADistanceNode&) = default;
  virtual ~UPGMADistanceNode() = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  virtual void writeNode(std::ofstream& outfile) const = 0;
  // Pure Virtual calculates the distance between nodes.
  virtual DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const = 0;
  // Pure Virtual calculates the zero distance between nodes.
  // This function is only re-defined and used if the distance metric needs to set a particular
  // condition for a zero distance. Most distance metrics will not need to re-define this function.
  virtual DistanceType_t zeroDistance(std::shared_ptr<const UPGMADistanceNode>) const { return 1.0; }

private:

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used by the classification functions.
using OutNodes = std::multimap<DistanceType_t , std::shared_ptr<PhyloNode>>;
class PhyloNode {

public:

  explicit PhyloNode(std::shared_ptr<const UPGMADistanceNode> leaf) : leaf_(leaf), distance_(0) {}
  ~PhyloNode() = default;

  void addOutNode(std::shared_ptr<PhyloNode> node) {
    out_nodes_.insert(std::pair<DistanceType_t , std::shared_ptr<PhyloNode>>(node->distance(), node));
  }

  DistanceType_t distance() const { return distance_; }
  void distance(DistanceType_t update) { distance_ = update; }

  std::shared_ptr<const UPGMADistanceNode> leaf() const { return leaf_; }
  const OutNodes& getMap() const { return out_nodes_; }

  // Recursively counts the total number of leaf nodes.
  bool isLeaf() const { return getMap().empty(); }
  size_t leafNodeCount() const;

private:

  std::shared_ptr<const UPGMADistanceNode> leaf_;
  DistanceType_t distance_;
  OutNodes out_nodes_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance matrix. Implements the PIMPL pattern to isolate Boost functionality.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DistanceMatrix {

public:

  explicit DistanceMatrix(size_t matrix_size);
  explicit DistanceMatrix(const DistanceMatrix& copy);
  virtual ~DistanceMatrix();  // Do not use the default destructor, see PIMPL fwd decl below.


  DistanceType_t minimum(size_t& i, size_t& j) const;
  DistanceType_t maximum(size_t& i, size_t& j) const;
  void reduce(size_t i, size_t j);
  void setDistance(size_t i, size_t j, DistanceType_t distance);
  DistanceType_t getDistance(size_t i, size_t j) const;

  virtual size_t getLeafCount(size_t) const { return 1; }

private:

  class BoostDistanceMatrix;       // Forward declaration of the boost strict diagonal implementation class
  std::unique_ptr<BoostDistanceMatrix> diagonal_impl_ptr_;    // PIMPL

  size_t size() const;
  void resize(size_t new_size);

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//using MatrixNode = PhyloNode;

class UPGMAMatrix : DistanceMatrix {

public:

  explicit UPGMAMatrix(std::shared_ptr<PhyloNodeVector> node_vector_ptr) : DistanceMatrix(node_vector_ptr->size()),
                                                                           node_vector_ptr_(node_vector_ptr) {
    initializeDistance();

  }
  ~UPGMAMatrix() override = default;


  void calculateReduce();

  bool writeNewick(const std::string& file_name) const;

private:

  DistanceType_t distance(std::shared_ptr<PhyloNode> row_node, std::shared_ptr<PhyloNode> column_node) const;
  DistanceType_t zeroDistance(std::shared_ptr<PhyloNode> row_node, std::shared_ptr<PhyloNode> column_node) const;
  void initializeDistance();
  virtual void normalizeDistance();
  void rescaleDistance();
  void identityZeroDistance();
  bool reduceNode(size_t row, size_t column, DistanceType_t minimum);
  void writeNode(std::shared_ptr<PhyloNode> node, std::ofstream& newick_file) const;
  size_t getLeafCount(size_t leaf_idx) const override;

  std::shared_ptr<PhyloNodeVector> node_vector_ptr_;


};


// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a population tree.
// We are comparing contigs across genomes. The distance metric will be based on the difference in
// the unphased variants held for each contig.
template<typename T, typename... Args>
void UPGMAUnphasedTree(const std::string& newick_file,
                       std::shared_ptr<const UnphasedPopulation> pop_unphased_ptr,
                       Args... args) {

  std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());

  for (auto genome : pop_unphased_ptr->getMap()) {

    std::shared_ptr<T> distance_ptr(std::make_shared<T>(genome.second, args...));
    std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
    node_vector_ptr->push_back(phylo_node_ptr);

  }

  UPGMAMatrix upgma_matrix(node_vector_ptr);

  upgma_matrix.calculateReduce();

  upgma_matrix.writeNewick(newick_file);

}


// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a population tree.
// We are comparing contigs across genomes, so a global dna distance metric is required.
template<typename T, typename... Args>
void UPGMAPopulationTree(const std::string& newick_file,
                         std::shared_ptr<const GlobalDNASequenceDistance> sequence_distance,
                         std::shared_ptr<const PhasedPopulation> pop_variant_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         Args... args) {

  std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());

  for (auto genome : pop_variant_ptr->getMap()) {

    std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance, genome.second, genome_db_ptr, args...));
    std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
    node_vector_ptr->push_back(phylo_node_ptr);

  }

  UPGMAMatrix upgma_matrix(node_vector_ptr);

  upgma_matrix.calculateReduce();

  upgma_matrix.writeNewick(newick_file);

}


// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a series of Gene trees.
// We are comparing dissimilar gene types, so only a local Amino distance metric should be used
template<typename T, typename... Args>
void UPGMAGeneTree(const std::string& path,
                   const std::string& newick_file,
                   std::shared_ptr<const LocalAminoSequenceDistance> sequence_distance,
                   std::shared_ptr<const PhasedPopulation> pop_variant_ptr,
                   std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                   const std::string& protein_family,
                   Args... args) {


  for (auto genome_variant : pop_variant_ptr->getMap()) {

    std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());

    for (auto contig : genome_db_ptr->getMap()) {

      for (auto gene : contig.second->getGeneMap()) {

        if (T::geneFamily(gene.second, genome_db_ptr, protein_family)) {

          std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance,
                                                              genome_variant.second,
                                                              genome_db_ptr,
                                                              gene.second,
                                                              protein_family,
                                                              args...));
          std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
          node_vector_ptr->push_back(phylo_node_ptr);

        }

      }

    }

    UPGMAMatrix upgma_matrix(node_vector_ptr);

    upgma_matrix.calculateReduce();

    std::string genome_newick = genome_variant.second->genomeId() + "_" + newick_file;
    genome_newick = Utility::filePath(genome_newick, path);
    upgma_matrix.writeNewick(genome_newick);

  }


}


// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a series of Gene trees.
// We are comparing between genes of the same type so we can use a local and global Amino distance classes
template<typename T, typename... Args>
void UPGMAGenePhyloTree(const std::string& path,
                        const std::string& newick_file,
                        std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                        std::shared_ptr<const PhasedPopulation> pop_variant_ptr,
                        std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                        const std::string& protein_family,
                        Args... args) {


    for (auto contig : genome_db_ptr->getMap()) {

      for (auto gene : contig.second->getGeneMap()) {

        if (T::geneFamily(gene.second, genome_db_ptr, protein_family)) {

          std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());

          for (auto genome_variant : pop_variant_ptr->getMap()) {

          std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance,
                                                              genome_variant.second,
                                                              genome_db_ptr,
                                                              gene.second,
                                                              protein_family,
                                                              args...));
          std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
          node_vector_ptr->push_back(phylo_node_ptr);

        }

        UPGMAMatrix upgma_matrix(node_vector_ptr);

        upgma_matrix.calculateReduce();

        std::string genome_newick = gene.second->id() + "_" + newick_file;
        genome_newick = Utility::filePath(genome_newick, path);
        upgma_matrix.writeNewick(genome_newick);

      }

    }

  }

}



// Function (not variadic) to combine the UPGMAMatrix and UPGMADistanceNode to compare a family of reference genes (unmutated)
// We are comparing between genes of the same type so we can use both local and global Amino distance classes
template<typename T>
void UPGMAGeneFamilyTree(const std::string& newick_file,
                          std::shared_ptr<const DNASequenceDistance> sequence_distance,
                          std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                          const std::string& protein_family) {


  std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());

  for (auto contig : genome_db_ptr->getMap()) {

    for (auto gene : contig.second->getGeneMap()) {

      // Is this gene a member of the requested family.
      if (T::geneFamily(gene.second, genome_db_ptr, protein_family)) {

        const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = GeneFeature::getCodingSequences(gene.second);

        if (coding_seq_ptr->empty()) {

          ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Gene contains no coding sequence : gene: {}", gene.second->id());

        }

        std::shared_ptr<const CodingSequence> coding_sequence = coding_seq_ptr->getFirst();
        std::shared_ptr<DNA5SequenceCoding> coding_dna_sequence;

        // Only add genes with valid coding sequences (no pseudo genes).
        if (contig.second->getDNA5SequenceCoding(coding_sequence, coding_dna_sequence)) {


          if (contig.second->verifyDNACodingSequence(coding_dna_sequence)) {

#define MIN_INTRON_LENGTH 10
            // Do we have a valid intron (VAR only)?
            std::shared_ptr<const DNA5SequenceCoding> intron_sequence = contig.second->sequence().intronSequence(coding_sequence);
            if (intron_sequence->length() > MIN_INTRON_LENGTH) {

              std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance,
                                                                  genome_db_ptr,
                                                                  gene.second,
                                                                  protein_family));

              std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
              node_vector_ptr->push_back(phylo_node_ptr);


              ExecEnv::log().info("ReferenceGeneDistance::getSequence(); Gene: {} intron ({}):\n {}-{}",
                                  gene.second->id(), intron_sequence->length(), gene.second->id(),
                                  intron_sequence->getSequenceAsString());

            }

          } // Valid Gene.

        } // Get coding sequence

      } // Is family member.

    } // All genes.

  } // All contigs.


  // Assemble the UPGMA structure.
  UPGMAMatrix upgma_matrix(node_vector_ptr);
  // Calculate.
  upgma_matrix.calculateReduce();
  // Report Results.
  upgma_matrix.writeNewick(newick_file);

}



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_STATISTICS_UPGMA_H
