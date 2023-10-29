//
// Created by kellerberrin on 26/11/18.
//




#ifndef KGL_UPGMA_H
#define KGL_UPGMA_H


#include "kgl_distance_sequence.h"


namespace kellerberrin::genome {   //  organization::project level namespace



// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a population tree.
// We are comparing contigs across genomes. The distance metric will be based on the difference in
// the unphased variants held for each contig_ref_ptr.
template<typename T, typename... Args>
void UnphasedDistanceTree(DistanceTreeBase& distance_tree,
                          const std::string& newick_file,
                          std::shared_ptr<const PopulationDB> pop_unphased_ptr,
                          Args... args) {

  std::shared_ptr<DistanceNodeVector> node_vector_ptr(std::make_shared<DistanceNodeVector>());

  for (auto genome : pop_unphased_ptr->getMap()) {

    std::shared_ptr<T> distance_ptr(std::make_shared<T>(genome.second, std::forward(args)...));
    std::shared_ptr<TreeDistanceNode> phylo_node_ptr(std::make_shared<TreeDistanceNode>(distance_ptr));
    node_vector_ptr->push_back(phylo_node_ptr);

  }

  distance_tree.calculateTree(node_vector_ptr);

  distance_tree.writeNewick(newick_file);

}


// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a population tree.
// We are comparing contigs across genomes, so a global dna distance metric is required.
template<typename T, typename... Args>
void PopulationDistanceTree(DistanceTreeBase& distance_tree,
                            const std::string& newick_file,
                            std::shared_ptr<const DNASequenceDistance> sequence_distance,
                            std::shared_ptr<const PopulationDB> pop_variant_ptr,
                            std::shared_ptr<const GenomeReference> genome_db_ptr,
                            Args... args) {

  std::shared_ptr<DistanceNodeVector> node_vector_ptr(std::make_shared<DistanceNodeVector>());

  for (auto genome : pop_variant_ptr->getMap()) {

    std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance, genome.second, genome_db_ptr, args...));
    std::shared_ptr<TreeDistanceNode> phylo_node_ptr(std::make_shared<TreeDistanceNode>(distance_ptr));
    node_vector_ptr->push_back(phylo_node_ptr);

  }

  distance_tree.calculateTree(node_vector_ptr);

  distance_tree.writeNewick(newick_file);

}


// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a series of Gene trees.
// We are comparing dissimilar gene types, so only a local Amino distance metric should be used
template<typename T, typename... Args>
void GeneDistanceTree(DistanceTreeBase& distance_tree,
                      const std::string& newick_file,
                      std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                      std::shared_ptr<const PopulationDB> pop_variant_ptr,
                      std::shared_ptr<const GenomeReference> genome_db_ptr,
                      const std::string& protein_family,
                      Args... args) {


  for (auto genome_variant : pop_variant_ptr->getMap()) {

    std::shared_ptr<DistanceNodeVector> node_vector_ptr(std::make_shared<DistanceNodeVector>());

    for (auto contig : genome_db_ptr->getMap()) {

      for (auto gene : contig.second->getGeneMap()) {

        if (T::geneFamily(gene.second, genome_db_ptr, protein_family)) {

          std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance,
                                                              genome_variant.second,
                                                              genome_db_ptr,
                                                              gene.second,
                                                              protein_family,
                                                              args...));
          std::shared_ptr<TreeDistanceNode> phylo_node_ptr(std::make_shared<TreeDistanceNode>(distance_ptr));
          node_vector_ptr->push_back(phylo_node_ptr);

        }

      }

    }

    distance_tree.calculateTree(node_vector_ptr);

    distance_tree.writeNewick(newick_file);

  }


}


// Variadic function to combine the UPGMAMatrix and UPGMADistanceNode to produce a series of Gene trees.
// We are comparing between genes of the same type so we can use a local and global Amino distance classes
template<typename T, typename... Args>
void GenePhyloTree(DistanceTreeBase& distance_tree,
                   const std::string& newick_file,
                   std::shared_ptr<const AminoSequenceDistance> sequence_distance,
                   std::shared_ptr<const PopulationDB> pop_variant_ptr,
                   std::shared_ptr<const GenomeReference> genome_db_ptr,
                   const std::string& protein_family,
                   Args... args) {


  for (auto contig : genome_db_ptr->getMap()) {

    for (auto gene : contig.second->getGeneMap()) {

      if (T::geneFamily(gene.second, genome_db_ptr, protein_family)) {

        std::shared_ptr<DistanceNodeVector> node_vector_ptr(std::make_shared<DistanceNodeVector>());

        for (auto genome_variant : pop_variant_ptr->getMap()) {

          std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance,
                                                              genome_variant.second,
                                                              genome_db_ptr,
                                                              gene.second,
                                                              protein_family,
                                                              args...));
          std::shared_ptr<TreeDistanceNode> phylo_node_ptr(std::make_shared<TreeDistanceNode>(distance_ptr));
          node_vector_ptr->push_back(phylo_node_ptr);

        }

        distance_tree.calculateTree(node_vector_ptr);

        distance_tree.writeNewick(newick_file);

      }

    }

  }

}



}   // end namespace genome



#endif //KGL_UPGMA_H
