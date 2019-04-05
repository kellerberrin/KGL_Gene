//
// Created by kellerberrin on 26/11/18.
//




#ifndef KGL_UPGMA_H
#define KGL_UPGMA_H


#include "kgl_upgma_node.h"
#include "kgl_sequence_offset.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




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
                         const std::string& intron_file,
                         std::shared_ptr<const DNASequenceDistance> sequence_distance,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         const std::string& protein_family) {


  // Open a file to receive csv delimited intron information.
  const char delimiter = ',';
  std::ofstream intron(intron_file);

  if (not intron.good()) {

    ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Unable to open Intron data file: {}", intron_file);

  }

#define I_PROMOTER "TGTATGTG"
#define I_COMPLEMENT_PROMOTER "ACATACAC"
#define I_5_PROMOTER "TCATA"
  // Header.
  intron << "Gene" << delimiter
         << "description" << delimiter
         << "Symbolic" << delimiter
         << "AltSymbolic" << delimiter
         << "IStart" << delimiter
         << "IEnd" << delimiter
         << "IStrand" << delimiter
         << I_PROMOTER << delimiter
         << I_COMPLEMENT_PROMOTER << delimiter
         << I_5_PROMOTER << delimiter
         << "Size" << delimiter
         << "ISequence" << '\n';

  std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());

  for (auto contig : genome_db_ptr->getMap()) {

    for (auto gene : contig.second->getGeneMap()) {

      std::string alt_symbolic;
      std::string symbolic;
      std::string description;
      std::shared_ptr<const OntologyRecord> ontology_record_ptr;
      if (genome_db_ptr->geneOntology().getGafFeatureVector(gene.second->id(), ontology_record_ptr)) {

        if (ontology_record_ptr) {

          description = Utility::findAndReplaceAll(ontology_record_ptr->description(), ",", ";");
          symbolic = Utility::findAndReplaceAll(ontology_record_ptr->symbolicReference(), ",", ";");
          alt_symbolic = Utility::findAndReplaceAll(ontology_record_ptr->altSymbolicReference(), ",", ";");

        } else {

          continue;

        }

      } else {

        continue;

      }

      // Is this gene a member of the requested family.
      if (Utility::toupper(description).find(protein_family) != std::string::npos) {

        const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = GeneFeature::getCodingSequences(gene.second);

        if (coding_seq_ptr->empty()) {

          ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Gene contains no coding sequence : gene: {}", gene.second->id());

        }

        std::shared_ptr<const CodingSequence> coding_sequence = coding_seq_ptr->getFirst();
        std::shared_ptr<DNA5SequenceCoding> coding_dna_sequence;


        if (contig.second->getDNA5SequenceCoding(coding_sequence, coding_dna_sequence)) {

          // Do we have a valid intron (VAR only)?
          std::vector<std::shared_ptr<DNA5SequenceCoding>> intron_sequence_array = contig.second->sequence().intronArraySequence(coding_sequence);
          std::shared_ptr<T> distance_ptr(std::make_shared<T>(sequence_distance,
                                                              genome_db_ptr,
                                                              gene.second,
                                                              protein_family));

          StrandSense strand;
          IntronOffsetMap intron_offset_map;
          SequenceOffset::intronOffsetAdapter(coding_sequence, strand, intron_offset_map);

          if (intron_offset_map.size() != intron_sequence_array.size()) {

            ExecEnv::log().error("UPGMAGeneFamilyTree, Intron map size: {} different from Inron sequence size: {}",
                                 intron_offset_map.size(), intron_sequence_array.size());
            continue;

          }

          auto intron_sequence_iter = intron_sequence_array.begin();


          for (auto intron_seq : intron_offset_map) {

            std::stringstream ss;
            distance_ptr->writeNode(ss);

            intron << gene.second->id() << delimiter
                   << description << delimiter
                   << symbolic << delimiter
                   << alt_symbolic << delimiter
                   << intron_seq.first << delimiter
                   << intron_seq.second << delimiter
                   << static_cast<char>(strand) << delimiter
                   << DNA5SequenceLinear::linearSequence((*intron_sequence_iter))->findAll(DNA5SequenceLinear(StringDNA5(I_PROMOTER))).size()
                   << delimiter
                   << DNA5SequenceLinear::linearSequence((*intron_sequence_iter))->findAll(DNA5SequenceLinear(StringDNA5(I_COMPLEMENT_PROMOTER))).size()
                   << delimiter
                   << DNA5SequenceLinear::linearSequence((*intron_sequence_iter))->findAll(DNA5SequenceLinear(StringDNA5(I_5_PROMOTER))).size()
                   << delimiter
                   << (*intron_sequence_iter)->length()
                   << delimiter
                   << (*intron_sequence_iter)->getSequenceAsString() << '\n';

            ++intron_sequence_iter;

          }

          // Only add genes with valid coding sequences (no pseudo genes).
          if (contig.second->verifyDNACodingSequence(coding_dna_sequence)) {

            // Only 1 intron (var genes)
            if (intron_sequence_array.size() == 1) {

#define MIN_INTRON_LENGTH 10
              if (intron_sequence_array.front()->length() >= MIN_INTRON_LENGTH) {

                std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
                node_vector_ptr->push_back(phylo_node_ptr);

              } // min length

            } // 1 intron

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



#endif //KGL_UPGMA_H
