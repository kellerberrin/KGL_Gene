//
// Created by kellerberrin on 30/01/18.
//


#include "kgl_upgma.h"


namespace kgl = kellerberrin::genome;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates the contigs and then compares the mutated contigs using the Myer Hirschberg sequence comparison
// algorthim. This is linear in space - but quadratic in time. Need to find a faster comparison algorithm.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DistanceType_t kgl::UPGMAContigDistance::distance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAContigDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAContigDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not up-cast node pointer");
    return 1.0;

  }

  DistanceType_t total_distance = 0;
  for (auto contig : mutated_contigs_) {

    auto result = node_ptr->mutated_contigs_.find(contig.first);

    if (result != node_ptr->mutated_contigs_.end()) {

      ExecEnv::log().info("distance(), Comparing Genome: {}, Contig: {} with Genome: {}",
                          genome_variant_ptr_->genomeId(), contig.first, node_ptr->genome_variant_ptr_->genomeId());
      CompareScore_t contig_score = contig.second->compareLevenshtein(result->second);
      total_distance += static_cast<DistanceType_t>(contig_score);
      ExecEnv::log().info("distance(), Calculated distance: {}", contig_score);

    } else {

      ExecEnv::log().error("distance(), Unexpected error, could not find contig: {}", contig.first);

    }

  }

  return total_distance;

}



void kgl::UPGMAContigDistance::mutateContigs() {

  std::shared_ptr<const DNA5SequenceContig> reference_contig_ptr;
  std::shared_ptr<DNA5SequenceContig> mutant_contig_ptr;
  for (auto contig : genome_db_ptr_->getMap()) {

    if (not genome_variant_ptr_->mutantContig(contig.first, genome_db_ptr_, reference_contig_ptr, mutant_contig_ptr)) {

      ExecEnv::log().error("Unexpected error mutating contig; genome: {} , contig: {}", genome_variant_ptr_->genomeId(), contig.first);
      // Fail gracefully, just insert the reference contig.
      mutated_contigs_.insert(std::pair<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>(contig.first, contig.second->sequence_ptr()));

    } else {

      mutated_contigs_.insert(std::pair<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>(contig.first, mutant_contig_ptr));

    }

  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates genome proteins and then compares the mutated proteins using the Myer Hirschberg sequence comparison
// algorthim. Myer Hirschberg is linear in space - but quadratic in time. Need to find a faster comparison algorithm.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::DistanceType_t kgl::UPGMAProteinDistance::distance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAProteinDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAProteinDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not up-cast node pointer");
    return 1.0;

  }

  DistanceType_t total_distance = 0;

  for (auto protein : getMap()) {

    auto result = node_ptr->getMap().find(protein.first);

    if (result != node_ptr->getMap().end()) {

      CompareScore_t contig_score = protein.second->compareLevenshtein(result->second);
      total_distance += std::fabs(static_cast<DistanceType_t>(contig_score));

    } else {

      ExecEnv::log().error("distance(), Unexpected error, could not find Genome: {}, sequence: {}",
                           node_ptr->genomeId(), protein.first);

    }

  }

  ExecEnv::log().info("distance(), Genome: {}, Genome: {}; Calculated distance: {}",
                      genomeId(), node_ptr->genomeId(), total_distance);
  return total_distance;

}



void kgl::UPGMAProteinDistance::mutateProteins() {

  mutated_proteins_.clear();
  for (auto contig : genome_db_ptr_->getMap()) {

    for (auto gene : contig.second->getGeneMap()) {

      const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene.second);
      for (auto sequence : coding_seq_ptr->getMap()) {

        std::shared_ptr<const ContigFeatures> contig_ptr = sequence.second->getGene()->contig();
        std::string gene_id = sequence.second->getGene()->id();
        std::string sequence_id = sequence.second->getCDSParent()->id();
        std::vector<std::shared_ptr<DNA5SequenceCoding>> mutant_sequence_vector;
        std::shared_ptr<DNA5SequenceCoding> reference_sequence;
        OffsetVariantMap variant_map;

        if (genome_variant_ptr_->mutantCodingDNA(contig_ptr->contigId(),
                                                 gene_id,
                                                 sequence_id,
                                                 genome_db_ptr_,
                                                 variant_map,
                                                 reference_sequence,
                                                 mutant_sequence_vector)) {

          if (mutant_sequence_vector.empty()) {

            ExecEnv::log().warn("mutateProteins(), No mutant proteins generated for : genome: {} sequence: {}",
                                genome_variant_ptr_->genomeId(), sequence_id);

          } else if (mutant_sequence_vector.size() == 1) {

            mutated_proteins_[gene_id] = mutant_sequence_vector.front();

          } else {

            ExecEnv::log().warn("mutateProteins(), {} mutant proteins generated for : genome: {} sequence: {} only the first mutant is analyzed",
                                mutant_sequence_vector.size(), genome_variant_ptr_->genomeId(), sequence_id);

            mutated_proteins_[gene_id] = mutant_sequence_vector.front();

          }

        }

      }

    }

  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance using protein families
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DistanceType_t kgl::UPGMAFamilyDistance::distance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAProteinDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAFamilyDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to UPGMAFamilyDistance");
    return 1.0;

  }

  DistanceType_t total_distance = 0;
  size_t gene_count = 0;

  for (auto protein : getMap()) {

    auto result = node_ptr->getMap().find(protein.first);

    if (result != node_ptr->getMap().end()) {

      std::shared_ptr<const OntologyRecord> ontology_record_ptr;
      if (ontology().getGafFeatureVector(protein.first, ontology_record_ptr)) {

        if (ontology_record_ptr->symbolicReference() == family_code_) {

          CompareScore_t contig_score = protein.second->compareLevenshtein(result->second);
          total_distance += std::fabs(static_cast<DistanceType_t>(contig_score));
          ++gene_count;

        }
        else {

//          ExecEnv::log().info("Record Family code: {}, Specified Family Code: {}", ontology_record_ptr->symbolicReference(), family_code_);

        }

      } else {

//        ExecEnv::log().info("Not found; Ontology Gene code: {}", protein.first);

      }

    } else {

      ExecEnv::log().error("distance(), Unexpected error, could not find Genome: {}, sequence: {}",
                           node_ptr->genomeId(), protein.first);

    }

  }

  ExecEnv::log().info("distance(), Genome: {}, Genome: {}; Calculated distance: {}, Gene Family: {}, Gene Count: {}",
                      genomeId(), node_ptr->genomeId(), total_distance, family_code_, gene_count);
  return total_distance;

}
