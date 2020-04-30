//
// Created by kellerberrin on 30/01/18.
//


#include <iomanip>
#include <kgl_sequence_compare_impl.h>
#include "kgl_upgma_node.h"


namespace kgl = kellerberrin::genome;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates genome proteins and then compares the mutated proteins.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::DistanceType_t kgl::UPGMAProteinDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAProteinDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAProteinDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to UPGMAProteinDistance");
    return 1.0;

  }

  CompareDistance_t total_distance = 0;
  size_t gene_count = 0;

  for (auto const& [protein_id, protein_ptr] : getMap()) {

    auto result = node_ptr->getMap().find(protein_id);

    if (result != node_ptr->getMap().end()) {

      CompareDistance_t contig_distance = sequence_distance_->amino_distance(*protein_ptr, *result->second);
      total_distance += contig_distance;
      ++gene_count;

    } else {

      ExecEnv::log().error("distance(), Unexpected error, could not find Genome: {}, sequence: {}",
                           node_ptr->genomeId(), protein_id);

    }

  }

  ExecEnv::log().info("distance(), Genome: {}, Genome: {}; {} distance: {}, Gene Family: {}, Gene Count: {}",
                      genomeId(), node_ptr->genomeId(), sequence_distance_->distanceType(), total_distance, protein_family_, gene_count);
  return total_distance;

}


void kgl::UPGMAProteinDistance::mutateProteins() {

  mutated_proteins_.clear();
  for (auto contig : genome_db_ptr_->getMap()) {

    for (auto gene : contig.second->getGeneMap()) {

      if (protein_family_ != PROTEIN_FAMILY_WILDCARD) {

        std::shared_ptr<const OntologyRecord> ontology_record_ptr;
        if (ontology().getGafFeatureVector(gene.second->id(), ontology_record_ptr)) {

          if (ontology_record_ptr->symbolicReference() == protein_family_) {

            getProtein(gene.second);

          }
        }
      }
      else {

        getProtein(gene.second);

      }

    }

  }

}


void kgl::UPGMAProteinDistance::getProtein(std::shared_ptr<const GeneFeature> gene_ptr) {

  const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr);
  for (auto const& sequence : coding_seq_ptr->getMap()) {

    std::shared_ptr<const ContigFeatures> contig_ptr = sequence.second->getGene()->contig();
    std::string gene_id = sequence.second->getGene()->id();
    std::string sequence_id = sequence.second->getCDSParent()->id();
    AminoSequence mutant_sequence;
    AminoSequence reference_sequence;
    OffsetVariantMap variant_map;

    if (genome_variant_ptr_->mutantProteins(contig_ptr->contigId(),
                                            ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                            gene_id,
                                            sequence_id,
                                            genome_db_ptr_,
                                            variant_map,
                                            reference_sequence,
                                            mutant_sequence)) {

      mutated_proteins_[gene_id] = std::make_shared<AminoSequence>(std::move(mutant_sequence));

    }

  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates a single gene and compares to other selected genes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::UPGMAGeneDistance::geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                                        std::shared_ptr<const RuntimeGenomeDatabase> genome_db_ptr,
                                        const std::string& protein_family) {

  std::shared_ptr<const OntologyRecord> ontology_record_ptr;
  if (genome_db_ptr->geneOntology().getGafFeatureVector(gene_ptr->id(), ontology_record_ptr)) {

    if (ontology_record_ptr->symbolicReference() == protein_family) {

      return true;

    }

  }

  return false;

}


void kgl::UPGMAGeneDistance::mutateProtein() {

  const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr_);

  if (coding_seq_ptr->empty()) {

    ExecEnv::log().critical("mutateProtein(), Gene contains no coding sequence : genome: {} gene: {}",
                            genome_variant_ptr_->genomeId(), gene_ptr_->id());

  }

  std::shared_ptr<const CodingSequence> sequence = coding_seq_ptr->getFirst();
  std::shared_ptr<const ContigFeatures> contig_ptr = sequence->getGene()->contig();
  std::string gene_id = sequence->getGene()->id();
  std::string sequence_id = sequence->getCDSParent()->id();

  if (coding_seq_ptr->size() > 1) {

    ExecEnv::log().warn("mutateProtein(),  Genome: {} gene: {} contains: {} sequences using sequence: {}",
                        genome_variant_ptr_->genomeId(), gene_ptr_->id(), coding_seq_ptr->size(), sequence_id);

  }


  AminoSequence mutant_sequence;
  AminoSequence reference_sequence;
  OffsetVariantMap variant_map;

  if (genome_variant_ptr_->mutantProteins(contig_ptr->contigId(),
                                          ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                          gene_id,
                                          sequence_id,
                                          genome_db_ptr_,
                                          variant_map,
                                          reference_sequence,
                                          mutated_protein_)) {

  }

}

kgl::DistanceType_t kgl::UPGMAGeneDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAGeneDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAGeneDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to UPGMAGeneDistance");
    return 1.0;

  }


  CompareDistance_t contig_score = sequence_distance_->amino_distance(mutated_protein_, node_ptr->mutated_protein_);
  DistanceType_t total_distance = static_cast<DistanceType_t>(contig_score);

  ExecEnv::log().info("distance(), Genome: {}, Gene: {}, Gene: {}; {} distance: {}, Gene Family: {}",
                      genome_variant_ptr_->genomeId(), gene_ptr_->id(), node_ptr->gene_ptr_->id(),
                      sequence_distance_->distanceType(), total_distance, protein_family_);

  return total_distance;

}


void kgl::UPGMAGeneDistance::writeNode(std::ostream& outfile) const {

  std::stringstream ss;

  double contig_proportion = static_cast<double>(gene_ptr_->sequence().begin()) / static_cast<double>(gene_ptr_->contig()->contigSize());
  contig_proportion = contig_proportion * 100.0;
  ss << gene_ptr_->id() << "_" << std::setprecision(3) << contig_proportion;

  std::shared_ptr<const OntologyRecord> ontology_record_ptr;
  if (genome_db_ptr_->geneOntology().getGafFeatureVector(gene_ptr_->id(), ontology_record_ptr)) {

    ss << "_" << ontology_record_ptr->altSymbolicReference();

  }

  outfile << ss.str();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares a single gene between isolates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::UPGMAATP4Distance::writeNode(std::ostream& outfile) const {

  const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr_);

  if (coding_seq_ptr->size() == 0) {

    ExecEnv::log().critical("write_node(), Gene contains no coding sequence : genome: {} gene: {}",
                            genome_variant_ptr_->genomeId(), gene_ptr_->id());

  }

  std::shared_ptr<const CodingSequence> sequence = coding_seq_ptr->getFirst();
  std::shared_ptr<const ContigFeatures> contig_ptr = sequence->getGene()->contig();
  std::string gene_id = sequence->getGene()->id();
  std::string sequence_id = sequence->getCDSParent()->id();

  if (coding_seq_ptr->size() > 1) {

    ExecEnv::log().warn("write_node(),  Genome: {} gene: {} contains: {} sequences using sequence: {}",
                        genome_variant_ptr_->genomeId(), gene_ptr_->id(), coding_seq_ptr->size(), sequence_id);

  }


  AminoSequence mutant_sequence;
  AminoSequence reference_sequence;
  OffsetVariantMap variant_map;

  if (not genome_variant_ptr_->mutantProteins(contig_ptr->contigId(),
                                              ContigVariant::HAPLOID_HOMOLOGOUS_INDEX,
                                              gene_id,
                                              sequence_id,
                                              genome_db_ptr_,
                                              variant_map,
                                              reference_sequence,
                                              mutant_sequence)) {

    ExecEnv::log().critical("write_node(), Cannot mutate sequence for : genome: {} gene: {}",
                            genome_variant_ptr_->genomeId(), gene_ptr_->id());

  }

  std::stringstream ss;

  ss << genome_variant_ptr_->genomeId();

  outfile << ss.str();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares reference (unmutated) genes between gene familes (VAR, RIF, STEVOR, MAURER).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<const kgl::CodingSequence>  kgl::ReferenceGeneDistance::getCodingSequence() {

  const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = GeneFeature::getCodingSequences(gene_ptr_);

  if (coding_seq_ptr->empty()) {

    ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Gene contains no coding sequence : gene: {}", gene_ptr_->id());

  }

  std::shared_ptr<const CodingSequence> coding_sequence = coding_seq_ptr->getFirst();

  if (coding_seq_ptr->size() > 1) {

    std::string gene_id = coding_sequence->getGene()->id();
    std::string sequence_id = coding_sequence->getCDSParent()->id();

    ExecEnv::log().warn("ReferenceGeneDistance::getSequence;  Gene: {} contains: {} CDS sequences. Using sequence: {}",
                        gene_ptr_->id(), coding_seq_ptr->size(), sequence_id);

  }

  return coding_sequence;

}


void kgl::ReferenceGeneDistance::writeNode(std::ostream& outfile) const {

  std::string alt_symbolic;
  std::shared_ptr<const OntologyRecord> ontology_record_ptr;
  if (genome_db_ptr_->geneOntology().getGafFeatureVector(gene_ptr_->id(), ontology_record_ptr)) {

    if (ontology_record_ptr) {

      alt_symbolic = "-";
      alt_symbolic += ontology_record_ptr->altSymbolicReference();

    } else {

      ExecEnv::log().error("ReferenceGeneDistance::writeNode; NULL OntologyRecord pointer returned");
      alt_symbolic = "-*";

    }


  } else {

    alt_symbolic = "-*";

  }


  outfile << gene_ptr_->id();

}


void  kgl::DNAGeneDistance::getExonSequence() {

  std::shared_ptr<const ContigFeatures> contig_ptr = gene_ptr_->contig();

  DNA5SequenceCoding dna_coding_sequence = contig_ptr->sequence().codingSequence(getCodingSequence());

  linear_sequence_ = DNA5SequenceLinear::downConvertToLinear(dna_coding_sequence);


}


void  kgl::DNAGeneDistance::getIntronSequence() {

  std::shared_ptr<const ContigFeatures> contig_ptr = gene_ptr_->contig();

  DNA5SequenceCoding dna_coding_sequence = contig_ptr->sequence().intronSequence(getCodingSequence());

  linear_sequence_ = DNA5SequenceLinear::downConvertToLinear(dna_coding_sequence);

}


kgl::DistanceType_t kgl::DNAGeneDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const DNAGeneDistance> node_ptr = std::dynamic_pointer_cast<const DNAGeneDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to DNAGeneDistance");
    return 1.0;

  }


  ExecEnv::log().info("distance();  {} Comparing | {}({}), {}({}) |; Gene Family: {}",
                      sequence_distance_->distanceType(), gene_ptr_->id(), linear_sequence_.length(),
                      node_ptr->gene_ptr_->id(), node_ptr->linear_sequence_.length(), protein_family_);


  CompareDistance_t contig_score = sequence_distance_->linear_distance(linear_sequence_, node_ptr->linear_sequence_);

  DistanceType_t total_distance = static_cast<DistanceType_t>(contig_score);

  ExecEnv::log().info("distance();  {} | {}({}), {}({}) |  =  {}; Gene Family: {}",
                      sequence_distance_->distanceType(), gene_ptr_->id(), linear_sequence_.length(),
                      node_ptr->gene_ptr_->id(), node_ptr->linear_sequence_.length(),
                      total_distance, protein_family_);

  return total_distance;

}


bool kgl::ReferenceGeneDistance::geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                                            std::shared_ptr<const RuntimeGenomeDatabase> genome_db_ptr,
                                            const std::string& protein_family) {

  std::shared_ptr<const OntologyRecord> ontology_record_ptr;
  if (genome_db_ptr->geneOntology().getGafFeatureVector(gene_ptr->id(), ontology_record_ptr)) {

    if (ontology_record_ptr->symbolicReference() == protein_family) {

      return true;

    }

  }

  return false;

}



void  kgl::AminoGeneDistance::getAminoSequence() {

  std::shared_ptr<const ContigFeatures> contig_ptr = gene_ptr_->contig();

  DNA5SequenceCoding dna_coding_sequence = contig_ptr->sequence().codingSequence(getCodingSequence());

  amino_sequence_ = contig_ptr->getAminoSequence(dna_coding_sequence);

}


kgl::DistanceType_t kgl::AminoGeneDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const AminoGeneDistance> node_ptr = std::dynamic_pointer_cast<const AminoGeneDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to DNAGeneDistance");
    return 1.0;

  }


  ExecEnv::log().info("distance();  {} Comparing | {}({}), {}({}) |; Gene Family: {}",
                      sequence_distance_->distanceType(), gene_ptr_->id(), amino_sequence_.length(),
                      node_ptr->gene_ptr_->id(), node_ptr->amino_sequence_.length(), protein_family_);


  CompareDistance_t contig_score = sequence_distance_->amino_distance(amino_sequence_, node_ptr->amino_sequence_);

  DistanceType_t total_distance = static_cast<DistanceType_t>(contig_score);

  ExecEnv::log().info("distance();  {} | {}({}), {}({}) |  =  {}; Gene Family: {}",
                      sequence_distance_->distanceType(), gene_ptr_->id(), amino_sequence_.length(),
                      node_ptr->gene_ptr_->id(), node_ptr->amino_sequence_.length(),
                      total_distance, protein_family_);

  return total_distance;

}
