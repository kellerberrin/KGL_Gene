//
// Created by kellerberrin on 11/11/23.
//

#include "kgl_distance_sequence.h"


#include "kgl_sequence_compare_impl.h"
#include "kgl_distance_sequence.h"
#include "kgl_mutation_transcript.h"


namespace kgl = kellerberrin::genome;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates genome proteins and then compares the mutated proteins.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::DistanceType_t kgl::ProteinDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const ProteinDistance> node_ptr = std::dynamic_pointer_cast<const ProteinDistance>(distance_node);
  if (not node_ptr) {

    ExecEnv::log().error("UPGMAProteinDistance::distance; Unexpected error, could not down-cast node pointer to UPGMAProteinDistance");
    return 1.0;

  }

  CompareDistance_t total_distance = 0;
  size_t gene_count = 0;

  for (auto const& [protein_id, protein_ptr] : getMap()) {

    auto find_iter = node_ptr->getMap().find(protein_id);
    if (find_iter != node_ptr->getMap().end()) {

      auto const& [feature_id, amino_sequence_ptr] = *find_iter;
      CompareDistance_t contig_distance = sequence_distance_(*protein_ptr, *amino_sequence_ptr);
      total_distance += contig_distance;
      ++gene_count;

    } else {

      ExecEnv::log().error("UPGMAProteinDistance::distance; Unexpected error, could not find Genome: {}, sequence: {}",
                           node_ptr->genomeId(), protein_id);

    }

  }

  ExecEnv::log().info("UPGMAProteinDistance::distance; Genome: {}, Genome: {}; distance: {}, Gene Family: {}, Gene Count: {}",
                      genomeId(), node_ptr->genomeId(), total_distance, protein_family_, gene_count);
  return total_distance;

}


void kgl::ProteinDistance::mutateProteins() {

  mutated_proteins_.clear();
  for (auto const& [contig_id, contig_ptr] : genome_ref_ptr_->getMap()) {

    for (auto const& [gene_id, gene_ptr] : contig_ptr->getGeneMap()) {


      getProtein(gene_ptr);


    }

  }

}


void kgl::ProteinDistance::getProtein(std::shared_ptr<const GeneFeature> gene_ptr) {

  auto transcript_array_ptr = kgl::GeneFeature::getTranscriptionSequences(gene_ptr);
  for (auto const& [transcript_id, transcript_ptr] : transcript_array_ptr->getMap()) {

    auto& contig_ref_ptr = transcript_ptr->getGene()->contig_ref_ptr();
    const std::string& contig_id = contig_ref_ptr->contigId();
    const std::string& gene_id = transcript_ptr->getGene()->id();

    auto contig_db_opt = genome_db_ptr_->getContig(contig_id);
    if (not contig_db_opt) {

      ExecEnv::log().warn("Contig: {} not found for Genome: {}", contig_id, genome_db_ptr_->genomeId());
      return;

    }
    const auto& contig_db_ptr = contig_db_opt.value();

    SequenceTranscript modified_transcript(contig_db_ptr, transcript_ptr);

    auto reference_sequence_opt = modified_transcript.getOriginalCoding();
    auto mutant_sequence_opt = modified_transcript.getModifiedCoding();

    if (modified_transcript.sequenceStatus() and reference_sequence_opt and mutant_sequence_opt) {

      auto reference_sequence = contig_ref_ptr->getAminoSequence(reference_sequence_opt.value());
      auto mutant_sequence = contig_ref_ptr->getAminoSequence(mutant_sequence_opt.value());


      mutated_proteins_[gene_id] = std::make_shared<AminoSequence>(std::move(mutant_sequence));

    }

  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates a single gene and compares to other selected genes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GeneDistance::geneFamily(std::shared_ptr<const GeneFeature> ,
                                   std::shared_ptr<const GenomeReference> ,
                                   const std::string& ) {

  return false;

}


void kgl::GeneDistance::mutateProtein() {

  auto transcipt_array_ptr = kgl::GeneFeature::getTranscriptionSequences(gene_ptr_);

  if (transcipt_array_ptr->empty()) {

    ExecEnv::log().critical("Gene contains no coding transcript_ptr : genome: {} gene: {}",
                            genome_db_ptr_->genomeId(), gene_ptr_->id());

  }

  const auto& transcript_ptr = transcipt_array_ptr->getFirst();
  const auto& contig_ref_ptr = transcript_ptr->getGene()->contig_ref_ptr();
  const std::string& contig_id = contig_ref_ptr->contigId();
  const std::string& transcript_id = transcript_ptr->getParent()->id();

  if (transcipt_array_ptr->size() > 1) {

    ExecEnv::log().warn("Genome: {} gene: {} contains: {} sequences using transcript_ptr: {}",
                        genome_db_ptr_->genomeId(), gene_ptr_->id(), transcipt_array_ptr->size(), transcript_id);

  }

  auto contig_db_opt = genome_db_ptr_->getContig(contig_id);
  if (not contig_db_opt) {

    ExecEnv::log().warn("Contig: {} not found for Genome: {}", contig_id, genome_db_ptr_->genomeId());
    return;

  }
  const auto& contig_db_ptr = contig_db_opt.value();

  const SequenceTranscript modified_transcript(contig_db_ptr, transcript_ptr);

  auto reference_sequence_opt = modified_transcript.getOriginalCoding();
  auto mutant_sequence_opt = modified_transcript.getModifiedCoding();

  if (modified_transcript.sequenceStatus() and reference_sequence_opt and mutant_sequence_opt) {

    auto reference_sequence = contig_ref_ptr->getAminoSequence(reference_sequence_opt.value());
    auto mutant_sequence = contig_ref_ptr->getAminoSequence(mutant_sequence_opt.value());

  }

}

kgl::DistanceType_t kgl::GeneDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const GeneDistance> node_ptr = std::dynamic_pointer_cast<const GeneDistance>(distance_node);
  if (not node_ptr) {

    ExecEnv::log().error("GeneDistance::distance; Unexpected error, could not down-cast node pointer to GeneDistance");
    return 1.0;

  }

  CompareDistance_t contig_score = sequence_distance_(mutated_protein_, node_ptr->mutated_protein_);
  DistanceType_t total_distance = static_cast<DistanceType_t>(contig_score);

  ExecEnv::log().info("GeneDistance::distance; Genome: {}, Gene: {}, Gene: {}; {}, Gene Family: {}",
                      genome_db_ptr_->genomeId(), gene_ptr_->id(), node_ptr->gene_ptr_->id(), total_distance, protein_family_);

  return total_distance;

}


void kgl::GeneDistance::writeNode(std::ostream& outfile) const {

  outfile << std::format("{}_{}", genome_db_ptr_->genomeId(), gene_ptr_->id());

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares a single gene between isolates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::UPGMAATP4Distance::writeNode(std::ostream& outfile) const {

  auto transcript_array_ptr = kgl::GeneFeature::getTranscriptionSequences(gene_ptr_);
  if (transcript_array_ptr->empty()) {

    ExecEnv::log().critical("write_node(), Gene contains no coding transcript_ptr : genome: {} gene: {}",
                            genome_db_ptr_->genomeId(), gene_ptr_->id());

  }

  const auto& transcript_ptr = transcript_array_ptr->getFirst();
  const auto& contig_ref_ptr = transcript_ptr->getGene()->contig_ref_ptr();
  const std::string& contig_id = contig_ref_ptr->contigId();
  const std::string& gene_id = transcript_ptr->getGene()->id();
  const std::string& transcript_id = transcript_ptr->getParent()->id();

  if (transcript_array_ptr->size() > 1) {

    ExecEnv::log().warn("Genome: {} gene: {} contains: {} sequences using transcript_ptr: {}",
                        genome_db_ptr_->genomeId(), gene_id, transcript_array_ptr->size(), transcript_id);

  }

  AminoSequence mutant_sequence;
  AminoSequence reference_sequence;

  auto contig_db_opt = genome_db_ptr_->getContig(contig_id);
  if (not contig_db_opt) {

    ExecEnv::log().warn("Contig: {} not found for Genome: {}", contig_id, genome_db_ptr_->genomeId());
    return;

  }
  const auto& contig_db_ptr = contig_db_opt.value();

  const SequenceTranscript modified_transcript(contig_db_ptr, transcript_ptr);

  auto reference_sequence_opt = modified_transcript.getOriginalCoding();
  auto mutant_sequence_opt = modified_transcript.getModifiedCoding();

  if (modified_transcript.sequenceStatus() and reference_sequence_opt and mutant_sequence_opt) {

    auto protein_reference = contig_ref_ptr->getAminoSequence(reference_sequence_opt.value());
    auto protein_mutant = contig_ref_ptr->getAminoSequence(mutant_sequence_opt.value());

  } else {

    ExecEnv::log().critical("Cannot mutate transcript_ptr for : genome: {} gene: {}",
                            genome_db_ptr_->genomeId(), gene_ptr_->id());

  }

  std::stringstream ss;

  ss << genome_db_ptr_->genomeId();

  outfile << ss.str();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares reference (unmutated) genes between gene familes (VAR, RIF, STEVOR, MAURER).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<const kgl::TranscriptionSequence>  kgl::ReferenceGeneDistance::getCodingSequence() {

  auto coding_seq_ptr = GeneFeature::getTranscriptionSequences(gene_ptr_);

  if (coding_seq_ptr->empty()) {

    ExecEnv::log().critical("Gene contains no coding sequence : gene: {}", gene_ptr_->id());

  }

  std::shared_ptr<const TranscriptionSequence> coding_sequence = coding_seq_ptr->getFirst();

  if (coding_seq_ptr->size() > 1) {

    std::string gene_id = coding_sequence->getGene()->id();
    std::string sequence_id = coding_sequence->getParent()->id();

    ExecEnv::log().warn("Gene: {} contains: {} CDS sequences. Using sequence: {}",
                        gene_ptr_->id(), coding_seq_ptr->size(), sequence_id);

  }

  return coding_sequence;

}


void kgl::ReferenceGeneDistance::writeNode(std::ostream& outfile) const {

  std::string alt_symbolic;

  alt_symbolic = "-*";


  outfile << gene_ptr_->id();

}


void  kgl::DNAGeneDistance::getExonSequence() {

  const std::shared_ptr<const ContigReference>& contig_ref_ptr = gene_ptr_->contig_ref_ptr();

  // Just get the first transcript of the gene.
  auto transcript_ptr = getCodingSequence();

  auto coding_seq_opt = contig_ref_ptr->codingSequence(transcript_ptr);
  if (not coding_seq_opt) {

    ExecEnv::log().warn("Unable create coding sequence from Contig: {},  Gene: {}, Transcript: {}",
                        contig_ref_ptr->contigId(),
                        transcript_ptr->getGene()->id(),
                        transcript_ptr->getParent()->id());
    return;
  }
  auto& dna_coding_sequence = coding_seq_opt.value();

  linear_sequence_ = DNA5SequenceLinear::downConvertToLinear(dna_coding_sequence);

}

kgl::DistanceType_t kgl::DNAGeneDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const DNAGeneDistance> node_ptr = std::dynamic_pointer_cast<const DNAGeneDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("DNAGeneDistance::distance; Unexpected error, could not down-cast node pointer to DNAGeneDistance");
    return 1.0;

  }

  bool verbose = false;
  if (verbose) {

    ExecEnv::log().info("calculateDistance(); Comparing | {}({}), {}({}) |; Gene Family: {}",
                        gene_ptr_->id(), linear_sequence_.length(),
                        node_ptr->gene_ptr_->id(), node_ptr->linear_sequence_.length(), protein_family_);
  }


  CompareDistance_t contig_score = sequence_distance_(linear_sequence_, node_ptr->linear_sequence_);

  DistanceType_t total_distance = static_cast<DistanceType_t>(contig_score);

  if (verbose) {

    ExecEnv::log().info("calculateDistance(); {}({}), {}({}) |  =  {}; Gene Family: {}",
                        gene_ptr_->id(), linear_sequence_.length(),
                        node_ptr->gene_ptr_->id(), node_ptr->linear_sequence_.length(),
                        total_distance, protein_family_);

  }


  return total_distance;

}



void  kgl::AminoGeneDistance::getAminoSequence() {

  std::shared_ptr<const ContigReference> contig_ref_ptr = gene_ptr_->contig_ref_ptr();

  // Just get the first transcript of the gene.
  auto transcript_ptr = getCodingSequence();

  auto coding_seq_opt = contig_ref_ptr->codingSequence(transcript_ptr);

  if (not coding_seq_opt) {

    ExecEnv::log().warn("Unable create coding sequence from Contig: {},  Gene: {}, Transcript: {}",
                        contig_ref_ptr->contigId(),
                        transcript_ptr->getGene()->id(),
                        transcript_ptr->getParent()->id());
    return;
  }
  auto& dna_coding_sequence = coding_seq_opt.value();

  amino_sequence_ = contig_ref_ptr->getAminoSequence(dna_coding_sequence);

}


kgl::DistanceType_t kgl::AminoGeneDistance::distance(std::shared_ptr<const VirtualDistanceNode>  distance_node) const {

  std::shared_ptr<const AminoGeneDistance> node_ptr = std::dynamic_pointer_cast<const AminoGeneDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("AminoGeneDistance::distance; Unexpected error, could not down-cast node pointer to DNAGeneDistance");
    return 1.0;

  }

  CompareDistance_t contig_score = sequence_distance_(amino_sequence_, node_ptr->amino_sequence_);

  DistanceType_t total_distance = static_cast<DistanceType_t>(contig_score);


  return total_distance;

}
