//
// Created by kellerberrin on 15/07/23.
//

#include "kgl_mutation_analysis.h"

#include <fstream>


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Mutation analysis and statistics particularly with respect to offsets with multiple minor alleles.
// Statistics are by Gene Transcription and Genome/Contig.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::MutateAnalysis::addTranscriptRecord(const TranscriptMutateRecord& record) {

  std::string map_key = record.genePtr()->id() + '|' + record.transcriptionPtr()->getParent()->id();

  auto [insert_iter, result] = transcript_map_.try_emplace(map_key, record);
  if (not result) {

    ExecEnv::log().error("MutateAnalysis::addTranscriptRecord; could not insert transcript record: {} (duplicate)", map_key);

  }

}


void kgl::MutateAnalysis::addGenomeRecords(const GenomeContigMutate& genome_mutate) {

  auto find_iter = genome_map_.find(genome_mutate.genomeId());
  if (find_iter == genome_map_.end()) {

    auto [insert_iter, result] = genome_map_.try_emplace(genome_mutate.genomeId(), genome_mutate);
    if (not result) {

      ExecEnv::log().error("MutateAnalysis::addGenomeRecords; cannot add genome: {} (duplicate)", genome_mutate.genomeId());

    }


  }

  auto& [map_genome, genome_contig_mutate] = *find_iter;
  genome_contig_mutate.addRecord(genome_mutate);


}


void kgl::MutateAnalysis::printMutationTranscript(const std::string& file_name) const {

  std::ofstream out_file(file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("utateAnalysis::printMutationTranscript; unable to open output file: {}", file_name);
    return;

  }

  out_file << "Contig" << DELIMITER_
           << "Gene" << DELIMITER_
           << "Transcript" << DELIMITER_
           << "Gene Type" << DELIMITER_
           << "Description" << DELIMITER_
           << "Transcript Begin" << DELIMITER_
           << "Transcript End" << DELIMITER_
           << "Transcript CDS" << DELIMITER_
           << "Transcript Size" << DELIMITER_
           << "Total Variants" << DELIMITER_
           << "Duplicate Variants"<< DELIMITER_
           << "Duplicate Genomes"<< DELIMITER_
           << "Mutated_Genomes" << DELIMITER_
           << "Total_Geromes" << '\n';

  for (auto const& [key, transcript_record] : transcript_map_) {

    out_file << transcript_record.genePtr()->contig()->contigId() << DELIMITER_
             << transcript_record.genePtr()->id() << DELIMITER_
             << transcript_record.transcriptionPtr()->getParent()->id() << DELIMITER_
             << (GeneFeature::proteinCoding(transcript_record.genePtr()) ? GeneFeature::PROTEIN_CODING_GENE_ : GeneFeature::NCRNA_GENE_)
             << DELIMITER_
             << transcript_record.genePtr()->descriptionText() << DELIMITER_
             << transcript_record.transcriptionPtr()->start() << DELIMITER_
             << transcript_record.transcriptionPtr()->end() << DELIMITER_
             << transcript_record.transcriptionPtr()->codingFeatures() << DELIMITER_
             << transcript_record.transcriptionPtr()->codingNucleotides() << DELIMITER_
             << transcript_record.totalVariants() << DELIMITER_
             << transcript_record.multipleVariants() << DELIMITER_
             << transcript_record.duplicateGenomes() << DELIMITER_
             << transcript_record.mutatedGenomes() << DELIMITER_
             << transcript_record.totalGenomes() << '\n';

  }

}


void kgl::MutateAnalysis::printGenomeContig(const std::string& file_name) const {

  std::ofstream write_file(file_name);

  if (not write_file.good()) {

    ExecEnv::log().error("MutateAnalysis::printGenomeContig; unable to open output file: {}", file_name);
    return;

  }

  // Write header.
  write_file << "Genome"  << DELIMITER_
             << "Total Variants" << DELIMITER_
             << "Multiple Variants" << '\n';

  // Write Genome Data.
  for (auto const& [genome_id, genome_record] : genome_map_) {

      write_file << genome_record.genomeId() << DELIMITER_
                 << genome_record.totalVariants() << DELIMITER_
                 << genome_record.multipleVariants() << '\n';

  }

}


