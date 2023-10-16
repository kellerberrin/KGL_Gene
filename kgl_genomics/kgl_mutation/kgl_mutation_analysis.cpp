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

    ExecEnv::log().error("Could not insert transcript record: {} (duplicate)", map_key);

  }

}


void kgl::MutateAnalysis::addGenomeRecords(const GenomeId_t& genome_id,
                                           size_t total_variants,
                                           size_t multiple_variants,
                                           CodingSequenceValidity modified,
                                           CodingSequenceValidity original) {

  auto find_iter = genome_map_.find(genome_id);
  if (find_iter == genome_map_.end()) {

    auto [insert_iter, result] = genome_map_.try_emplace(genome_id, GenomeContigMutate(genome_id));
    if (not result) {

      ExecEnv::log().error("Cannot add genome: {} (duplicate)", genome_id);
      return;

    }
    find_iter = insert_iter;

  }

  auto& [map_genome, genome_contig_mutate] = *find_iter;
  genome_contig_mutate.addRecord(total_variants, multiple_variants, modified, original);

}


void kgl::MutateAnalysis::printMutationTranscript(const std::string& file_name) const {

  std::ofstream out_file(file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("Unable to open output file: {}", file_name);
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
           << "Total SNP" << DELIMITER_
           << "Total Frameshift" << DELIMITER_
           << "Duplicate Variants"<< DELIMITER_
           << "Duplicate Genomes"<< DELIMITER_
           << "Upstream Variants"<< DELIMITER_
           << "Upstream Genomes"<< DELIMITER_
           << "Mutated_Genomes" << DELIMITER_
           << "Total_Geromes" << '\n';

  for (auto const& [key, transcript_record] : transcript_map_) {

    out_file << transcript_record.genePtr()->contig_ref_ptr()->contigId() << DELIMITER_
             << transcript_record.genePtr()->id() << DELIMITER_
             << transcript_record.transcriptionPtr()->getParent()->id() << DELIMITER_
        << (GeneFeature::proteinCoding(transcript_record.genePtr()) ? GeneFeature::PROTEIN_CODING_GENE_ : GeneFeature::NCRNA_GENE_)
        << DELIMITER_
        << transcript_record.genePtr()->descriptionText() << DELIMITER_
        << transcript_record.transcriptionPtr()->start() << DELIMITER_
        << transcript_record.transcriptionPtr()->end() << DELIMITER_
        << transcript_record.transcriptionPtr()->getFeatureMap().size() << DELIMITER_
        << transcript_record.transcriptionPtr()->codingNucleotides() << DELIMITER_
        << transcript_record.mutateStats().total_variants_ << DELIMITER_
        << transcript_record.mutateStats().total_snp_variants_ << DELIMITER_
        << transcript_record.mutateStats().total_frameshift_variants_ << DELIMITER_
        << transcript_record.mutateStats().duplicate_variants_ << DELIMITER_
        << transcript_record.mutateStats().duplicate_genomes_ << DELIMITER_
        << transcript_record.mutateStats().upstream_delete_variants_ << DELIMITER_
        << transcript_record.mutateStats().upstream_delete_genomes_ << DELIMITER_
             << transcript_record.mutateStats().mutant_genomes_ << DELIMITER_
             << transcript_record.mutateStats().total_genomes_ << '\n';

  }

}


void kgl::MutateAnalysis::printMutationValidity(const std::string& file_name) const {

  std::ofstream out_file(file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("Unable to open output file: {}", file_name);
    return;

  }

  out_file << "Contig" << DELIMITER_
           << "Gene" << DELIMITER_
           << "Transcript" << DELIMITER_
           << "Gene Type" << DELIMITER_
           << "Description" << DELIMITER_
           << "Transcript Begin" << DELIMITER_
           << "Transcript Size" << DELIMITER_
           << "Total Variants" << DELIMITER_
           << "Total Sequences" << DELIMITER_
           << "NCRNA Sequences" << DELIMITER_
           << "Valid Original" << DELIMITER_
           << "Valid Modified" << DELIMITER_
           << "Modified Not Mod3" << DELIMITER_
           << "Modified No Start Codon" << DELIMITER_
           << "Modified No Stop Codon" << DELIMITER_
           << "Modified Nonsense Mutation" << DELIMITER_
           << "Total_Geromes" << '\n';

  for (auto const& [key, transcript_record] : transcript_map_) {

    out_file << transcript_record.genePtr()->contig_ref_ptr()->contigId() << DELIMITER_
             << transcript_record.genePtr()->id() << DELIMITER_
    << transcript_record.transcriptionPtr()->getParent()->id() << DELIMITER_
    << (GeneFeature::proteinCoding(transcript_record.genePtr()) ? GeneFeature::PROTEIN_CODING_GENE_ : GeneFeature::NCRNA_GENE_)
    << DELIMITER_
    << transcript_record.genePtr()->descriptionText() << DELIMITER_
    << transcript_record.transcriptionPtr()->start() << DELIMITER_
    << transcript_record.transcriptionPtr()->codingNucleotides() << DELIMITER_
    << transcript_record.mutateStats().total_variants_ << DELIMITER_
    << transcript_record.mutateStats().modified_validity_.totalSequence() << DELIMITER_
    << transcript_record.mutateStats().modified_validity_.ncRNA() << DELIMITER_
    << transcript_record.mutateStats().original_validity_.validProtein() << DELIMITER_
    << transcript_record.mutateStats().modified_validity_.validProtein()  << DELIMITER_
    << transcript_record.mutateStats().modified_validity_.notMod3() << DELIMITER_
    << transcript_record.mutateStats().modified_validity_.noStartCodon() << DELIMITER_
    << transcript_record.mutateStats().modified_validity_.noStopCodon() << DELIMITER_
             << transcript_record.mutateStats().modified_validity_.nonsenseMutation() << DELIMITER_
         << transcript_record.mutateStats().total_genomes_ << '\n';

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


