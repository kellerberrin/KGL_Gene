//
// Created by kellerberrin on 15/07/23.
//

#ifndef KGL_MUTATION_ANALYSIS_H
#define KGL_MUTATION_ANALYSIS_H


#include "kgl_genome_genome.h"


namespace kellerberrin::genome {   //  organization::project level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Mutation analysis and statistics particularly with respect to offsets with multiple minor alleles.
// Statistics are by Gene/Transcription and Genome/Contig.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Simple object to pass mutation statistics
struct MutateStats {

  size_t total_variants_{0};
  size_t total_snp_variants_{0};
  size_t total_frameshift_variants_{0};
  size_t total_genomes_{0};
  size_t mutant_genomes_{0};
  size_t duplicate_variants_{0};
  size_t duplicate_genomes_{0};
  size_t upstream_delete_variants_{0};
  size_t upstream_delete_genomes_{0};
  SequenceValidityStatistics modified_validity_;
  SequenceValidityStatistics original_validity_;

};



// Object to return and hold multi-threaded mutation results.
class TranscriptMutateRecord {

public:

  TranscriptMutateRecord(std::shared_ptr<const GeneFeature> gene_ptr,
                         std::shared_ptr<const TranscriptionSequence> transcription_ptr,
                         MutateStats mutate_stats) : gene_ptr_(std::move(gene_ptr)),
                                                     transcription_ptr_(std::move(transcription_ptr)),
                                                     mutate_stats_(mutate_stats) {}

  TranscriptMutateRecord(const TranscriptMutateRecord &) = default;
  ~TranscriptMutateRecord() = default;

  [[nodiscard]] const std::shared_ptr<const GeneFeature> &genePtr() const { return gene_ptr_; }
  [[nodiscard]] const std::shared_ptr<const TranscriptionSequence> &transcriptionPtr() const { return transcription_ptr_; }
  [[nodiscard]] const MutateStats& mutateStats() const { return mutate_stats_; }

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const TranscriptionSequence> transcription_ptr_;
  MutateStats mutate_stats_;

};


class GenomeContigMutate {

public:

  explicit GenomeContigMutate(GenomeId_t genome_id) : genome_id_(std::move(genome_id)) {}
  GenomeContigMutate(const GenomeContigMutate &) = default;
  ~GenomeContigMutate() = default;

  [[nodiscard]] const GenomeId_t &genomeId() const { return genome_id_; }
  [[nodiscard]] size_t totalVariants() const { return total_variants_; }
  [[nodiscard]] size_t multipleVariants() const { return multiple_variants_; }
  [[nodiscard]] const SequenceValidityStatistics& modifiedValidity() const { return modified_validity_; }
  [[nodiscard]] const SequenceValidityStatistics& originalValidity() const { return original_validity_; }

  void addRecord(size_t total_variants,
                 size_t multiple_variants,
                 CodingSequenceValidity modified,
                 CodingSequenceValidity original) {

    total_variants_ += total_variants;
    multiple_variants_ += multiple_variants;
    modified_validity_.updateValidity(modified);
    original_validity_.updateValidity(original);

  }

private:

  GenomeId_t genome_id_;
  size_t total_variants_{0};
  size_t multiple_variants_{0};
  SequenceValidityStatistics modified_validity_;
  SequenceValidityStatistics original_validity_;

};


using TranscriptRecordMap = std::map<std::string, TranscriptMutateRecord>;
using GenomeRecordMap = std::map<GenomeId_t , GenomeContigMutate>;

class MutateAnalysis {

public:

  MutateAnalysis() = default;
  ~MutateAnalysis() = default;

  void addTranscriptRecord(const TranscriptMutateRecord &record);
  void addGenomeRecords(const GenomeId_t& genome_id,
                        size_t total_variants,
                        size_t multiple_variants,
                        CodingSequenceValidity modified,
                        CodingSequenceValidity original);

  [[nodiscard]] const TranscriptRecordMap &getTranscriptMap() const { return transcript_map_; }
  [[nodiscard]] const GenomeRecordMap &getGenomeMap() const { return genome_map_; }

  void printMutationTranscript(const std::string &file_name) const;
  void printMutationValidity(const std::string& file_name) const;
  void printGenomeContig(const std::string &file_name) const;

private:

  constexpr static const char DELIMITER_ = ',';
  TranscriptRecordMap transcript_map_;
  GenomeRecordMap genome_map_;

};


} // Namespace



#endif //KGL_MUTATION_ANALYSIS_H
