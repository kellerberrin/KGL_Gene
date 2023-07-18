//
// Created by kellerberrin on 15/07/23.
//

#ifndef KGL_MUTATION_ANALYSIS_H
#define KGL_MUTATION_ANALYSIS_H


#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Mutation analysis and statistics particularly with respect to offsets with multiple minor alleles.
// Statistics are by Gene/Transcription and Genome/Contig.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Object to return and hold multi-threaded mutation results.
class TranscriptMutateRecord {

public:

  TranscriptMutateRecord(std::shared_ptr<const GeneFeature> gene_ptr,
                         std::shared_ptr<const TranscriptionSequence> transcription_ptr,
                         size_t total_variants,
                         size_t multiple_variants,
                         size_t mutated_genomes,
                         size_t duplicate_genomes,
                         size_t total_genomes) : gene_ptr_(std::move(gene_ptr)),
                                                 transcription_ptr_(std::move(transcription_ptr)),
                                                 total_variants_(total_variants),
                                                 multiple_variants_(multiple_variants),
                                                 mutated_genomes_(mutated_genomes),
                                                 duplicate_genomes_(duplicate_genomes),
                                                 total_genomes_(total_genomes) {}

  TranscriptMutateRecord(const TranscriptMutateRecord &) = default;
  ~TranscriptMutateRecord() = default;

  [[nodiscard]] const std::shared_ptr<const GeneFeature> &genePtr() const { return gene_ptr_; }
  [[nodiscard]] const std::shared_ptr<const TranscriptionSequence> &transcriptionPtr() const { return transcription_ptr_; }
  [[nodiscard]] size_t totalVariants() const { return total_variants_; }
  [[nodiscard]] size_t multipleVariants() const { return multiple_variants_; }
  [[nodiscard]] size_t mutatedGenomes() const { return mutated_genomes_; }
  [[nodiscard]] size_t duplicateGenomes() const { return duplicate_genomes_; }
  [[nodiscard]] size_t totalGenomes() const { return total_genomes_; }

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const TranscriptionSequence> transcription_ptr_;
  size_t total_variants_{0};
  size_t multiple_variants_{0};
  size_t mutated_genomes_{0};
  size_t duplicate_genomes_{0};
  size_t total_genomes_{0};

};


class GenomeContigMutate {

public:

  GenomeContigMutate(GenomeId_t genome_id,
                     size_t total_variants,
                     size_t multiple_variants) : genome_id_(std::move(genome_id)),
                                                 total_variants_(total_variants),
                                                 multiple_variants_(multiple_variants){}
  GenomeContigMutate(const GenomeContigMutate &) = default;
  ~GenomeContigMutate() = default;

  [[nodiscard]] const GenomeId_t &genomeId() const { return genome_id_; }
  [[nodiscard]] size_t totalVariants() const { return total_variants_; }
  [[nodiscard]] size_t multipleVariants() const { return multiple_variants_; }
  void addRecord(const GenomeContigMutate& record) {

    total_variants_ += record.total_variants_;
    multiple_variants_ += record.multiple_variants_;

  }

private:

  GenomeId_t genome_id_;
  size_t total_variants_;
  size_t multiple_variants_;

};


using TranscriptRecordMap = std::map<std::string, TranscriptMutateRecord>;
using GenomeRecordMap = std::map<GenomeId_t , GenomeContigMutate>;

class MutateAnalysis {

public:

  MutateAnalysis() = default;
  ~MutateAnalysis() = default;

  void addTranscriptRecord(const TranscriptMutateRecord &record);
  void addGenomeRecords(const GenomeContigMutate &record);

  [[nodiscard]] const TranscriptRecordMap &getTranscriptMap() const { return transcript_map_; }
  [[nodiscard]] const GenomeRecordMap &getGenomeMap() const { return genome_map_; }

  void printMutationTranscript(const std::string &file_name) const;
  void printGenomeContig(const std::string &file_name) const;

private:

  constexpr static const char DELIMITER_ = ',';
  TranscriptRecordMap transcript_map_;
  GenomeRecordMap genome_map_;

};


} // Namespace



#endif //KGL_MUTATION_ANALYSIS_H
