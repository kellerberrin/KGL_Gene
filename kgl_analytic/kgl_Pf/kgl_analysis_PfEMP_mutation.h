//
// Created by kellerberrin on 6/07/23.
//

#ifndef KGL_ANALYSIS_PFEMP_MUTATION_H
#define KGL_ANALYSIS_PFEMP_MUTATION_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace

// Object to return and hold multi-threaded mutation results.
class TranscriptMutateRecord {

public:

  TranscriptMutateRecord(std::shared_ptr<const GeneFeature> gene_ptr,
                         std::shared_ptr<const TranscriptionSequence> transcription_ptr,
                         size_t total_variants,
                         size_t multiple_variants) : gene_ptr_(std::move(gene_ptr)),
                                                    transcription_ptr_(std::move(transcription_ptr)),
                                                    total_variants_(total_variants),
                                                    multiple_variants_(multiple_variants) {}
  TranscriptMutateRecord(const TranscriptMutateRecord&) = default;
  ~TranscriptMutateRecord() = default;

  [[nodiscard]] const std::shared_ptr<const GeneFeature>& genePtr() const { return gene_ptr_; }
  [[nodiscard]] const std::shared_ptr<const TranscriptionSequence>& transcriptionPtr() const { return transcription_ptr_; }
  [[nodiscard]] size_t totalVariants() const { return total_variants_; }
  [[nodiscard]] size_t multipleVariants() const { return multiple_variants_; }

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const TranscriptionSequence> transcription_ptr_;
  size_t total_variants_{0};
  size_t multiple_variants_{0};

};


using TranscriptRecordMap = std::map<std::string, TranscriptMutateRecord>;
class MutateAnalysis {

public:

  MutateAnalysis() = default;
  ~MutateAnalysis() = default;

  void addTranscriptRecord(const TranscriptMutateRecord& record);

  [[nodiscard]] const TranscriptRecordMap& getTranscriptMap() const { return transcript_map_; }

  void printMutationAnalysis(const std::string& file_name) const;

private:

  constexpr static const char DELIMITER_ = ',';
  TranscriptRecordMap transcript_map_;

};



// Used to find multiple different variants at a particular Genome/Contig/Offset.
// Multiple variants indicate different minor alleles at the same offset.
// P.Falciparum is haploid at the blood stage, so this indicates complexity of infection
// or an issue with the generation of the VCF file.
using MultipleVariantMap = std::map<std::string, std::pair<std::vector<std::shared_ptr<const Variant>>, std::vector<GenomeId_t>>>;

using GeneContigMap = std::map<std::shared_ptr<const GeneFeature>, ContigId_t>;
class MutateGenes {

public:

  explicit MutateGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) : genome_ptr_(genome_ptr) {

    initializeGeneContigMap(genome_ptr_);

  }
  ~MutateGenes() = default;
  // Return the id of all genes from a particular contig.
  std::vector<std::shared_ptr<const GeneFeature>> contigGenes(const ContigId_t& contig_id) const;
  // Process a gene and transcript.
  // .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
  std::pair<size_t, size_t> mutateTranscript(const std::shared_ptr<const GeneFeature>& gene_ptr,
                                             const FeatureIdent_t& transcript,
                                             const std::shared_ptr<const PopulationDB>& population,
                                             const std::shared_ptr<const GenomeReference>& reference_genome) const;

  static MultipleVariantMap multipleOffsetVariants(const std::shared_ptr<const PopulationDB>& population);
  static std::multimap<double, std::tuple<GenomeId_t, size_t, size_t>>
    genomeSummary( const MultipleVariantMap& multiple_variant_map,
                   const std::shared_ptr<const PopulationDB>& population_ptr);

  static void printMultipleAllele(const std::shared_ptr<const PopulationDB>& population_ptr, const std::string& file_name);

private:

  std::shared_ptr<const GenomeReference> genome_ptr_;
  GeneContigMap gene_contig_map_;

  void initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr);

  // .first total variants, .second multiple (duplicate) variants per offset.
  static std::pair<size_t, size_t> threadMutation(std::shared_ptr<const GenomeDB> genome_ptr,
                                           const std::shared_ptr<const GeneFeature>& gene_ptr,
                                           const FeatureIdent_t& transcript_id);
  // .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
  std::pair<size_t, size_t> mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                           const FeatureIdent_t& transcript_id,
                                           const std::shared_ptr<const PopulationDB>& gene_population_ptr) const;

};


} // Namespace.

#endif //KGL_KGL_ANALYSIS_PFEMP_MUTATION_H
