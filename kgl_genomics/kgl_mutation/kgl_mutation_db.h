//
// Created by kellerberrin on 6/07/23.
//

#ifndef KGL_DB_MUTATION_H
#define KGL_DB_MUTATION_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"
#include "kgl_mutation_analysis.h"


namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Map holds unique canonical variants for a specified region or transcript.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::map<ContigOffset_t, std::shared_ptr<const Variant>>;
// This object contains a map of unique canonical variants ready to mutate the specified transcript.
class UniqueMutationTranscript {

public:

  UniqueMutationTranscript(GenomeId_t genome_id,
                           ContigId_t contig_id,
                           FeatureIdent_t gene_id,
                           FeatureIdent_t transcript_id,
                           OffsetVariantMap unique_mutation_map) : genome_id_(std::move(genome_id)),
                                                                   contig_id_(std::move(contig_id)),
                                                                   gene_id_(std::move(gene_id)),
                                                                   transcript_id_(std::move(transcript_id)),
                                                                   unique_mutation_map_(std::move(unique_mutation_map)) {}
  UniqueMutationTranscript(const UniqueMutationTranscript&) = default;
  ~UniqueMutationTranscript() = default;

  [[nodiscard]] const GenomeId_t& genomeId() const { return genome_id_; }
  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  [[nodiscard]] const FeatureIdent_t& geneId() const { return gene_id_; }
  [[nodiscard]] const FeatureIdent_t& transcriptId() const { return transcript_id_; }
  [[nodiscard]] const OffsetVariantMap& uniqueMutationMap() const { return unique_mutation_map_; }

private:

  GenomeId_t genome_id_;
  ContigId_t contig_id_;
  FeatureIdent_t gene_id_;
  FeatureIdent_t transcript_id_;
  OffsetVariantMap unique_mutation_map_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Used to find multiple different variants at a particular Genome/Contig/Offset.
// Multiple variants indicate different minor alleles at the same offset.
// P.Falciparum is haploid at the blood stage, so this indicates complexity of infection
// or an issue with the generation of the VCF file.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



using GeneContigMap = std::map<std::shared_ptr<const GeneFeature>, ContigId_t>;
class MutateGenes {

public:

  explicit MutateGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) : genome_ptr_(genome_ptr) {

    initializeGeneContigMap(genome_ptr_);

  }
  ~MutateGenes() = default;
  // Return the id of all genes from a particular contig.
  std::vector<std::shared_ptr<const GeneFeature>> contigGenes(const ContigId_t& contig_id) const;
  // Mutate the population.
  void mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr);
  // Access the mutation analysis object.
  const MutateAnalysis& mutateAnalysis() const { return mutate_analysis_; }

private:

  std::shared_ptr<const GenomeReference> genome_ptr_;
  GeneContigMap gene_contig_map_;
  mutable MutateAnalysis mutate_analysis_;

  void initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr);

  // Process a gene and transcript.
  void mutateTranscript(const std::shared_ptr<const GeneFeature>& gene_ptr,
                        const FeatureIdent_t& transcript_id,
                        const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                        const std::shared_ptr<const PopulationDB>& population,
                        const std::shared_ptr<const GenomeReference>& reference_genome) const;

  // .first total variants, .second multiple (duplicate) variants per offset.
  static std::tuple<UniqueMutationTranscript, size_t, size_t>
  genomeTranscriptMutation( const std::shared_ptr<const GenomeDB>& genome_ptr,
                            const std::shared_ptr<const GeneFeature>& gene_ptr,
                            const FeatureIdent_t& transcript_id);

  // .first total variants, .second multiple (duplicate) variants per offset.
  static std::tuple<UniqueMutationTranscript, size_t, size_t>
  altGenomeTranscriptMutation( const std::shared_ptr<const GenomeDB>& genome_ptr,
                            const std::shared_ptr<const GeneFeature>& gene_ptr,
                            const FeatureIdent_t& transcript_id);

  // .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
  std::tuple<size_t, size_t, size_t, size_t> mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                    const FeatureIdent_t& transcript_id,
                                                    const std::shared_ptr<const PopulationDB>& gene_population_ptr) const;

};


} // Namespace.

#endif //KGL_DB_MUTATION_H
