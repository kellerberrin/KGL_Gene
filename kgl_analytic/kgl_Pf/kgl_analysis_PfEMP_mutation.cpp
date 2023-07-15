//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_analysis_PfEMP_mutation.h"
#include "kgl_variant_filter_features.h"
#include "kgl_variant_filter_db.h"

#include <ranges>
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
           << "Duplicate Variants"<< '\n';

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
             << transcript_record.multipleVariants() << '\n';

  }

}


kgl::MultipleVariantMap kgl::MutateAnalysis::multipleOffsetVariants(const std::shared_ptr<const PopulationDB>& population_ptr) {

  MultipleVariantMap multiple_variant_map;

  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      for (auto const& [offset, offset_ptr] : contig_ptr->getMap()) {

        if (offset_ptr->getVariantArray().size() > 1) {

          std::map<std::string, std::shared_ptr<const Variant>> unique_variant_map;
          for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

            // Overwrites identical (homozygous) variants.
            unique_variant_map[variant_ptr->HGVS()] = variant_ptr;

          }

          // If we have different unique variants.
          if (unique_variant_map.size() > 1) {

            // Create a lexicographic unique key using the map ordering.
            std::string unique_key;
            for (auto const& [hgvs, variant_ptr] : unique_variant_map) {

              unique_key += hgvs;
              unique_key += '|';

            }

            auto find_iter = multiple_variant_map.find(unique_key);
            // Insert new record.
            if (find_iter == multiple_variant_map.end()) {

              std::vector<std::shared_ptr<const Variant>> variant_vector;
              for (auto const& [unique_hgvs, unique_variant_ptr] : unique_variant_map) {

                variant_vector.push_back(unique_variant_ptr);

              }

              std::pair<std::vector<std::shared_ptr<const Variant>>, std::vector<GenomeId_t>> value_pair{ variant_vector, std::vector<GenomeId_t>()};
              auto [insert_iter, result] = multiple_variant_map.try_emplace(unique_key, value_pair);
              if (not result) {

                ExecEnv::log().error("MutateGenes::multipleOffsetVariants; Unexpected, unable to insert key: {}", unique_key);
                continue;

              }

              find_iter = insert_iter;

            } // Insert New Record.

            // Add Genome.
            auto& [hgvs_key, value_pair] = *find_iter;
            auto& [variant_vector, genome_vector] = value_pair;
            genome_vector.push_back(genome_id);

          } // Multiple Unique Variants?

        } // Offset Size > 1

      } // All Offsets

    } // All Contigs

  } // All Genomes

  return multiple_variant_map;

}

//
std::multimap<double, std::tuple<kgl::GenomeId_t, size_t, size_t>>
kgl::MutateAnalysis::genomeSummary( const MultipleVariantMap& multiple_variant_map,
                                    const std::shared_ptr<const PopulationDB>& population_ptr) {

  // .first is total genome variant .second is the number of multiple minor alleles.
  std::map<kgl::GenomeId_t, std::pair<size_t, std::string>> genome_summary;

  for (auto const& [hgvs_key, value_pair] : multiple_variant_map) {

    auto const& [variant_vector, genome_vector] = value_pair;
    for (auto const& genome_id : genome_vector) {

      auto find_iter = genome_summary.find(genome_id);
      if (find_iter == genome_summary.end()) {

        genome_summary[genome_id] = {1, hgvs_key};

      } else {

        auto& [genome_key, genome_pair] = *find_iter;
        auto& [genome_count, key] = genome_pair;
        ++genome_count;

      }

    }

  }

  std::multimap<double, std::tuple<kgl::GenomeId_t, size_t, size_t>> genome_count_multimap;
  for (auto const& [genome_id, genome_count] : genome_summary) {


    auto genome_opt = population_ptr->getGenome(genome_id);
    if (not genome_opt) {

      ExecEnv::log().error("MutateGenes::genomeSummary; cannor fund genome: {}", genome_id);
      continue;

    }

    auto const& [count, key] = genome_count;
    size_t variant_count = genome_opt.value()->variantCount();
    std::tuple<kgl::GenomeId_t, size_t, size_t> value{genome_id, count, variant_count};
    double proportion{0.0};
    if (variant_count > 0) {

      proportion = static_cast<double>(count) / static_cast<double>(variant_count);

    }

    genome_count_multimap.insert({proportion, value});

  }

  return genome_count_multimap;

}

void kgl::MutateAnalysis::printGenomeContig(const std::string& file_name) const {

  std::ofstream write_file(file_name);

  if (not write_file.good()) {

    ExecEnv::log().error("MutateAnalysis::printGenomeContig; unable to open output file: {}", file_name);
    return;

  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Used to find multiple different variants at a particular Genome/Contig/Offset.
// Multiple variants indicate different minor alleles at the same offset.
// P.Falciparum is haploid at the blood stage, so this indicates complexity of infection
// or an issue with the generation of the VCF file.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::MutateGenes::initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [gene_offset, gene_ptr] : contig_ptr->getGeneMap()) {

      auto [iter, result] = gene_contig_map_.try_emplace(gene_ptr, contig_id);
      if (not result) {

        ExecEnv::log().warn("MutateGenes::initializeGeneContigMap; Unable to insert duplicate gene: {}", gene_ptr->id());

      }

    }

  }

}


std::vector<std::shared_ptr<const kgl::GeneFeature>> kgl::MutateGenes::contigGenes(const ContigId_t& contig_id) const {

  std::vector<std::shared_ptr<const GeneFeature>> gene_vector;

  for (auto const& [gene_ptr, gene_contig_id] : gene_contig_map_) {

    if (gene_contig_id == contig_id) {

      gene_vector.push_back(gene_ptr);

    }

  }

  return gene_vector;

}

// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
std::pair<size_t, size_t> kgl::MutateGenes::mutateTranscript( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                              const FeatureIdent_t& transcript_id,
                                                              const std::shared_ptr<const PopulationDB>& population_ptr,
                                                              const std::shared_ptr<const GenomeReference>& reference_genome) const {

  // Remove homozygous variants
  std::shared_ptr<const PopulationDB> unique_population_ptr = population_ptr->viewFilter(UniqueUnphasedFilter());
  // Only variants for the gene and transcript.
  std::shared_ptr<const PopulationDB> gene_population_ptr = unique_population_ptr->viewFilter(FilterGeneTranscriptVariants(gene_ptr, transcript_id));
  // Convert to a population with canonical variants.
  std::shared_ptr<const PopulationDB> mut_population_ptr = gene_population_ptr->canonicalPopulation();

  // Checked the structure of the canonical gene population.
  auto const [variant_count, validated_variants] = mut_population_ptr->validate(reference_genome);
  if (variant_count != validated_variants) {

    ExecEnv::log().error("MutateGenes::mutateTranscript; Failed Reference Genome Validation, Canonical Population Variant: {}, Canonical Validated: {}",
                         variant_count, validated_variants);

  }

  size_t mut_active_genome{0};
  for (auto const& [genome_id, genome_ptr] : mut_population_ptr->getMap()) {

    if (genome_ptr->variantCount() > 0) {

      ++mut_active_genome;

    }

  }

  size_t gene_active_genome{0};
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    if (genome_ptr->variantCount() > 0) {

      ++gene_active_genome;

    }

  }

  auto const [total_variants, duplicate_variants] = mutateGenomes(gene_ptr, transcript_id, mut_population_ptr);

  ExecEnv::log().info("MutateGenes::mutateTranscript; Filtered Gene: {}, Transcript: {}, Description: {}",
                      gene_ptr->id(), transcript_id, gene_ptr->descriptionText());
  ExecEnv::log().info("MutateGenes::mutateTranscript; Total Variants: {}, Processed Variants: {}, Duplicate Variants: {}, Total Genomes: {}, Canonical Mutant Genomes: {}, Mutant Genomes: {}",
                      gene_population_ptr->variantCount(),
                      total_variants,
                      duplicate_variants,
                      gene_population_ptr->getMap().size(),
                      mut_active_genome,
                      gene_active_genome);

  return {total_variants, duplicate_variants};

}


// Multi-tasked filtering for large populations.
// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
std::pair<size_t, size_t> kgl::MutateGenes::mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                           const FeatureIdent_t& transcript_id,
                                                           const std::shared_ptr<const PopulationDB>& gene_population_ptr) const {

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  // the tuple is bool, genome_id, genome_id, transcript_id,
  std::vector<std::future<std::pair<size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<std::pair<size_t, size_t>> future = thread_pool.enqueueFuture(&MutateGenes::threadMutation, genome_ptr, gene_ptr, transcript_id);
    future_vector.push_back(std::move(future));

  }


  size_t total_variants{0};
  size_t duplicate_variants{0};
  // Wait for all the threads to return.
  for (auto& future : future_vector) {

    auto const [genome_total, genome_duplicate] = future.get();
    total_variants += genome_total;
    duplicate_variants += genome_duplicate;

  }

  return {total_variants, duplicate_variants};

}

// .first total variants, .second multiple (duplicate) variants per offset.
std::pair<size_t, size_t> kgl::MutateGenes::threadMutation( std::shared_ptr<const GenomeDB> genome_ptr,
                                                            const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                            const FeatureIdent_t& transcript_id) {

  class ProcessVariants {

  public:

    ProcessVariants(std::shared_ptr<const GenomeDB> genome_ptr,
                    std::shared_ptr<const GeneFeature> gene_ptr) : genome_ptr_(std::move(genome_ptr)), gene_ptr_(std::move(gene_ptr)) {}
    ~ProcessVariants() = default;

    bool processVariants(const std::shared_ptr<const Variant>& variant_ptr) {

      auto [ref_seq, alt_seq, canonical_offset] = variant_ptr->canonicalSequences();

      auto [insert_iter, result] = ordered_variant_map_.try_emplace(canonical_offset, variant_ptr);
      if (not result) {
        // Duplicate variant at offset.

        ++multiple_variants_;

      }

      return true;

    }

    size_t multipleVariants() const { return multiple_variants_; }
    std::map<ContigOffset_t, std::shared_ptr<const Variant>> ordered_variant_map_;

  private:

    std::shared_ptr<const GenomeDB> genome_ptr_;
    std::shared_ptr<const GeneFeature> gene_ptr_;
    size_t multiple_variants_{0};

  };

  ProcessVariants process_object(genome_ptr, gene_ptr);
  genome_ptr->processAll(process_object, &ProcessVariants::processVariants);

  return {genome_ptr->variantCount(), process_object.multipleVariants() };

}

