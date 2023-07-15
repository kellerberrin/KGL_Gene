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

  if (genome_map_.empty()) {

    return;

  }

  // Write header.
  auto first_iter = genome_map_.begin();
  auto const [map_key, map_record] = *first_iter;

  write_file << "Genome";
  for (auto const& [contig_key, record] : map_record.contigMap()) {

    write_file << DELIMITER_ << contig_key
               << DELIMITER_ << "Total Variants"
               << DELIMITER_ << "Multiple Variants";

  }
  write_file << '\n';

  // Write Genome Data.
  for (auto const& [genome_id, genome_record] : genome_map_) {

    write_file << genome_id;
    for (auto const& [contig_id, contig_record] : genome_record.contigMap()) {

      write_file << DELIMITER_ << contig_id
                 << DELIMITER_ << contig_record.totalVariants()
                 << DELIMITER_ << contig_record.multipleVariants();

    }
    write_file << '\n';

  }

}

