//
// Created by kellerberrin on 9/04/23.
//


#include "kgl_analysis_PfEMP_variant.h"
#include "kgl_variant_filter.h"

#include <fstream>



void kgl::GenomeGeneVariantAnalysis::setGeneVector(const GeneVector& gene_vector) {

  gene_vector_ = gene_vector;

}


// Get variants only occurring within codingFeatures for all mRNA sequences.
void kgl::GenomeGeneVariantAnalysis::getGeneVariants(const std::shared_ptr<const PopulationDB>& population_ptr) {

  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    // Get/Create the matching genome from the gene population

    auto gene_genome_opt = gene_population_ptr->getCreateGenome(genome_id);

    if (not gene_genome_opt) {

      ExecEnv::log().error("GenomeGeneVariantAnalysis::genePopulationVariants; Unable to get/create genome: {}", genome_id);
      continue;

    }

    auto gene_genome_ptr = gene_genome_opt.value();

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      if (contig_ptr->getMap().empty()) {

        continue;

      }

      for (auto const& gene_ptr : gene_vector_) {

        if (gene_ptr->contig()->contigId() == contig_id) {

          // Get/Create the gene contig from the gene population.
          auto gene_contig_opt = gene_genome_ptr->getCreateContig(gene_ptr->id());

          if (not gene_contig_opt) {

            ExecEnv::log().error("GenomeGeneVariantAnalysis::genePopulationVariants; Unable to get/create gene Contig: {}, Genome: {}",
                                 gene_ptr->id(), gene_genome_ptr->genomeId());
            continue;

          }

          auto gene_contig_ptr = gene_contig_opt.value();

          auto coding_sequence_array = GeneFeature::getTranscriptionSequences(gene_ptr);

          // Only analyze the first sequence.
          if (not coding_sequence_array->empty()) {

            auto sequence_ptr = coding_sequence_array->getFirst();

            for (const auto &[feature_id, feature_ptr]: sequence_ptr->getFeatureMap()) {

              // Retrieve a gene contig containing all variants within the gene coding sequence.
              gene_contig_ptr->merge(contig_ptr->subset(feature_ptr->sequence().begin(), feature_ptr->sequence().end()));

            } // For all cds

          } // if not empty

        } // if same contig.

      } // genemap

    } // contig

  } // genome

}


void kgl::GenomeGeneVariantAnalysis::writeGeneResults(const std::string& variant_file_name) {

  std::ofstream variant_file(variant_file_name);

  if (not variant_file.good()) {

    ExecEnv::log().error("GenomeGeneVariantAnalysis::setGeneVector; Unable to open gene variant results file: {}", variant_file_name);
    return;

  }

  variant_file << "Gene_ID"
               << CSV_DELIMITER_
               << "Coding Length"
               << CSV_DELIMITER_
               << "Variant Count"
               << CSV_DELIMITER_
               << "Variant Rate"
               << CSV_DELIMITER_
               << "Unique Variants"
               << CSV_DELIMITER_
               << "Unique Variant Rate"
               << CSV_DELIMITER_
               << "Variant 0 Genomes"
               << CSV_DELIMITER_
               << "Variant 1 Genomes"
               << CSV_DELIMITER_
               << "Variant Top1 Genomes"
               << CSV_DELIMITER_
               << "Variant Top2 Genomes"
               << CSV_DELIMITER_
               << "Variant Top3 Genomes"
               << CSV_DELIMITER_
               << "Variant Top4 Genomes"
               << CSV_DELIMITER_
               << "Variant Top5 Genomes"
               << CSV_DELIMITER_
               << "Gene Detail"
               << '\n';

  for (auto const& gene_ptr : gene_vector_) {

    auto sequence_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);

    size_t coding_length{0};
    if (not sequence_array_ptr->empty()) {

      coding_length = sequence_array_ptr->getFirst()->codingNucleotides();

    }

    auto aggregated_gene_contig_ptr = std::make_unique<ContigDB>(gene_ptr->id());
    for (auto const& [gene_genome_id, gene_genome_ptr] : gene_population_ptr->getMap()) {

      auto gene_contig_opt = gene_genome_ptr->getContig(gene_ptr->id());
      if (gene_contig_opt) {

        aggregated_gene_contig_ptr->merge(gene_contig_opt.value());

      }

    }
    // Actual variants
    size_t variant_count = aggregated_gene_contig_ptr->variantCount();
    // Determine the number of unique variants.
    auto unique_variant_contig = aggregated_gene_contig_ptr->filterVariants(UniqueUnphasedFilter());
    size_t unique_variant_count = unique_variant_contig->variantCount();


    double variant_rate{0.0};
    double unique_variant_rate{0.0};
    if (coding_length > 0) {

      variant_rate = static_cast<double>(variant_count) / static_cast<double>(coding_length);
      unique_variant_rate = static_cast<double>(unique_variant_count) / static_cast<double>(coding_length);

    }

    // Generate genome per variant statistics.
    GeneGenomeAnalysis gene_genome_analysis(gene_ptr, unique_variant_contig);
    gene_genome_analysis.analyzeGenePopulation(gene_population_ptr);

    const size_t top_count_size{5};
    std::array<size_t, top_count_size> top_count_variants{0};
    size_t index = 0;
    for (auto const& [count, genome_count] : gene_genome_analysis.getCountSorted()) {

      if (index >= top_count_variants.size()) {

        break;

      }

      top_count_variants[index] = count;

      ++index;

    }

    variant_file << gene_ptr->id()
                 << CSV_DELIMITER_
                 << coding_length
                 << CSV_DELIMITER_
                 << variant_count
                 << CSV_DELIMITER_
                 << variant_rate
                 << CSV_DELIMITER_
                 << unique_variant_count
                 << CSV_DELIMITER_
                 << unique_variant_rate
                 << CSV_DELIMITER_
                 << gene_genome_analysis.zeroVariants().size()
                 << CSV_DELIMITER_
                 << gene_genome_analysis.getSingletons().size()
                 << CSV_DELIMITER_
                 << top_count_variants[0]
                 << CSV_DELIMITER_
                 << top_count_variants[1]
                 << CSV_DELIMITER_
                 << top_count_variants[2]
                 << CSV_DELIMITER_
                 << top_count_variants[3]
                 << CSV_DELIMITER_
                 << top_count_variants[4]
                 << CSV_DELIMITER_
                 << gene_ptr->featureText(CSV_DELIMITER_)
                 << '\n';
  } // Per gene.

}


kgl::GeneGenomeAnalysis::GeneGenomeAnalysis(std::shared_ptr<const GeneFeature> gene_ptr,
                                            const std::shared_ptr<const ContigDB>& gene_unique_variants)
                                            : gene_ptr_(std::move(gene_ptr)),
                                              gene_genome_analysis_ptr_(std::make_shared<GenomeCountMap>()) {


  ;
  gene_genome_analysis_ptr_->clear();
  gene_unique_variants->processAll(*this, &GeneGenomeAnalysis::addVariant);

}

bool kgl::GeneGenomeAnalysis::addVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  auto variant_hash = variant_ptr->HGVS();
  auto genome_count_ptr = std::make_shared<GenomeCount>();
  genome_count_ptr->variant_ptr_ = variant_ptr;
  auto [iter, result] = gene_genome_analysis_ptr_->insert({variant_hash, genome_count_ptr});
  if (not result) {

    ExecEnv::log().warn("GeneGenomeAnalysis::addVariant; attempt to add duplicate variant (should be unique)");

  }

  return true;

}

void kgl::GeneGenomeAnalysis::analyzeGenePopulation(const std::shared_ptr<const PopulationDB>& gene_population_ptr) {

  for (auto const& [gene_genome_id, gene_genome_ptr] : gene_population_ptr->getMap()) {

    // Gene contig should exist for each genome even if empty (no variants)
    auto gene_contig_opt = gene_genome_ptr->getContig(gene_ptr_->id());
    if (not gene_contig_opt) {

//      ExecEnv::log().error("GeneGenomeAnalysis::analyzeGenePopulation; Genome: {}, No Gene Contig: {}", gene_genome_id, gene_ptr_->id());
      zero_variants_.push_back(gene_genome_id);
      continue;

    }

    auto gene_contig = gene_contig_opt.value();
    if (gene_contig->variantCount() == 0) {

      zero_variants_.push_back(gene_genome_id);

    } else {

      // Create an object to process all variants.
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class ProcessVariants {

      public:

        ProcessVariants(GenomeId_t genome_id, std::shared_ptr<GenomeCountMap> variant_genome_count_ptr)
        : genome_id_(std::move(genome_id)), variant_genome_count_ptr_(std::move(variant_genome_count_ptr)) {}
        ~ProcessVariants() = default;

        bool processAllVariants(const std::shared_ptr<const Variant>& variant_ptr) {

          auto result = variant_genome_count_ptr_->find(variant_ptr->HGVS());
          if (result == variant_genome_count_ptr_->end()) {

            ExecEnv::log().error("GeneGenomeAnalysis::analyzeGenePopulation;; cannot find variant HGVS entry: {}", variant_ptr->HGVS());
            return true;

          }

          auto& [HGVS_str, genome_count_ptr] = *result;
          genome_count_ptr->genome_vector_.push_back(genome_id_);

          return true;

        }

      private:

        GenomeId_t genome_id_;
        std::shared_ptr<GenomeCountMap> variant_genome_count_ptr_;


      };
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

      ProcessVariants process_variants(gene_genome_id, gene_genome_analysis_ptr_);
      gene_contig->processAll(process_variants, &ProcessVariants::processAllVariants);

    }


  } // For all genomes.

}

kgl::GenomeCountSorted kgl::GeneGenomeAnalysis::getCountSorted() const {

  GenomeCountSorted count_sorted_map;
  for (auto& [variant_str, count_ptr] : *gene_genome_analysis_ptr_) {

      count_sorted_map.insert({count_ptr->genome_vector_.size(), count_ptr});

  }

  return count_sorted_map;

}

std::vector<std::shared_ptr<const kgl::GenomeCount>> kgl::GeneGenomeAnalysis::getSingletons() const {

  std::vector<std::shared_ptr<const kgl::GenomeCount>> singletons;
  for (auto& [variant_str, count_ptr] : *gene_genome_analysis_ptr_) {

    if (count_ptr->genome_vector_.size() == 1) {

      singletons.push_back(count_ptr);

    }

  }

  return singletons;

}
