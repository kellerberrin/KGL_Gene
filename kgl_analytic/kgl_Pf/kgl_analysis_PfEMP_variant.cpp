//
// Created by kellerberrin on 9/04/23.
//


#include "kgl_analysis_PfEMP_variant.h"
#include "kgl_variant_filter_db_contig.h"
#include "kgl_variant_filter_db_offset.h"
#include "kel_utility.h"

#include <fstream>


namespace kgl = kellerberrin::genome;



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

        if (gene_ptr->contig_ref_ptr()->contigId() == contig_id) {

          // Get/Create the gene contig_ref_ptr from the gene population.
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

              // Retrieve a gene contig_ref_ptr containing all variants within the gene coding sequence.
              std::shared_ptr<const ContigDB> feature_contig_ptr = contig_ptr->viewFilter(ContigModifyFilter(feature_ptr->sequence().begin(),
                                                                                                             feature_ptr->sequence().end()));
              gene_contig_ptr->merge(feature_contig_ptr);

            } // For all cds

          } // if not empty

        } // if same contig_ref_ptr.

      } // genemap

    } // contig_ref_ptr

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
               << "Unique filter"
               << CSV_DELIMITER_
               << "Unique Variant Rate"
               << CSV_DELIMITER_
               << "Single Genome filter"
               << CSV_DELIMITER_
               << "Single Variant Genomes"
               << CSV_DELIMITER_
               << "No Variant Genomes"
               << CSV_DELIMITER_
               << "Genotypes"
               << CSV_DELIMITER_
               << "Gini"
               << CSV_DELIMITER_
               << "HRel"
               << CSV_DELIMITER_
               << "2-IQV"
               << CSV_DELIMITER_
               << "Top1 Genotype"
               << CSV_DELIMITER_
               << "Top2 Genotype"
               << CSV_DELIMITER_
               << "Top3 Genotype"
               << CSV_DELIMITER_
               << "Top4 Genotype"
               << CSV_DELIMITER_
               << "Top5 Genotype"
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
    std::shared_ptr<const ContigDB> unique_variant_contig = aggregated_gene_contig_ptr->viewFilter(UniqueUnphasedFilter());
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
//    auto top_count_variants = gene_genome_analysis.getTopCounts<5>();

    // Generate unique genotypes.
    GenotypeAnalysis genotype_analysis(gene_ptr);
    genotype_analysis.analyzeGenePopulation(gene_population_ptr);
    auto top_count_genotypes = genotype_analysis.getTopCounts<5>();

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
                 << gene_genome_analysis.getSingletonVariants().size()
                 << CSV_DELIMITER_
                 << gene_genome_analysis.getSingletonGenomes().size()
                 << CSV_DELIMITER_
                 << gene_genome_analysis.zeroVariants().size()
                 << CSV_DELIMITER_
                 << genotype_analysis.getGenoTypeMap().size()
                 << CSV_DELIMITER_
                 << genotype_analysis.Gini()
                 << CSV_DELIMITER_
                 << genotype_analysis.HRel()
                 << CSV_DELIMITER_
                 << genotype_analysis.IQVAdj()
                 << CSV_DELIMITER_
                 << top_count_genotypes[0]
                 << CSV_DELIMITER_
                 << top_count_genotypes[1]
                 << CSV_DELIMITER_
                 << top_count_genotypes[2]
                 << CSV_DELIMITER_
                 << top_count_genotypes[3]
                 << CSV_DELIMITER_
                 << top_count_genotypes[4]
                 << CSV_DELIMITER_
                 << gene_ptr->featureText(CSV_DELIMITER_)
                 << '\n';
  } // Per gene.

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::GeneGenomeAnalysis::GeneGenomeAnalysis(std::shared_ptr<const GeneFeature> gene_ptr,
                                            const std::shared_ptr<const ContigDB>& gene_unique_variants)
                                            : gene_ptr_(std::move(gene_ptr)),
                                              gene_genome_analysis_ptr_(std::make_shared<GenomeCountMap>()) {

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

    auto gene_contig_opt = gene_genome_ptr->getContig(gene_ptr_->id());
    if (not gene_contig_opt) {

      zero_variants_.push_back(gene_genome_id);
      continue;

    }

    auto gene_contig = gene_contig_opt.value();
    if (gene_contig->variantCount() == 0) {

      zero_variants_.push_back(gene_genome_id);

    } else {

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // A local object to process all variants.
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class ProcessVariants {

      public:

        ProcessVariants(GenomeId_t genome_id, std::shared_ptr<GenomeCountMap> variant_genome_count_ptr)
        : genome_id_(std::move(genome_id)), variant_genome_count_ptr_(std::move(variant_genome_count_ptr)) {}
        ~ProcessVariants() = default;

        bool processAllVariants(const std::shared_ptr<const Variant>& variant_ptr) {

          auto result = variant_genome_count_ptr_->find(variant_ptr->HGVS());
          if (result == variant_genome_count_ptr_->end()) {

            ExecEnv::log().error("GeneGenomeAnalysis::analyzeGenePopulation; cannot find variant HGVS entry: {}", variant_ptr->HGVS());
            return true;

          }

          auto& [HGVS_str, genome_count_ptr] = *result;
          auto [iter, insert_result] = genome_count_ptr->genome_set_.insert(genome_id_);
          if (not insert_result) {

            auto [homo_iter, homo_result] = genome_count_ptr->homozygous_set_.insert(genome_id_);
            if (not homo_result) {

              ExecEnv::log().error("GeneGenomeAnalysis::analyzeGenePopulation; cannot insert genome: {} into homozygous set", genome_id_);

            }

          }

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

      count_sorted_map.insert({count_ptr->genome_set_.size(), count_ptr});

  }

  return count_sorted_map;

}

std::vector<std::shared_ptr<const kgl::GenomeCount>> kgl::GeneGenomeAnalysis::getSingletonVariants() const {

  std::vector<std::shared_ptr<const kgl::GenomeCount>> singletons;
  for (auto& [variant_str, count_ptr] : *gene_genome_analysis_ptr_) {

    if (count_ptr->genome_set_.size() == 1) {

      singletons.push_back(count_ptr);

    }

  }

  return singletons;

}

std::set<kgl::GenomeId_t> kgl::GeneGenomeAnalysis::getSingletonGenomes() const {

  std::set<GenomeId_t> singleton_genomes;
  for (auto& [variant_str, count_ptr] : *gene_genome_analysis_ptr_) {

    if (count_ptr->genome_set_.size() == 1) {

      auto genome = *(count_ptr->genome_set_.begin());
      singleton_genomes.insert(genome);

    }

  }

  return singleton_genomes;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::GenotypeAnalysis::analyzeGenePopulation(const std::shared_ptr<const PopulationDB>& gene_population_ptr) {

  for (auto const& [genome_id, gene_genome_ptr] : gene_population_ptr->getMap()) {

  // Find the gene contig_ref_ptr (or not).
    auto gene_contig_opt = gene_genome_ptr->getContig(gene_ptr_->id());
    if (not gene_contig_opt) {

      zero_variants_.push_back(genome_id);
      continue;

    }

    auto gene_contig_ptr = gene_contig_opt.value();
    // Don't record no variant genotypes.
    if (gene_contig_ptr->variantCount() == 0) {

      zero_variants_.push_back(genome_id);
      continue;

    }

    size_t genotype_hash = genotypeHash(gene_contig_ptr);

    auto find_iter = genotype_map_.find(genotype_hash);
    if (find_iter == genotype_map_.end()) {

      auto genotype_ptr = std::make_shared<GenotypeCount>();
      genotype_ptr->genotype_hash_ = genotype_hash;
      genotype_ptr->genotype_ = gene_contig_ptr;
      genotype_ptr->genomes_.push_back(genome_id);

      auto [insert_iter, insert_result] = genotype_map_.insert({genotype_hash, genotype_ptr});
      if (not insert_result) {

        ExecEnv::log().error("GenotypeAnalysis::analyzeGenePopulation; could not insert (duplicate) genotype for gene: {}, genome: {}",
                             gene_ptr_->id(), genome_id);

      }

    } else {

      auto& [geno_hash, genotype_ptr] = *find_iter;
      genotype_ptr->genomes_.push_back(genome_id);

    }

  } // For all genomes.

}


size_t kgl::GenotypeAnalysis::genotypeHash(const std::shared_ptr<const ContigDB>& genotype) {

  struct GenotypeVariants {

    // The Pf7 data is presented as diploid to assess COI. Only generate the hash from unique variants (ignore homozygous variants).
    std::map<std::string, std::shared_ptr<const Variant>> unique_genotype_variants_;

    bool processAllVariants(const std::shared_ptr<const Variant>& variant_ptr) {

      auto hash_str = variant_ptr->HGVS();
      unique_genotype_variants_.insert({hash_str, variant_ptr});
      return true;

    }

  };

  GenotypeVariants genotype_variants;
  genotype->processAll(genotype_variants, &GenotypeVariants::processAllVariants);

  size_t geno_hash{0};

  // A bit dodgy, there may be hash collisions. This should be tested.
  for (auto const& [hash_str, variant_ptr] : genotype_variants.unique_genotype_variants_) {

    geno_hash += Utility::hash(hash_str);

  }

  return geno_hash;

}

kgl::GenotypeCountSorted kgl::GenotypeAnalysis::getCountSorted() const {

  GenotypeCountSorted  count_sorted_map;
  for (auto const& [geno_hash, genotype_ptr] : genotype_map_) {

    count_sorted_map.insert({genotype_ptr->genomes_.size(), genotype_ptr});

  }

  return count_sorted_map;

}

// The gini coefficient
double kgl::GenotypeAnalysis::Gini() const {

  const size_t genotype_count = genotype_map_.size() + 1; // categories.
  if (genotype_count <= 1) {

    return 1.0;

  }

  size_t total_genomes = zero_variants_.size();
  for (auto const& [geno_hash, genotype_ptr] : genotype_map_) {

    total_genomes += genotype_ptr->genomes_.size();

  }
  // Generate the frequencies.
  std::vector<double> genotype_frequencies;
  genotype_frequencies.push_back(zero_variants_.size());
  for (auto const& [geno_hash, genotype_ptr] : genotype_map_) {

    double frequency = static_cast<double>(genotype_ptr->genomes_.size());
    genotype_frequencies.push_back(frequency);

  }

  double sum{0.0};
  for (size_t i = 0; i < (genotype_count - 1); ++i) {

    for (size_t j = i+ 1; j < genotype_count; ++j) {

        sum += std::fabs(genotype_frequencies[i] - genotype_frequencies[j]);

    }

  }

  double scaling = 1.0 / (static_cast<double>(total_genomes) * static_cast<double>(genotype_count - 1));

  return 1.0 - (scaling * sum);

}

// A modified Shannon entropy measure for different category (genotype) counts.
double kgl::GenotypeAnalysis::HRel() const {

  const size_t genotype_count = genotype_map_.size() + 1; // categories.
  if (genotype_count <= 1) {

    return 0.0;

  }

  size_t total_genomes = zero_variants_.size();
  for (auto const& [geno_hash, genotype_ptr] : genotype_map_) {

    total_genomes += genotype_ptr->genomes_.size();

  }

  if (total_genomes == 0) {

    return 0.0;

  }
  // Generate the frequencies.
  double sum{0.0};
  if (not zero_variants_.empty()) {

    double proportion = static_cast<double>(zero_variants_.size()) / static_cast<double>(total_genomes);
    sum = std::log2(proportion) * proportion;

  }

  for (auto const& [geno_hash, genotype_ptr] : genotype_map_) {

    if (genotype_ptr->genomes_.empty()) {

      ExecEnv::log().warn("GenotypeAnalysis::HRel; Unexpected empty genotype");
      continue;

    }

    double proportion = static_cast<double>(genotype_ptr->genomes_.size()) / static_cast<double>(total_genomes);
    sum += std::log2(proportion) * proportion;

  }

  return (-1.0 * sum) / std::log2(genotype_count);

}

// Index of qualitative variation.
double kgl::GenotypeAnalysis::IQV() const {

  const size_t genotype_count = genotype_map_.size() + 1; // categories.
  if (genotype_count <= 1) {

    return 1.0;

  }

  size_t total_genomes = zero_variants_.size();
  for (auto const& [geno_hash, genotype_ptr] : genotype_map_) {

    total_genomes += genotype_ptr->genomes_.size();

  }
  
  if (total_genomes == 0) {

    return 1.0;

  }

  // Generate the frequencies.
  double proportion = static_cast<double>(zero_variants_.size()) / (static_cast<double>(total_genomes) * 100.0);
  double sum = (proportion * proportion);
  for (auto const& [geno_hash, genotype_ptr] : genotype_map_) {

    proportion = static_cast<double>(genotype_ptr->genomes_.size()) / (static_cast<double>(total_genomes) * 100.0);
    sum += (proportion * proportion);

  }

  double scale = static_cast<double>(genotype_count) / static_cast<double>(genotype_count - 1);

  return scale * (1.0 - sum);

}
