//
// Created by kellerberrin on 29/09/18.
//


#include "kgl_finestructure_analysis.h"

#include <fstream>
#include <kgl_filter.h>


namespace kgl = kellerberrin::genome;


bool kgl::FineStructureAnalysis::generateFiles(const std::string& Filename,
                                               std::shared_ptr<const PhasedPopulation> population_ptr,
                                               double bases_per_centimorgan) {

  if (not generateIDFile(Filename, population_ptr)) {

    ExecEnv::log().error("FineStructureAnalysis::generateFiles; problem generating ID file: '{}.ids'", Filename);
    return false;

  }

  if (not generatePhaseFile(Filename, population_ptr, bases_per_centimorgan)) {

    ExecEnv::log().error("FineStructureAnalysis::generateFiles; problem generating Phase files: {}", Filename);
    return false;

  }


  return true;

}


bool kgl::FineStructureAnalysis::generateIDFile(const std::string& Filename,
                                                std::shared_ptr<const PhasedPopulation> population_ptr) {

  std::string full_filename = Filename + ".ids";

  std::ofstream id_file(full_filename);

  if (not id_file.good()) {

    ExecEnv::log().error("FineStructureAnalysis::generateFiles; another to open output file: {}", full_filename);

  }

  for (auto genome : population_ptr->getMap()) {

    std::string id_line = genome.first + DELIMITER_ + population_ptr->populationId() + DELIMITER_ + "1" + "\n";
    id_file << id_line;

  }

  return true;

}



bool kgl::FineStructureAnalysis::generatePhaseFile(const std::string& Filename,
                                                   std::shared_ptr<const PhasedPopulation> population_ptr,
                                                   double bases_per_centimorgan) {

// Restrict this to SNPs.
  std::shared_ptr<const PhasedPopulation> snp_population_ptr = population_ptr->filterVariants(SNPFilter());

  using UniqueContigs = std::set<ContigId_t>;
  UniqueContigs unique_contigs;
// Generate a list of all non empty contigs.
  for (auto genome : snp_population_ptr->getMap()) {

    for (auto contig : genome.second->getMap()) {

      if (contig.second->variantCount() > 0) {

        unique_contigs.insert(contig.first);

      }

    }

  }

  ExecEnv::log().info("FineStructureAnalysis::generatePhaseFile; population : {} has: {} non-empty contigs",
                      population_ptr->populationId(), unique_contigs.size());

// For all contigs and genomes calculate how many variants (SNP & INDEL) in each contig

  size_t contig_index = 0;
  for (auto contig : unique_contigs) {

    ++contig_index;

    // Each element of the vector is an homologous chromosome/contig.
    using VariantMapVector = std::map<ContigOffset_t, std::vector<std::shared_ptr<const Variant>>>;
    VariantMapVector variant_map_vector;
    size_t heterozygous_count = 0;
    size_t ploidy = 0;

    for (auto genome : snp_population_ptr->getMap()) {

      std::shared_ptr<ContigVariant> contig_variant_ptr;
      genome.second->getContigVariant(contig, contig_variant_ptr);

      if (not contig_variant_ptr) continue; // continue on to next genome.

      ploidy = contig_variant_ptr->getVector().size();

      size_t homologous_index = 0;
      for (auto homologous : contig_variant_ptr->getVector()) {

        for (auto variant : homologous->getMap()) {

          auto find_variant = variant_map_vector.find(variant.second->offset());

          if (find_variant == variant_map_vector.end()) {

            // create a vector of nulls for all homologous contigs.
            std::vector<std::shared_ptr<const Variant>> variant_vector(ploidy, nullptr);
            // add the variant to the vector at the correct homologous offset.
            variant_vector[homologous_index] = variant.second;
            // add to the map.
            std::pair<ContigOffset_t, std::vector<std::shared_ptr<const Variant>>> insert_pair(variant.second->offset(), variant_vector);
            auto result = variant_map_vector.insert(insert_pair);
            if (not result.second) {

              ExecEnv::log().error("FineStructureAnalysis::generatePhaseFile; Unable to add variant: {}",
                                   variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));

            }

          } else { // variant offset already added

            if (find_variant->second[homologous_index]) {

              if (not variant.second->equivalent(*(find_variant->second[homologous_index]))) {

                ExecEnv::log().vinfo("FineStructureAnalysis::generatePhaseFile; heterozygous variants: \n{}\n{}",
                                    variant.second->output(' ', VariantOutputIndex::START_0_BASED, true),
                                    find_variant->second[homologous_index]->output(' ', VariantOutputIndex::START_0_BASED, true));
                ++heterozygous_count;

              }

            }

          }

        }

        ++homologous_index;

      }

    }  // create variant map.

    ExecEnv::log().info("FineStructureAnalysis::generatePhaseFile; Contig: {} has: {} variants, with: {} heterzygous",
                        contig, variant_map_vector.size(), heterozygous_count);

    // Open the file and add header info.
    std::stringstream ss;
    ss << Filename << "_" << contig_index << ".phase";
    std::ofstream phase_file(ss.str());

    if (not phase_file.good()) {

      ExecEnv::log().error("FineStructureAnalysis::generatePhaseFile; Could not open phase file: {} for output", ss.str());
      return false;

    }

    // Header info - number of haplotypes.
    phase_file << population_ptr->getMap().size() * ploidy << '\n';
    // Header info - number of variants.
    phase_file << variant_map_vector.size() << '\n';

    // All the variant offsets
    phase_file << "P";

    for (auto variant : variant_map_vector) {

      phase_file << DELIMITER_ << variant.first;

    }

    phase_file << '\n';

    // The SNP data for the individual genomes.
    for (auto genome : snp_population_ptr->getMap()) {

      std::shared_ptr<ContigVariant> contig_variant_ptr;
      genome.second->getContigVariant(contig, contig_variant_ptr);
      if (not contig_variant_ptr) {

        ExecEnv::log().error("FineStructureAnalysis::generatePhaseFile; Genome: {} does not contain contig: {}",
                             genome.first, contig);

        return false;

      }

      for (auto variant : variant_map_vector) {

        for (size_t idx = 0; idx < variant.second.size(); ++idx) {

          if (variant.second.size() != contig_variant_ptr->ploidy()) {

            ExecEnv::log().error("FineStructureAnalysis::generatePhaseFile; Genome: {} ploidy mismatch for contig: {}",
                                 genome.first, contig);

            return false;

          }

          auto result = contig_variant_ptr->getHomologous(idx)->getMap().find(variant.first);

          if (result == contig_variant_ptr->getHomologous(idx)->getMap().end()) {

            phase_file << '0';

          } else {

            auto vcf_ptr = std::dynamic_pointer_cast<const VCFVariant>(result->second);

            if (not vcf_ptr) {

              ExecEnv::log().error("FineStructureAnalysis::generatePhaseFile; Unexpected variant: {}",
                                   result->second->output(' ', VariantOutputIndex::START_0_BASED, true));

              return false;

            }

            phase_file << vcf_ptr->alternate().getSequenceAsString();

          }

        }

      }

      phase_file << '\n';


      // Open the file and add header info.
      std::stringstream ss;
      ss << Filename << "_" << contig_index << ".recombfile";
      std::ofstream recomb_file(ss.str());

      if (not phase_file.good()) {

        ExecEnv::log().error("FineStructureAnalysis::generatePhaseFile; Could not open recomb file: {} for output", ss.str());
        return false;

      }

      // Generate a recomb file (assumes constant recomb rate).
      recomb_file << "start.pos" << DELIMITER_ << "recom.rate.perbp" << '\n';

      for (auto iter = variant_map_vector.begin(); iter != variant_map_vector.end(); ++iter) {

        auto next_iter = iter;
        ++next_iter;

        if (next_iter == variant_map_vector.end()) {

          recomb_file << iter->first << DELIMITER_ << '0' << '\n';

        } else {

          double recomb_prob =  static_cast<double>(next_iter->first - iter->first) / (bases_per_centimorgan * 100.0);
          recomb_file << iter->first << DELIMITER_ << recomb_prob << '\n';

        }

      }

    }

  } // contig.

// Generates a recomb file for every non-empty contig.

  return true;

}

