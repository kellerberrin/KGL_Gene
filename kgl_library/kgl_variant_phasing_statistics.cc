//
// Created by kellerberrin on 3/04/18.
//

#include "kgl_variant_phasing_statistics.h"


namespace kgl = kellerberrin::genome;



bool kgl::ContigPhasingStatistics::insertPhasingStatistic(ContigOffset_t  offset,
                                                          const PhasedSNPVector& phased_snp_vector,
                                                          ContigStatsMap& stats_map) {

  std::shared_ptr<OffsetPhasingStatistic> phasing_statisic_ptr(std::make_shared<OffsetPhasingStatistic>(offset, phased_snp_vector));
  auto result = stats_map.insert(ContigStatsMap::value_type(offset, phasing_statisic_ptr));
  return result.second;

}


bool kgl::GenomePhasingStatistics::phasedSNPs(const VCFGenome& vcf_genome) {

  bool return_result = true;

  for (auto contig : vcf_genome.getMap()) {

    std::shared_ptr<ContigPhasingStatistics> contig_stats_map_ptr(std::make_shared<ContigPhasingStatistics>());
    auto current_variant = contig.second->begin();
    while (current_variant != contig.second->end()) {

      // Take a peek at the next variant to see if it is an SNP with the same contig offset.
      auto next_variant = current_variant;
      ++next_variant;
      if (next_variant == contig.second->end()) {

        break;

      }

      if (current_variant->second->isSNP()) {

        PhasedSNPVector different_snp_vector;
        PhasedSNPVector identical_snp_vector;
        different_snp_vector.push_back(current_variant->second);
        identical_snp_vector.push_back(current_variant->second);

        while (current_variant->second->offset() == next_variant->second->offset()) {

          if (next_variant->second->isSNP()) {

            if (current_variant->second->equivalent(*(next_variant->second))) {

              identical_snp_vector.push_back(next_variant->second);

            } else {

              different_snp_vector.push_back(next_variant->second);

            }

          }

          ++next_variant;
          if (next_variant == contig.second->end()) {

            break;

          }

        }

        if (different_snp_vector.size() >= 2) {

          if (not contig_stats_map_ptr->differentSNP(current_variant->second->offset(), different_snp_vector)) {

            ExecEnv::log().error("Unable to add phased SNP vector contig: {}, offset: {}", contig.first, current_variant->second->offset());
            return_result = false;

          }

        }
        if (identical_snp_vector.size() >= 2) {

          if (not contig_stats_map_ptr->identicalSNP(current_variant->second->offset(), identical_snp_vector)) {

            ExecEnv::log().error("Unable to add phased SNP vector contig: {}, offset: {}", contig.first, current_variant->second->offset());
            return_result = false;

          }

        }
        if (identical_snp_vector.size() == 1 and identical_snp_vector.size()) {

          if (not contig_stats_map_ptr->singleSNP(current_variant->second->offset(), identical_snp_vector)) {

            ExecEnv::log().error("Unable to add phased SNP vector contig: {}, offset: {}", contig.first, current_variant->second->offset());
            return_result = false;

          }

        }

      } // current == SNP


      current_variant = next_variant;

    }

    auto result = contig_map_.insert(StatContigMap::value_type(contig.first, contig_stats_map_ptr));
    if (not result.second) {

      ExecEnv::log().error("Unable to add phased SNP vector map to contig: {}", contig.first);
      return_result = false;

    }

  }

  return return_result;

}



size_t kgl::GenomePhasingStatistics::differentSNPCount() const {

  size_t differentSNPs = 0;
  for (auto contig : getMap()) {

    differentSNPs += contig.second->differentSNP().size();

  }

  return differentSNPs;

}



size_t kgl::GenomePhasingStatistics::identicalSNPCount() const {

  size_t differentSNPs = 0;
  for (auto contig : getMap()) {

    differentSNPs += contig.second->identicalSNP().size();

  }

  return differentSNPs;

}


size_t kgl::GenomePhasingStatistics::singleSNPCount() const {

  size_t differentSNPs = 0;
  for (auto contig : getMap()) {

    differentSNPs += contig.second->singleSNP().size();

  }

  return differentSNPs;

}



bool kgl::PopulationPhasingStatistics::phasedSNPs(const VCFPopulation& vcf_population) {

  bool return_result = true;

  for (auto genome : vcf_population.getMap()) {

    std::shared_ptr<GenomePhasingStatistics> genome_stats_ptr(std::make_shared<GenomePhasingStatistics>());

    if (not genome_stats_ptr->phasedSNPs(*genome.second)) {

      ExecEnv::log().error("Problem calculating phased SNP for genome: {}", genome.first);

    }

    auto result = genome_map_.insert(GenomeStatMap::value_type(genome.first, genome_stats_ptr));
    if (not result.second) {

      ExecEnv::log().error("Problem insert phased SNP statistics for genome: {}", genome.first);
      return_result = false;

    }

  }

  return return_result;

}


void kgl::PopulationPhasingStatistics::outputPopulation() const {

  for (auto genome : getMap()) {

    size_t different = genome.second->differentSNPCount();
    size_t identical = genome.second->identicalSNPCount();
    size_t single = genome.second->singleSNPCount();
    size_t allSNP = different + identical + single;

    double percent_different = (static_cast<double>(different) * 100.0) / static_cast<double>(allSNP);
    double percent_identical = (static_cast<double>(identical) * 100.0) / static_cast<double>(allSNP);
    double percent_single = (static_cast<double>(single) * 100.0) / static_cast<double>(allSNP);

    ExecEnv::log().info("Diploid SNP types; Genome: {} Different: {}({}%), Identical: {}({}%), Single: {}({}%)",
                        genome.first,
                        different, percent_different,
                        identical, percent_identical,
                        single, percent_single);

  }

}

