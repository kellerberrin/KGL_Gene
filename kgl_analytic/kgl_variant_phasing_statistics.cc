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


bool kgl::DiploidPhasingStatistics::phasedSNPs(const VCFGenome& vcf_genome) {

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

          if (not contig_stats_map_ptr->heterozygousSNP(current_variant->second->offset(), different_snp_vector)) {

            ExecEnv::log().error("Unable to add phased SNP vector contig: {}, offset: {}", contig.first, current_variant->second->offset());
            return_result = false;

          }

        }
        if (identical_snp_vector.size() >= 2) {

          if (not contig_stats_map_ptr->homozygousSNP(current_variant->second->offset(), identical_snp_vector)) {

            ExecEnv::log().error("Unable to add phased SNP vector contig: {}, offset: {}", contig.first, current_variant->second->offset());
            return_result = false;

          }

        }
        if (identical_snp_vector.size() == 1 and identical_snp_vector.size()) {

          if (not contig_stats_map_ptr->singleHeterozygousSNP(current_variant->second->offset(), identical_snp_vector)) {

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



size_t kgl::DiploidPhasingStatistics::heterozygousSNPCount() const {

  size_t differentSNPs = 0;
  for (auto contig : getMap()) {

    differentSNPs += contig.second->heterozygousSNP().size();

  }

  return differentSNPs;

}



size_t kgl::DiploidPhasingStatistics::homozygousSNPCount() const {

  size_t differentSNPs = 0;
  for (auto contig : getMap()) {

    differentSNPs += contig.second->homozygousSNP().size();

  }

  return differentSNPs;

}


size_t kgl::DiploidPhasingStatistics::singleHeterozygousSNPCount() const {

  size_t differentSNPs = 0;
  for (auto contig : getMap()) {

    differentSNPs += contig.second->singleHeterozygousSNP().size();

  }

  return differentSNPs;

}



bool kgl::PopulationPhasingStatistics::phasedSNPs(const VCFPopulation& vcf_population) {

  bool return_result = true;

  for (auto genome : vcf_population.getMap()) {

    std::shared_ptr<DiploidPhasingStatistics> genome_stats_ptr(std::make_shared<DiploidPhasingStatistics>());

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


bool kgl::PopulationPhasingStatistics::getPhasing(const GenomeId_t& genome_id,
                                                  const ContigId_t& contig_id,
                                                  ContigOffset_t offset,
                                                  std::shared_ptr<const OffsetPhasingStatistic>& snp_phasing) const {

  snp_phasing = nullptr;
  auto genome_result = getMap().find(genome_id);

  if (genome_result == getMap().end()) {

    return false;

  }

  auto contig_result = genome_result->second->getMap().find(contig_id);

  if (contig_result == genome_result->second->getMap().end()) {

    return false;

  }

  bool found = false;

  auto offset_result = contig_result->second->heterozygousSNP().find(offset);

  if (offset_result != contig_result->second->heterozygousSNP().end()) {

    found = true;
    snp_phasing = offset_result->second;

  }

  offset_result = contig_result->second->homozygousSNP().find(offset);

  if (offset_result != contig_result->second->homozygousSNP().end()) {

    if (found) {

      ExecEnv::log().error("SNP at Genome: {}, Contig: {}, Offset: {}, recorded as both Heterozygous AND Homozygous");

    }

    found = true;
    snp_phasing = offset_result->second;

  }

  offset_result = contig_result->second->singleHeterozygousSNP().find(offset);

  if (offset_result != contig_result->second->singleHeterozygousSNP().end()) {

    if (found) {

      ExecEnv::log().error("SNP at Genome: {}, Contig: {}, Offset: {}, recorded as both Homozygous AND Singular");

    }

    found = true;
    snp_phasing = offset_result->second;

  }

  return found;

}


void kgl::PopulationPhasingStatistics::outputPopulation() const {

  for (auto genome : getMap()) {

    size_t heterozygous = genome.second->heterozygousSNPCount();
    size_t homozygous = genome.second->homozygousSNPCount();
    size_t single = genome.second->singleHeterozygousSNPCount();
    size_t allSNP = heterozygous + homozygous + single;

    double percent_heterozygous = (static_cast<double>(heterozygous) * 100.0) / static_cast<double>(allSNP);
    double percent_homozygous = (static_cast<double>(homozygous) * 100.0) / static_cast<double>(allSNP);
    double percent_single = (static_cast<double>(single) * 100.0) / static_cast<double>(allSNP);

    ExecEnv::log().info("Diploid SNP types; Genome: {} Heterozygous: {}({}%), Homozygous: {}({}%), Singular: {}({}%)",
                        genome.first,
                        heterozygous, percent_heterozygous,
                        homozygous, percent_homozygous,
                        single, percent_single);

  }

}


bool kgl::PopulationPhasingStatistics::genomePhasing(const GenomeId_t& genome_id,
                                                     size_t& heterozygous,
                                                     size_t& homozygous,
                                                     size_t& singleheterozygous) const {

  auto genome_result = getMap().find(genome_id);

  if (genome_result == getMap().end()) {

    heterozygous = 0;
    homozygous = 0;
    singleheterozygous =0;
    return false;

  }

  heterozygous = genome_result->second->heterozygousSNPCount();
  homozygous = genome_result->second->homozygousSNPCount();
  singleheterozygous = genome_result->second->singleHeterozygousSNPCount();

  return true;

}
