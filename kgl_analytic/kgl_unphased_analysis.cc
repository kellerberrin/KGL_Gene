//
// Created by kellerberrin on 11/08/18.
//

#include "kgl_unphased_analysis.h"

#include <fstream>


namespace kgl = kellerberrin::genome;



bool kgl::AggregateVariantDistribution::variantDistribution(std::shared_ptr<const UnphasedPopulation> unphased_population) {


  for (auto genome : unphased_population->getMap()) {

    for (auto contig : genome.second->getMap()) {

      for (auto offset : contig.second->getMap()) {

        for (auto variant : offset.second) {

          analysis_genome_.addVariant(variant);

        } // offset vector

      } // offset

    } // contig

  } // genome

  ExecEnv::log().info("AggregatVariantDistribution; variants processed: {}", analysis_genome_.variantCount());

  return true;

}


bool kgl::AggregateVariantDistribution::writeDistribution(std::shared_ptr<const GenomeDatabase> genome_db,
                                                         const std::string& filename,
                                                         const char delimiter) const {

  std::ofstream outfile;
  outfile.open(filename, std::ofstream::out | std::ofstream::app);

  if (not outfile.good()) {

    ExecEnv::log().error("AggregatVariantDistribution; could open results file: {}", filename);
    return false;

  }

  writeHeader(outfile, delimiter);
  writeData(genome_db, outfile, delimiter);

  return outfile.good();

}



bool kgl::AggregateVariantDistribution::writeHeader(std::ostream& output, const char delimiter) const {

  output << "contig" << delimiter;
  output << "interval_start" << delimiter;
  output << "interval_size" << delimiter;
  output << "snp_count" << delimiter;
  output << "delete_count" << delimiter;
  output << "delete_bases" << delimiter;
  output << "insert_count" << delimiter;
  output << "insert_bases" << '\n';

  return output.good();

}


bool kgl::AggregateVariantDistribution::writeData(std::shared_ptr<const GenomeDatabase> genome_db,
                                                 std::ostream& output,
                                                 const char delimiter) const {
  for (auto contig : genome_db->getMap()) {

    ContigSize_t contig_size = contig.second->contigSize();
    ContigId_t contig_id = contig.second->contigId();

    ContigOffset_t contig_from = 0;
    ContigOffset_t contig_to = interval_;

    auto contig_iter  = analysis_genome_.getMap().find(contig_id);

    if (contig_iter != analysis_genome_.getMap().end()) {

      while (contig_from < contig_size) {

        output << contig_id << delimiter;
        output << contig_from << delimiter;
        output << (contig_to - contig_from) << delimiter;


        size_t snp_count = 0;
        size_t delete_count = 0;
        size_t delete_base_count = 0;
        size_t insert_count = 0;
        size_t insert_base_count = 0;

        auto begin_offset_iter = contig_iter->second->getMap().lower_bound(contig_from);
        auto end_offset_iter = contig_iter->second->getMap().upper_bound(contig_to);

        for (auto offset_iter = begin_offset_iter; offset_iter != end_offset_iter; ++offset_iter) {

          for (auto variant : offset_iter->second) {

            if (variant->isSNP()) {

              ++snp_count;

            } else if (variant->isDelete()) {

              ++delete_count;
              delete_base_count += variant->size();

            } else if (variant->isInsert()) {

              ++insert_count;
              insert_base_count += variant->size();

            } else {

              ExecEnv::log().critical("AggregatVariantDistribution; encounted unknown variant type");

            }

          }

        }

        output << snp_count << delimiter;
        output << delete_count << delimiter;
        output << delete_base_count << delimiter;
        output << insert_count << delimiter;
        output << insert_base_count << '\n';

        contig_from = contig_to;
        if (contig_to + interval_ < contig_size) {

          contig_to += interval_;

        } else {

          contig_to = contig_size;

        }

      }

    }

  }

  return output.good();

}



bool kgl::HeterozygousStatistics::heterozygousStatistics(std::shared_ptr<const UnphasedPopulation> unphased_population) {
// Write data.

  for (auto genome : unphased_population->getMap()) {

    size_t heterozygous = 0;
    size_t homozygous = 0;
    size_t singleheterozygous = 0;

    if (not unphased_population->genomePhasingStats(genome.first, false, heterozygous, homozygous, singleheterozygous)) {

      ExecEnv::log().error("No phasing statistics found for Genome: {}", genome.first);

    }

    auto result = hetero_results_.find(genome.first);

    if (result == hetero_results_.end()) {

      HeteroResultArray result_array;
      result_array[HOMOZYGOUS_] = homozygous;
      result_array[HETEROZYGOUS_] = heterozygous;
      result_array[SINGLEHETEROZYGOUS_] = singleheterozygous;

      std::pair<GenomeId_t , HeteroResultArray> add_genome(genome.first, result_array);

      auto insert_result = hetero_results_.insert(add_genome);

      if (not insert_result.second) {

        ExecEnv::log().error("HeterozygousStatistics; error inserting genome: {}", genome.first);
        return false;

      }

    } else {

      result->second[HOMOZYGOUS_] += homozygous;
      result->second[HETEROZYGOUS_] += heterozygous;
      result->second[SINGLEHETEROZYGOUS_] += singleheterozygous;

    }

  }

  return true;

}



bool kgl::HeterozygousStatistics::writeHeterozygousStatistics(const std::string& file_name, const char delimiter) const {
// Write data.

  std::ofstream out_file(file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("heterozygousStatistics; Unable to open output file: {}", file_name);
    return false;

  }

  // Write header.
  out_file << "genome" << delimiter;
  out_file << "homozygous" << delimiter;
  out_file << "heterozygous" << delimiter;
  out_file << "singleheterozygous" << delimiter;
  out_file << std::endl;

  for (auto genome : getMap()) {

    out_file << genome.first << delimiter;
    out_file << genome.second[HOMOZYGOUS_] << delimiter;
    out_file << genome.second[HETEROZYGOUS_] << delimiter;
    out_file << genome.second[SINGLEHETEROZYGOUS_] << delimiter;
    out_file << std::endl;

  }

  return out_file.good();

}
