//
// Created by kellerberrin on 11/08/18.
//

#include "kgl_distribution_analysis.h"
#include "kgl_sequence_complexity.h"

#include <fstream>


namespace kgl = kellerberrin::genome;



bool kgl::AggregateVariantDistribution::variantDistribution(std::shared_ptr<const UnphasedPopulation> unphased_population) {


  for (auto const& genome : unphased_population->getMap()) {

    for (auto const& contig : genome.second->getMap()) {

      for (auto const& offset : contig.second->getMap()) {

        for (auto const& variant_ptr : offset.second) {

            if (not addVariant(variant_ptr)) {

              ExecEnv::log().error( "AggregatVariantDistribution; Cannot add variant: {}",
                                    variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));

            } // if

          } // variants.

        } // offset vector

      } // contig

    } // genome

  ExecEnv::log().info("AggregatVariantDistribution; Unphased variants processed: {}", unphased_population->variantCount());

  return true;

}


bool kgl::AggregateVariantDistribution::variantDistribution(std::shared_ptr<const PhasedPopulation> population_ptr) {


  for (auto const& genome : population_ptr->getMap()) {

    for (auto const& contig : genome.second->getMap()) {

      for (auto const& homologous : contig.second->getVector()) {

        for (auto const& variant : homologous->getMap()) {

          if (not addVariant(variant.second)) {

            ExecEnv::log().error( "AggregatVariantDistribution; Cannot add variant: {}",
                                 variant.second->output(' ', VariantOutputIndex::START_0_BASED, false));

          }

        } // offset vector

      } // offset

    } // contig

  } // genome

  ExecEnv::log().info("AggregatVariantDistribution; Phased variants processed: {}", population_ptr->variantCount());

  return true;

}


bool kgl::AggregateVariantDistribution::addVariant(std::shared_ptr<const Variant> variant) {

  auto contig_map = interval_contig_map_.find(variant->contigId());

  if (contig_map == interval_contig_map_.end()) {

    std::pair<ContigId_t, IntervalVariantMap> insert_contig(variant->contigId(), IntervalVariantMap());
    auto insert_result = interval_contig_map_.insert(insert_contig);

    if (not insert_result.second) {

      ExecEnv::log().error("AggregateVariantDistribution::addVariant; unable to add contig: {}", variant->contigId());
      return false;

    }

    contig_map = insert_result.first;

  }

  std::pair<ContigOffset_t, std::shared_ptr<const Variant>> insert_offset(variant->offset(), variant);
  contig_map->second.insert(insert_offset);

  return true;

}





bool kgl::AggregateVariantDistribution::writeDistribution(std::shared_ptr<const GenomeReference> genome_db,
                                                          size_t interval_size,
                                                          const ContigId_t analysis_contig,
                                                          ContigOffset_t start_offset,
                                                          ContigOffset_t end_offset,
                                                          bool display_sequence,
                                                          const std::string& filename,
                                                          char delimiter) const {

  std::ofstream outfile(filename);

  if (not outfile.good()) {

    ExecEnv::log().error("AggregateVariantDistribution; could not open results file: {}", filename);
    return false;

  }

  if (not writeHeader(outfile, delimiter, display_sequence)) {

    ExecEnv::log().error("AggregateVariantDistribution; could not write header to file: {}", filename);
    return false;

  }
  if (not writeData(genome_db, interval_size, analysis_contig, start_offset, end_offset, display_sequence, outfile, delimiter)) {

    ExecEnv::log().error("AggregateVariantDistribution; could not write results to file: {}", filename);
    return false;

  }

  return outfile.good();

}



bool kgl::AggregateVariantDistribution::writeHeader(std::ostream& output, char delimiter, bool display_sequence) const {

  output << "Contig" << delimiter;
  output << "Interval_start" << delimiter;
  output << "Interval_size" << delimiter;
  output << "Snp_count" << delimiter;
  output << "Variant_count" << delimiter;
  output << "LempelZiv" << delimiter;
  output << "ShannonEntropy" << delimiter;
  output << "CpG" << delimiter;
  output << "Prop_A" << delimiter;
  output << "Prop_C" << delimiter;
  output << "Prop_G" << delimiter;
  output << "Prop_T";

  if (display_sequence) {

    output << delimiter << "Sequence" << '\n';

  } else {

    output << '\n';

  }

  return output.good();

}


bool kgl::AggregateVariantDistribution::writeData(std::shared_ptr<const GenomeReference> genome_db,
                                                  size_t interval_size,
                                                  const ContigId_t analysis_contig,
                                                  ContigOffset_t start_offset,
                                                  ContigOffset_t end_offset,
                                                  bool display_sequence,
                                                  std::ostream& output,
                                                  char delimiter) const {
  for (auto [contig_id, contig_ptr] : genome_db->getMap()) {

    ContigOffset_t contig_from = 0;
    ContigSize_t contig_size = contig_ptr->contigSize();

    if (not analysis_contig.empty()) {

      if (analysis_contig == contig_id) {

        contig_from = start_offset;
        if (end_offset > contig_from) {

          contig_size = end_offset;

        }

        ExecEnv::log().info("AggregateVariantDistribution; Contig: {}, From: {}, To: {}", analysis_contig, contig_from, contig_size);

      } else {

        continue;

      }

    }

    ContigOffset_t contig_to = interval_size + contig_from;

    auto contig_iter  = interval_contig_map_.find(contig_id);

    if (contig_iter != interval_contig_map_.end()) {

      while (contig_from < contig_size) {

        output << contig_id << delimiter;
        output << contig_from << delimiter;
        output << (contig_to - contig_from) << delimiter;

        size_t snp_count = 0;
        size_t variant_count = 0;

        auto begin_offset_iter = contig_iter->second.lower_bound(contig_from);
        auto end_offset_iter = contig_iter->second.upper_bound(contig_to);

        for (auto offset_iter = begin_offset_iter; offset_iter != end_offset_iter; ++offset_iter) {

          if (offset_iter->second->isSNP()) {

            ++snp_count;

          }

          ++variant_count;

        }

        output << snp_count << delimiter;
        output << variant_count << delimiter;


        if (contig_to >= contig_ptr->contigSize()) {

          ExecEnv::log().info("AggregateVariantDistribution::writeData; processed contig: {}", contig_id);

        }

        size_t interval_size = contig_to - contig_from;

        DNA5SequenceLinear sequence = contig_ptr->sequence_ptr()->subSequence(contig_from, interval_size);
        if (sequence.length() == 0) {

          ExecEnv::log().info("AggregateVariantDistribution::writeData; zero sized sequence, offset: {}, size: {}, contig: {} contig size: {}",
                              contig_from, interval_size, contig_id, contig_ptr->contigSize());

        }



        if (sequence.length() != interval_size) {

          ExecEnv::log().error("AggregateVariantDistribution::writeData; unexpected sequence size: {} returned from contig: {}, offset: {}, size: {}",
                                 sequence.length(), contig_id, contig_from, interval_size);
          break;

        }

        output << SequenceComplexity::complexityLempelZiv(sequence) << delimiter;
        output << SequenceComplexity::alphabetEntropy<DNA5>(sequence) << delimiter;
        output << SequenceComplexity::relativeCpGIslands(sequence) << delimiter;
        double A_prop;
        double C_prop;
        double G_prop;
        double T_prop;
        SequenceComplexity::proportionNucleotides(sequence, A_prop, C_prop, G_prop, T_prop);
        output << A_prop << delimiter;
        output << C_prop << delimiter;
        output << G_prop << delimiter;
        output << T_prop;

        if (display_sequence) {

          output << delimiter << sequence.getSequenceAsString() << '\n';

        } else {

          output << '\n';

        }

        contig_from = contig_to;
        if (contig_to + interval_size < contig_size) {

          contig_to += interval_size;

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

    if (not unphased_population->genomePhasingStats(genome.first, heterozygous, homozygous, singleheterozygous)) {

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
