//
// Created by kellerberrin on 5/5/20.
//

#include "kgl_analysis_interval.h"
#include "kgl_sequence_complexity.h"

#include <fstream>


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::IntervalAnalysis::initializeAnalysis( const std::string& work_directory,
                                                const RuntimeParameterMap& named_parameters,
                                                std::shared_ptr<const GenomeCollection> reference_genomes) {

  if (reference_genomes->getMap().size() > 1 or reference_genomes->getMap().empty()) {

    ExecEnv::log().error("Analytic: {} called with {} genomes.  Only 1 genome can be analysed at a time.",
                    ident(), reference_genomes->getMap().size());
    return false;

  }

  if (not getParameters(work_directory, named_parameters)) {

    return false;

  }

  // The first genome is the only genome.
  std::shared_ptr<const GenomeReference> genome = reference_genomes->getMap().begin()->second;
  setupIntervalStructure(genome);

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::IntervalAnalysis::iterateAnalysis( std::shared_ptr<const GenomeCollection>, std::shared_ptr<const UnphasedPopulation> population) {

  return variantIntervalCount(population);

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::IntervalAnalysis::finalizeAnalysis(std::shared_ptr<const GenomeCollection> reference_genomes) {

  // The first genome is the only genome.
  std::shared_ptr<const GenomeReference> genome = reference_genomes->getMap().begin()->second;

  return writeResults(genome, false, OUTPUT_DELIMITER_);

}

bool kgl::IntervalAnalysis::getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters) {

  // Get the analysis interval
  auto result = named_parameters.find(INTERVAL_SIZE_);
  if (result == named_parameters.end()) {

    ExecEnv::log().error("Analytic: {}; Expected Parameter: {} to be defined. {} is deactivated. Available named Parameters:", ident(), INTERVAL_SIZE_, ident());
    for (auto const& [parameter_ident, parameter_value] : named_parameters) {

      ExecEnv::log().info("initializeAnalysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

    }

    return false;

  }

  try {

    interval_size_ = std::stoll(result->second);

  }
  catch(...) {

    ExecEnv::log().error("Analytic: {}; Expected Parameter: {} Value: {} is invalid, expected unsigned integer. {} is deactivated.", ident(), INTERVAL_SIZE_, result->second, ident());
    return false;

  }
  // Get the output filename
  result = named_parameters.find(OUTPUT_FILE_);
  if (result == named_parameters.end()) {
    ExecEnv::log().error("Analytic: {}; Expected Parameter: {} to be defined. {} is deactivated. Available named Parameters:", ident(), OUTPUT_FILE_, ident());
    for (auto const& [parameter_ident, parameter_value] : named_parameters) {

      ExecEnv::log().info("initializeAnalysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

    }
    return false;
  }
  output_file_name_ = Utility::filePath(result->second, work_directory);

  ExecEnv::log().info("Analysis: {}, initialized with interval size: {}, output file: {}", ident(), interval_size_, output_file_name_);

  return true;

}

void kgl::IntervalAnalysis::setupIntervalStructure(std::shared_ptr<const GenomeReference> genome) {

  for (auto const& [contig_id, contig_ptr] : genome->getMap()) {

    ContigSize_t contig_size = contig_ptr->contigSize();
    size_t vector_size = (contig_size / interval_size_) + 1;
    auto result = interval_map_.insert(std::pair<ContigId_t, IntervalVector>(contig_id, IntervalVector(vector_size, {0,0})));
    if (not result.second) {

      ExecEnv::log().error("IntervalAnalysis::setupIntervalStructure, Genome: {} Duplicate contig: {}", genome->genomeId(), contig_id);

    }

  }

}


bool kgl::IntervalAnalysis::variantIntervalCount(std::shared_ptr<const UnphasedPopulation> population_ptr) {

  size_t variant_count{0};
  for (auto const& genome : population_ptr->getMap()) {

    for (auto const& [contig_id, contig_ptr] : genome.second->getMap()) {

      auto result = interval_map_.find(contig_id);
      if (result == interval_map_.end()) {

        ExecEnv::log().error("IntervalAnalysis::variantIntervalCount; Cannot find contig: {} mismatch between Reference Genome and Variant Population");
        return false;

      }

      IntervalVector& interval_vector = result->second;
      for (auto const& [offset, array] : contig_ptr->getMap()) {

        size_t count_index = offset /  interval_size_;
        if (count_index >= interval_vector.size()) {

          ExecEnv::log().error("IntervalAnalysis::variantIntervalCount; Calculated Offset: {} exceeds count array size: {}", count_index, result->second.size());
          return false;

        }


        interval_vector[count_index].first += array.size();
        variant_count += array.size();
        for (auto const& variant : array) {

          if (variant->isSNP()) {

            ++interval_vector[count_index].second;

          }

        } // array

      } // offset

    } // contig

  } // genome

  ExecEnv::log().info("IntervalAnalysis::variantIntervalCount; Variants processed: {}", variant_count);

  return true;

}


bool kgl::IntervalAnalysis::writeResults( std::shared_ptr<const GenomeReference> genome_db,
                                          bool display_sequence,
                                          char delimiter) const {

  std::ofstream outfile(output_file_name_);

  if (not outfile.good()) {

    ExecEnv::log().error("AggregateVariantDistribution; could not open results file: {}", output_file_name_);
    return false;

  }

  if (not writeHeader(outfile, delimiter, display_sequence)) {

    ExecEnv::log().error("AggregateVariantDistribution; could not write header to file: {}", output_file_name_);
    return false;

  }
  if (not writeData(genome_db, display_sequence, outfile, delimiter)) {

    ExecEnv::log().error("AggregateVariantDistribution; could not write results to file: {}", output_file_name_);
    return false;

  }

  return outfile.good();

}


bool kgl::IntervalAnalysis::writeHeader(std::ostream& output, char delimiter, bool display_sequence) const {

  output << "Contig" << delimiter;
  output << "Interval_start" << delimiter;
  output << "Interval_size" << delimiter;
  output << "Snp_count" << delimiter;
  output << "Variant_count" << delimiter;
  output << "LempelZiv" << delimiter;
  output << "ShannonEntropy" << delimiter;
  output << "CpG" << delimiter;
  output << "A%" << delimiter;
  output << "C%" << delimiter;
  output << "G%" << delimiter;
  output << "T%";

  if (display_sequence) {

    output << delimiter << "Sequence" << '\n';

  } else {

    output << '\n';

  }

  return output.good();

}


bool kgl::IntervalAnalysis::writeData( std::shared_ptr<const GenomeReference> genome_db,
                                       bool display_sequence,
                                       std::ostream& output,
                                       char delimiter) const {

  for (auto [contig_id, contig_ptr] : genome_db->getMap()) {

    ContigOffset_t contig_from = 0;
    ContigSize_t contig_size = contig_ptr->contigSize();
    ContigOffset_t contig_to = interval_size_ + contig_from;
    auto contig_iter  = interval_map_.find(contig_id);
    size_t count_index = 0;

    if (contig_iter != interval_map_.end()) {

      const IntervalVector& interval_vector = contig_iter->second;
      while (contig_from < contig_size) {

        output << contig_id << delimiter;
        output << contig_from << delimiter;
        output << (contig_to - contig_from) << delimiter;

        output << interval_vector[count_index].second << delimiter; // snp
        output << interval_vector[count_index].first << delimiter; // variant


        if (contig_to >= contig_ptr->contigSize()) {

          ExecEnv::log().info("IntervalAnalysis::writeData; processed contig: {}", contig_id);

        }

        size_t interval_size = contig_to - contig_from;

        DNA5SequenceLinear sequence = contig_ptr->sequence_ptr()->subSequence(contig_from, interval_size);
        if (sequence.length() == 0) {

          ExecEnv::log().info("IntervalAnalysis::writeData; zero sized sequence, offset: {}, size: {}, contig: {} contig size: {}",
                              contig_from, interval_size, contig_id, contig_ptr->contigSize());

        }

        if (sequence.length() != interval_size) {

          ExecEnv::log().error("IntervalAnalysis::writeData; unexpected sequence size: {} returned from contig: {}, offset: {}, size: {}",
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
        output << (A_prop * 100.0) << delimiter;
        output << (C_prop * 100.0) << delimiter;
        output << (G_prop * 100.0) << delimiter;
        output << (T_prop * 100.0);

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

        ++count_index;

      }

    }

  }

  return output.good();

}
