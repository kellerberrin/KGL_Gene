//
// Created by kellerberrin on 5/5/20.
//

#include "kgl_analysis_interval.h"
#include "kgl_sequence_complexity.h"

#include <fstream>


namespace kgl = kellerberrin::genome;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IntervalData members.

void kgl::IntervalData::addArrayVariantCount(size_t size) {

  if (size < ARRAY_VARIANT_COUNT_) {

    ++array_variant_count_[size - 1];

  } else {

    ++array_variant_count_[ARRAY_VARIANT_COUNT_ - 1];

  }

  ++variant_offset_count_;

}

void kgl::IntervalData::emptyIntervalOffset(const ContigOffset_t& previous_variant_offset, const ContigOffset_t& variant_offset) {

  if (variant_offset < offset() or variant_offset >= (offset() + interval())) {

    ExecEnv::log().error("IntervalData::emptyIntervalOffset; Contig: {}, Offset: {} out of Range, Interval Offset: {}, Interval size: {}",
                         contigId(), variant_offset, offset(), interval());

  }

  ContigSize_t relative_offset = variant_offset - offset();
  SignedOffset_t empty_interval = variant_offset - previous_variant_offset;

  if (empty_interval < 0) {

    ExecEnv::log().error("IntervalData::emptyIntervalOffset; Contig: {}, Negative Interval : {}, Offset : {}, Previous Offset : {}",
                        contigId(), empty_interval, variant_offset, previous_variant_offset);
    return;

  }

  if (empty_interval > max_empty_interval_.second) {

    max_empty_interval_.first = relative_offset;
    max_empty_interval_.second = empty_interval;

  }

  sum_empty_interval_ += empty_interval;

}

double kgl::IntervalData::meanEmptyInterval() const {

  return static_cast<double>(sum_empty_interval_) / static_cast<double>(variant_offset_count_ + 1);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IntervalAnalysis members.


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

      ExecEnv::log().info("Analysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

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

      ExecEnv::log().info("Analysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

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
    IntervalVector interval_vector;
    for (size_t index = 0; index < vector_size; ++index) {

      ContigOffset_t offset = index *  interval_size_;
      ContigSize_t size = (contig_size - offset) >= interval_size_ ? interval_size_ : (contig_size - offset);
      interval_vector.emplace_back(contig_id, offset, size);

    }
    auto result = interval_map_.insert(std::pair<ContigId_t, IntervalVector>(contig_id, std::move(interval_vector)));
    if (not result.second) {

      ExecEnv::log().error("IntervalAnalysis::setupIntervalStructure, Genome: {} Duplicate contig: {}", genome->genomeId(), contig_id);

    }

  }

}


bool kgl::IntervalAnalysis::variantIntervalCount(std::shared_ptr<const UnphasedPopulation> population_ptr) {

// We are profiling variants against a reference genome. Therefore we need to compress the population of variants
// into a single genome.
  std::shared_ptr<const UnphasedGenome> compressed_genome = population_ptr->compressPopulation();

  size_t variant_count{0};
  // For contigs in the compressed genome.
  for (auto const& [contig_id, contig_ptr] : compressed_genome->getMap()) {

    auto result = interval_map_.find(contig_id);
    if (result == interval_map_.end()) {

      ExecEnv::log().error("IntervalAnalysis::variantIntervalCount; Cannot find contig: {} mismatch between Reference Genome and Variant Population", contig_id);
      return false;

    }

    IntervalVector& interval_vector = result->second;
    // For all intervals.
    for (auto& interval_data : interval_vector) {

      auto lower_bound = contig_ptr->getMap().lower_bound(interval_data.offset());
      ContigOffset_t upperbound_offset = interval_data.offset() + interval_data.interval() - 1;
      auto upper_bound = contig_ptr->getMap().upper_bound(upperbound_offset);

      // For all variant array within the interval.
      ContigOffset_t previous_offset = interval_data.offset() - 1;
      for (auto array_ptr = lower_bound; array_ptr != upper_bound; ++array_ptr) {

        interval_data.emptyIntervalOffset(previous_offset, array_ptr->first); // Variant empty interval calculated from this.
        interval_data.addVariantCount(array_ptr->second.size());
        interval_data.addArrayVariantCount(array_ptr->second.size());
        variant_count += array_ptr->second.size();
        size_t snp_count{0};

        // Count SNP.
        for (auto const& variant : array_ptr->second) {

          if (variant->isSNP()) {

            ++snp_count;

          }

        } // count.

        interval_data.addSNPCount(snp_count);
        previous_offset = array_ptr->first;

      } // variant array.

      // The empty interval to the end of the data interval.
      interval_data.emptyIntervalOffset(previous_offset, upperbound_offset);

    } // interval

  } // contig

  ExecEnv::log().info("Analysis: {},  Variants processed: {}", ident(), variant_count);

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
  output << "Variant_count" << delimiter;
  output << "Snp_count" << delimiter;
  output << "Variant_offset" << delimiter;
//  output << "Variant=1" << delimiter;
//  output << "Variant=2" << delimiter;
//  output << "Variant=3" << delimiter;
//  output << "Variant=4" << delimiter;
//  output << "Variant>=5" << delimiter;
  output << "MaxEmptyOffset" << delimiter;
  output << "MaxEmptyInterval" << delimiter;
  output << "AvEmptyInterval" << delimiter;
  std::vector<DNA5::Alphabet> symbols = DNA5::enumerateAlphabet();
  for (auto const symbol : symbols) {

    output << static_cast<char>(symbol) << delimiter;

  }
//  output << "ZivLempel" << delimiter;
  output << "ShannonEntropy" << delimiter;
  output << "CpG";

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


    auto contig_iter  = interval_map_.find(contig_id);
    if (contig_iter == interval_map_.end()) {

      ExecEnv::log().error("IntervalAnalysis::writeData; could not find variant interval vector for contig: {}", contig_id);
      continue; // next contig.

    }

    const IntervalVector& interval_vector = contig_iter->second;

    if (interval_vector.empty()) {

      ExecEnv::log().warn("IntervalAnalysis::writeData; zero sized interval vector for contig: {}", contig_id);
      continue; // next contig.

    }

    ExecEnv::log().info("IntervalAnalysis::writeData; processing contig: {}", contig_id);
    ContigOffset_t contig_offset = 0;
    ContigSize_t contig_size = contig_ptr->contigSize();

    for (size_t count_index = 0; count_index <  interval_vector.size(); ++count_index) {

      if (contig_offset > contig_size) {

        ExecEnv::log().error("IntervalAnalysis::writeData; calculated offset; {} exceeds contig size: {}", contig_offset, contig_size);
        break;

      }

      size_t interval_size = 0;
      if ((contig_offset + interval_size_) <  contig_size) {

        interval_size =  interval_size_;

      } else {

        interval_size = contig_size - contig_offset;

      }

      DNA5SequenceLinear sequence = contig_ptr->sequence_ptr()->subSequence(contig_offset, interval_size);
      if (sequence.length() == 0) {

        ExecEnv::log().warn("IntervalAnalysis::writeData; zero sized sequence, offset: {}, size: {}, contig: {} contig size: {}",
                            contig_offset, interval_size, contig_id, contig_ptr->contigSize());
        break;

      }

      if (sequence.length() != interval_size) {

        ExecEnv::log().error("IntervalAnalysis::writeData; unexpected sequence size: {} returned from contig: {}, offset: {}, size: {}",
                             sequence.length(), contig_id, contig_offset, interval_size);
        break;

      }

      output << interval_vector[count_index].contigId() << delimiter;
      output << interval_vector[count_index].offset() << delimiter;
      output << interval_vector[count_index].interval() << delimiter;

      output << interval_vector[count_index].variantCount() << delimiter; // variant
      output << interval_vector[count_index].SNPCount() << delimiter; // snp
      output << interval_vector[count_index].variantOffsetCount() << delimiter; // variant offsets.
//      output << interval_vector[count_index].arrayVariantCount()[0] << delimiter;
//      output << interval_vector[count_index].arrayVariantCount()[1] << delimiter;
//      output << interval_vector[count_index].arrayVariantCount()[2] << delimiter;
//      output << interval_vector[count_index].arrayVariantCount()[3] << delimiter;
//      output << interval_vector[count_index].arrayVariantCount()[4] << delimiter;
      output << interval_vector[count_index].maxEmptyInterval().first << delimiter;
      output << interval_vector[count_index].maxEmptyInterval().second << delimiter;
      output << interval_vector[count_index].meanEmptyInterval() << delimiter;

      // count the symbols in the sequence
      const std::vector<std::pair<DNA5::Alphabet, size_t>>& symbol_vector = sequence.countSymbols();

      for (auto const& count : symbol_vector) {

        output << (static_cast<double>(count.second) * 100.0) / static_cast<double>(sequence.length()) << delimiter;

      }

//      output << SequenceComplexity::complexityLempelZiv(sequence) << delimiter;
      output << SequenceComplexity::alphabetEntropy<DNA5>(sequence, symbol_vector) << delimiter;
      output << SequenceComplexity::relativeCpGIslands(sequence);
      // Output the relative symbol proportions.

      if (display_sequence) {

        output << delimiter << sequence.getSequenceAsString() << '\n';

      } else {

        output << '\n';

      }

      contig_offset += interval_size; // next interval;

    } // offset within contig.

    ExecEnv::log().info("IntervalAnalysis::writeData; processed contig: {}", contig_id);

  } // contig

  return output.good();

}
