//
// Created by kellerberrin on 15/06/18.
//


#include "kgl_variant_classify.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant_single.h"
#include <fstream>
#include <algorithm>


namespace kgl = kellerberrin::genome;



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects accept unphased variants and produces ref/alt read statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Classify the ref/alt counts by variant (row) and then genome (column)
kgl::VariantClassifier::VariantClassifier(std::shared_ptr<const UnphasedPopulation> vcf_population_ptr) {

  // Get a list of all genomes.
  genomes_ = vcf_population_ptr->genomeList();

  // for all variants index by contig/offset
  for (auto genome : vcf_population_ptr->getMap()) {

    for(auto contig : genome.second->getMap()) {

      for(auto offset : contig.second->getMap()) {

        for(auto variant : offset.second) {

          ContigOffsetVariants ordered_variant(variant);
          auto result = variant_map_.find(ordered_variant);

          if (result != variant_map_.end()) {


            // We ignore duplicate inserts, since this is just the variant caller calling
            // 2 identical variants under the ploidy = 2 setting. Both variants have identical ref/alt ratios.
            // with a small (generally zero) ref count and a much larger alt count.
            result->second.insert(GenomeOffsetMap::value_type(variant->genomeId(), variant));

          } else {

            GenomeOffsetMap genome_map;

            auto insert_result = variant_map_.insert(VariantMap::value_type(ordered_variant, genome_map));

            if (not insert_result.second) {

              ExecEnv::log().error("Attempt to insert duplicate variant: {}",
                                   variant->output(' ',VariantOutputIndex::START_0_BASED, false));

            } else {

              // We ignore duplicate inserts, since this is just the variant caller calling
              // 2 identical variants under the ploidy = 2 setting. Both variants have identical ref/alt ratios.
              // with a small (generally zero) ref count and a much larger alt count.
              insert_result.first->second.insert(GenomeOffsetMap::value_type(variant->genomeId(), variant));

            }

          } // if found

        } // for variant

      } // for offset

    } // for contig

  } // for genome

}

// Write ref/alt counts into separate files by variant order.
bool kgl::VariantClassifier::writeVariants(char delimiter,
                                           const std::string& ref_file_name,
                                           const std::string& alt_file_name,
                                           size_t min_count) const {

  std::ofstream ref_count_file(ref_file_name);

  if (not ref_count_file.good()) {

    ExecEnv::log().error("VariantClassifier; Unable to open ref count file: {}", ref_file_name);
    return false;

  }

  std::ofstream alt_count_file(alt_file_name);

  if (not alt_count_file.good()) {

    ExecEnv::log().error("VariantClassifier; Unable to open alt count file: {}", alt_file_name);
    return false;

  }

  // write the header.
  ref_count_file << "Contig" << delimiter << "Offset" << delimiter;
  alt_count_file << "Contig" << delimiter << "Offset" << delimiter;

  for (auto genome : getGenomes()) {

    ref_count_file << genome << delimiter;
    alt_count_file << genome << delimiter;

  }

  ref_count_file << std::endl;
  alt_count_file << std::endl;

  // for all variant offsets.
  for (auto variant_offset : getMap()) {

    ref_count_file << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset() << delimiter;
    alt_count_file << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset() << delimiter;

    for (auto genome : getGenomes()) {

      auto find_result = variant_offset.second.find(genome);

      if (find_result != variant_offset.second.end()) {

        std::shared_ptr<const CountEvidence> count_evidence_ptr = std::dynamic_pointer_cast<const CountEvidence>(find_result->second->evidence());

        if (count_evidence_ptr) {

          size_t total_count = count_evidence_ptr->refCount() + count_evidence_ptr->altCount();

          if (total_count >= min_count) {

            ref_count_file << count_evidence_ptr->refCount() << delimiter;
            alt_count_file << count_evidence_ptr->altCount() << delimiter;

          } else {

            ref_count_file << '0' << delimiter;
            alt_count_file << '0' << delimiter;

          }

        } else {

          ExecEnv::log().error("Variant without count evidence: {}", find_result->second->output(' ', VariantOutputIndex::START_0_BASED, false));

        }

      } else { // genome not found

        ref_count_file << '0' << delimiter;
        alt_count_file << '0' << delimiter;

      }

    }

    ref_count_file << std::endl;
    alt_count_file << std::endl;

  }

  return true;

}


// Write ref/alt counts into separate files by variant order.
bool kgl::VariantClassifier::writePlaf(char delimiter,
                                       const std::string& plaf_file_name,
                                       const std::string& ref_file_prefix,
                                       const std::string& alt_file_prefix,
                                       size_t min_count) const {

  std::ofstream plaf_file(plaf_file_name);

  if (not plaf_file.good()) {

    ExecEnv::log().error("VariantClassifier; Unable to open plaf file: {}", plaf_file_name);
    return false;

  }

  plaf_file << "#Contig" <<	delimiter << "Offset"	<< delimiter << "VariantFreq" << '\n';

  std::vector<std::shared_ptr<std::ofstream>> genome_ref_files;
  std::vector<std::shared_ptr<std::ofstream>> genome_alt_files;
  // Open all the genome files.
  for (auto genome : getGenomes()) {

    std::string ref_file_name = ref_file_prefix + genome + ".tab";

    std::shared_ptr<std::ofstream> ref_file_ptr(std::make_shared<std::ofstream>(ref_file_name));
    if (not ref_file_ptr->good()) {

      ExecEnv::log().error("VariantClassifier; Unable to open ref file: {}", ref_file_name);
      return false;

    }

    std::string alt_file_name = alt_file_prefix + genome + ".tab";

    std::shared_ptr<std::ofstream> alt_file_ptr(std::make_shared<std::ofstream>(alt_file_name));
    if (not alt_file_ptr->good()) {

      ExecEnv::log().error("VariantClassifier; Unable to open ref file: {}", alt_file_name);
      return false;

    }

    genome_ref_files.push_back(ref_file_ptr);
    genome_alt_files.push_back(alt_file_ptr);

  }

  if (genome_ref_files.size() != getGenomes().size()) {

    ExecEnv::log().error("VariantClassifier; Mismatch between genomes size: {} and number of genome ref files: {}",
                         getGenomes().size(), genome_ref_files.size());
    return false;


  }

  if (genome_alt_files.size() != getGenomes().size()) {

    ExecEnv::log().error("VariantClassifier; Mismatch between genomes size: {} and number of genome alt files: {}",
                         getGenomes().size(), genome_alt_files.size());
    return false;


  }

  // write the header to each file.
  for (auto file_ptr : genome_ref_files) {

    (*file_ptr) << "#Contig" <<	delimiter << "Offset"	<< delimiter << "Variant" << '\n';

  }

  for (auto file_ptr : genome_alt_files) {

    (*file_ptr) << "#Contig" <<	delimiter << "Offset"	<< delimiter << "Variant" << '\n';

  }

  // for all variant offsets.
  for (auto variant_offset : getMap()) {

    size_t allele_count = 0;
    size_t genome_offset = 0;
    for (auto genome : getGenomes()) {

      auto find_result = variant_offset.second.find(genome);

      if (find_result != variant_offset.second.end()) {

        std::shared_ptr<const CountEvidence> count_evidence_ptr = std::dynamic_pointer_cast<const CountEvidence>(
        find_result->second->evidence());

        if (count_evidence_ptr) {

          size_t total_count = count_evidence_ptr->refCount() + count_evidence_ptr->altCount();

          if (total_count >= min_count) {

            ++allele_count;
            *(genome_ref_files[genome_offset]) << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                                               << delimiter << count_evidence_ptr->refCount() << '\n';

            *(genome_alt_files[genome_offset]) << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                                               << delimiter << count_evidence_ptr->altCount() << '\n';


          } else {

            *(genome_ref_files[genome_offset]) << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                                               << delimiter << static_cast<size_t>(0) << '\n';

            *(genome_alt_files[genome_offset]) << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                                               << delimiter << static_cast<size_t>(0) << '\n';

          }

        } else {

          ExecEnv::log().error("Variant without count evidence: {}",
                               find_result->second->output(' ', VariantOutputIndex::START_0_BASED, false));

        }

      } else {

        *(genome_ref_files[genome_offset]) << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                                           << delimiter << static_cast<size_t>(0) << '\n';

        *(genome_alt_files[genome_offset]) << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                                           << delimiter << static_cast<size_t>(0) << '\n';

      }

      ++genome_offset;

    } // all genomes;

    if (allele_count > 0) {

      double allele_population_freq = static_cast<double>(allele_count) / static_cast<double>(getGenomes().size());

      plaf_file << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                << delimiter << allele_population_freq << std::endl;

    } else {


      plaf_file << variant_offset.first.variant()->contigId() << delimiter << variant_offset.first.variant()->offset()
                << delimiter << static_cast<double>(0.0) << std::endl;

    }

  }  // all variants.

  return true;

}


bool kgl::VariantClassifier::writeOrderedVariants(char delimiter,
                                                  const std::string& ref_file_name,
                                                  const std::string& alt_file_name,
                                                  size_t min_count,
                                                  size_t max_lines) const {

  OrderedAltCount ordered_count;
  orderVariantCount(ordered_count, min_count);
  return writeOrderedCount(ordered_count, delimiter, ref_file_name, alt_file_name, max_lines);

}


// Write ref/alt counts simultaneously into separate files by total genome alt count order.
// The maximum number of genomes with non-zero alt counts.
// This is used in initial 'R' testing of the the J. O'brien MCMC mixture algorithm
bool kgl::VariantClassifier::orderVariantCount(OrderedAltCount& ordered_count, size_t min_count) const {

  // for all variant offsets.
  for (auto variant_offset : getMap()) {

    std::vector<size_t> ref_count;
    std::vector<size_t> alt_count;

    for (auto genome : getGenomes()) {

      auto find_result = variant_offset.second.find(genome);

      if (find_result != variant_offset.second.end()) {

        std::shared_ptr<const CountEvidence> count_evidence_ptr = std::dynamic_pointer_cast<const CountEvidence>(find_result->second->evidence());

        if (count_evidence_ptr) {

          size_t total_count = count_evidence_ptr->refCount() + count_evidence_ptr->altCount();

          if (total_count >= min_count) {

            ref_count.push_back(count_evidence_ptr->refCount());
            alt_count.push_back(count_evidence_ptr->altCount());

          } else {

            ref_count.push_back(0);
            alt_count.push_back(0);

          }

        } else {

          ExecEnv::log().error("Variant without count evidence: {}", find_result->second->output(' ', VariantOutputIndex::START_0_BASED, false));

        }

      } else { // genome not found

        ref_count.push_back(0);
        alt_count.push_back(0);

      }

    } // for genomes.

    auto non_zero_lambda = [](size_t count){ return count > 0; };
    size_t non_zero = std::count_if (alt_count.begin(), alt_count.end(), non_zero_lambda);
    non_zero += std::count_if (ref_count.begin(), ref_count.end(), non_zero_lambda);
    std::pair<std::vector<size_t>, std::vector<size_t>> count_pair(ref_count, alt_count);
    ordered_count.insert(OrderedAltCount::value_type(non_zero, count_pair));

  } // for variants.

  return true;

}



bool kgl::VariantClassifier::writeOrderedCount(OrderedAltCount& ordered_count,
                                               char delimiter,
                                               const std::string& ref_file_name,
                                               const std::string& alt_file_name,
                                               size_t max_lines) const {

  std::ofstream ref_count_file(ref_file_name);

  if (not ref_count_file.good()) {

    ExecEnv::log().error("VariantClassifier; Unable to open ref count file: {}", ref_file_name);
    return false;

  }

  std::ofstream alt_count_file(alt_file_name);

  if (not alt_count_file.good()) {

    ExecEnv::log().error("VariantClassifier; Unable to open alt count file: {}", alt_file_name);
    return false;

  }
/*
  for (auto genome : getGenomes()) {

    ref_count_file << genome << delimiter;
    alt_count_file << genome << delimiter;

  }

  ref_count_file << std::endl;
  alt_count_file << std::endl;
*/
  size_t line_count = 0;
  for (auto count_rit = ordered_count.rbegin();  count_rit != ordered_count.rend(); ++count_rit) {

    bool first = true;
    for (auto ref_count : count_rit->second.first) {

      if (not first) ref_count_file << delimiter;
      ref_count_file << ref_count;
      if (first) first = false;

    }

    ref_count_file << std::endl;

    first = true;
    for (auto alt_count : count_rit->second.second) {

      if (not first) alt_count_file << delimiter;
      alt_count_file << alt_count;
      if (first) first = false;

    }

    alt_count_file << std::endl;

    ++line_count;
    if (line_count >= max_lines) break;

  }

  return true;

}

