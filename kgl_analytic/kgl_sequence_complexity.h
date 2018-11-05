//
// Created by kellerberrin on 26/10/18.
//

#ifndef KGL_SEQUENCE_COMPLEXITY_H
#define KGL_SEQUENCE_COMPLEXITY_H

#include "kgl_sequence_base.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class SequenceComplexity {

public:

  SequenceComplexity() = delete;
  ~SequenceComplexity() = delete;

  template<typename Alphabet>
  static size_t kmerCount(std::shared_ptr<const AlphabetSequence<Alphabet>> sequence,
                          std::shared_ptr<const AlphabetSequence<Alphabet>> kmer);

  static double relativeCpGIslands(std::shared_ptr<const DNA5SequenceLinear> sequence);
  // Calculate Entropy
  template<typename Alphabet>
  static double alphabetEntropy(std::shared_ptr<const AlphabetSequence<Alphabet>> sequence);
  // Different for different sequence lengths.
  static size_t complexityLempelZiv(std::shared_ptr<const DNA5SequenceLinear> sequence);
  // Calculate Nucleotide content.
  static void proportionNucleotides(std::shared_ptr<const DNA5SequenceLinear> sequence,
                                    double& A_prop,
                                    double& C_prop,
                                    double& G_prop,
                                    double& T_prop);

private:

  // Calculate GC content.
  static void countNucleotides(std::shared_ptr<const DNA5SequenceLinear> sequence,
                               size_t& A_count,
                               size_t& C_count,
                               size_t& G_count,
                               size_t& T_count);

};


// Calculate Entropy
template<typename Alphabet>
double SequenceComplexity::alphabetEntropy(std::shared_ptr<const AlphabetSequence<Alphabet>> sequence) {

  if (sequence->length() == 0) {

    ExecEnv::log().error("SequenceComplexity::alphabetEntropy; zero sized sequence");
    return 0.0;

  }

  auto enumerated_alphabet = Alphabet::enumerateAlphabet();

  std::map<char, size_t> lookup_map;

  // Create a lookup map for efficiency
  for (auto letter : enumerated_alphabet) {

    std::pair<char, size_t> insert_pair(Alphabet::convertToChar(letter), 0);

    auto result = lookup_map.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("SequenceComplexity::alphabetEntropy; duplicate alphabet item: {}, for sequence: {}",
                           Alphabet::convertToChar(letter), sequence->getSequenceAsString());

    }

  }

  // Count the sequence into the lookup map.
  for (size_t idx = 0; idx < sequence->length(); ++idx) {

    auto result = lookup_map.find(Alphabet::convertToChar(sequence->at(idx)));

    if (result == lookup_map.end()) {

      ExecEnv::log().error("SequenceComplexity::alphabetEntropy; missing alphabet item: {}, for sequence: {}",
                           Alphabet::convertToChar(sequence->at(idx)), sequence->getSequenceAsString());

    } else {

      ++(result->second);

    }

  }

  double sequence_size = static_cast<double>(sequence->length());

  // Calculate entropy.
  double entropy = 0.0;

  for (auto letter : lookup_map) {

    if (letter.second > 0) {

      double ratio = static_cast<double>(letter.second) / sequence_size;

      entropy += ratio * std::log(ratio);

    }

  }

  // Adjust for alphabet size
  entropy = entropy * (-1.0 / std::log(static_cast<double>(enumerated_alphabet.size())));

  return entropy;

}



template<typename Alphabet>
size_t SequenceComplexity::kmerCount(std::shared_ptr<const AlphabetSequence<Alphabet>> sequence,
                                     std::shared_ptr<const AlphabetSequence<Alphabet>> kmer) {

  std::map<std::string, size_t> word_map;

  std::string sequence_string = sequence->getSequenceAsString();

  std::string kmer_string = kmer->getSequenceAsString();

  size_t kmer_size = kmer->length();

  std::pair<std::string, size_t> insert_kmer(kmer_string, 0);

  auto result = word_map.insert(insert_kmer);

  if (not result.second) {

    ExecEnv::log().error("SequenceComplexity::kmerCount; unable to insert kmer: {} for sequence: {}", kmer_string, sequence_string);
    return 0;

  }

  for (size_t i = 0; i < sequence->length(); i++) {

    if (i + kmer_size >= sequence->length()) break;

    std::string word = sequence_string.substr(i, kmer_size);

    auto result = word_map.find(word);

    if (result != word_map.end()) {

      ++(result->second);

    }

  }

  return word_map.begin()->second;

}





}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_SEQUENCE_COMPLEXITY_H
