//
// Created by kellerberrin on 26/10/18.
//

#ifndef KGL_SEQUENCE_COMPLEXITY_H
#define KGL_SEQUENCE_COMPLEXITY_H

#include "kgl_sequence/kgl_sequence_base.h"

namespace kellerberrin::genome {   //  organization level namespace


class SequenceComplexity {

public:

  SequenceComplexity() = delete;
  ~SequenceComplexity() = delete;

  template<typename Alphabet>
  [[nodiscard]] static size_t kmerCount( const AlphabetSequence<Alphabet>& sequence,
                                         const AlphabetSequence<Alphabet>& kmer);

  [[nodiscard]] static double relativeCpGIslands(const DNA5SequenceCoding& sequence);
  [[nodiscard]] static double relativeCpGIslands(const DNA5SequenceLinear& sequence);
  // Calculate Entropy
  template<typename Alphabet>
  [[nodiscard]] static double alphabetEntropy( const AlphabetSequence<Alphabet>& sequence,
                                               const std::vector<std::pair<typename Alphabet::Alphabet, size_t>>& symbol_vector);
  // Different for different sequence lengths.
  template<typename Alphabet>
  [[nodiscard]] static size_t complexityLempelZiv(const AlphabetSequence<Alphabet>& sequence);
  // Calculate Nucleotide content.

private:

};

inline double SequenceComplexity::relativeCpGIslands(const DNA5SequenceCoding& sequence) {

  size_t count = sequence.AlphabetSequence<CodingDNA5>::countTwoSymbols(CodingDNA5::Alphabet::C, CodingDNA5::Alphabet::G);

  // Combinatorial theory tells us there is an expected CpG every 32 random nucleotides.
  return (static_cast<double>(count) * 32.0) / static_cast<double>(sequence.length());

}

inline double SequenceComplexity::relativeCpGIslands(const DNA5SequenceLinear& sequence) {

  size_t count = sequence.AlphabetSequence<DNA5>::countTwoSymbols(DNA5::Alphabet::C, DNA5::Alphabet::G);

  // Combinatorial theory tells us there is an expected CpG every 32 random nucleotides.
  return (static_cast<double>(count) * 32.0) / static_cast<double>(sequence.length());

}


// Calculate Entropy
template<typename Alphabet>
double SequenceComplexity::alphabetEntropy( const AlphabetSequence<Alphabet>& sequence,
                                            const std::vector<std::pair<typename Alphabet::Alphabet, size_t>>& symbol_vector) {

  if (sequence.length() == 0) {

    ExecEnv::log().error("SequenceComplexity::alphabetEntropy; zero sized sequence");
    return 0.0;

  }

  auto sequence_size = static_cast<double>(sequence.length());

  // Calculate entropy.
  double entropy = 0.0;

  for (auto symbol : symbol_vector) {

    if (symbol.second > 0) {

      double ratio = static_cast<double>(symbol.second) / sequence_size;

      entropy += ratio * std::log(ratio);

    }

  }

  // Adjust for alphabet size
  entropy = entropy * (-1.0 / std::log(static_cast<double>(symbol_vector.size())));

  return entropy;

}



template<typename Alphabet>
size_t SequenceComplexity::kmerCount(const AlphabetSequence<Alphabet>& sequence, const AlphabetSequence<Alphabet>& kmer) {

  std::map<std::string, size_t> word_map;

  std::string sequence_string = sequence.getSequenceAsString();

  std::string kmer_string = kmer.getSequenceAsString();

  size_t kmer_size = kmer.length();

  std::pair<std::string, size_t> insert_kmer(kmer_string, 0);

  auto result = word_map.insert(insert_kmer);

  if (not result.second) {

    ExecEnv::log().error("SequenceComplexity::kmerCount; unable to insert kmer: {} for sequence: {}", kmer_string, sequence_string);
    return 0;

  }

  for (size_t i = 0; i < sequence.length(); i++) {

    if (i + kmer_size >= sequence.length()) break;

    std::string word = sequence_string.substr(i, kmer_size);

    auto result = word_map.find(word);

    if (result != word_map.end()) {

      ++(result->second);

    }

  }

  return word_map.begin()->second;

}

template<typename Alphabet>
inline size_t SequenceComplexity::complexityLempelZiv(const AlphabetSequence<Alphabet>& sequence) {

  size_t u = 0;
  size_t v = 1;
  size_t w = 1;
  size_t v_max = 1;
  size_t length = sequence.length();
  size_t complexity = 1;


  while (true) {

    if (sequence.at(u + v - 1) == sequence.at(w + v - 1)) {

      v += 1;

      if (w + v >= length) {

        complexity += 1;
        break;

      }

    } else {

      if (v > v_max) {

        v_max = v;

      }

      u += 1;

      if (u == w) {

        complexity += 1;
        w += v_max;

        if (w > length) {

          break;

        } else {

          u = 0;
          v = 1;
          v_max = 1;

        }

      } else {

        v = 1;

      }

    }

  }

  return complexity;

}


}   // end namespace genome



#endif //KGL_SEQUENCE_COMPLEXITY_H
