//
// Created by kellerberrin on 26/10/18.
//

#include "kgl_sequence_complexity.h"
#include <map>


namespace kgl = kellerberrin::genome;


double kgl::SequenceComplexity::sequenceEntropy(std::shared_ptr<const DNA5SequenceLinear> sequence, size_t word_size) {

  std::map<std::string, size_t> word_map;

  std::string sequence_string = sequence->getSequenceAsString();

  for (size_t i = 0; i < sequence->length(); i++) {

    if (i + word_size >= sequence->length()) break;

    std::string word = sequence_string.substr(i, word_size);

    auto result = word_map.find(word);

    if (result == word_map.end()) {

      std::pair<std::string, size_t> insert_word(word, 1);

      auto insert_result = word_map.insert(insert_word);

      if (not insert_result.second) {

        ExecEnv::log().warn("sequenceEntropy; Unable to insert word: {}", word);

      }

    } else {


      ++result->second;

    }

  }

  double entropy = 0;

  for (auto word : word_map) {

    double word_ratio = static_cast<double>(word.second) / static_cast<double>(sequence->length());

    double log_word_ratio = std::log(word_ratio) / word_size;

    entropy += word_ratio * log_word_ratio * -1.0;

  }

  return entropy;

}



double kgl::SequenceComplexity::propGC(std::shared_ptr<const DNA5SequenceLinear> sequence) {

  size_t count = 0;
  for(ContigSize_t index = 0; index < sequence->length(); ++index) {

    if (sequence->at(index) == DNA5::Alphabet::C or sequence->at(index) == DNA5::Alphabet::G) {

      ++count;

    }

  }

  if (sequence->length() > 0) {

    return static_cast<double>(count) / static_cast<double>(sequence->length());

  } else {

    ExecEnv::log().warn("SequenceComplexity::propGC; analysis of zero sized sequence");
    return 0.0;

  }

}





size_t kgl::SequenceComplexity::complexityLempelZiv(std::shared_ptr<const DNA5SequenceLinear> sequence) {

  size_t u = 0;
  size_t v = 1;
  size_t w = 1;
  size_t v_max = 1;
  size_t length = sequence->length();
  size_t complexity = 1;


  while (true) {

    if (sequence->at(u + v - 1) == sequence->at(w + v - 1)) {

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

