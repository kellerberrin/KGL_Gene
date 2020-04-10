//
// Created by kellerberrin on 26/10/18.
//

#include "kgl_sequence_complexity.h"
#include <map>


namespace kgl = kellerberrin::genome;


void kgl::SequenceComplexity::countNucleotides(const DNA5SequenceLinear& sequence,
                                               size_t& A_count,
                                               size_t& C_count,
                                               size_t& G_count,
                                               size_t& T_count) {

  A_count = 0;
  C_count = 0;
  G_count = 0;
  T_count = 0;

  for(ContigSize_t index = 0; index < sequence.length(); ++index) {

    switch(sequence.at(index)) {

      case DNA5::Alphabet::A:
        ++A_count;
        break;

      case DNA5::Alphabet::C:
        ++C_count;
        break;

      case DNA5::Alphabet::G:
        ++G_count;
        break;

      case DNA5::Alphabet::T:
        ++T_count;
        break;

      case DNA5::Alphabet::N:
        break;

      default:
      ExecEnv::log().error("SequenceComplexity::countNucleotides; bad nucleotide in sequence: {}", sequence.getSequenceAsString());
      break;

    }

  }

}


double  kgl::SequenceComplexity::relativeCpGIslands(const DNA5SequenceLinear& sequence) {

  if (sequence.length() == 0) {

    ExecEnv::log().error("SequenceComplexity::relativeCpGIslands; zero sized sequence");
    return 0.0;

  }

  size_t A_count = 0;
  size_t C_count = 0;
  size_t G_count = 0;
  size_t T_count = 0;

  countNucleotides(sequence, A_count, C_count, G_count,  T_count);

  if (C_count + G_count == 0) {

    return 0.0;

  }

  const DNA5SequenceLinear CpG(StringDNA5("CG"));

  size_t CpG_count = kmerCount<DNA5>(sequence, CpG);

  double expected_CpG = static_cast<double>(sequence.length()) / 32.0;

  return static_cast<double>(CpG_count) / expected_CpG;

}



void kgl::SequenceComplexity::proportionNucleotides(const DNA5SequenceLinear& sequence,
                                                    double& A_prop,
                                                    double& C_prop,
                                                    double& G_prop,
                                                    double& T_prop) {


  A_prop = 0.0;
  C_prop = 0.0;
  G_prop = 0.0;
  T_prop = 0.0;

  if (sequence.length() == 0) {

    ExecEnv::log().error("SequenceComplexity::proportionNucleotides; zero sized sequence");
    return;

  }

  size_t A_count = 0;
  size_t C_count = 0;
  size_t G_count = 0;
  size_t T_count = 0;

  countNucleotides(sequence, A_count, C_count, G_count,  T_count);

  A_prop = static_cast<double>(A_count) / static_cast<double>(sequence.length());
  C_prop = static_cast<double>(C_count) / static_cast<double>(sequence.length());
  G_prop = static_cast<double>(G_count) / static_cast<double>(sequence.length());
  T_prop = static_cast<double>(T_count) / static_cast<double>(sequence.length());

}



size_t kgl::SequenceComplexity::complexityLempelZiv(const DNA5SequenceLinear& sequence) {

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

