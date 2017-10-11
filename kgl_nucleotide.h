// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
// Created by kellerberrin on 20/09/17.
//

#ifndef KGL_NUCLEOTIDE_H
#define KGL_NUCLEOTIDE_H

#include <cstdint>
#include <memory>
#include <string>
#include <queue>
#include "kgl_logging.h"
#include "kgl_genome_types.h"
#include "kgl_exec_env.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Standard contig nucleotide column layout.
// Implement the NUCLEOTIDE_COLUMNS as "A", "C", "G", "T"/"U", "-", "+" in that order.

class StandardNucleotideColumn {

public:

  explicit StandardNucleotideColumn() = default;
  ~StandardNucleotideColumn() = default;

  static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 7;
  static constexpr Nucleotide_t DELETE_NUCLEOTIDE = '-';
  static constexpr Nucleotide_t INSERT_SEQUENCE = '+';

  ContigOffset_t nucleotideToColumn(const Nucleotide_t nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case 'A':
      case 'a': return 0;

      case 'C':
      case 'c': return 1;

      case 'G':
      case 'g': return 2;

      case 'U':
      case 'u':
      case 'T':
      case 't': return 3;

      case 'N':
      case 'n': return 4;

      case DELETE_NUCLEOTIDE: return 5;

      case INSERT_SEQUENCE: return 6;

      default:
        ExecEnv::log().critical("nucleotideToColumn(), Count data array accessed with unknown nucleotide: {}",
                                nucleotide);
        return 0; // Never reached, to keep the compiler happy.

    }

  }

private:

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_NUCLEOTIDE_H
