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

class NucleotideColumn_DNA5 {

public:

  using NucleotideType = Nucleotide_DNA5_t;

  explicit NucleotideColumn_DNA5() = default;
  ~NucleotideColumn_DNA5() = default;

  static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 7;
  static constexpr Nucleotide_DNA5_t A_NUCLEOTIDE = 'A';
  static constexpr ContigOffset_t A_NUCLEOTIDE_OFFSET = 0;
  static constexpr Nucleotide_DNA5_t C_NUCLEOTIDE = 'C';
  static constexpr ContigOffset_t C_NUCLEOTIDE_OFFSET = 1;
  static constexpr Nucleotide_DNA5_t G_NUCLEOTIDE = 'G';
  static constexpr ContigOffset_t G_NUCLEOTIDE_OFFSET = 2;
  static constexpr Nucleotide_DNA5_t T_NUCLEOTIDE = 'T';
  static constexpr ContigOffset_t T_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_DNA5_t N_NUCLEOTIDE = 'N';
  static constexpr ContigOffset_t N_NUCLEOTIDE_OFFSET = 4;
  static constexpr Nucleotide_DNA5_t DELETE_NUCLEOTIDE = '-';
  static constexpr ContigOffset_t DELETE_NUCLEOTIDE_OFFSET = 5;
  static constexpr Nucleotide_DNA5_t INSERT_SEQUENCE = '+';
  static constexpr ContigOffset_t INSERT_NUCLEOTIDE_OFFSET = 6;

  static ContigOffset_t nucleotideToColumn(const Nucleotide_DNA5_t nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case A_NUCLEOTIDE:
      case 'a': return A_NUCLEOTIDE_OFFSET;

      case C_NUCLEOTIDE:
      case 'c': return C_NUCLEOTIDE_OFFSET;

      case G_NUCLEOTIDE:
      case 'g': return G_NUCLEOTIDE_OFFSET;

      case T_NUCLEOTIDE:
      case 't':
      case 'U':
      case 'u': return T_NUCLEOTIDE_OFFSET;

      case N_NUCLEOTIDE:
      case 'n': return N_NUCLEOTIDE_OFFSET;

      case DELETE_NUCLEOTIDE: return DELETE_NUCLEOTIDE_OFFSET;

      case INSERT_SEQUENCE: return INSERT_NUCLEOTIDE_OFFSET;

      default:
        ExecEnv::log().critical("nucleotideToColumn(), Nucleotide array accessed with unknown nucleotide: {}",
                                nucleotide);
        return 0; // Never reached, to keep the compiler happy.

    }

  }

  static Nucleotide_DNA5_t offsetToNucleotide(ContigOffset_t offset) {

    // Translate the nucleotide to an array column
    switch (offset) {

      case A_NUCLEOTIDE_OFFSET: return A_NUCLEOTIDE;

      case C_NUCLEOTIDE_OFFSET: return C_NUCLEOTIDE;

      case G_NUCLEOTIDE_OFFSET: return G_NUCLEOTIDE;

      case T_NUCLEOTIDE_OFFSET: return T_NUCLEOTIDE;

      case N_NUCLEOTIDE_OFFSET: return N_NUCLEOTIDE;

      case DELETE_NUCLEOTIDE_OFFSET: return DELETE_NUCLEOTIDE;

      case INSERT_NUCLEOTIDE_OFFSET: return INSERT_SEQUENCE;

      default:
        ExecEnv::log().critical("indexToNucleotide(), unknown nucleotide: {}", offset);
        return N_NUCLEOTIDE; // Never reached, to keep the compiler happy.

    }

  }


private:

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_NUCLEOTIDE_H
