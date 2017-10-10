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
// Created by kellerberrin on 11/09/17.
//

#ifndef KGL_GENOME_TYPES_H
#define KGL_GENOME_TYPES_H

#include <string>
#include <cstdint>
#include <cassert>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Types used to represent genome data.

using Nucleotide_t = char;
using Sequence_t = std::string;
using ContigId_t = std::string;
using ContigFeatureId_t = std::string;
using ContigOffset_t = uint64_t;              // Paris Japonica has 150 billion base pairs, use 64 bit integer.
using ContigSize_t = ContigOffset_t;
using NucleotideReadCount_t = uint32_t;
using CDSPhaseType_t = unsigned char;

// Usage and semantics : ContigOffset_t should be used when referring to the Genome and std::size_t should be used
// when referring to the underlying data structure. In reality both are 64 bit integers (Paris Japonica). This is
// asserted below.

static_assert( sizeof(ContigSize_t) == sizeof(std::size_t)
             , "kgl::ContigSize_t and std::size_t should both be 64bit integers");

}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_GENOME_TYPES_H
