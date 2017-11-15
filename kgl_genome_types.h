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
using Amino_t = char;
using CharSequence_t = std::basic_string<Nucleotide_t>;
using ContigId_t = std::string;
using GenomeId_t = std::string;
using ContigFeatureId_t = std::string;
using ContigOffset_t = uint64_t;              // Paris Japonica has 150 billion base pairs, use 64 bit integer.
using ContigSize_t = ContigOffset_t;
using NucleotideReadCount_t = uint32_t;
using CDSPhaseType_t = unsigned char;
using FeatureIdent_t = std::string;
using FeatureType_t = std::string;
using VariantType_t = std::string;


// Usage and semantics : ContigOffset_t should be used when referring to the Genome and std::size_t should be used
// when referring to the underlying data structure. In reality both are 64 bit integers (Paris Japonica). This is
// asserted below.

static_assert( sizeof(ContigSize_t) == sizeof(std::size_t)
             , "kgl::ContigSize_t and std::size_t should both be 64bit integers");

}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_GENOME_TYPES_H
