//
// Created by kellerberrin on 11/09/17.
//

#ifndef KGL_GENOME_TYPES_H
#define KGL_GENOME_TYPES_H

#include <string>
#include <cstdint>
#include <cassert>


namespace kellerberrin::genome {   //  organization::project level namespace

// Types used to represent genome data.

using Alphabet_t = char;   // Storage type for Amino, DNA5 or ReadCountColumns alphabets
using Nucleotide_t = Alphabet_t;  // Semantic alias for nucleotides.
using Amino_t = Alphabet_t; // Semantic alias for Amino Acids.
using AlphabetSequence_t = std::basic_string<Alphabet_t>;  // Just a std::string for now.
using ContigId_t = std::string;
using GenomeId_t = std::string;
using PopulationId_t = std::string;
using HaplotypeId_t = std::string;
using GenotypeId_t = std::string;
using CompareScore_t = long;  // Signed integer for sequence comparison.
using CompareDistance_t = double;  // double for sequence parentDistance.
using ContigFeatureId_t = std::string;
using ContigOffset_t = uint64_t;              // Paris Japonica has 150 billion base pairs, use 64 bit integers.
using AlleleOffset_t = uint64_t;              // Difference between the reference offset and the start of the allele (1 for Indel, 0 for SNP).
using SignedOffset_t = int64_t;
using ContigSize_t = ContigOffset_t;
using NucleotideReadCount_t = uint32_t;
using CDSPhaseType_t = unsigned char;
using FeatureIdent_t = std::string;
using FeatureType_t = std::string;
using VariantType_t = std::string;
using OntologyIdent_t = std::string;
using Phred_t = double;


// Usage and semantics : ContigOffset_t should be used when referring to the Genome and std::size_t should be used
// when referring to the underlying data structure. In reality both are 64 bit integers (Paris Japonica). This is
// asserted below.

static_assert( sizeof(ContigSize_t) == sizeof(std::size_t)
             , "kgl::ContigSize_t and std::size_t should both be 64bit integers");



}   // end namespace

#endif // KGL_GENOME_TYPES_H
