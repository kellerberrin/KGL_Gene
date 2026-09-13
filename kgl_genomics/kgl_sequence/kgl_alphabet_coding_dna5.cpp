//
// Created by kellerberrin on 25/11/17.
//


#include "kgl_alphabet_coding_dna5.h"
#include "kgl_alphabet_dna5_tables.h"


namespace kgl = kellerberrin::genome;


bool kgl::CodingDNA5::validAlphabet(Alphabet nucleotide) {

  // We DO NOT use a switch here.
  // Because the switch always assumes we can only have 5 base types (and a memory corrupted sequence may not).
  return detail::NucleotideTables<CodingDNA5>::VALID_TABLE[static_cast<unsigned char>(nucleotide)];

}


// Convert a base to an array offset.
kgl::ContigOffset_t kgl::CodingDNA5::symbolToColumn(Alphabet nucleotide) {

  // The table lookup is total so a corrupted raw value cannot index out of bounds.
  return detail::NucleotideTables<CodingDNA5>::COLUMN_TABLE[static_cast<unsigned char>(nucleotide)];

}


const std::vector<kgl::CodingDNA5::Alphabet>& kgl::CodingDNA5::enumerateAlphabet() {

  return detail::NucleotideTables<CodingDNA5>::ALPHABET_VECTOR;

}