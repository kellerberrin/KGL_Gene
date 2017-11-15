//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_alphabet_base.h"

namespace kgl = kellerberrin::genome;

// Instantiate the NucleotideColumn_DNA5 constants for the edification of the linker.


const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::A_NUCLEOTIDE;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::A_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::BaseDNA5::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::C_NUCLEOTIDE;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::C_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::BaseDNA5::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::G_NUCLEOTIDE;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::G_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::BaseDNA5::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::U_NUCLEOTIDE;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::U_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::BaseDNA5::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::T_NUCLEOTIDE;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::T_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::BaseDNA5::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::N_NUCLEOTIDE;
const kgl::Nucleotide_ExtendedDNA5 kgl::BaseDNA5::N_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::BaseDNA5::N_NUCLEOTIDE_OFFSET;


// Convert char to Alphabet enum type.
kgl::BaseDNA5::Alphabet kgl::BaseDNA5::convertChar(char chr_base) {

  switch (std::toupper(chr_base)) {

    case 'A':return Alphabet::A;

    case 'C': return Alphabet::C;

    case 'G': return Alphabet::G;

    case 'U':
    case 'T': return Alphabet::T;

    case 'N': return Alphabet::N;

    default:
      ExecEnv::log().error("BaseDNA5::convertchar(), Invalid nucleotide: {}", chr_base);
      return Alphabet::N;

  }

}


const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::NUCLEOTIDE_COLUMNS;
const kgl::Nucleotide_ExtendedDNA5 kgl::NucleotideColumn_DNA5::DELETE_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::DELETE_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::NucleotideColumn_DNA5::INSERT_A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::INSERT_A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::NucleotideColumn_DNA5::INSERT_C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::INSERT_C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::NucleotideColumn_DNA5::INSERT_G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::INSERT_G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::NucleotideColumn_DNA5::INSERT_U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::INSERT_U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::NucleotideColumn_DNA5::INSERT_T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::INSERT_T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_ExtendedDNA5 kgl::NucleotideColumn_DNA5::INSERT_N_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::NucleotideColumn_DNA5::INSERT_N_NUCLEOTIDE_OFFSET;

