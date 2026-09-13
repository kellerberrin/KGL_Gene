//
// Created by kellerberrin on 10/01/18.
//

#ifndef KGL_TABLE_REGISTRY_H
#define KGL_TABLE_REGISTRY_H


#include "kgl_table_impl.h"
#include "kgl_table_ncbi.h"
#include "kgl_table_organism.h"

#include <array>
#include <string_view>
#include <stdexcept>


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compile-time decoding of the compact NCBI table encodings (T2) and the table registry (T1).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace Tables {

  // Decodes a compact table encoding into an AminoTableColumn array at compile time.
  // The encoding is 64 groups of 5 characters "base1 base2 base3 amino start" separated by single spaces.
  // A malformed encoding is a compile-time error (the throw below rejects the constant evaluation).
  consteval AminoTableColumn* DecodeTable(std::string_view encoded_table, AminoTableColumn* decoded) {

    if (encoded_table.size() != (AMINO_TABLE_SIZE * 6) - 1) {

      throw std::logic_error("Malformed amino table encoding");

    }

    for (size_t index = 0; index < static_cast<size_t>(AMINO_TABLE_SIZE); ++index) {

      const size_t group_offset = index * 6;   // 5 data characters plus a space separator.
      decoded[index] = AminoTableColumn{ .amino_acid = encoded_table[group_offset + 3],
                                         .start = encoded_table[group_offset + 4],
                                         .base1 = encoded_table[group_offset],
                                         .base2 = encoded_table[group_offset + 1],
                                         .base3 = encoded_table[group_offset + 2] };

    }

    return decoded;

  }

  // The decoded NCBI amino tables.
  inline constexpr auto STANDARD_TABLE = []() consteval {

    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    DecodeTable(NCBITables::ENCODED_TABLE_1, columns.data());
    return columns;

  }();

  inline constexpr auto TABLE_2_COLUMNS = []() consteval {

    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    DecodeTable(NCBITables::ENCODED_TABLE_2, columns.data());
    return columns;

  }();

  inline constexpr auto TABLE_3_COLUMNS = []() consteval {

    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    DecodeTable(NCBITables::ENCODED_TABLE_3, columns.data());
    return columns;

  }();

  inline constexpr auto TABLE_4_COLUMNS = []() consteval {

    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    DecodeTable(NCBITables::ENCODED_TABLE_4, columns.data());
    return columns;

  }();

  inline constexpr auto TABLE_5_COLUMNS = []() consteval {

    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    DecodeTable(NCBITables::ENCODED_TABLE_5, columns.data());
    return columns;

  }();

  inline constexpr auto P_FALCIPARUM_COLUMNS = []() consteval {

    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    DecodeTable(OrganismTables::ENCODED_PF_TABLE, columns.data());
    return columns;

  }();

  // The decoded translation table descriptors.
  inline constexpr TranslationTable TABLE_1{ STANDARD_TABLE.data(), NCBITables::AMINO_TABLE_NAME.data(), NCBITables::AMINO_TABLE_DESC.data() };
  inline constexpr TranslationTable TABLE_2{ TABLE_2_COLUMNS.data(), NCBITables::AMINO_TABLE_2_NAME.data(), NCBITables::AMINO_TABLE_2_DESC.data() };
  inline constexpr TranslationTable TABLE_3{ TABLE_3_COLUMNS.data(), NCBITables::AMINO_TABLE_3_NAME.data(), NCBITables::AMINO_TABLE_3_DESC.data() };
  inline constexpr TranslationTable TABLE_4{ TABLE_4_COLUMNS.data(), NCBITables::AMINO_TABLE_4_NAME.data(), NCBITables::AMINO_TABLE_4_DESC.data() };
  inline constexpr TranslationTable TABLE_5{ TABLE_5_COLUMNS.data(), NCBITables::AMINO_TABLE_5_NAME.data(), NCBITables::AMINO_TABLE_5_DESC.data() };
  inline constexpr TranslationTable P_FALCIPARUM{ P_FALCIPARUM_COLUMNS.data(), OrganismTables::PF_TABLE_NAME.data(), OrganismTables::PF_TABLE_DESC.data() };

  // Transcription spot checks (the T2 verification protocol; mandatory before deleting the braced lists).
  static_assert(STANDARD_TABLE[14].amino_acid == 'M' and STANDARD_TABLE[14].start == 'M', "Standard table: ATG must code M start");
  static_assert(STANDARD_TABLE[48].start == '*' and STANDARD_TABLE[50].start == '*' and STANDARD_TABLE[56].start == '*', "Standard table: TAA/TAG/TGA must be stops");
  static_assert(TABLE_2_COLUMNS[8].start == '*' and TABLE_2_COLUMNS[10].start == '*', "Table 2: AGA/AGG must be stops");
  static_assert(TABLE_2_COLUMNS[48].start == '*' and TABLE_2_COLUMNS[50].start == '*', "Table 2: TAA/TAG must be stops");
  static_assert(TABLE_2_COLUMNS[12].start == 'M', "Table 2: ATA must be a start (M)");
  static_assert(TABLE_3_COLUMNS[31].amino_acid == 'T', "Table 3: CTT must code T (yeast)");
  static_assert(TABLE_4_COLUMNS[62].start == 'M', "Table 4: TTG must be a start (M)");
  static_assert(TABLE_5_COLUMNS[46].start == 'M', "Table 5: GTG must be a start (M)");
  static_assert(P_FALCIPARUM_COLUMNS[14].amino_acid == 'M' and P_FALCIPARUM_COLUMNS[14].start == 'M', "P. falciparum: ATG must code M start");

  // NCBI_TABLE_1 is the default translation table.
  inline constexpr const TranslationTable* STANDARDTABLE = &TABLE_1;

  // The registry. TABLE_5 replaces the duplicate TABLE_2 entry (a copy-paste slip, B1).
  inline constexpr std::array TABLEARRAY{ STANDARDTABLE,
                                          &TABLE_2,
                                          &TABLE_3,
                                          &TABLE_4,
                                          &TABLE_5,
                                          &P_FALCIPARUM };

  inline constexpr size_t TABLEARRAYSIZE = std::size(TABLEARRAY);

}   // namespace Tables


}   // end namespace


#endif //KGL_TABLE_REGISTRY_H