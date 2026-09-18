//
// kgl_genetic_code.h — NCBI genetic code tables and the merged translation table/facade.
//
// The proven compact NCBI encodings are retained (diffable against the NCBI page) and decoded
// at compile time. The overloaded 'M'/'*'/'-' start field becomes a CodonRole enum, and the
// CONTAINS_BASE_N sentinel is deleted (N is checked up front).
//

#ifndef KGL_GENETIC_CODE_H
#define KGL_GENETIC_CODE_H


#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "kgl_genome_types.h"
#include "kgl_alphabet.h"
#include "kgl_sequence.h"
#include "kgl_sequence_dna.h"
#include "kgl_sequence_codon.h"


namespace kellerberrin::genome {   //  organization::project level namespace


enum class CodonRole : std::uint8_t { Normal, Start, Stop };

struct AminoTableColumn { Amino amino; CodonRole role; };

struct TranslationTable {

  std::string_view table_name;
  std::string_view table_description;
  std::span<const AminoTableColumn, 64> amino_table;

};


namespace Tables {

  inline constexpr std::size_t AMINO_TABLE_SIZE = 64;

  // Decode the compact 64x"b1 b2 b3 amino start" encoding at compile time; a malformed
  // encoding (or an unexpected 5th character) is a compile-time error.
  consteval AminoTableColumn* decodeTable(std::string_view encoded_table, AminoTableColumn* decoded) {

    if (encoded_table.size() != (AMINO_TABLE_SIZE * 6) - 1) {
      throw std::logic_error("Malformed amino table encoding");
    }

    for (std::size_t index = 0; index < AMINO_TABLE_SIZE; ++index) {

      const std::size_t group_offset = index * 6;
      const char role_char = encoded_table[group_offset + 4];
      CodonRole role = CodonRole::Normal;
      if (role_char == 'M') {
        role = CodonRole::Start;
      } else if (role_char == '*') {
        role = CodonRole::Stop;
      } else if (role_char != '-') {
        throw std::logic_error("Malformed amino table role character");
      }

      decoded[index] = AminoTableColumn{ .amino = aminoFromChar(encoded_table[group_offset + 3]),
                                         .role = role };

    }

    return decoded;

  }

  // The compact NCBI encodings (canonical order: index = 16*b1 + 4*b2 + b3, A=0 C=1 G=2 T=3).
  inline constexpr std::string_view ENCODED_TABLE_1 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAI- ATCI- ATGMM ATTI- "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGLM CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGV- GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGA** TGCC- TGGW- TGTC- TTAL- TTCF- TTGLM TTTF-";

  inline constexpr std::string_view ENCODED_TABLE_2 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGA** AGCS- AGG** AGTS- ATAMM ATCIM ATGMM ATTIM "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGL- CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGVM GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTAL- TTCF- TTGL- TTTF-";

  inline constexpr std::string_view ENCODED_TABLE_3 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAMM ATCI- ATGMM ATTI- "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAT- CTCT- CTGT- CTTT- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGV- GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTAL- TTCF- TTGL- TTTF-";

  inline constexpr std::string_view ENCODED_TABLE_4 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAIM ATCIM ATGMM ATTIM "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGLM CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGVM GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTALM TTCF- TTGLM TTTF-";

  inline constexpr std::string_view ENCODED_TABLE_5 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAS- AGCS- AGGS- AGTS- ATAMM ATCIM ATGMM ATTIM "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGL- CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGVM GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTAL- TTCF- TTGLM TTTF-";

  inline constexpr std::string_view ENCODED_PF_TABLE =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAI- ATCI- ATGMM ATTI- "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGL- CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGV- GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGA** TGCC- TGGW- TGTC- TTAL- TTCF- TTGL- TTTF-";

  inline constexpr auto TABLE_1_COLUMNS = []() consteval {
    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    decodeTable(ENCODED_TABLE_1, columns.data());
    return columns;
  }();

  inline constexpr auto TABLE_2_COLUMNS = []() consteval {
    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    decodeTable(ENCODED_TABLE_2, columns.data());
    return columns;
  }();

  inline constexpr auto TABLE_3_COLUMNS = []() consteval {
    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    decodeTable(ENCODED_TABLE_3, columns.data());
    return columns;
  }();

  inline constexpr auto TABLE_4_COLUMNS = []() consteval {
    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    decodeTable(ENCODED_TABLE_4, columns.data());
    return columns;
  }();

  inline constexpr auto TABLE_5_COLUMNS = []() consteval {
    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    decodeTable(ENCODED_TABLE_5, columns.data());
    return columns;
  }();

  inline constexpr auto P_FALCIPARUM_COLUMNS = []() consteval {
    std::array<AminoTableColumn, AMINO_TABLE_SIZE> columns{};
    decodeTable(ENCODED_PF_TABLE, columns.data());
    return columns;
  }();

  // The decoded translation table descriptors.
  inline constexpr TranslationTable TABLE_1{ "NCBI_TABLE_1", "The Standard Amino Code", TABLE_1_COLUMNS };
  inline constexpr TranslationTable TABLE_2{ "NCBI_TABLE_2", "The Vertebrate Mitochondrial Code", TABLE_2_COLUMNS };
  inline constexpr TranslationTable TABLE_3{ "NCBI_TABLE_3", "The Yeast Mitochondrial Code", TABLE_3_COLUMNS };
  inline constexpr TranslationTable TABLE_4{ "NCBI_TABLE_4", "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code", TABLE_4_COLUMNS };
  inline constexpr TranslationTable TABLE_5{ "NCBI_TABLE_5", "The Invertebrate Mitochondrial Code", TABLE_5_COLUMNS };
  inline constexpr TranslationTable P_FALCIPARUM{ "P_FALCIPARUM", "Organism Plasmodium Falciparum codon table downloaded from http://plasmodb.org", P_FALCIPARUM_COLUMNS };

  // Spot checks pin the transcription. Codon index = 16*b1 + 4*b2 + b3.
  static_assert(TABLE_1_COLUMNS[14].amino == Amino::M and TABLE_1_COLUMNS[14].role == CodonRole::Start, "ATG must code M start");
  static_assert(TABLE_1_COLUMNS[30].amino == Amino::L and TABLE_1_COLUMNS[30].role == CodonRole::Start, "Table 1: CTG must code L start");
  static_assert(TABLE_1_COLUMNS[48].role == CodonRole::Stop and TABLE_1_COLUMNS[50].role == CodonRole::Stop and TABLE_1_COLUMNS[56].role == CodonRole::Stop, "TAA/TAG/TGA must be stops");
  static_assert(TABLE_1_COLUMNS[62].amino == Amino::L and TABLE_1_COLUMNS[62].role == CodonRole::Start, "Table 1: TTG must code L start");
  static_assert(TABLE_2_COLUMNS[8].role == CodonRole::Stop and TABLE_2_COLUMNS[10].role == CodonRole::Stop, "Table 2: AGA/AGG must be stops");
  static_assert(TABLE_2_COLUMNS[12].role == CodonRole::Start, "Table 2: ATA must be a start");
  static_assert(TABLE_2_COLUMNS[56].amino == Amino::W and TABLE_2_COLUMNS[56].role == CodonRole::Normal, "Table 2: TGA must code W (not a stop)");
  static_assert(TABLE_3_COLUMNS[31].amino == Amino::T, "Table 3: CTT must code T");
  static_assert(TABLE_4_COLUMNS[12].amino == Amino::I and TABLE_4_COLUMNS[12].role == CodonRole::Start, "Table 4: ATA must code I start");
  static_assert(TABLE_4_COLUMNS[56].amino == Amino::W and TABLE_4_COLUMNS[56].role == CodonRole::Normal, "Table 4: TGA must code W (not a stop) - regression pin");
  static_assert(TABLE_4_COLUMNS[62].role == CodonRole::Start, "Table 4: TTG must be a start");
  static_assert(TABLE_5_COLUMNS[12].amino == Amino::M and TABLE_5_COLUMNS[12].role == CodonRole::Start, "Table 5: ATA must code M start");
  static_assert(TABLE_5_COLUMNS[46].role == CodonRole::Start, "Table 5: GTG must be a start");
  static_assert(P_FALCIPARUM_COLUMNS[14].amino == Amino::M and P_FALCIPARUM_COLUMNS[14].role == CodonRole::Start, "P. falciparum: ATG must code M start");
  static_assert(P_FALCIPARUM_COLUMNS[30].amino == Amino::L and P_FALCIPARUM_COLUMNS[30].role == CodonRole::Normal, "P. falciparum: CTG codes L but is not a start");
  static_assert(P_FALCIPARUM_COLUMNS[56].amino == Amino::Stop, "P. falciparum: TGA must be a stop");

  inline constexpr const TranslationTable* STANDARDTABLE = &TABLE_1;

  // The registry. TABLE_5 replaces the duplicate TABLE_2 entry (a copy-paste slip, B1).
  inline constexpr std::array TABLEARRAY{ STANDARDTABLE,
                                          &TABLE_2,
                                          &TABLE_3,
                                          &TABLE_4,
                                          &TABLE_5,
                                          &P_FALCIPARUM };

  inline constexpr std::size_t TABLEARRAYSIZE = std::size(TABLEARRAY);

  // Per-table start-amino mask, derived consteval in the same pass
  // isStartAmino becomes O(1) and the mask is
  // pinned by static_asserts below. Indexed by (amino_char - 'A') over [0, 26).
  template<const TranslationTable* TablePtr>
  inline constexpr auto START_AMINO_MASK = []() consteval {

    std::array<bool, 26> mask{};
    for (const auto& column : TablePtr->amino_table) {

      if (column.role == CodonRole::Start) {

        const char amino_char = AminoAcid::convertToChar(column.amino);
        if (amino_char >= 'A' and amino_char <= 'Z') {
          mask[static_cast<std::size_t>(amino_char - 'A')] = true;
        }

      }

    }
    return mask;

  }();

  // The per-registry-index mask array. AminoTranslationTable stores its registry index, so
  // isStartAmino() is a single indexed load — no pointer scan and no switch (v4 improvement).
  using StartAminoMask = std::array<bool, 26>;
  inline constexpr std::array<StartAminoMask, TABLEARRAYSIZE> START_MASKS{ START_AMINO_MASK<&TABLE_1>,
                                                                           START_AMINO_MASK<&TABLE_2>,
                                                                           START_AMINO_MASK<&TABLE_3>,
                                                                           START_AMINO_MASK<&TABLE_4>,
                                                                           START_AMINO_MASK<&TABLE_5>,
                                                                           START_AMINO_MASK<&P_FALCIPARUM> };

  // The mask for a valid (non-corrupted) index; corrupted indices have no start aminos.
  [[nodiscard]] constexpr bool startAmino(std::size_t mask_index, char amino_char) noexcept {
    return mask_index < TABLEARRAYSIZE
           and START_MASKS[mask_index][static_cast<std::size_t>(amino_char - 'A')];
  }

  static_assert(START_AMINO_MASK<&TABLE_1>['M' - 'A'], "Table 1: M must be a start amino");
  static_assert(START_AMINO_MASK<&TABLE_1>['L' - 'A'], "Table 1: L must be a start amino (TTG/CTG)");
  static_assert(not START_AMINO_MASK<&TABLE_2>['L' - 'A'], "Table 2: L must not be a start amino (TTG/CTG are not starts)");
  static_assert(START_AMINO_MASK<&TABLE_2>['I' - 'A'], "Table 2: I must be a start amino (ATA)");
  static_assert(START_AMINO_MASK<&TABLE_4>['I' - 'A'], "Table 4: I must be a start amino (ATA)");

}   // namespace Tables


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AminoTranslationTable — merges the reference's AminoTranslationTable and TranslateToAmino.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AminoTranslationTable {

public:

  AminoTranslationTable() = default;
  ~AminoTranslationTable() = default;

  /// ReturnType the translation table name.
  [[nodiscard]] std::string translationTableName() const { return std::string(tableNameView()); }
  [[nodiscard]] std::string tableName() const { return std::string(tableNameView()); }

  /// Set the translation table by name (case-insensitive; standard on miss).
  [[nodiscard]] bool setTranslationTable(std::string_view table_name) {

    std::string upper_name(table_name);
    std::ranges::transform(upper_name, upper_name.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    const auto table_iter = std::ranges::find_if(Tables::TABLEARRAY,
                                                 [&upper_name](const TranslationTable* table) { return table->table_name == upper_name; });
    if (table_iter == Tables::TABLEARRAY.end()) {

      ExecEnv::log().warn("Amino translation table: {} not found, Standard table 1 used", table_name);
      table_ = Tables::STANDARDTABLE;
      index_ = 0;
      return false;

    }

    table_ = *table_iter;
    index_ = static_cast<std::size_t>(std::distance(Tables::TABLEARRAY.begin(), table_iter));
    return true;

  }

  /// Deprecated misspelling retained for downstream compatibility.
  [[nodiscard]] bool settranslationTable(const std::string& table_name) { return setTranslationTable(table_name); }

  /// Returns an amino acid for the codon. Amino::Z if any base is 'N'.
  [[nodiscard]] AminoAcid::Alphabet getAmino(const Codon& codon) const noexcept {

    if (codon.containsBaseN()) {
      return AminoAcid::AMINO_UNKNOWN;
    }
    return table_->amino_table[index(codon)].amino;

  }

  /// Returns true if the amino acid is a stop amino acid.
  [[nodiscard]] bool isStopAmino(AminoAcid::Alphabet amino) const noexcept {
    return amino == AminoAcid::AMINO_STOP;
  }

  /// Returns true if any codon in the table that codes the amino acid is a start codon.
  /// O(1): the stored registry index selects the consteval mask directly (v4: no pointer scan).
  [[nodiscard]] bool isStartAmino(AminoAcid::Alphabet amino) const noexcept {

    const char amino_char = AminoAcid::convertToChar(amino);
    if (amino_char < 'A' or amino_char > 'Z') {
      return false;
    }
    return Tables::startAmino(index_, amino_char);

  }

  /// Returns true if the codon is a stop codon.
  [[nodiscard]] bool isStopCodon(const Codon& codon) const noexcept {

    if (codon.containsBaseN()) {
      return false;
    }
    return table_->amino_table[index(codon)].role == CodonRole::Stop;

  }

  /// Returns true if the codon is a start codon.
  [[nodiscard]] bool isStartCodon(const Codon& codon) const noexcept {

    if (codon.containsBaseN()) {
      return false;
    }
    return table_->amino_table[index(codon)].role == CodonRole::Start;

  }

  /// ReturnType the amino acid sequence for the coding sequence.
  [[nodiscard]] AminoSequence getAminoSequence(const DNA5SequenceCoding& coding_sequence) const;

  [[nodiscard]] bool checkStartCodon(const AminoSequence& amino_sequence) const;
  [[nodiscard]] bool checkStopCodon(const AminoSequence& amino_sequence) const;
  [[nodiscard]] std::pair<std::size_t, bool> firstStopSequenceSize(const AminoSequence& amino_sequence) const;

private:

  const TranslationTable* table_{ Tables::STANDARDTABLE };
  std::size_t index_{0};   // The registry index of `table_` (names are unique; O(1) mask dispatch).

  /// The table name as a view (avoids a temporary std::string in the accessors).
  [[nodiscard]] std::string_view tableNameView() const noexcept { return table_->table_name; }

  // The table index for an N-free codon. Valid (non-N) bases have columns <= 3, so the index
  // is statically proven <= 63; the shared helper keeps the "N-check then index" invariant in
  // one place (improvement: deepseek repeated this 16/4/1 formula in three members).
  [[nodiscard]] static std::size_t index(const Codon& codon) noexcept {

    return 16 * DNA5::symbolToColumn(codon[0])
         + 4 * DNA5::symbolToColumn(codon[1])
         + DNA5::symbolToColumn(codon[2]);

  }

};


using TranslateToAmino = AminoTranslationTable;


}   // end namespace


#endif //KGL_GENETIC_CODE_H
