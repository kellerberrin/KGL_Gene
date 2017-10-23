//
// Created by kellerberrin on 22/10/17.
//

#include "kgl_table.h"

namespace kgl = kellerberrin::genome;

// Link definitions for translation tables.
constexpr const kgl::AminoTableColumn kgl::Tables::StandardTranslationTable[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::Tables::TABLE_1;

constexpr const kgl::AminoTableColumn kgl::Tables::TranslationTable_2[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::Tables::TABLE_2;

constexpr const kgl::AminoTableColumn kgl::Tables::TranslationTable_3[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::Tables::TABLE_3;

constexpr const kgl::AminoTableColumn kgl::Tables::TranslationTable_4[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::Tables::TABLE_4;

constexpr const kgl::AminoTableColumn kgl::Tables::TranslationTable_5[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::Tables::TABLE_5;

// Link definitions for translation table array.

constexpr const kgl::TranslationTable* kgl::Tables::TableArray[];
