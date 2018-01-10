//
// Created by kellerberrin on 10/01/18.
//

#include "kgl_table_impl.h"
#include "kgl_table_ncbi.h"
#include "kgl_table_organism.h"

namespace kgl = kellerberrin::genome;


// NBCI tables from the NCBI website: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

constexpr const kgl::AminoTableColumn kgl::NCBITables::StandardTranslationTable[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::NCBITables::TABLE_1;

constexpr const kgl::AminoTableColumn kgl::NCBITables::TranslationTable_2[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::NCBITables::TABLE_2;

constexpr const kgl::AminoTableColumn kgl::NCBITables::TranslationTable_3[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::NCBITables::TABLE_3;

constexpr const kgl::AminoTableColumn kgl::NCBITables::TranslationTable_4[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::NCBITables::TABLE_4;

constexpr const kgl::AminoTableColumn kgl::NCBITables::TranslationTable_5[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::NCBITables::TABLE_5;

// Organism specific tables.

constexpr const kgl::AminoTableColumn kgl::OrganismTables::PFTranslationTable[Tables::AMINO_TABLE_SIZE];

constexpr const kgl::TranslationTable kgl::OrganismTables::P_FALCIPARUM;

// Public data items.
// NCBI_TABLE_1 is the default translation table.
const kgl::TranslationTable* kgl::Tables::STANDARDTABLE = &kgl::NCBITables::TABLE_1;

const kgl::TranslationTable* kgl::Tables::TABLEARRAY[] = { kgl::Tables::STANDARDTABLE,
                                                           &kgl::NCBITables::TABLE_2,
                                                           &kgl::NCBITables::TABLE_3,
                                                           &kgl::NCBITables::TABLE_4,
                                                           &kgl::NCBITables::TABLE_2,
                                                           &kgl::OrganismTables::P_FALCIPARUM };

const size_t kgl::Tables::TABLEARRAYSIZE = sizeof(kgl::Tables::TABLEARRAY) / sizeof(kgl::TranslationTable*);