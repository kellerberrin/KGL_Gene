//
// Created by kellerberrin on 10/01/18.
//

#ifndef KGL_TABLE_ORGANISM_H
#define KGL_TABLE_ORGANISM_H


#include <array>
#include <string_view>
#include "kgl_genome_types.h"
#include "kgl_table_impl.h"


namespace kellerberrin::genome {   //  organization level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines Organism Specific DNA/RNA to Amino Acid translation tables.
// Plasmodium Falciparum codon table downloaded from http://plasmodb.org (same as the standard table)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace OrganismTables {

  inline constexpr std::string_view PF_TABLE_NAME = "P_FALCIPARUM";
  inline constexpr std::string_view PF_TABLE_DESC = "Organism Plasmodium Falciparum codon table downloaded from http://plasmodb.org";

  // The compact encoding (see kgl_table_ncbi.h for the format).
  // Note - CTG and TTG code L but are NOT start codons in this table (unlike the standard table).
  inline constexpr std::string_view ENCODED_PF_TABLE =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAI- ATCI- ATGMM ATTI- "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGL- CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGV- GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGA** TGCC- TGGW- TGTC- TTAL- TTCF- TTGL- TTTF-";

}   // namespace OrganismTables


}   // end namespace


#endif //KGL_TABLE_ORGANISM_H