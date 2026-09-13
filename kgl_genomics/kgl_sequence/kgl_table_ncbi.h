//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_TABLE_NCBI_H
#define KGL_TABLE_NCBI_H


#include <array>
#include <string_view>
#include "kgl_genome_types.h"
#include "kgl_table_impl.h"


namespace kellerberrin::genome {   //  organization level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines DNA/RNA to Amino Acid translation tables.
// These are found at the NCBI website: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// At the current time only 5 tables are implemented (there are 31) as it is a somewhat tedious task.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The NCBI tables, re-encoded as compact compile-time data (T2).
// Each table is a string of 64 groups of 5 characters: "base1 base2 base3 amino start".
// The rows are in canonical NCBI order: codon index = 16*base1 + 4*base2 + base3 with A=0, C=1, G=2, T=3.
// The decoding into AminoTableColumn arrays is performed at compile time (see DecodeTable below);
// the spot-check static_asserts below pin the transcription (ATG start, TAA/TAG/TGA stop, table 2 AGA/AGG stop).

namespace NCBITables {

  inline constexpr std::string_view AMINO_TABLE_NAME = "NCBI_TABLE_1";
  inline constexpr std::string_view AMINO_TABLE_DESC = "The Standard Amino Code";

  inline constexpr std::string_view ENCODED_TABLE_1 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAI- ATCI- ATGMM ATTI- "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGLM CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGV- GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGA** TGCC- TGGW- TGTC- TTAL- TTCF- TTGLM TTTF-";

  // The Vertebrate Mitochondrial Code (transl_table=2)
  inline constexpr std::string_view AMINO_TABLE_2_NAME = "NCBI_TABLE_2";
  inline constexpr std::string_view AMINO_TABLE_2_DESC = "The Vertebrate Mitochondrial Code";

  inline constexpr std::string_view ENCODED_TABLE_2 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGA** AGCS- AGG** AGTS- ATAMM ATCIM ATGMM ATTIM "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGL- CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGVM GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTAL- TTCF- TTGL- TTTF-";

  // The Yeast Mitochondrial Code (transl_table=3)
  inline constexpr std::string_view AMINO_TABLE_3_NAME = "NCBI_TABLE_3";
  inline constexpr std::string_view AMINO_TABLE_3_DESC = "The Yeast Mitochondrial Code";

  inline constexpr std::string_view ENCODED_TABLE_3 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAMM ATCI- ATGMM ATTI- "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAT- CTCT- CTGT- CTTT- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGV- GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTAL- TTCF- TTGL- TTTF-";

  // The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
  inline constexpr std::string_view AMINO_TABLE_4_NAME = "NCBI_TABLE_4";
  inline constexpr std::string_view AMINO_TABLE_4_DESC =
  "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code";

  inline constexpr std::string_view ENCODED_TABLE_4 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAR- AGCS- AGGR- AGTS- ATAIM ATCIM ATGMM ATTIM "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGLM CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGVM GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTALM TTCF- TTGLM TTTF-";

  // The Invertebrate Mitochondrial Code (transl_table=5)
  inline constexpr std::string_view AMINO_TABLE_5_NAME = "NCBI_TABLE_5";
  inline constexpr std::string_view AMINO_TABLE_5_DESC = "The Invertebrate Mitochondrial Code";

  inline constexpr std::string_view ENCODED_TABLE_5 =
    "AAAK- AACN- AAGK- AATN- ACAT- ACCT- ACGT- ACTT- "
    "AGAS- AGCS- AGGS- AGTS- ATAMM ATCIM ATGMM ATTIM "
    "CAAQ- CACH- CAGQ- CATH- CCAP- CCCP- CCGP- CCTP- "
    "CGAR- CGCR- CGGR- CGTR- CTAL- CTCL- CTGL- CTTL- "
    "GAAE- GACD- GAGE- GATD- GCAA- GCCA- GCGA- GCTA- "
    "GGAG- GGCG- GGGG- GGTG- GTAV- GTCV- GTGVM GTTV- "
    "TAA** TACY- TAG** TATY- TCAS- TCCS- TCGS- TCTS- "
    "TGAW- TGCC- TGGW- TGTC- TTAL- TTCF- TTGLM TTTF-";

}   // namespace NCBITables


}   // end namespace


#endif //KGL_TABLE_NCBI_H