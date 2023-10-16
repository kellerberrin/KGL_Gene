//
// Created by kellerberrin on 30/6/21.
//

#include "kga_analysis_mutation_data.h"


namespace kgl = kellerberrin::genome;

const std::vector<std::string> kgl::MutationAnalysisData::OMIMGeneSymbol() {

  std::vector<std::string> gene_vector;
  for (auto const&[ensembl, symbol] : omim_ensembl_symbol_) {

    gene_vector.push_back(symbol);

  }

  return gene_vector;

}

const std::vector<std::string> kgl::MutationAnalysisData::OMIMGeneEnsembl() {

  std::vector<std::string> gene_vector;
  for (auto const&[ensembl, symbol] : omim_ensembl_symbol_) {

    gene_vector.push_back(ensembl);

  }

  return gene_vector;

}

// Malaria active genes harvested from the Uniprot website.
const std::vector<std::string> kgl::MutationAnalysisData::UniprotGeneSymbol() {

  std::vector<std::string> gene_vector;
  for (auto const&[ensembl, symbol] : uniprot_ensembl_symbol_) {

    gene_vector.push_back(symbol);

  }

  return gene_vector;

}

// Malaria active genes harvested from the Uniprot website.
const std::vector<std::string> kgl::MutationAnalysisData::UniprotGeneEnsembl() {

  std::vector<std::string> gene_vector;
  for (auto const&[ensembl, symbol] : uniprot_ensembl_symbol_) {

    gene_vector.push_back(ensembl);

  }

  return gene_vector;

}



// From the OMIM entry #611162 available at https://www.omim.org/entry/611162
const std::vector<std::pair<std::string, std::string>> kgl::MutationAnalysisData::omim_ensembl_symbol_{

    {"ENSG00000162692", "VCAM1"},
    {"ENSG00000162706", "CADM3"},
    {"ENSG00000213088", "ACKR1"},
    {"ENSG00000143226", "FCGR2A"},
    {"ENSG00000072694", "FCGR2B"},
    {"ENSG00000058668", "ATP2B4"},
    {"ENSG00000196352", "CD55"},
    {"ENSG00000203710", "CR1"},
    {"ENSG00000136732", "GYPC"},
    {"ENSG00000114737", "CISH"},
    {"ENSG00000250361", "GYPB"},
    {"ENSG00000170180", "GYPA"},
    {"ENSG00000228978", "TNF"},
    {"ENSG00000236315", "NCR3"},
    {"ENSG00000172461", "FUT9"},
    {"ENSG00000135218", "CD36"},
    {"ENSG00000281879", "ABO"},
    {"ENSG00000244734", "HBB"},
    {"ENSG00000150455", "TIRAP"},
    {"ENSG00000206172", "HBA1"},
    {"ENSG00000140832", "MARVELD3"},
    {"ENSG00000007171", "NOS2"},
    {"ENSG00000004939", "SLC4A1"},
    {"ENSG00000172270", "BSG"},
    {"ENSG00000090339", "ICAM1"},
    {"ENSG00000276163", "LAIR1"},
    {"ENSG00000276452", "LILRB1"},
    {"ENSG00000101000", "PROCR"},
    {"ENSG00000160211", "G6PD"}

};

// From the Uniprot website.
const std::vector<std::pair<std::string, std::string>> kgl::MutationAnalysisData::uniprot_ensembl_symbol_ = {

    {"ENSG00000213088", "ACKR1"},
    {"ENSG00000072694", "FCGR2B"},
    {"ENSG00000203710", "CR1"},
    {"ENSG00000136732", "GYPC"},
    {"ENSG00000114737", "CISH"},
    {"ENSG00000170180", "GYPA"},
    {"ENSG00000228978", "TNF"},
    {"ENSG00000236315", "NCR3"},
    {"ENSG00000135218", "CD36"},
    {"ENSG00000158856", "DMTN"},
    {"ENSG00000029534", "ANK1"},
    {"ENSG00000244734", "HBB"},
    {"ENSG00000150455", "TIRAP"},
    {"ENSG00000157837", "SPPL3"},
    {"ENSG00000070182", "SPTB"},
    {"ENSG00000007171", "NOS2"},
    {"ENSG00000004939", "SLC4A1"},
    {"ENSG00000261371", "PECAM1"},
    {"ENSG00000172270", "BSG"},
    {"ENSG00000005206", "SPPL2B"},
    {"ENSG00000090339", "ICAM1"},
    {"ENSG00000160211", "G6PD"},
    {"ENSG00000250361", "GYPB"}

};


// Malaria active genes harvested from the Uniprot website.
const std::vector<std::string> kgl::MutationAnalysisData::adHocLILRB1GeneSymbol() {

  std::vector<std::string> gene_vector;
  for (auto const&[ensembl, symbol] : adhoc_LILRB1_ensembl_symbol_) {

    gene_vector.push_back(symbol);

  }

  return gene_vector;

}

// Malaria active genes harvested from the Uniprot website.
const std::vector<std::string> kgl::MutationAnalysisData::adHocLILRB1GenesEnsembl() {

  std::vector<std::string> gene_vector;
  for (auto const&[ensembl, symbol] : adhoc_LILRB1_ensembl_symbol_) {

    gene_vector.push_back(ensembl);

  }

  return gene_vector;

}


// Used in adhoc LILRB1 gene  polymorphism analysis (chromosome 19).
const std::vector<std::pair<std::string, std::string>> kgl::MutationAnalysisData::adhoc_LILRB1_ensembl_symbol_ = {

    { "ENSG00000276452", "LILRB1"},  // LILRB1
    { "ENSG00000276163", "LAIR1" },  //  LAIR1
    { "ENSG00000277816", "LILRB3" }, //   LILRB3
    { "ENSG00000275584", "LILRA6"},  // LILRA6	leukocyte immunoglobulin like receptor A6
    { "ENSG00000277414", "LILRB5"},  // LILRB5	leukocyte immunoglobulin like receptor B5
    { "ENSG00000274513", "LILRB2" }, //	LILRB2	leukocyte immunoglobulin like receptor B2
    { "ENSG00000278355", "LILRA5" }, //	LILRA5	leukocyte immunoglobulin like receptor A5
    { "ENSG00000276798", "LILRA4"},  //	LILRA4	leukocyte immunoglobulin like receptor A4
    { "ENSG00000274084", "LAIR2" },  //	LAIR2	leukocyte associated immunoglobulin like receptor 2
    { "ENSG00000274000", "LILRA2" }, //	LILRA2	leukocyte immunoglobulin like receptor A2
    { "ENSG00000274935", "LILRA1" }, //	LILRA1	leukocyte immunoglobulin like receptor A1
    { "ENSG00000278555", "LILRB4" },  // 	LILRB4	leukocyte immunoglobulin like receptor B4
    { "ENSG00000276433", "KIR3DL3" }, //	KIR3DL3	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail_ 3
    { "ENSG00000273947", "KIR2DL3" }, //	KIR2DL3	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail_ 3
    { "ENSG00000276820", "KIR2DL1" }, //	KIR2DL1	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail_ 1
    { "ENSG00000276779", "KIR2DL4" }, //	KIR2DL4	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail_ 4
    { "ENSG00000273775", "KIR3DL1" }, //	KIR3DL1	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail_ 1
    {"ENSG00000274324", "KIR2DS4" }, //	KIR2DS4	killer cell immunoglobulin like receptor%2C two Ig domains and short cytoplasmic tail_ 4
    { "ENSG00000273735",	"KIR3DL2"}  // KIR3DL2	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail_ 2

};

// Used in adhoc polymorphism analysis (chromosome 19).
const std::vector<std::string> kgl::MutationAnalysisData::alt_LILRB1_ensembl_symbol_ = {

    "ENSG00000167613",   // Alt  LAIR1
    "ENSG00000104972", // Alt LILRB1
    "ENSG00000204577", // Alt LILRB3
    "ENSG00000244482", // Alt LILRA6
    "ENSG00000105609", // Alt LILRB5	leukocyte immunoglobulin like receptor B5
    "ENSG00000131042", // Alt LILRB2	leukocyte immunoglobulin like receptor B2
    "ENSG00000187116", //	Alt LILRA5	leukocyte immunoglobulin like receptor A5
    "ENSG00000239961", // Alt LILRA4	leukocyte immunoglobulin like receptor A4
    "ENSG00000167618", //	Alt LAIR2	leukocyte associated immunoglobulin like receptor 2
    "ENSG00000239998", //	Alt LILRA2	leukocyte immunoglobulin like receptor A2
    "ENSG00000104974", //	Alt LILRA1	leukocyte immunoglobulin like receptor A1
    "ENSG00000186818",  // 	Alt LILRB4	leukocyte immunoglobulin like receptor B4
    "ENSG00000242019", //	Alt KIR3DL3	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail_ 3
    "ENSG00000243772", //	Alt KIR2DL3	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail_ 3
    "ENSG00000125498", //	Alt KIR2DL1	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail_ 1
    "ENSG00000189013", //	Alt KIR2DL4	killer cell immunoglobulin like receptor%2C two Ig domains and long cytoplasmic tail_ 4
    "ENSG00000167633", //	Alt KIR3DL1	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail_ 1
    "ENSG00000221957", //	Alt KIR2DS4	killer cell immunoglobulin like receptor%2C two Ig domains and short cytoplasmic tail_ 4
    "ENSG00000240403"	// Alt KIR3DL2	killer cell immunoglobulin like receptor%2C three Ig domains and long cytoplasmic tail_ 2

};





// Used in adhoc FCGR gene  polymorphism analysis (chromosome 1).
const std::vector<std::string> kgl::MutationAnalysisData::adhoc_FCGR_ensembl_symbol_ = {

    "ENSG00000150337",   // FCGR1A
    "ENSG00000143226",   // FCGR2A
    "ENSG00000072694",   // FCGR2B
    "ENSG00000162747",  // FCGR3B
    "ENSG00000203747"   // FCGR3A

};

// The list of genes to be analyzed variant by variant. Must be ensembl codes (for now).
const std::vector<std::string> kgl::MutationAnalysisData::ontology_derived_gene_list_{

    "ENSG00000284690", "ENSG00000282992", "ENSG00000275019", "ENSG00000262576", "ENSG00000256797",
    "ENSG00000254521", "ENSG00000213402", "ENSG00000204345", "ENSG00000198178", "ENSG00000196371",
    "ENSG00000189184", "ENSG00000188211", "ENSG00000186407", "ENSG00000185475", "ENSG00000185187",
    "ENSG00000183840", "ENSG00000183019", "ENSG00000180549", "ENSG00000179213", "ENSG00000172794",
    "ENSG00000172322", "ENSG00000171840", "ENSG00000170956", "ENSG00000170425", "ENSG00000169704",
    "ENSG00000167850", "ENSG00000167123", "ENSG00000166589", "ENSG00000165682", "ENSG00000164713",
    "ENSG00000163600", "ENSG00000163485", "ENSG00000162897", "ENSG00000161649", "ENSG00000159674"

};

const std::vector<std::string> kgl::MutationAnalysisData::malaria_MeSH_list_{
/*
  "MESH:D008288",  // Malaria (Main Subject).
  "MESH:D016779",  // Cerebral Malaria
*/
  "MESH:D016778",   // Malaria Falciparum
  "MESH:D010963",    // Plasmodium falciparum
/*
  "MESH:D016780",  // Malaria Vivax
  "MESH:D008289",  // Malaria Avian
  "MESH:D010965",  // Plasmodium Malariae
  "MESH:D017780",  // Malaria Vaccines
  "MESH:D001742",  // Blackwater Fever.
  "MESH:D010922",  // Placental Infection
  "MESH:D007239",  // Malarial infection
  "MESH:D018512",  // Malaria parasitemias|asexual parasitemia versus gametocytemia|parasitemia
  "MESH:D000740",  // Parasitaemia and anaemia|anemia|anaemia|malarial anaemia|anaemia - Anaemia
  "MESH:D031261",	 // Malarial mortality
  "MESH:D010272",	 // Malarial parasites
  "MESH:D007022",	 // Blood-stage malarial infection
  "MESH:D011488",	 // Malarial proteins
  "MESH:D014397",	 // malarial drug consumption|malarial consumption
  "MESH:C531736"	 // acute malaria coinfection|acute malaria
*/
};