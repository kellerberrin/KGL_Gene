//
// Created by kellerberrin on 30/6/21.
//

#include "kgl_analysis_mutation_data.h"


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

