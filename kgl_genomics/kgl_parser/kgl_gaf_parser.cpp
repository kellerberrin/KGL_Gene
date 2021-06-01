//
// Created by kellerberrin on 26/01/18.
//

#include <fstream>

#include "kel_utility.h"
#include "kel_exec_env.h"
#include "kgl_gaf_parser.h"


namespace kgl = kellerberrin::genome;



void kgl::OntologyRecord::addGafRecord(const std::shared_ptr<const kol::GAFRecord>& gaf_record_ptr) {

  go_records_.emplace(gaf_record_ptr->goIdent(), gaf_record_ptr);

}



bool kgl::GeneOntology::readGafFile(const std::string &file_name) {


  ExecEnv::log().info("Reading Gaf file: {}", file_name);

  std::ifstream gaf_file;

  // Open input file.

  gaf_file.open(file_name);

  if (not gaf_file.good()) {

    ExecEnv::log().critical("I/O error; could not open Gaf file: {}", file_name);

  }


  size_t counter{0};
  std::string record_str;
  while (not std::getline(gaf_file, record_str).eof()) {

    if ((record_str)[0] == HEADER_CHAR_) {

      continue;   // ignore header records.

    }

    std::shared_ptr<kol::GAFRecord> gaf_record_ptr(std::make_shared<kol::GAFRecord>());
    if (gaf_record_ptr->parseGafRecord(record_str)) {

      semanticGafParse(gaf_record_ptr);

    }

    ++counter;

    if (counter % REPORT_INCREMENT_ == 0) {

      ExecEnv::log().info("Read: {} Gaf records", counter);

    }

  }

  gaf_file.close();

  ExecEnv::log().info("Processed: {} Gaf records", counter);


  return true;

}


void kgl::GeneOntology::semanticGafParse(const std::shared_ptr<const kol::GAFRecord>& gaf_record_ptr) {



  auto result = gaf_record_map_.find(gaf_record_ptr->geneUniprotId());
  if (result == gaf_record_map_.end()) {

  // create an ontology object and insert the gaf record.
    std::shared_ptr<OntologyRecord> ontology_ptr(std::make_shared<OntologyRecord>(gaf_record_ptr->geneUniprotId(),
                                                                                gaf_record_ptr->geneSymbolicId(),
                                                                                gaf_record_ptr->altSymbolicRef(),
                                                                                gaf_record_ptr->description()));
    ontology_ptr->addGafRecord(gaf_record_ptr);

    auto [iter, result] = gaf_record_map_.try_emplace(gaf_record_ptr->geneUniprotId(), ontology_ptr);
    if (not result) {

      ExecEnv::log().error("GeneOntology::semanticGafParse; Unexpected, Could not insert gaf (GO:) record for (duplicate) gene: {}", gaf_record_ptr->geneUniprotId());

    }

  } else {
  // Just insert

    auto& [gene_id, ontology_ptr] = *result;
    ontology_ptr->addGafRecord(gaf_record_ptr);

  }

}


std::optional<std::shared_ptr<const kgl::OntologyRecord>> kgl::GeneOntology::getGafFeatureVector(const FeatureIdent_t& gene_id) const {


  auto result = gaf_record_map_.find(gene_id);
  if (result == gaf_record_map_.end()) {

    return std::nullopt;

  }

  const auto& [gene_ident, ontology_ptr] = *result;
  return ontology_ptr;

}



bool kgl::GeneOntology::readIdFile(const std::string& filename) {

  ParseGeneIdents parse_gene_idents;
  if (not parse_gene_idents.parseIdentFile(filename)) {

    return false;

  }

  synonym_vector_ = parse_gene_idents.getSynonymVector();

  return true;

}


// The Gaf data structure with records sorted by the Symbolic Reference field.
void kgl::ResortGaf::sortBySymbolic(const GafRecordMap& gaf_map) {

  static bool warning{false};
  gaf_record_map_.clear();
  for (auto const& [gaf_id, gaf_record_ptr] : gaf_map) {

    auto [iter, result] = gaf_record_map_.try_emplace(gaf_record_ptr->symbolicReference(), gaf_record_ptr);

    if (not result and not warning) {

      ExecEnv::log().warn( "ResortGaf::sortBySymbolic; Cannot insert (duplicate) Gaf record: {}, Symbolic Ref: {}",
                           gaf_id, gaf_record_ptr->symbolicReference());
      warning = true;

    }

  }

}


// The Gaf data structure with records sorted by the Gene id field.
void kgl::ResortGaf::sortByGeneId(const GafRecordMap& gaf_map) {

  static bool warning{false};
  gaf_record_map_.clear();
  for (auto const& [gaf_id, gaf_record_ptr] : gaf_map) {

    auto [iter, result] = gaf_record_map_.try_emplace(gaf_record_ptr->gene_uniprot_id(), gaf_record_ptr);

    if (not result and not warning) {

      ExecEnv::log().warn( "ResortGaf::sortByGeneId; Cannot insert (duplicate) Gaf record: {}, Gene Id: {}",
                            gaf_id, gaf_record_ptr->gene_uniprot_id());
      warning = true;
    }

  }

}


void kgl::ResortIds::sortByHGNC(const GeneSynonymVector& ident_vector) {

  static bool warning{false};
  gene_id_map_.clear();
  for (auto const& [HGNC_id, ensembl_id] :  ident_vector) {

    if (not HGNC_id.empty() and not ensembl_id.empty()) {

      auto [iter, result] = gene_id_map_.try_emplace(HGNC_id, ensembl_id);
      if (not result and not warning) {

        ExecEnv::log().warn("ResortIds::sortByHGNC; Cannot insert (duplicate) HGNC id: {}, ensembl id: {}", HGNC_id, ensembl_id);
        warning = true;

      }

    }

  }

}


void kgl::ResortIds::sortByEnsembl(const GeneSynonymVector& ident_vector) {

  static bool warning{false};
  gene_id_map_.clear();
  for (auto const& [HGNC_id, ensembl_id] :  ident_vector) {

    if (not HGNC_id.empty() and not ensembl_id.empty()) {

      auto [iter, result] = gene_id_map_.try_emplace(ensembl_id, HGNC_id);
      if (not result and not warning) {

        ExecEnv::log().warn("ResortIds::sortByEnsembl; Cannot insert (duplicate) ensembl id: {}, HGNC id: {}", ensembl_id, HGNC_id);
        warning = true;

      }

    }

  }

}
