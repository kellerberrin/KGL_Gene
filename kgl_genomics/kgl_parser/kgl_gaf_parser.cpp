//
// Created by kellerberrin on 26/01/18.
//

#include <fstream>

#include "kel_exec_env.h"
#include "kgl_gaf_parser.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace bt = boost;
namespace kgl = kellerberrin::genome;



void kgl::OntologyRecord::addGafRecord(GAFRecord gaf_record) {

  go_records_.emplace(gaf_record.ontolotgy_id, std::move(gaf_record));

}



bool kgl::GeneOntology::readGafFile(const std::string &file_name) {


  ExecEnv::log().info("Reading Gaf file: {}", file_name);

  std::ifstream gaf_file;

  // Open input file.

  gaf_file.open(file_name);

  if (not gaf_file.good()) {

    ExecEnv::log().critical("I/O error; could not open Gaf file: {}", file_name);

  }

  try {

    long counter = 0;

    while (true) {

      std::string record_str;

      if (std::getline(gaf_file, record_str).eof()) break;

      if ((record_str)[0] == '!') continue;   // ignore header records.

      parseGafRecord(record_str);

      ++counter;

      if (counter % REPORT_INCREMENT_ == 0) {

        ExecEnv::log().info("Read: {} Gaf records", counter);

      }

    }

    gaf_file.close();

    ExecEnv::log().info("Processed: {} Gaf records", counter);

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("Gaf file: {}, unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}


void kgl::GeneOntology::parseGafRecord(const std::string& record_str) {

  std::vector<std::string> field_vec;
  bt::char_separator<char> item_key_sep("\t","", bt::keep_empty_tokens);
  bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);
  for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

    field_vec.push_back(*iter_item);

  }

  if (field_vec.size() != GO_FIELD_COUNT_) {

    ExecEnv::log().error("Unexpected; Gaf record should have: {} fields, parsed: {} fields", GO_FIELD_COUNT_, field_vec.size());
    ExecEnv::log().error("Gaf record : {}", record_str);

  } else {

    semanticGafParse(field_vec);

  }


}


void kgl::GeneOntology::semanticGafParse(const std::vector<std::string> &field_vec) {


  GAFRecord gaf_record;

  gaf_record.DB_ID = field_vec[0];              // required (1)
  gaf_record.gene_id = field_vec[1];            // required (1) key
  gaf_record.symbolic_ref = field_vec[2];	      // required (1)
  gaf_record.Qualifier = field_vec[3];          // optional (n)
  gaf_record.ontolotgy_id = field_vec[4];       // required (n)
  gaf_record.DB_Reference = field_vec[5];       // required (n)
  gaf_record.Evidence_Code = field_vec[6];	    // required (n)
  gaf_record.With_From = field_vec[7];	        // optional (n)
  gaf_record.Aspect = field_vec[8];             // required (n) biological (P), cellular (C) or molecular (F)
  gaf_record.description = field_vec[9];	      // optional (1)
  gaf_record.alt_symbolic_ref = field_vec[10];       // required (1)
  gaf_record.DB_Object_Type = field_vec[11];	  // required (n)
  gaf_record.Taxon = field_vec[12];             // required (1)
  gaf_record.Date = field_vec[13];              // required (1)
  gaf_record.Assigned_By = field_vec[14];       // required (n)
  gaf_record.Annotation = field_vec[15];        // optional (n)
  gaf_record.Gene_Product = field_vec[16];      // optional (n)

  auto result = gaf_record_map_.find(gaf_record.gene_id);
  if (result == gaf_record_map_.end()) {

  // create an ontology object and insert the gaf record.
    std::shared_ptr<OntologyRecord> ontology_ptr(std::make_shared<OntologyRecord>(gaf_record.gene_id,
                                                                                gaf_record.symbolic_ref,
                                                                                gaf_record.alt_symbolic_ref,
                                                                                gaf_record.description));
    ontology_ptr->addGafRecord(gaf_record);

    std::pair<FeatureIdent_t, std::shared_ptr<OntologyRecord>> insert_pair(gaf_record.gene_id, ontology_ptr);
    auto insert_result = gaf_record_map_.insert(insert_pair);

    if (not insert_result.second) {

      ExecEnv::log().error("Unexpected, Could not insert gaf (GO:) record for gene: {}", gaf_record.gene_id);

    }

  } else {
  // Just insert

    result->second->addGafRecord(gaf_record);

  }

}


std::optional<std::shared_ptr<const kgl::OntologyRecord>> kgl::GeneOntology::getGafFeatureVector(const FeatureIdent_t& gene_id) const {


  auto result = gaf_record_map_.find(gene_id);

  if (result == gaf_record_map_.end()) {

    return std::nullopt;

  }

  const auto [gene_ident, ontology_ptr] = *result;
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

  gaf_record_map_.clear();
  for (auto const& [gaf_id, gaf_record_ptr] : gaf_map) {

    auto [iter, result] = gaf_record_map_.try_emplace(gaf_record_ptr->symbolicReference(), gaf_record_ptr);

    if (not result) {

      ExecEnv::log().warn( "ResortGaf::sortBySymbolic; Cannot insert (duplicate) Gaf record: {}, Symbolic Ref: {}",
                           gaf_id, gaf_record_ptr->symbolicReference());

    }

  }

}


// The Gaf data structure with records sorted by the Gene id field.
void kgl::ResortGaf::sortByGeneId(const GafRecordMap& gaf_map) {

  gaf_record_map_.clear();
  for (auto const& [gaf_id, gaf_record_ptr] : gaf_map) {

    auto [iter, result] = gaf_record_map_.try_emplace(gaf_record_ptr->gene_id(), gaf_record_ptr);

    if (not result) {

      ExecEnv::log().warn( "ResortGaf::sortByGeneId; Cannot insert (duplicate) Gaf record: {}, Gene Id: {}",
                            gaf_id, gaf_record_ptr->gene_id());

    }

  }

}


void kgl::ResortIds::sortByHGNC(const GeneSynonymVector& ident_vector) {

  gene_id_map_.clear();
  for (auto const& [HGNC_id, ensembl_id] :  ident_vector) {

    if (not HGNC_id.empty() and not ensembl_id.empty()) {

      auto [iter, result] = gene_id_map_.try_emplace(HGNC_id, ensembl_id);
      if (not result) {

        ExecEnv::log().warn("ResortIds::sortByHGNC; Cannot insert (duplicate) HGNC id: {}, ensembl id: {}", HGNC_id, ensembl_id);

      }

    }

  }

}


void kgl::ResortIds::sortByEnsembl(const GeneSynonymVector& ident_vector) {

  gene_id_map_.clear();
  for (auto const& [HGNC_id, ensembl_id] :  ident_vector) {

    if (not HGNC_id.empty() and not ensembl_id.empty()) {

      auto [iter, result] = gene_id_map_.try_emplace(ensembl_id, HGNC_id);
      if (not result) {

        ExecEnv::log().warn("ResortIds::sortByEnsembl; Cannot insert (duplicate) ensembl id: {}, HGNC id: {}", ensembl_id, HGNC_id);

      }

    }

  }

}
