//
// Created by kellerberrin on 26/01/18.
//

#include <fstream>

#include "kel_utility.h"

#include "kol_ParserAnnotationGaf.h"

#include "kel_exec_env.h"
#include "kgl_gaf_parser.h"


namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;



bool kgl::GeneOntology::readGafFile(const std::string &file_name) {

  ExecEnv::log().info("Reading Gaf file: {}", file_name);

  gaf_record_vector_ = kol::ParserAnnotationGaf::readAnnotationFile(file_name);

  ExecEnv::log().info("Processed: {} Gaf records", gaf_record_vector_.size());

  return true;

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
