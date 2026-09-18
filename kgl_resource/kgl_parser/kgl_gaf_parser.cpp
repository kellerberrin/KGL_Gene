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

