//
// Created by kellerberrin on 26/5/21.
//
#include "kol_ParserAnnotationGaf.h"

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <fstream>


namespace kol = kellerberrin::ontology;


std::shared_ptr<const kol::TermAnnotation> kol::ParserAnnotationGaf::parseAnnotationFile( const PolicyEvidence &policy,
                                                                                          const std::string &file_name,
                                                                                          AnnotationGeneName gene_id_type) {


  std::vector<std::shared_ptr<const GAFRecord>> gaf_records;

  // Check that the supplied policy is valid.
  if (not policy.isValid()) {

    ExecEnv::log().error("ParserAnnotationGaf::parseAnnotationFile; invalid annotation policy for file: {}", file_name);
    return std::make_shared<TermAnnotation>(gaf_records);

  }

  gaf_records = readAnnotationFile(file_name);

  std::shared_ptr<const TermAnnotation> annoData(std::make_shared<TermAnnotation>(policy, gaf_records, gene_id_type));

  return annoData;

}


std::vector<std::shared_ptr<const kol::GAFRecord>> kol::ParserAnnotationGaf::readAnnotationFile(const std::string &file_name) {


  std::vector<std::shared_ptr<const GAFRecord>> gaf_records;

  std::ifstream anno_file(file_name);
  if (not anno_file.good()) {

    ExecEnv::log().error("ParserAnnotationGaf::parseAnnotationFile; problem opening annotation file: {}", file_name);
    return  std::vector<std::shared_ptr<const GAFRecord>>();

  }

  std::string line;
  size_t line_count{0};
  while (not std::getline(anno_file, line).eof()) {

    ++line_count;
    if (line[0] == COMMENT_CHAR_) {

      continue;

    }

    auto view_vector = Utility::view_tokenizer(line, TAB_FIELD_DELIMITER_);
    if (view_vector.size() != EXPECTED_FIELD_COUNT_) {

      ExecEnv::log().error("ParserAnnotationGaf::readAnnotationFile; line: {}, expected fields: {}, found: {}, file: {}",
                           line_count, EXPECTED_FIELD_COUNT_, view_vector.size(), file_name);
      continue;

    }

    if (not view_vector[FIELD_OFFSET_QUALIFIER_].empty()) {

      auto qualifier = Utility::toupper(std::string(view_vector[FIELD_OFFSET_QUALIFIER_]));
      if (qualifier.find(NOT_QUALIFIER_) != std::string::npos) {

        continue;

      }

    }

    std::shared_ptr<GAFRecord> record_ptr(std::make_shared<GAFRecord>());
    if (not record_ptr->parseGafRecord(line)) {

      ExecEnv::log().error("ParserAnnotationGaf::readAnnotationFile; error parsing Gaf record from line: {}", line);

    } else {

      gaf_records.push_back(record_ptr);

    }

  }

  return gaf_records;

}
