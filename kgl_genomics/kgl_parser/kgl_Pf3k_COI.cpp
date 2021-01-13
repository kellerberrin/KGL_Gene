//
// Created by kellerberrin on 9/1/21.
//

#include "kgl_Pf3k_COI.h"


namespace kgl = kellerberrin::genome;


bool kgl::Pf3kCOIParser::parseCOIPf3k(const std::string& file_name) {

  std::shared_ptr<SquareTextRows> parsed_COI_file = flat_file_parser_.parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);
  if (not pf3k_coi_ptr_->parseFlatFile(*parsed_COI_file)) {

    ExecEnv::log().error("Pf3kCOIParser::parseCOIPf3k; failed to parse flat file format into column indexed format");
    return false;

  }

  if (pf3k_coi_ptr_->indexedFile().getHeaderMap().size() != Pf3kCOIDB::FIELD_COUNT) {

    ExecEnv::log().error( "Pf3kCOIParser::parseCOIPf3k; expected fields: {}, actual fields: {}",
                          Pf3kCOIDB::FIELD_COUNT, pf3k_coi_ptr_->indexedFile().getHeaderMap().size());
    return false;

  }

  return true;

}
