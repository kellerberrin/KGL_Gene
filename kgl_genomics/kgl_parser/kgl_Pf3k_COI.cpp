//
// Created by kellerberrin on 9/1/21.
//

#include "kgl_Pf3k_COI.h"


namespace kgl = kellerberrin::genome;


bool kgl::Pf3kCOIParser::parseCOIPf3k(const std::string& file_name) {

  std::shared_ptr<SquareTextRows> parsed_COI_file = flat_file_parser_.parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);

  return true;

}
