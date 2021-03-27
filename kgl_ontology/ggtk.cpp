//
// Created by kellerberrin on 27/3/21.
//

#include "ggtk.hpp"

namespace kellerberrin::ontology {   //  organization::project level namespace




auto test_ggtk(const std::string& goa_file) {

  GoaAnnotationParser annotation_parser;

  return annotation_parser.parseAnnotationFile(goa_file);

}


}