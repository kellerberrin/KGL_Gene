//
// Created by kellerberrin on 19/05/18.
//

#ifndef KGL_KGD_ALLELEREAD_H
#define KGL_KGD_ALLELEREAD_H

#include "kgd_allelefreq.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


using GenomeAlleleVector = std::vector<GenomeAlleles>;
class AlleleReader {


public: // move the following to private

  AlleleReader() = default;
  ~AlleleReader() = default;

  // Methods

  bool parseFile(const std::string& file_name, const std::string& delimiters = DELIMITER_CHARACTERS_);


private:

  GenomeAlleleVector genome_allele_vector_;

  constexpr static const char COMMENT_CHARACTER_ = '#';
  constexpr static const char* DELIMITER_CHARACTERS_ = "\t,"; // accept tabbed and csv files.

  // Methods

  bool parseDataLine(const std::vector<std::string>& text_fields);


};


}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_ALLELEREAD_H
