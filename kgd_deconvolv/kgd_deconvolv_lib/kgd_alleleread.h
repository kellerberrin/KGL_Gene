//
// Created by kellerberrin on 19/05/18.
//

#ifndef KGL_KGD_ALLELEREAD_H
#define KGL_KGD_ALLELEREAD_H

#include "kgd_allelefreq.h"



namespace kellerberrin::deconvolv {    // organization level namespace


using GenomeAlleleVector = std::vector<GenomeAlleles>;
class AlleleReader {


public: // move the following to private

  AlleleReader() = default;
  ~AlleleReader() = default;

  // Methods

  bool parseFile(const std::string& file_name, const std::string& delimiters = DELIMITER_CHARACTERS_);


private:

  GenomeAlleleVector genome_allele_vector_;
  std::vector<std::string> headers_;

  constexpr static const char COMMENT_CHARACTER_ = '!';
  constexpr static const char HEADER_CHARACTER_ = '#';
  constexpr static const char* DELIMITER_CHARACTERS_ = "\t,"; // accept tabbed and csv files.
  constexpr static const size_t MINIMUM_FIELD_COUNT_ = 3;  // #chrom offset value

  // Methods

  bool parseDataLine(const std::vector<std::string>& text_fields);
  bool parseHeaderLine(const std::vector<std::string>& text_fields);

};


}   // end  namespace


#endif //KGL_KGD_ALLELEREAD_H
