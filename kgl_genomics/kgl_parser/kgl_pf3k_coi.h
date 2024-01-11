//
// Created by kellerberrin on 9/1/21.
//

#ifndef KGL_COI_PF3K_H
#define KGL_COI_PF3K_H

#include "kgl_properties_resource.h"
#include "kgl_square_parser.h"
#include "kgl_genome_types.h"


#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Read the Pf3k Complexity of Infection database file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Pf3kCOIResource : public ResourceBase {

public:

  Pf3kCOIResource(std::string identifier, const SquareTextRows& square_text);
  ~Pf3kCOIResource() override = default;

  [[nodiscard]] const SquareTextIndexed& indexedFile() const { return indexed_file_; }

  [[nodiscard]] std::optional<size_t> genomeCOI(const GenomeId_t& genome) const;
  [[nodiscard]] std::optional<std::string> genomeSite(const GenomeId_t& genome) const { return genomeField(genome, SITE_FIELD); }
  [[nodiscard]] std::optional<std::string> genomeCountry(const GenomeId_t& genome) const { return genomeField(genome, COUNTRY_FIELD); }
  [[nodiscard]] std::optional<std::string> genomeContinent(const GenomeId_t& genome) const { return genomeField(genome, CONTINENT_FIELD); }

  [[nodiscard]] std::optional<std::string> genomeField(const GenomeId_t& genome, size_t field_index) const;

  // Annotate the genome map with site and country information.
  [[nodiscard]] std::map<GenomeId_t, std::string> annotatedGenomeMap() const;

private:

  SquareTextIndexed indexed_file_;

  // Field offsets.
  constexpr static const size_t FIELD_COUNT = 10;
  constexpr static const size_t SAMPLE_FIELD = 0;
  constexpr static const size_t P1_FIELD = 1;
  constexpr static const size_t P2_FIELD = 2;
  constexpr static const size_t P3_FIELD = 3;
  constexpr static const size_t P4_FIELD = 4;
  constexpr static const size_t P5_FIELD = 5;
  constexpr static const size_t COI_FIELD = 6;
  constexpr static const size_t SITE_FIELD = 7;
  constexpr static const size_t COUNTRY_FIELD = 8;
  constexpr static const size_t CONTINENT_FIELD = 9;

  [[nodiscard]] bool parseFlatFile(const SquareTextRows& square_text) { return indexed_file_.parseTextRows(square_text); }

};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Pf3kCOIParser  {

public:

  Pf3kCOIParser() = default;
  ~Pf3kCOIParser() = default;

  [[nodiscard]] bool parseCOIPf3k(const std::string& file_name);
  [[nodiscard]] const SquareTextRows& parsedText() const { return *parsed_text_ptr_; }

public:

  SquareTextParser flat_file_parser_;
  std::shared_ptr<SquareTextRows> parsed_text_ptr_;

};


} // namespace


#endif //KGL_KGL_COI_PF3K_H
