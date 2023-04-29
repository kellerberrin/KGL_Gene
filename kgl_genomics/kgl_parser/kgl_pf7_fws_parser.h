//
// Created by kellerberrin on 3/04/23.
//

#ifndef KGL_PF7_FWS_PARSER_H
#define KGL_PF7_FWS_PARSER_H



#include "kgl_properties_resource.h"
#include "kgl_square_parser.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Gene nomenclature resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Pf7FwsRecord {

  std::string Pf7Sample_id;
  double FWS_value;

};

using Pf7FwsVector = std::vector<Pf7FwsRecord>;
using Pf7FwsMap = std::map<std::string, Pf7FwsRecord>;
enum class FwsFilterType { GREATER_EQUAL, LESS_EQUAL };

class Pf7FwsResource : public ResourceBase {

public:

  Pf7FwsResource(std::string identifier, Pf7FwsVector Pf7Fws_vector)
      : ResourceBase(ResourceProperties::PF7FWS_RESOURCE_ID_, std::move(identifier)),
        Pf7Fws_vector_(std::move( Pf7Fws_vector)) {

    indexPf7FwsData();

  }
  ~Pf7FwsResource() override = default;

  [[nodiscard]] const Pf7FwsMap& getMap() const { return Pf7Fws_map_; }

  [[nodiscard]] double getFWS(const GenomeId_t& genome) const;

  // Population must be PF7
  [[nodiscard]] std::shared_ptr<PopulationDB> filterFWS( FwsFilterType filter_type,
                                                         double fws_threshold,
                                                         const std::shared_ptr<const PopulationDB>& Pf7_unfiltered_ptr) const;

private:

  const Pf7FwsVector Pf7Fws_vector_;
  Pf7FwsMap Pf7Fws_map_;

  void indexPf7FwsData();

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParsePf7Fws : public SquareTextParser {

public:

  ParsePf7Fws() = default;
  ~ParsePf7Fws() = default;

  [[nodiscard]] bool parsePf7FwsFile(const std::string& file_name);
  [[nodiscard]] const Pf7FwsVector& getPf7FwsVector() const { return Pf7Fws_vector_; }

private:

  Pf7FwsVector  Pf7Fws_vector_;

  constexpr static const size_t MINIMUM_ROW_COUNT_{1};
  constexpr static const size_t COLUMN_COUNT_{2};

  // Field record offsets.
  constexpr static const size_t SAMPLE_ID_OFFSET_{0};
  constexpr static const size_t FWS_OFFSET_{1};

};


} // namespace


#endif //KGL_PF7_FWS_PARSER_H
