//
// Created by kellerberrin on 31/03/23.
//

#ifndef KGL_PF7_SAMPLE_PARSER_H
#define KGL_PF7_SAMPLE_PARSER_H


#include "kgl_resource_db.h"
#include "kgl_square_parser.h"


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Gene nomenclature resource object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Pf7SampleRecord {

  std::string Pf7Sample_id;
  std::string study_;
  std::string country_;
  std::string location1_;
  std::string country_latitude_;
  std::string country_longitude_;
  std::string location1_latitude_;
  std::string location1_longitude_;
  std::string year_;
  std::string ena_;
  std::string all_samples_;
  std::string population_;
  std::string callable_;
  std::string qc_pass_;
  std::string qc_fail_reason_;
  std::string sample_type_;
  std::string sample_in_pf6_;

};

using Pf7SampleVector = std::vector<Pf7SampleRecord>;
using Pf7SampleMap = std::map<std::string, Pf7SampleRecord>;

class Pf7SampleResource : public ResourceBase {

public:

  Pf7SampleResource(std::string identifier,
                    Pf7SampleVector Pf7sample_vector) : ResourceBase(std::move(identifier)),
                                                        Pf7Sample_vector_(std::move(Pf7sample_vector)) {

    IndexPf7SampleData();

  }
  ~Pf7SampleResource() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::PF7_SAMPLE_DATA; }
  [[nodiscard]] const Pf7SampleMap& getMap() const { return Pf7Sample_map_; }

private:

  const Pf7SampleVector Pf7Sample_vector_;
  Pf7SampleMap Pf7Sample_map_;

  void IndexPf7SampleData();

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse a gene ident file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParsePf7Sample : public SquareTextParser {

public:

  ParsePf7Sample() = default;
  ~ParsePf7Sample() = default;

  [[nodiscard]] bool parsePf7SampleFile(const std::string& file_name);
  [[nodiscard]] const Pf7SampleVector& getPf7SampleVector() const { return Pf7Sample_vector_; }

private:

  Pf7SampleVector  Pf7Sample_vector_;

  constexpr static const size_t MINIMUM_ROW_COUNT_{1};
  constexpr static const size_t COLUMN_COUNT_{17};

  // Field record offsets.
  constexpr static const size_t SAMPLE_ID_OFFSET_{0};
  constexpr static const size_t STUDY_OFFSET_{1};
  constexpr static const size_t COUNTRY_OFFSET_{2};
  constexpr static const size_t LOCATION1_OFFSET_{3};
  constexpr static const size_t COUNTRY_LATITUDE_OFFSET_{4};
  constexpr static const size_t COUNTRY_LONGITUDE_OFFSET_{5};
  constexpr static const size_t LOCATION1_LATITUDE_OFFSET_{6};
  constexpr static const size_t LOCATION1_LONGITUDE_OFFSET_{7};
  constexpr static const size_t YEAR_OFFSET_{8};
  constexpr static const size_t ENA_OFFSET_{9};
  constexpr static const size_t ALL_SAMPLES_OFFSET_{10};
  constexpr static const size_t POPULATION_OFFSET_{11};
  constexpr static const size_t CALLABLE_OFFSET_{12};
  constexpr static const size_t QC_PASS_OFFSET_{13};
  constexpr static const size_t QC_FAIL_REASON_OFFSET_{14};
  constexpr static const size_t SAMPLE_TYPE_OFFSET_{15};
  constexpr static const size_t SAMPLE_IN_PF6_OFFSET_{16};

};


} // namespace



#endif //KGL_KGL_PF7_SAMPLE_PARSER_H
