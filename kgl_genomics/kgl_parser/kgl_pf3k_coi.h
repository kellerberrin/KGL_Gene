//
// Created by kellerberrin on 9/1/21.
//

#ifndef KGL_COI_PF3K_H
#define KGL_COI_PF3K_H


#include "kgl_data_file_type.h"
#include "kgl_data_file_impl.h"
#include "kgl_square_parser.h"


#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Read the Pf3k Complexity of Infection database file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Pf3kCOIDB : public DataDB {

public:

  Pf3kCOIDB(std::string file_ident, DataSourceEnum data_source) : DataDB(data_source), file_ident_(std::move(file_ident)) {}
  ~Pf3kCOIDB() override = default;

  [[nodiscard]] const std::string& fileId() const override { return file_ident_; }
  [[nodiscard]] const SquareTextIndexed& indexedFile() const { return indexed_file_; }
  bool parseFlatFile(const SquareTextRows& square_text) { return indexed_file_.parseTextRows(square_text); }

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

private:

  const std::string file_ident_;
  SquareTextIndexed indexed_file_;

};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Pf3kCOIParser  {

public:

  explicit Pf3kCOIParser(std::shared_ptr<Pf3kCOIDB> pf3k_coi_ptr) : pf3k_coi_ptr_(std::move(pf3k_coi_ptr)) {}
  ~Pf3kCOIParser() = default;

  bool parseCOIPf3k(const std::string& file_name);

public:



  std::shared_ptr<Pf3kCOIDB> pf3k_coi_ptr_;
  SquareTextParser flat_file_parser_;

};


} // namespace


#endif //KGL_KGL_COI_PF3K_H
