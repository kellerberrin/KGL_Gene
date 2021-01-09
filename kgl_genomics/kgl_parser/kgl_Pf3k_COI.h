//
// Created by kellerberrin on 9/1/21.
//

#ifndef KGL_COI_PF3K_H
#define KGL_COI_PF3K_H


#include "kgl_variant_db_type.h"
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

private:

  const std::string file_ident_;

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
