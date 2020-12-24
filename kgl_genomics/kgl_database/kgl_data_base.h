//
// Created by kellerberrin on 31/7/20.
//

#ifndef KGL_DATA_BASE_H
#define KGL_DATA_BASE_H

#include <string>

namespace kellerberrin::genome {   //  organization::project


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A do nothing base class to attach generic pointers to parsed data files.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// enum class DataTypeEnum


enum class DataTypeEnum { DiploidPopulation, HaploidPopulation, UnphasedPopulation, PedAncestor, NotImplemented};


class DataObjectBase {

public:

  DataObjectBase(std::string identifier) : data_identifier_(std::move(identifier)) {}
  virtual ~DataObjectBase() = default;

  [[nodiscard]] const std::string& Id() const { return data_identifier_; }
  void setId(const std::string& data_id) { data_identifier_ = data_id; }

  virtual DataTypeEnum dataType() const = 0;


private:

   std::string data_identifier_;

};


} // namespace


#endif //KGL_KGL_DATA_FILE_H
