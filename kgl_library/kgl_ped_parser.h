//
// Created by kellerberrin on 1/8/20.
//

#ifndef KGL_PED_PARSER_H
#define KGL_PED_PARSER_H


#include "kgl_data_base.h"
#include "kgl_data_file_impl.h"

namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file data structure.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GenomePEDData : public DataObjectBase {

public:

  explicit GenomePEDData(std::string ped_ident) : DataObjectBase(std::move(ped_ident)) {}
  ~GenomePEDData() override = default;

  DataTypeEnum dataType() override { return DataTypeEnum::PedAncestor; }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file parser.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ParsePedFile : public FileDataIO {

public:

  explicit ParsePedFile(std::shared_ptr<GenomePEDData> ped_data,
                        std::string file_name) : FileDataIO(std::move(file_name)), ped_data_(ped_data) {}
  ~ParsePedFile() = default;

  void readParsePEDImpl() {}

private:

  std::shared_ptr<GenomePEDData> ped_data_;

};



} // namespace

#endif //KGL_KGL_PED_PARSER_H
