//
// Created by kellerberrin on 12/03/18.
//

#ifndef KGL_PF3K_AUX_CSV_H
#define KGL_PF3K_AUX_CSV_H

#include <string>
#include <vector>
#include <map>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reads up the Pf3k auxillary data description file in csv format, 1 line for each sample.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using AuxAttributeVector = std::vector<std::string>;
using Pf3kSampleMap = std::map<std::string, AuxAttributeVector>;


class Pf3kAuxData {

public:

  Pf3kAuxData() = default;
  ~Pf3kAuxData() = default;

  bool readParseAuxData(const std::string& aux_file_name);

private:

  AuxAttributeVector aux_data_header_;  // Always assumed to be the first line (uppercase, query case conversion automatic).
  Pf3kSampleMap aux_sample_infromation_;

};






}   // namespace genome
}   // namespace kellerberrin









#endif // KGL_PF3K_AUX_CSV_H
