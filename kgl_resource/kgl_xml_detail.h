//
// Created by kellerberrin on 27/8/26.
//

#ifndef KGL_XML_DETAIL_H
#define KGL_XML_DETAIL_H


#include <string_view>
#include <initializer_list>
#include <string>


namespace kellerberrin::genome {   //  organization::project level namespace


/// Build a dotted runtime key such as "runTime.executeList.active".
/// Shared helper used by all XML section parsers to avoid duplicating
/// the same concatenation logic in every translation unit.
[[nodiscard]] inline std::string runtimeKey(std::initializer_list<std::string_view> parts) {

  std::string key("runTime");
  for (auto part : parts) {

    key.push_back('.');
    key.append(part);

  }
  return key;

}

/// Common XML tag constants that appear in multiple section parsers.
constexpr std::string_view RUNTIME_ROOT_{"runTime"};
constexpr std::string_view ACTIVE_{"active"};
constexpr std::string_view VALUE_{"value"};


} // namespace


#endif //KGL_XML_DETAIL_H