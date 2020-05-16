//
// Created by kellerberrin on 10/5/20.
//

#ifndef KGL_VARIANT_FACTORY_PARSE_INFO_H
#define KGL_VARIANT_FACTORY_PARSE_INFO_H

#include "kel_exec_env.h"

#include <string>
#include <array>


namespace kellerberrin::genome {   //  organization level namespace




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An efficient single pass parser that uses std::string_view for CPU efficiency.
// In general std::string_view is work of the devil and a seg fault waiting to happen.
// But if the underlying string has the same guaranteed lifetime as the associated std::string_view(s)
// then a seg fault may not be inevitable.
// This is guaranteed by enclosing the parsed std::string_views and source INFO string in the same object.
using InfoParserToken = std::pair<std::string_view, size_t>;
using InfoParserMap = std::map<std::string_view, InfoParserToken>;
class VCFInfoParser {

public:

  // std::move the info string into this object for efficiency.
  explicit VCFInfoParser(std::string&& info) : info_(std::move(info)), info_view_(info_) {

    if (not parseInfo()) {

      ExecEnv::log().error("VCFInfoParser::VCFInfoParser, Problem parsing info field");

    }

  }
  ~VCFInfoParser() = default;

  [[nodiscard]] const InfoParserMap& getMap() const { return parsed_token_map_; }
  [[nodiscard]] const std::string& info() const { return info_; }
  [[nodiscard]] std::optional<std::string> getInfoString(const std::string& key) const;
  [[nodiscard]] std::optional<std::vector<std::string>> getInfoStringArray(const std::string& key) const;
  [[nodiscard]] std::optional<int64_t> getInfoInteger(const std::string& key) const;
  [[nodiscard]] std::optional<std::vector<std::optional<int64_t>>> getInfoIntegerArray(const std::string& key) const;
  [[nodiscard]] std::optional<double> getInfoFloat(const std::string& key) const;
  [[nodiscard]] std::optional<std::vector<std::optional<double>>> getInfoFloatArray(const std::string& key) const;


private:

  std::string info_;
  std::string_view info_view_;
  InfoParserMap parsed_token_map_;

  constexpr static const char INFO_FIELD_DELIMITER_{';'};
  constexpr static const char INFO_VALUE_DELIMITER_{'='};
  constexpr static const char INFO_VECTOR_DELIMITER_{','};
  constexpr static const char* INFO_VECTOR_DELIMITER_STR_{","};
  constexpr static const char* INFO_VECTOR_MISSING_VALUE_STR_{"."};

  [[nodiscard]] bool parseInfo();
  [[nodiscard]] static std::optional<int64_t> convertToInteger(const std::string& key, const std::string& value);
  [[nodiscard]] static std::optional<double> convertToFloat(const std::string& key, const std::string& value);


};



} // namespace




#endif //KGL_KGL_VARIANT_FACTORY_GRCH_INFO_H
