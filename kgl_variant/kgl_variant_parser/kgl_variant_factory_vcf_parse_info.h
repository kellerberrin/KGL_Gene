//
// Created by kellerberrin on 10/5/20.
//

#ifndef KGL_VARIANT_FACTORY_PARSE_INFO_H
#define KGL_VARIANT_FACTORY_PARSE_INFO_H

#include "kel_exec_env.h"


#include <map>
#include <string>
#include <optional>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An efficient single pass parser that uses std::string_view for CPU efficiency.
// In general std::string_view is work of the devil and a seg fault waiting to happen.
// But if the underlying string has the same guaranteed lifetime as the associated std::string_view(s)
// then a seg fault may not be inevitable.
// This is guaranteed by enclosing the parsed std::string_views and source INFO string in the same object below.
// For CPU efficiency reasons, there are two implementations of this object, see below.

//#define USE_64BIT_TYPES 1

#ifdef USE_64BIT_TYPES

using InfoIntegerType = int64_t;
using InfoFloatType = double;

#else

using InfoIntegerType = int32_t;
using InfoFloatType = float;

#endif


using InfoParserToken = std::pair<std::string_view, size_t>;      // Tokens are a string_view and a field number (number of ',')
using InfoParserMap = std::map<std::string_view, InfoParserToken>; // The token parser data structure.


class VCFInfoParser {

public:

  // std::move the info string into this object for efficiency.
  explicit VCFInfoParser(std::string&& info) : info_(std::move(info)), info_view_(info_) {

    if (not infoTokenParser()) {

      ExecEnv::log().error("VCFInfoParser::VCFInfoParser, Token parser had a problem parsing info field");

    }

  }
  ~VCFInfoParser() = default;

  [[nodiscard]] std::optional<InfoParserToken> getToken(const std::string& key) const;

  [[nodiscard]] static InfoIntegerType convertToInteger(const std::string& value);
  [[nodiscard]] static InfoFloatType convertToFloat(const std::string& value);

  constexpr static const char INFO_VECTOR_DELIMITER_{','};

  // Largest negative values are interpreted as missing values.
  constexpr static const InfoFloatType MISSING_VALUE_FLOAT_ = std::numeric_limits<InfoFloatType>::lowest();
  constexpr static const InfoIntegerType MISSING_VALUE_INTEGER_ = std::numeric_limits<InfoIntegerType>::lowest();

private:

  const std::string info_; // The unparsed 'raw' VCF info record.
  const std::string_view info_view_;  // This a string_view of the info_ string above
  InfoParserMap parsed_token_map_;  // The parsed token map.

  constexpr static const char INFO_FIELD_DELIMITER_{';'};
  constexpr static const char INFO_VALUE_DELIMITER_{'='};
  constexpr static const char* INFO_VECTOR_MISSING_VALUE_STR_{"."};

  [[nodiscard]] bool infoTokenParser();

};


} // namespace


#endif //KGL_KGL_VARIANT_FACTORY_GRCH_INFO_H
