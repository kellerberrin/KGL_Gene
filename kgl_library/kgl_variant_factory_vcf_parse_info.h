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

using InfoParserBoolean = char;
using InfoParserString = std::string;
using InfoParserStringArray = std::vector<std::string>;
using InfoParserInteger = InfoIntegerType;
using InfoParserIntegerArray = std::vector<InfoIntegerType>;
using InfoParserFloat = InfoFloatType;
using InfoParserFloatArray = std::vector<InfoFloatType>;

using InfoParserToken = std::pair<std::string_view, size_t>;
using InfoParserMap = std::map<std::string_view, InfoParserToken>; // The token parser data structure.

using InfoParserArrayToken = std::vector<std::string_view>;
using InfoParserArrayMap = std::map<std::string_view, InfoParserArrayToken>;   // The array parser data structure.

/// This object is CPU time sensitive so two implementations (parsers) have been devised and tested for CPU efficiency.
/// Normally these would have their interfaces defined in a base class but we wish to avoid the overhead of
/// abstract function calls. So we switch between these parser implementation using a compile time definition.
/// This will allow for time benchmarking and for further refinement of the parser object.
/// The 'token' parser appears to be most time efficient because it makes an optimization based the fact the most
/// vector fields are unitary. This is based on testing of the Homo Sapien data. However the alternative implementation,
/// although slightly slower is completely general.

// Comment for the 'array' parser implementation, for the the slightly faster 'token' parser.
#define USE_TOKEN_PARSER 1

class VCFInfoParser {

public:

  // std::move the info string into this object for efficiency.
  explicit VCFInfoParser(std::string&& info) : info_(std::move(info)), info_view_(info_) {

#ifdef USE_TOKEN_PARSER

    if (not infoTokenParser()) {

      ExecEnv::log().error("VCFInfoParser::VCFInfoParser, Token parser had a problem parsing info field");

    }

#else

    if (not infoArrayParser()) {

      ExecEnv::log().error("VCFInfoParser::VCFInfoParser, Array parser had a problem parsing info field");

    }

#endif


  }
  ~VCFInfoParser() = default;

  [[nodiscard]] const std::string& info() const { return info_; }
  [[nodiscard]] std::optional<InfoParserToken> getToken(const std::string& key) const;

  [[nodiscard]] bool getInfoBoolean(const std::string& key) const;
  [[nodiscard]] InfoParserString getInfoString(const std::string& key) const;
  [[nodiscard]] InfoParserStringArray getInfoStringArray(const std::string& key) const;
  [[nodiscard]] InfoParserInteger getInfoInteger(const std::string& key) const;
  [[nodiscard]] InfoParserIntegerArray getInfoIntegerArray(const std::string& key) const;
  [[nodiscard]] InfoParserFloat getInfoFloat(const std::string& key) const;
  [[nodiscard]] InfoParserFloatArray getInfoFloatArray(const std::string& key) const;

  [[nodiscard]] static InfoParserInteger convertToInteger(const std::string& value);
  [[nodiscard]] static InfoParserFloat convertToFloat(const std::string& value);

private:

  const std::string info_; // The unparsed 'raw' VCF info record.
  const std::string_view info_view_;  // This a view of the info_ string above
  InfoParserMap parsed_token_map_;
  InfoParserArrayMap parsed_array_map_;

  constexpr static const char INFO_FIELD_DELIMITER_{';'};
  constexpr static const char INFO_VALUE_DELIMITER_{'='};
  constexpr static const char INFO_VECTOR_DELIMITER_{','};
  constexpr static const char* INFO_VECTOR_MISSING_VALUE_STR_{"."};

  [[nodiscard]] bool infoArrayParser();
  [[nodiscard]] bool infoTokenParser();

  bool compareParsers();

  // Largest negative values are interpreted as missing values.
  constexpr static const InfoParserFloat MISSING_VALUE_FLOAT = std::numeric_limits<InfoParserFloat>::lowest();
  constexpr static const InfoParserInteger MISSING_VALUE_INTEGER = std::numeric_limits<InfoParserInteger>::lowest();
  // Empty string is missing.
  constexpr static const char* MISSING_VALUE_STRING = "";
  // Missing arrays are just empty arrays.
  inline const static InfoParserStringArray MISSING_STRING_VECTOR;
  inline const static InfoParserIntegerArray MISSING_INTEGER_VECTOR;
  inline const static InfoParserFloatArray MISSING_FLOAT_VECTOR;
  // If we use std::optional.
  constexpr static const std::nullopt_t MISSING_VALUE_OPTIONAL = std::nullopt;


};



} // namespace




#endif //KGL_KGL_VARIANT_FACTORY_GRCH_INFO_H
