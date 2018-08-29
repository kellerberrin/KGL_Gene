//
// Created by kellerberrin on 28/08/18.
//

#ifndef KGL_PHYLOGENETIC_OPTION_H
#define KGL_PHYLOGENETIC_OPTION_H

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The requirement for program control and input options for KGL has exceeded the ability to be comfortably
// implemented as command line options.
//
// These are now implemented in an options file. The file is specified by the command line (mandatory) and is
// relative to the work directory (generally in a sub directory).
//
// The specification of the text options file is simple:
// On each line the first field is the option name and optional subsequent fields are the option values.
// Option Values are separated by commas ",", "=" or tabs "\t", other whitespace is considered to be part of the option.
// If there are no option values then the option is considered to be boolean (a flag).
// If the first character is '#' then the line is a comment and is ignored by the parser.
// It is good practise to document all the KGL options using comments.
// Undefined options will generate an error and KGL will terminate.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <map>
#include <array>



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object defines allowed option types. These are defined in the implementation file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class OptionOptional { REQUIRED, OPTIONAL };
enum class OptionArrayType { ARRAY, SINGULAR, FLAG };  // If FLAG is specified then no argument is expected.
enum class OptionArgumentType { STRING, FLOAT, INTEGER };

// If the option is an array then the "option_defaults" arguments must be separated using a valid field separator char
// defined below as one of the chars in RuntimeOptions::FIELD_DELIMITER_ (private).
// e.g. both "2,3,4,5" and "2=3=4=5" define the same array of 4 integers.

struct PredefinedOptionType {

  const char* name;
  OptionOptional optional;
  OptionArrayType option_array_type;
  OptionArgumentType option_arument_type;
  const char* option_defaults;
  const char* option_description;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Reads the actual option values from the specified option file, one option and associated arguments per line,
// and confirms that the options conform to the PredefinedOptionType specification defined above.
// All valid options must be predefined.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
using RuntimeOptionTag = std::string;
using RuntimeOptionVector = std::vector<std::string>;
using RuntimOptionsMap = std::map<RuntimeOptionTag , RuntimeOptionVector>; // All the options found in an options file.
class RuntimeOptions {

public:

  RuntimeOptions() = default;
  ~RuntimeOptions() = default;

  bool readRuntimeOptions(const std::string& file_name);
  bool getMixtureFile(std::string& file_name) const { return getRuntimeOption(MIXTURE_FILE_, file_name); }
  static void printHelp(std::ostream& stream);

private:

  RuntimOptionsMap runtime_option_map_;
  static PredefinedOptionType defined_options_[];  // defined in the implementation file.

  constexpr static const char* FIELD_DELIMITER_ = "=,\t"; // Valid separators
  constexpr static const char COMMENT_CHAR_ = '#';

  // Option names.
  constexpr static const char* MIXTURE_FILE_ = "mixtureFile";
  constexpr static const char* VCF_PLOIDY_ = "vcf_ploidy";

  const RuntimOptionsMap& getMap() const { return runtime_option_map_; }
  // Returns false if no genome found, mixture_statistics is zeroed.
  bool getRuntimeOptionArray(const RuntimeOptionTag& option, RuntimeOptionVector& option_values) const;
  bool getRuntimeOption(const RuntimeOptionTag& option, std::string& option_value) const;
  // Read the mixture file. Can be more than 1 char used as field delimiter.
  bool checkRequiredOptions();
  bool getDefinedOption(const RuntimeOptionTag& option, PredefinedOptionType& defined_option);
  bool parseOption(const RuntimeOptionTag& option, const RuntimeOptionVector& arguments);

};



}   // namespace genome
}   // namespace kellerberrin





#endif //KGL_PHYLOGENETIC_OPTION_H
