//
// Created by kellerberrin on 23/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_H



#include "kel_exec_env.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"

#include <string>
#include <array>


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The InfoEvidenceType object defines the Info data subscriber type (alternate allele etc), and the internal and external representation of the data.

// The data subscriber type.
enum class InfoEvidenceSubscriber {
  GeneralScalar,   // General single valued data.
  GeneralArray,    // General array valued data
  AlternateAllele,  // Applies to the Alternate Alleles. Can be a single value for all alt. alleles (assumed) or a vector of values for each alt. allele.
  AllAllele, // Applies to all Alleles including the reference allele. Can be a single value for all alleles (assumed) or a vector of values for each allele.
  Genotype,  // Applies to all defined genotypes in the VCF file (if any). Can be a single value (assumed) or a vector of values for each genotype.
  NotImplemented  // Unknown (number, type) tuple.
};

// Defines how the data is represented externally.
enum class InfoEvidenceExtern : size_t {
  Boolean = 0,
  Integer = 1,
  Float = 2,
  String = 3,
};

// Defines how the data is represented internally.
// For size and fragmentation efficiency there are only 3 raw vectors of chars, floats and integers available internally.
enum class InfoEvidenceIntern {
  intern_char, // Size is fixed (boolean variables) and known at Info data subscription time.
  intern_integer, // Size is fixed and known at Info data subscription time.
  intern_float, // Size is fixed and known at Info data subscription time.
  intern_string, // Size varies between records.
  intern_integer_array, // Size is fixed and known at Info data subscription time.
  intern_float_array, // Size is fixed and known at Info data subscription time.
  intern_string_array,  // Size varies between records.
  intern_unity_integer_array,   // Assumed to be single (scalar) value, but can be a vector.  Assumed size 1, actual size known after Info data parsed.
  intern_unity_float_array,   // Assumed to be single (scalar) value, but can be a vector.  Assumed size 1, actual size known after Info data parsed.
  intern_unity_string_array,    // Assumed to be single (scalar) value, but can be a vector. Assumed size 1, actual size known after Info data parsed.
  NotImplemented  // Unknown (number, type) tuple.
};



class InfoEvidenceType {

public:

  InfoEvidenceType(const InfoEvidenceSubscriber type_application,
                   const InfoEvidenceExtern extern_representation,
                   const InfoEvidenceIntern intern_representation)
  : subscriber_type_(type_application)
  , extern_representation_(extern_representation)
  , intern_representation_(intern_representation) {}
  InfoEvidenceType(const InfoEvidenceType&) = default;
  ~InfoEvidenceType() = default;


  [[nodiscard]] InfoEvidenceSubscriber subscriberType() const { return subscriber_type_; }  // The info number field
  [[nodiscard]] InfoEvidenceExtern ExternalInfoType() const { return extern_representation_; }     // The info type field
  [[nodiscard]] InfoEvidenceIntern InternalInfoType() const { return intern_representation_; }     // The info type field

  // The size of the data type does not change for different Info records
  // e.g. A scalar float (true) which is fixed, unlike a string (false) which will always vary between records.


private:

  const InfoEvidenceSubscriber subscriber_type_;       // The info number field
  const InfoEvidenceExtern extern_representation_;     // How the data is presented externally
  const InfoEvidenceIntern intern_representation_;     // How the data is represented internally

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object that defines Info data attributes by examining the "Number" and "Type" fields and doing a table lookup.
// returns an InfoEvidenceType object defined above.

// Never created.
class InfoTypeLookup {

public:

  static InfoEvidenceType evidenceType(const VCFInfoRecord &vcf_info_item);

private:

  InfoTypeLookup() = default;
  ~InfoTypeLookup() = default;

  // Text constants used in the VCF header. See the VCF 4.2 file specification.
  constexpr static const char *INTEGER_ = "Integer";
  constexpr static const char *FLOAT_ = "Float";
  constexpr static const char *FLAG_ = "Flag";
  constexpr static const char *CHAR_STRING_ = "Character";
  constexpr static const char *STRING_ = "String";
  constexpr static const char *SCALAR_ = "1";
  constexpr static const char *AlTERNATIVE_ALLELE_ = "A";
  constexpr static const char *ALL_ALLELE_ = "R";
  constexpr static const char *ALL_GENOTYPES_ = "G";
  constexpr static const char *FLAG_SCALAR_ = "0";
  constexpr static const char *INDETERMINATE_COUNT_ = ".";

  // The Info data type lookup table. THis takes the "Number" and "Type" fields from the VCF record and looks up
  // the data subscriber (all alleles, all genotypes, alternate alleles etc.), the external representation of the data
  // as a scalar or a vector. The internal data representation is designed for maximum data storage size efficiency. All
  // floats are indexed into a single float vector, bools and strings are allocated from a single char vector, and integers
  // are indexed from a single integer vector. See the InfoDataBlock for more information.
  // In addition Genotype ('G'), Alternate Allele ('A'), All Allele ('R') subscriber data is
  // assumed to be scalar for size efficiency, however they can also be stored as vectors if the parsed data returns a vector.
  using InfoTypeLookupTable = std::pair<InfoEvidenceType, std::function<bool(std::string number, std::string type)>>;
  static const std::vector<InfoTypeLookupTable> type_definitions_;
  static bool isVectorType(const std::string& type); // returns true if a number > 1 used in the lambdas below.

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Table to define the INFO field data types. An associated lambda does the actual lookup.
inline const std::vector<InfoTypeLookup::InfoTypeLookupTable> InfoTypeLookup::type_definitions_ = {

// Boolean types (allocated on the char heap). Assume no boolean vector type exists.
{         {InfoEvidenceSubscriber::GeneralScalar,   InfoEvidenceExtern::Boolean, InfoEvidenceIntern::intern_char},
[](const std::string& type, const std::string& number) { return type == FLAG_ and number == FLAG_SCALAR_; }    },

{         {InfoEvidenceSubscriber::AlternateAllele, InfoEvidenceExtern::Boolean, InfoEvidenceIntern::intern_char},
[](const std::string& type, const std::string& number) { return type == FLAG_ and number == AlTERNATIVE_ALLELE_; }    },

{        {InfoEvidenceSubscriber::AllAllele,       InfoEvidenceExtern::Boolean, InfoEvidenceIntern::intern_char},
[](const std::string& type, const std::string& number) { return type == FLAG_ and number == ALL_ALLELE_; }      },

{        {InfoEvidenceSubscriber::Genotype,        InfoEvidenceExtern::Boolean, InfoEvidenceIntern::intern_char},
[](const std::string& type, const std::string& number) { return type == FLAG_ and number == ALL_GENOTYPES_; }      },

// String types.
{        {InfoEvidenceSubscriber::GeneralScalar,   InfoEvidenceExtern::String, InfoEvidenceIntern::intern_string},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == SCALAR_; }       },

{        {InfoEvidenceSubscriber::AlternateAllele, InfoEvidenceExtern::String, InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == AlTERNATIVE_ALLELE_; }     },

{        {InfoEvidenceSubscriber::AllAllele,       InfoEvidenceExtern::String, InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == ALL_ALLELE_; }      },

{       {InfoEvidenceSubscriber::Genotype,        InfoEvidenceExtern::String,  InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == ALL_GENOTYPES_; }     },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::String,  InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == INDETERMINATE_COUNT_; }      },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::String,  InfoEvidenceIntern::intern_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and isVectorType(number); }      },

// Float types.
{       {InfoEvidenceSubscriber::GeneralScalar,   InfoEvidenceExtern::Float, InfoEvidenceIntern::intern_float},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == SCALAR_; }     },

{       {InfoEvidenceSubscriber::AlternateAllele, InfoEvidenceExtern::Float, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == AlTERNATIVE_ALLELE_; }    },

{       {InfoEvidenceSubscriber::AllAllele,       InfoEvidenceExtern::Float, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == ALL_ALLELE_; }     },

{       {InfoEvidenceSubscriber::Genotype,        InfoEvidenceExtern::Float, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == ALL_GENOTYPES_; }     },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::Float, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == INDETERMINATE_COUNT_; }      },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::Float, InfoEvidenceIntern::intern_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and isVectorType(number); }      },

// Integer types.
{      {InfoEvidenceSubscriber::GeneralScalar,   InfoEvidenceExtern::Integer, InfoEvidenceIntern::intern_integer},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == SCALAR_; }      },

{      {InfoEvidenceSubscriber::AlternateAllele, InfoEvidenceExtern::Integer, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == AlTERNATIVE_ALLELE_; }    },

{      {InfoEvidenceSubscriber::AllAllele,       InfoEvidenceExtern::Integer, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == ALL_ALLELE_; }     },

{      {InfoEvidenceSubscriber::Genotype,        InfoEvidenceExtern::Integer, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == ALL_GENOTYPES_; }     },

{      {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::Integer, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == INDETERMINATE_COUNT_; }     },

{      {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::Integer, InfoEvidenceIntern::intern_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and isVectorType(number); }      },

};


} // namespace



#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_H
