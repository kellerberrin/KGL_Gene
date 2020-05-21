//
// Created by kellerberrin on 16/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_H


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
enum class InfoEvidenceExtern {
  Boolean,
  Integer,
  Float,
  String,
  IntegerArray,
  FloatArray,
  StringArray,
  NotImplemented  // Unknown (number, type) tuple.
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
  bool fixedDataType() const;


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

{        {InfoEvidenceSubscriber::AlternateAllele, InfoEvidenceExtern::StringArray, InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == AlTERNATIVE_ALLELE_; }     },

{        {InfoEvidenceSubscriber::AllAllele,       InfoEvidenceExtern::StringArray, InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == ALL_ALLELE_; }      },

{       {InfoEvidenceSubscriber::Genotype,        InfoEvidenceExtern::StringArray, InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == ALL_GENOTYPES_; }     },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::StringArray, InfoEvidenceIntern::intern_unity_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and number == INDETERMINATE_COUNT_; }      },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::StringArray, InfoEvidenceIntern::intern_string_array},
[](const std::string& type, const std::string& number) { return type == STRING_ and isVectorType(number); }      },

// Float types.
{       {InfoEvidenceSubscriber::GeneralScalar,   InfoEvidenceExtern::Float, InfoEvidenceIntern::intern_float},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == SCALAR_; }     },

{       {InfoEvidenceSubscriber::AlternateAllele, InfoEvidenceExtern::FloatArray, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == AlTERNATIVE_ALLELE_; }    },

{       {InfoEvidenceSubscriber::AllAllele,       InfoEvidenceExtern::FloatArray, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == ALL_ALLELE_; }     },

{       {InfoEvidenceSubscriber::Genotype,        InfoEvidenceExtern::FloatArray, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == ALL_GENOTYPES_; }     },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::FloatArray, InfoEvidenceIntern::intern_unity_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and number == INDETERMINATE_COUNT_; }      },

{       {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::FloatArray, InfoEvidenceIntern::intern_float_array},
[](const std::string& type, const std::string& number) { return type == FLOAT_ and isVectorType(number); }      },

// Integer types.
{      {InfoEvidenceSubscriber::GeneralScalar,   InfoEvidenceExtern::Integer, InfoEvidenceIntern::intern_integer},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == SCALAR_; }      },

{      {InfoEvidenceSubscriber::AlternateAllele, InfoEvidenceExtern::IntegerArray, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == AlTERNATIVE_ALLELE_; }    },

{      {InfoEvidenceSubscriber::AllAllele,       InfoEvidenceExtern::IntegerArray, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == ALL_ALLELE_; }     },

{      {InfoEvidenceSubscriber::Genotype,        InfoEvidenceExtern::IntegerArray, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == ALL_GENOTYPES_; }     },

{      {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::IntegerArray, InfoEvidenceIntern::intern_unity_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and number == INDETERMINATE_COUNT_; }     },

{      {InfoEvidenceSubscriber::GeneralArray,    InfoEvidenceExtern::IntegerArray, InfoEvidenceIntern::intern_integer_array},
[](const std::string& type, const std::string& number) { return type == INTEGER_ and isVectorType(number); }      },

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The data object contains only 5 data structures
// 1. A vector of Floats.
// 2. A vector of ints.
// 3. A vector of chars.
// 4  A vector of std::string_view to index into the char vector to define strings.
// 5. A vector of { offset, size, type } to define integer or float vectors at run time when parsed info is presented.
// If we assume that AlternateAllele, AllAllele, Genotype, only contain scalars then 1 to 5 are known at subscribe time.
// These indexes are held for each subscribed variable below.
//
// If, however, the AlternateAllele, AllAllele, Genotype subscribers are vectors then 4. and 5. are calculated when the parsed INFO
// data is presented to the InfoDataBlock object and the data vectors for 1. to 5. are adjusted accordingly.

// To simplify initial coding, the AlternateAllele, AllAllele, Genotype are assumed to be scalars.
// The more complex vector case will be added later.

// Internal variable index structure (only uses 64 bits for size efficiency).
struct InfoDataIndex {

public:

  InfoDataIndex() = default;
  ~InfoDataIndex() = default;

  [[nodiscard]] size_t infoVariableIndex() const { return static_cast<size_t>(info_variable_index_); }
  [[nodiscard]] size_t infoOffset() const { return static_cast<size_t>(info_element_offset_); }
  [[nodiscard]] size_t infoSize() const { return static_cast<size_t>(info_element_count_); }

  void infoVariableIndex(size_t info_variable_index) { info_variable_index_ = static_cast<VariableIndexImpl>(info_variable_index); }
  void infoOffset(size_t info_element_offset)  { info_element_offset_ = static_cast<InfoOffsetImpl>(info_element_offset); }
  void infoSize(size_t info_element_count) { info_element_count_ = static_cast<InfoCountImpl>(info_element_count); }


private:

  // Minimize storage size.
  using VariableIndexImpl = uint16_t;
  using InfoOffsetImpl = uint32_t;
  using InfoCountImpl = uint16_t;

  VariableIndexImpl info_variable_index_{0};  // corresponds to the variable index in the header, sort key for the array block.
  InfoOffsetImpl info_element_offset_{0};
  InfoCountImpl info_element_count_{0};

};


class InfoEvidenceHeader; // forward.
struct ItemOffset {

  size_t offset{0};
  bool is_array{false};

};
class InfoSubscribedField {

public:

  InfoSubscribedField(VCFInfoRecord vcfInfoRecord, std::shared_ptr<const InfoEvidenceHeader> info_evidence_header)
  : vcfInfoRecord_(std::move(vcfInfoRecord)),
    type_(InfoTypeLookup::evidenceType(vcfInfoRecord_)),
    info_evidence_header_(std::move(info_evidence_header)) {}

//  InfoSubscribedField(const InfoSubscribedField &) = default;
  ~InfoSubscribedField() = default;

  [[nodiscard]] const VCFInfoRecord &infoRecord() const { return vcfInfoRecord_; }
  [[nodiscard]] InfoEvidenceType evidenceType() const { return type_; }
  [[nodiscard]] ItemOffset dataOffset() const { return data_offset_; }
  void dataOffset(const ItemOffset& item_offset) { data_offset_ = item_offset; }

private:

  const VCFInfoRecord vcfInfoRecord_;  // The original VCF Header Record.
  const InfoEvidenceType type_;  // THe inferred subscriber type, external type and internal type.
  std::shared_ptr<const InfoEvidenceHeader> info_evidence_header_; // Ensure the index knows which header it belongs to.
  ItemOffset data_offset_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the per variant INFO data block.
// Generated from the parser data block by specifying the Allele index.
class EvidenceFactory; // Fwd.
struct DataInfoTypeCount {

public:

  DataInfoTypeCount() = default;
  ~DataInfoTypeCount() = default;

  size_t arrayCount() const { return array_count_; }
  size_t floatCount() const { return float_count_; }
  size_t integerCount() const { return integer_count_; }
  size_t stringCount() const { return string_count_; }
  size_t charCount() const { return char_count_; }

  void arrayCount(size_t count) { array_count_ = count; }
  void floatCount(size_t count) { float_count_ = count; }
  void integerCount(size_t count) { integer_count_ = count; }
  void stringCount(size_t count) { string_count_ = count; }
  void charCount(size_t count) { char_count_ = count; }

  // Notionally allocates data on the 5 data arrays and returns a data block
  // The function is to be used sequentially on all subscribed variables.
  // The final values are used to actually allocate memory in the InfoDataBlock object.
  ItemOffset staticIncrementAndAllocate(InfoEvidenceIntern internal_type);  // Set up the indexes and pre-allocate fixed sized fields (run once).
  void dynamicIncrementAndAllocate(InfoEvidenceIntern internal_type, size_t data_item_count, size_t data_item_size); // Allocate additional vector space at runtime (each VCF record).

private:

  size_t array_count_{0};
  size_t float_count_{0};
  size_t integer_count_{0};
  size_t string_count_{0};
  size_t char_count_{0};

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoSubscribedFields. There is only one of these held by all variants with INFO evidence fields.


using InfoSubscribedMap = std::map<std::string, InfoSubscribedField>;
class InfoEvidenceHeader {

public:

  explicit InfoEvidenceHeader() {}

  InfoEvidenceHeader(const InfoEvidenceHeader &) = delete;

  ~InfoEvidenceHeader() = default;

  std::optional<const InfoSubscribedField> getSubscribedField(const std::string &field_id) const;
  [[nodiscard]] const InfoSubscribedMap &getMap() const { return info_subscribed_map_; }

  // Note that this routine takes a shared_ptr to itself obtained from the info data factory. This is passed onto subscribed field objects.
  [[nodiscard]] bool setupEvidenceHeader(const VCFInfoRecord& vcf_info_record, std::shared_ptr<const InfoEvidenceHeader> self_ptr);
  // Static storage and field offsets.
  void setupStaticStorage();
  const DataInfoTypeCount& staticStorage() const { return static_storage_; }

private:

  InfoSubscribedMap info_subscribed_map_;
  DataInfoTypeCount static_storage_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoDataBlock {

public:


  InfoDataBlock(std::shared_ptr<InfoEvidenceHeader> info_evidence_header) : info_evidence_header_(
  std::move(info_evidence_header)) {}
  InfoDataBlock(const InfoEvidenceHeader &) = delete;
  ~InfoDataBlock() = default;

  void allocateMemory(const DataInfoTypeCount& type_count) {

    char_memory_ = std::make_unique<char[]>(type_count.charCount());
    integer_memory_ = std::make_unique<InfoIntegerType[]>(type_count.integerCount());
    float_memory_ = std::make_unique<InfoFloatType[]>(type_count.floatCount());
    array_memory_ = std::make_unique<InfoDataIndex[]>(type_count.arrayCount());
    string_memory_ = std::make_unique<std::string_view[]>(type_count.stringCount());

  }

private:

  std::shared_ptr<InfoEvidenceHeader> info_evidence_header_; // The data header.

  std::unique_ptr<char[]> char_memory_;
  std::unique_ptr<InfoIntegerType []> integer_memory_;
  std::unique_ptr<InfoFloatType []> float_memory_;
  std::unique_ptr<InfoDataIndex[]> array_memory_;
  std::unique_ptr<std::string_view[]> string_memory_;


};


class InfoDataBlockNaive : public InfoDataBlock {

public:

  friend EvidenceFactory;

  InfoDataBlockNaive(std::shared_ptr<InfoEvidenceHeader> info_evidence_header) : InfoDataBlock(std::move(info_evidence_header)) {}
  InfoDataBlockNaive(const InfoEvidenceHeader &) = delete;
  ~InfoDataBlockNaive() = default;

  DataInfoTypeCount dataPayload();

private:

  // The stored data.
  std::vector<bool> bool_data_;
  std::vector<InfoParserFloat> float_data_;
  std::vector<InfoParserInteger> integer_data_;
  std::vector<InfoParserString> string_data_;
  std::vector<InfoParserFloatArray> float_array_data_;
  std::vector<InfoParserIntegerArray> integer_array_data_;
  std::vector<InfoParserStringArray> string_array_data_;

  void clearAll();

};






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant.

// createVariantEvidence() either returns a single data block, InfoDataBlock for a single alternate allele or
// a MultipleAlleleDataBlock data block for a multiple alternative allele VCF record.
using InfoDataEvidence = std::optional<std::shared_ptr<InfoDataBlock>>;
class EvidenceFactory {

public:

  explicit EvidenceFactory(const EvidenceInfoSet& evidence_set)  // An ordered set of Info field identifiers.
  : evidence_map_(evidence_set),
    info_evidence_header_(std::make_shared<InfoEvidenceHeader>()) {}
  ~EvidenceFactory() = default;

  // A map of all available Info fields parsed from the VCF header
  // This also initializes the InfoEvidenceHeader object.
  void availableInfoFields(const VCFInfoRecordMap& vcf_info_map);
  [[nodiscard]] const VCFInfoRecordMap& availableInfoFields() const { return all_available_map_; }
  // For each input VCFRecord info text field (std::moved for efficiency), create a parsed data object.
  [[nodiscard]] InfoDataEvidence createVariantEvidence(std::string&& info);


private:

  // The evidence fields specified in the runtime XML file.
  const EvidenceInfoSet evidence_map_;
  // All available info fields.
  VCFInfoRecordMap all_available_map_;
  // The Info header block, contains all the subscribed field information.
  std::shared_ptr<InfoEvidenceHeader> info_evidence_header_;

  // ******** Temp.
  std::unique_ptr<InfoDataBlockNaive> parseSubscribed(std::string&& info);
  DataInfoTypeCount parseSubscribed_alt(std::string&& info);

  // If the user specifies just specifies "None" (case insensitive) then no Info fields will be subscribed.
  constexpr static const char *NO_FIELD_SUBSCRIBED_ = "NONE";

};



} // namespace





#endif //KGL_VARIANT_FACTORY_VCF_EVIDENCE_H
