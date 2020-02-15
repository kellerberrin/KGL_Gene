//
// Created by kellerberrin on 28/02/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_PARSE_IMPL_H
#define KGL_VARIANT_FACTORY_VCF_PARSE_IMPL_H



#include <map>

#include "kgl_utility.h"
#include "kgl_genome_types.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db.h"


#include <seqan/vcf_io.h>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF parser. Miscellaneous parser functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returned from the cigar functions.
enum class CigarEditType : char { UNCHANGED = 'M', INSERT = 'I', DELETE = 'D', CHANGED = 'X'};
using CigarEditItem = std::pair<size_t, CigarEditType>; // Used to specify edit as vector '1M,1X,3D,3I'.
using CigarVector = std::vector<CigarEditItem>;

using ActiveContigMap = std::map<ContigId_t, ContigSize_t>;
using VcfInfoKeyMap = std::map<std::string, std::string>;

class ParseVCFMiscImpl {

public:

  ParseVCFMiscImpl() = default;
  ~ParseVCFMiscImpl() = default;


  static bool parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                             const seqan::VcfHeader& header,
                             ActiveContigMap& active_contig_map,
                             bool cigar_required);

// assumes input "key_1=value_1; ...;key_n=value_n"
  static bool tokenizeVcfInfoKeyValues(const std::string& key_value_text,
                                       std::map<std::string, std::string>& key_value_map);

  static bool parseCigar(const std::string& cigar,
                         size_t& check_reference_size,
                         size_t& check_alternate_size,
                         std::vector<CigarEditItem>& parsed_cigar);

// tokenize a string
  static bool tokenize(const std::string& parse_text, const std::string& separator_text, std::vector<std::string>& item_vector);

  // Generate a CIGAR from two sequences.
  static std::string generateCigar(const std::string& reference, const std::string& alternate);

  // Generate a CIGAR from a cigar vector
  static std::string generateCigar(const CigarVector& cigar_vector);

// Use edlib to generate a cigar vector.
  static CigarVector generateEditVector(const std::string& reference, const std::string& alternate);

  // Given a reference count and a cigar vector compute a number that calculates the equivalent
  // size of the alternate string.
  // For UNCHANGED = 'M' and CHANGED = 'X' cigar items the reference_count and alternate count are incremented.
  // For INSERT = 'I' the alternate is incremented and the reference_count is not.
  // For DELETE = 'D' the reference count is incremented and the alternate is not.
  static size_t alternateCount(size_t reference_count, const CigarVector& cigar_vector);


private:

  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  static bool tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                         VcfInfoKeyMap& key_value_pairs);

// Use edlib to generate a cigar string.
  static void generateEditString(const std::string& reference, const std::string& alternate, std::vector<CigarEditType>& edit_vector);

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* ID_CIGAR_VALUE_ = "CIGAR";
  constexpr static const char* ID_READ_DEPTH_ = "DPB";
  constexpr static const char* ID_PROPORTION_ = "AF";


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Info field parser.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VCFInfoField {


public:

  VCFInfoField(const std::string& key_value_text) {

    if (not ParseVCFMiscImpl::tokenizeVcfInfoKeyValues(key_value_text, info_key_value_map_)) {

      ExecEnv::log().error("Unable to parse VCF INFO: {}", key_value_text);

    }

  }
  virtual ~VCFInfoField() = default;

  bool getInfoField(const std::string& key, std::string& value) const {

    auto key_it = info_key_value_map_.find(key);

    if (key_it != info_key_value_map_.end()) {

      value = key_it->second;
      return true;

    }

    value = "";
    return false;

  }

  const VcfInfoKeyMap& getMap() const { return info_key_value_map_; }

private:

  VcfInfoKeyMap info_key_value_map_;

};


}   // end namespace








#endif //KGL_VARIANT_FACTORY_VCF_PARSE_IMPL_H
