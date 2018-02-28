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


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF parser. Miscellaneous parser functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ActiveContigMap = std::map<ContigId_t, ContigSize_t>;

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
                         std::vector<std::pair<char, size_t>>& parsed_cigar);

private:

  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  static bool tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                         std::map<std::string, std::string>& key_value_pairs);

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* ID_CIGAR_VALUE_ = "CIGAR";
  constexpr static const char* ID_READ_DEPTH_ = "DPB";
  constexpr static const char* ID_PROPORTION_ = "AF";


};


}   // namespace genome
}   // namespace kellerberrin








#endif //KGL_VARIANT_FACTORY_VCF_PARSE_IMPL_H
