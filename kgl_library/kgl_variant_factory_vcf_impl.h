//
// Created by kellerberrin on 22/01/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_VCF_IMPL_H


#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_single.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <seqan/vcf_io.h>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF (freebayes) parser. Low-level implementation. Do not include this file in any source files except the following:
// kgl_variant_factory_vcf_fbimpl.cc
// kgl_variant_factory_vcf_gatkimpl.cc
// kgl_variant_factory_vcf_impl.cc
// VcfFactory::FreeBayesVCFImpl does all the heavy lifting using the 3rd party libraries, seqan and boost.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ActiveContigMap = std::map<kgl::ContigId_t, kgl::ContigSize_t>;

class kgl::VcfFactory::ParseVCFImpl {

public:

  ParseVCFImpl() = default;
  virtual ~ParseVCFImpl() = default;


protected:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 5000;

  size_t vcf_record_count_;
  size_t vcf_record_ignored_;
  size_t vcf_record_error_;
  size_t vcf_record_rejected_;
  size_t vcf_variant_count_;

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* ID_CIGAR_VALUE_ = "CIGAR";
  constexpr static const char* ID_READ_DEPTH_ = "DPB";
  constexpr static const char* ID_PROPORTION_ = "AF";

  bool parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                      const seqan::VcfHeader& header,
                      ActiveContigMap& active_contig_map,
                      bool cigar_required) const;

  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  bool tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                  std::map<std::string, std::string>& key_value_pairs) const;

// assumes input "key_1=value_1; ...;key_n=value_n"
  bool tokenizeVcfInfoKeyValues(const std::string& key_value_text,
                                std::map<std::string, std::string>& key_value_map) const;

  bool parseCigar(const std::string& cigar,
                  size_t& check_reference_size,
                  size_t& check_alternate_size,
                  std::vector<std::pair<char, size_t>>& parsed_cigar) const;

};



#endif //KGL_KGL_VARIANT_FACTORY_VCF_IMPL_H
