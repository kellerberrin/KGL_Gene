//
// Created by kellerberrin on 22/01/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_VCF_IMPL_H


#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_readvcf_impl.h"

#include <seqan/vcf_io.h>

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF parser. Multi threaded implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ActiveContigMap = std::map<ContigId_t, ContigSize_t>;

class ParseVCFImpl {

public:

  ParseVCFImpl(std::shared_ptr<PopulationVariant> pop_variant_ptr,
               std::shared_ptr<const GenomeDatabase> genome_db_ptr,
               const std::string& vcf_file_name,
               Phred_t variant_quality) : pop_variant_ptr_(pop_variant_ptr),
                                          genome_db_ptr_(genome_db_ptr),
                                          vcf_file_name_(vcf_file_name),
                                          variant_quality_(variant_quality) {

    reader_ptr_ = std::make_shared<VCFReaderMT<ParseVCFImpl>>(vcf_file_name, this, &ParseVCFImpl::ProcessVCFRecord);

  }


  virtual ~ParseVCFImpl() = default;

  virtual void readParseVCFImpl();

  virtual void ProcessVCFRecord(const seqan::VcfRecord& record_ptr) = 0;

  const std::vector<std::string>& getGenomeNames() const { return reader_ptr_->getFieldNames(); }

protected:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 10000;

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

  size_t addThreadSafeGenomeVariant(std::shared_ptr<GenomeVariant> genome_variants,
                                    std::shared_ptr<const Variant> variant_ptr) const;

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

  std::shared_ptr<PopulationVariant> pop_variant_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  const std::string vcf_file_name_;
  Phred_t variant_quality_;

  std::shared_ptr<VCFReaderMT<ParseVCFImpl>> reader_ptr_;

private:

  mutable std::mutex mutex_;  // mutex to lock the GenomeVariant structure when inserting variants.

};




}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_KGL_VARIANT_FACTORY_VCF_IMPL_H
