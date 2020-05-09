//
// Created by kellerberrin on 20/4/20.
//

#ifndef KGL_VARIANT_FACTORY_GRCH_IMPL_H
#define KGL_VARIANT_FACTORY_GRCH_IMPL_H




#include "kel_utility.h"
#include "kgl_variant_file_impl.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_record_vcf_impl.h"



namespace kellerberrin::genome {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is an efficient single pass parser.
// In general std::string_view is work of the devil and a seg fault waiting to happen.
// But if the underlying string has the same guaranteed lifetime as the associated std::string_view then a seg fault
// may not be inevitable.
using InfoParserMap = std::map<std::string_view, std::string_view>;
class GrchInfoParser {

public:

  // std::move the info string into this object.
  explicit GrchInfoParser(std::string info) : info_(std::move(info)), info_view_(info_) {}
  ~GrchInfoParser() = default;

  [[nodiscard]] bool parseInfo();
  [[nodiscard]] const InfoParserMap& getMap() const { return parsed_token_map_; }
  [[nodiscard]] const std::string& info() const { return info_; }

private:

  std::string info_;
  std::string_view info_view_;

  InfoParserMap parsed_token_map_;

  constexpr static const char INFO_FIELD_DELIMITER_{';'};
  constexpr static const char INFO_VALUE_DELIMITER_{'='};

};



class GrchVCFImpl : public VCFReaderMT {

public:

  GrchVCFImpl(std::shared_ptr<UnphasedGenome> vcf_genome_ptr,
              std::shared_ptr<const GenomeReference> genome_db_ptr,
              const std::string &vcf_file_name,
              const ContigAliasMap& contig_alias_map) : VCFReaderMT(vcf_file_name),
                                                        vcf_genome_ptr_(std::move(vcf_genome_ptr)),
                                                        genome_db_ptr_(std::move(genome_db_ptr)),
                                                        contig_alias_map_(contig_alias_map) {}

  ~GrchVCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord &vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo &header_info) override;

  void readParseVCFImpl();


private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_{100000};
  constexpr static const char MULIPLE_ALT_SEPARATOR_{','};
  const std::string alt_separator_ = {MULIPLE_ALT_SEPARATOR_};

  std::shared_ptr<UnphasedGenome> vcf_genome_ptr_;
  std::shared_ptr<const GenomeReference> genome_db_ptr_;
  ContigAliasMap contig_alias_map_;

// Progress counters.

  size_t variant_count_{0};
  mutable std::mutex add_variant_mutex_;

  std::map<ContigId_t, std::pair<ContigSize_t, size_t>> contig_count_;

  bool createAddVariant( const GenomeId_t& genome_name,
                         const ContigId_t& contig_id,
                         ContigOffset_t contig_offset,
                         const std::string& reference_text,
                         const std::string& alternate_text,
                         const std::shared_ptr<const VariantEvidence> evidence_ptr);

  bool addThreadSafeVariant(std::shared_ptr<const Variant>& variant_ptr);

};

} // namespace

#endif //KGL_KGL_VARIANT_FACTORY_GRCH_IMPL_H
