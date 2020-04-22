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


class GrchVCFImpl : public VCFReaderMT {

public:

  GrchVCFImpl(std::shared_ptr<UnphasedGenome> vcf_genome_ptr,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string &vcf_file_name) : VCFReaderMT(vcf_file_name),
                                                  vcf_genome_ptr_(std::move(vcf_genome_ptr)),
                                                  genome_db_ptr_(std::move(genome_db_ptr)){

    for (auto const& [contig_id, contig_ptr] : genome_db_ptr_->getMap()) {

      contig_count_[contig_id]  = std::pair<ContigSize_t, size_t>(contig_ptr->contigSize(), 0);

    }

  }

  ~GrchVCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord &vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo &header_info) override;

  void readParseVCFImpl();


private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 10000;

  std::shared_ptr<UnphasedGenome> vcf_genome_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

// Progress counters.

  size_t vcf_variant_count_{0};
  std::atomic<uint64_t> record_count_{0};
  size_t variant_count_{0};

  std::map<ContigId_t, std::pair<ContigSize_t, size_t>> contig_count_;

};

} // namespace

#endif //KGL_KGL_VARIANT_FACTORY_GRCH_IMPL_H
