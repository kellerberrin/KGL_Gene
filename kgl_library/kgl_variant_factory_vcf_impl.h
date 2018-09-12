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

  ParseVCFImpl(std::shared_ptr<UnphasedPopulation> unphased_population_ptr,
               std::shared_ptr<const GenomeDatabase> genome_db_ptr,
               const std::string& vcf_file_name) : genome_db_ptr_(genome_db_ptr),
                                                   vcf_file_name_(vcf_file_name),
                                                   unphased_population_ptr_(unphased_population_ptr){

    reader_ptr_ = std::make_shared<VCFReaderMT<ParseVCFImpl>>(vcf_file_name, this, &ParseVCFImpl::ProcessVCFRecord);

    // Pre-initializes the UnphasedPopulation with a list of genomes/contigs.
    // Since some genomes may not have a variant (3D7).
    // The initialization of the reader above generates a list of genomes in the VCF file.
    setupPopulationStructure(genome_db_ptr_);

  }


  virtual ~ParseVCFImpl() = default;

  virtual void readParseVCFImpl();

  virtual void ProcessVCFRecord(size_t vcf_record_count, const seqan::VcfRecord& vcf_record) = 0;

  const std::vector<std::string>& getGenomeNames() const { return reader_ptr_->getFieldNames(); }

  const ContigId_t contigId(int32_t rId) const { return reader_ptr_->getContig(rId); }

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

  size_t addThreadSafeGenomeVariant(std::shared_ptr<const Variant>& variant_ptr);
  void setupPopulationStructure(std::shared_ptr<const GenomeDatabase> genome_db_ptr);

  std::shared_ptr<const GenomeDatabase> genome_db_ptr_; // read access only.
  const std::string vcf_file_name_;
  std::shared_ptr<VCFReaderMT<ParseVCFImpl>> reader_ptr_;

private:

  mutable std::mutex mutex_;  // mutex to lock the UnphasedPopulation structure when inserting variants.
  // This object is write accessed by multiple threads, it MUST BE mutex guarded for any access.
  std::shared_ptr<UnphasedPopulation> unphased_population_ptr_;   // Un-phased variants.

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_KGL_VARIANT_FACTORY_VCF_IMPL_H
