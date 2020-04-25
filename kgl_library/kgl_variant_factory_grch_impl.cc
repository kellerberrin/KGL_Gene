//
// Created by kellerberrin on 20/4/20.
//

#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant_vcf.h"


namespace kgl = kellerberrin::genome;


void kgl::GrchVCFImpl::processVCFHeader(const VcfHeaderInfo&) {


}


void kgl::GrchVCFImpl::readParseVCFImpl() {



  // multi-threaded
  readVCFFile();
  // single threaded

  for (auto const& [contig_id, count] : contig_count_) {

    ExecEnv::log().info("GrchVCFImpl; Contig id: {}, length: {}, variant count :{}", contig_id, count.first, count.second);

  }

}


void kgl::GrchVCFImpl::ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) {

  // Parse the info fields into an array and assign to a shared ptr.
//  auto mutable_info = const_cast<std::string&>(vcf_record.info);
//  std::shared_ptr<std::string> info_ptr = std::make_shared<std::string>(std::move(mutable_info));
  std::shared_ptr<std::string> null_str_ptr = std::make_shared<std::string>("");
  // Evidence object
  std::shared_ptr<VariantEvidence> evidence_ptr(std::make_shared<VariantEvidence>(null_str_ptr));

  // Convert VCF contig to genome contig.
  std::string contig = contig_alias_map_.lookupAlias(Utility::toupper(vcf_record.contig_id));

  // Add the variant.
  if (not createAddVariant( vcf_genome_ptr_->genomeId(),
                            contig,
                            vcf_record.offset,
                            vcf_record.ref,
                            vcf_record.alt,
                            evidence_ptr)) {

    ExecEnv::log().error("Parsing GRCh VCF, Problem parsing vcf_record");

  }

  ++variant_count_;

  if (vcf_record_count % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("Processed :{} records, total variants: {}", static_cast<size_t>(vcf_record_count), variant_count_);
    ExecEnv::log().info("Contig: {}, offset: {}", vcf_record.contig_id, vcf_record.offset);

  }

}


bool kgl::GrchVCFImpl::createAddVariant(const GenomeId_t& genome_name,
                                        const ContigId_t& contig_id,
                                        ContigOffset_t contig_offset,
                                        const std::string& reference_text,
                                        const std::string& alternate_text,
                                        const std::shared_ptr<const VariantEvidence> evidence_ptr)  {

  StringDNA5 reference_str(reference_text);
  StringDNA5 alternate_str(alternate_text);

  std::shared_ptr<const Variant> variant_ptr(std::make_shared<VCFVariant>(genome_name,
                                                                          contig_id,
                                                                          VariantSequence::UNPHASED,
                                                                          contig_offset,
                                                                          evidence_ptr,
                                                                          std::move(reference_str),
                                                                          std::move(alternate_str)));

  return addThreadSafeVariant(variant_ptr);

}


bool kgl::GrchVCFImpl::addThreadSafeVariant(std::shared_ptr<const Variant>& variant_ptr) {

  std::scoped_lock<std::mutex> lock(add_variant_mutex_); // Write Locked

  return vcf_genome_ptr_->addVariant(variant_ptr);

}
