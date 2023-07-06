//
// Created by kellerberrin on 20/4/20.
//

#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant.h"

#include "kel_mem_alloc.h"

#include <string>

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

void kgl::GrchVCFImpl::processVCFHeader(const VCFHeaderInfo& header_info) {

  // Investigate header.
  VCFContigMap vcf_contig_map;
  VCFInfoRecordMap vcf_info_map;
  if (not VCFParseHeader::parseVcfHeader( header_info, vcf_contig_map, vcf_info_map)) {

    ExecEnv::log().error("GrchVCFImpl::processVCFHeader, Problem parsing header information in VCF file. No variants processed.");

  }

  // Store all the available VCF info fields.
  evidence_factory_.availableInfoFields(vcf_info_map);

  if (VCFParseHeader::checkVCFReferenceContigs(vcf_contig_map, genome_db_ptr_)) {

    ExecEnv::log().info("VCF File and Reference Genome Contig Size and Name All Match");

  } else {

    ExecEnv::log().info("VCF File and Reference Genome Contig Size and Name Mis-Match.");

  }

}


void kgl::GrchVCFImpl::readParseVCFImpl(const std::string &vcf_file_name) {

   // multi-threaded
  readVCFFile(vcf_file_name);
  // single threaded

}


void kgl::GrchVCFImpl::ProcessVCFRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) {

  // Parse the info fields into a map..
  // For performance reasons the info field is std::moved - don't reference again.
  auto mutable_info = const_cast<std::string&>(vcf_record_ptr->info);
  std::shared_ptr<const DataMemoryBlock> info_evidence_ptr = evidence_factory_.createVariantEvidence(std::move(mutable_info));

  // Look at the filter field for "Pass"
  bool passed_filter = Utility::toupper(vcf_record_ptr->filter) == PASSED_FILTERS_;

  // Convert VCF contig to genome contig.
  std::string contig = contig_alias_map_.lookupAlias(vcf_record_ptr->contig_id);

  // Check for multiple alt sequences
  size_t position = vcf_record_ptr->alt.find_first_of(MULIPLE_ALT_SEPARATOR_);  // Check for ',' separators
  // The alt field can be blank (deletion).
  if (position == std::string::npos or vcf_record_ptr->alt.empty()) {

    // We have no format data and only 1 variant specified.
    // These variables declared to make this obvious.
    uint32_t variant_count = 1;
    uint32_t variant_index = 0;

    VariantEvidence evidence(vcf_record_ptr->line_number,
                             unphased_population_ptr_->dataSource(),
                             passed_filter,
                             info_evidence_ptr,
                             nullptr,
                             variant_index,
                             variant_count);
    VariantEvidence evidence1(evidence);
    // Add the variant.
    std::shared_ptr<const Variant> variant_ptr(std::make_shared<const Variant>( contig,
                                                                                vcf_record_ptr->offset,
                                                                                VariantPhase::UNPHASED,
                                                                                vcf_record_ptr->id,
                                                                                DNA5SequenceLinear(StringDNA5(vcf_record_ptr->ref)),
                                                                                DNA5SequenceLinear(StringDNA5(vcf_record_ptr->alt)),
                                                                                evidence));

    if (not addThreadSafeVariant(variant_ptr, genome_db_ptr_->genomeId())) {

      ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Unable to add variant to population");

    } else {

      ++variant_count_;

    }

  } else {

    std::vector<std::string> alt_vector = Utility::charTokenizer(vcf_record_ptr->alt, MULIPLE_ALT_SEPARATOR_);

    if (alt_vector.empty()) {

      ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Zero sized alt vector, alt: {}", vcf_record_ptr->alt);

    }

    uint32_t variant_count = alt_vector.size();
    uint32_t variant_index = 0;
    for (auto const& alternate : alt_vector) {

      // We have no format data and multiple alternate alleles.
      // Setup the evidence object.
      VariantEvidence evidence(vcf_record_ptr->line_number,
                               unphased_population_ptr_->dataSource(),
                               passed_filter,
                               info_evidence_ptr,
                               nullptr,
                               variant_index,
                               variant_count);
      VariantEvidence evidence1(evidence);
      // Add the variant.
      std::shared_ptr<const Variant> variant_ptr(std::make_shared<const Variant>( contig,
                                                                                  vcf_record_ptr->offset,
                                                                                  VariantPhase::UNPHASED,
                                                                                  vcf_record_ptr->id,
                                                                                  DNA5SequenceLinear(StringDNA5(vcf_record_ptr->ref)),
                                                                                  DNA5SequenceLinear(StringDNA5(alternate)),
                                                                                  evidence));

      if (not addThreadSafeVariant(variant_ptr, genome_db_ptr_->genomeId())) {

        ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Unable to add variant to population");

      } else {

        ++variant_count_;

      }

    }

  }

  if (vcf_record_ptr->line_number % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("Processed :{} records, total variants: {}", vcf_record_ptr->line_number, variant_count_);
    ExecEnv::log().info("Contig: {}, offset: {}", contig, vcf_record_ptr->offset);

  }

}


bool kgl::GrchVCFImpl::addThreadSafeVariant(const std::shared_ptr<const Variant>& variant_ptr, const GenomeId_t& genome) const {

  // The population structure can be updated concurrently (embedded mutexes).

  if (filterVariant(variant_ptr, genome)) {

    std::vector<GenomeId_t> genome_vector;
    genome_vector.push_back(genome);

    return unphased_population_ptr_->addVariant(variant_ptr, genome_vector);

  } else {

    return true;

  }

}


bool kgl::SNPdbVCFImpl::filterVariant(const std::shared_ptr<const Variant>& /* variant_ptr*/ ,  const GenomeId_t& /* genome */) const {

  return true;

}
