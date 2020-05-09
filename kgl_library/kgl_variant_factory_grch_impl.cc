//
// Created by kellerberrin on 20/4/20.
//

#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant_vcf.h"


namespace kgl = kellerberrin::genome;


// Efficient parser for the info field.
// Implemented as a simple finite state parser.
bool kgl::GrchInfoParser::parseInfo() {

  enum class ParserStates { KeyToken, ValueToken} parser_state = ParserStates::KeyToken;
  size_t info_length = info_view_.length();
  size_t key_token_count = 0;
  size_t key_token_offset = 0;
  size_t value_token_count = 0;
  size_t value_token_offset = 0;
  std::string_view key_view;
  std::string_view value_view;
  const std::string_view empty_value;
  for (size_t index = 0; index < info_length; ++index) {

    switch(parser_state) {

      case ParserStates::KeyToken:
        if (info_view_[index] == INFO_VALUE_DELIMITER_) {

          key_view = info_view_.substr(key_token_offset, key_token_count);
          parser_state = ParserStates::ValueToken;
          value_token_offset = index + 1;
          value_token_count = 0;

        } else if (info_view_[index] == INFO_FIELD_DELIMITER_) {

          key_view = info_view_.substr(key_token_offset, key_token_count);
          auto result = parsed_token_map_.emplace(key_view, empty_value);
          if (not result.second) {

            ExecEnv::log().warn("GrchInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

          }
          key_token_offset = index + 1;
          key_token_count = 0;

        } else {

          ++key_token_count;

        }
        break;

      case ParserStates::ValueToken:
        if (info_view_[index] == INFO_FIELD_DELIMITER_) {

          value_view = info_view_.substr(value_token_offset, value_token_count);
          auto result = parsed_token_map_.emplace(key_view, value_view);
          if (not result.second) {

            ExecEnv::log().warn("GrchInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

          }
          parser_state = ParserStates::KeyToken;
          key_token_offset = index + 1;
          key_token_count = 0;

        } else {

          ++value_token_count;

        }
        break;

    } // switch

  } // for loop

  // Check that we are in ParserStates::ValueToken
  // and value_token_offset + value_token_count does not exceed the length of the info field.
  if (parser_state == ParserStates::ValueToken) {

    if (value_token_offset + value_token_count > info_length) {

      ExecEnv::log().error("GrchInfoParser::parseInfo, Final Value Token Offset: {}, Size: {} exceeds the Info size :{}",
                            value_token_offset, value_token_count, info_length);
      return false;

    }

    value_view = info_view_.substr(value_token_offset, value_token_count);
    auto result = parsed_token_map_.emplace(key_view, value_view);
    if (not result.second) {

      ExecEnv::log().warn("GrchInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

    }

  } else {

    ExecEnv::log().error("GrchInfoParser::parseInfo, parser terminates in unexpected state");
    return false;

  }

 return true;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

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

//  std::map<std::string, std::string> info_key_value_map;
//  if (not ParseVCFMiscImpl::tokenizeVcfInfoKeyValues(vcf_record.info, info_key_value_map)) {

//    ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord; Unable to parse VCF record info field: {}", vcf_record.info);

//  }

  auto mutable_info = const_cast<std::string&>(vcf_record.info);
  GrchInfoParser info_parser(std::move(mutable_info));
  if (not info_parser.parseInfo()) {

    ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Problem parsing info field");

  }

  std::shared_ptr<std::string> null_str_ptr = std::make_shared<std::string>("");
  // Evidence object
  std::shared_ptr<VariantEvidence> evidence_ptr(std::make_shared<VariantEvidence>(null_str_ptr, vcf_record_count));

  // Convert VCF contig to genome contig.
  std::string contig = contig_alias_map_.lookupAlias(vcf_record.contig_id);

  // Check for multiple alt sequences
  size_t position = vcf_record.alt.find_first_of(MULIPLE_ALT_SEPARATOR_);  // Check for ',' separators
  // The alt field can be blank (deletion).
  if (position == std::string::npos or vcf_record.alt.empty()) {

    // Add the variant.
    if (not createAddVariant( vcf_genome_ptr_->genomeId(),
                              contig,
                              vcf_record.offset,
                              vcf_record.ref,
                              vcf_record.alt,
                              evidence_ptr)) {

      ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Problem parsing vcf_record.");

    }

    ++variant_count_;

  } else {

    std::vector<std::string> alt_vector;
    if (not ParseVCFMiscImpl::tokenize(vcf_record.alt, alt_separator_, alt_vector)) {

      ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Unable to parse VCF record alt field: {}", vcf_record.alt);

    } else {

      if (alt_vector.empty()) {

        ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Zero sized alt vector, alt: {}", vcf_record.alt);

      }

      for (auto const& alt : alt_vector) {

        // Add the variant.
        if (not createAddVariant( vcf_genome_ptr_->genomeId(),
                                  contig,
                                  vcf_record.offset,
                                  vcf_record.ref,
                                  alt,
                                  evidence_ptr)) {

          ExecEnv::log().error("GrchVCFImpl::ProcessVCFRecord, Parsing GRCh VCF, Problem parsing vcf_record");

        }

        ++variant_count_;

      }

    }

  }

  if (vcf_record_count % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("Processed :{} records, total variants: {}, info keys: {}", vcf_record_count, variant_count_, info_parser.getMap().size());
    ExecEnv::log().info("Contig: {}, offset: {}", contig, vcf_record.offset);

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
