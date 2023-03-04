//
// Created by kellerberrin on 3/03/23.
//

#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_variant_vcf_impl.h"

#include <string>
#include <vector>
#include <memory>

// Implementation file classes need to be defined within namespaces.
namespace kgl = kellerberrin::genome;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::ParseVCF::~ParseVCF() {

  enqueue_thread_.joinThreads();
  if (stream_ptr_) {

    stream_ptr_->close();

  }
  vcf_pipeline_.clear();

}

bool kgl::ParseVCF::open(const std::string& vcf_file_name, size_t decompression_threads, size_t vcf_parse_threads) {

  auto stream_opt = BaseStreamIO::getStreamIO(vcf_file_name, decompression_threads);
  if (not stream_opt) {

    // Enqueue an eof marker further up the pipeline.
    enqueueEOF();
    return false;

  }

  stream_ptr_ = std::move(stream_opt.value());

  // The stream is now open, start parsing VCF fields.
  vcf_pipeline_.activatePipeline(vcf_parse_threads, &kgl::ParseVCF::moveToVcfRecord, this);

  enqueue_thread_.enqueueVoid(&ParseVCF::enqueueLineRecord, this);

  return true;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parse the line record into a VCF record and enqueue the record.

void kgl::ParseVCF::enqueueLineRecord() {

  while (true) {

    IOLineRecord line_record = stream_ptr_->readLine();

    // Check for EOF condition.
    if (line_record.EOFRecord()) {

      // push the eof marker back on the queue and stop processing.
      enqueueEOF();
      break;

    }

    // Peek at the line record.
    auto line_view = line_record.getView();

    // Check for zero length lines.
    if (line_view.empty()) {

      ExecEnv::log().warn( "ParseVCF::enqueueLineRecord; file: {}, unexpected zero length line found at line: {}"
                          , getFileName(), line_record.lineCount());
      continue;

    }

    // Skip header records by examining the first char for '#'
    if (line_view[0] == HEADER_CHAR_) {

      continue;

    }

    vcf_pipeline_.push(std::move(line_record));

  }

}

// Parse the VCF line into basic VCF fields.
std::unique_ptr<kgl::VcfRecord> kgl::ParseVCF::moveToVcfRecord(IOLineRecord line) {


  auto [line_num, line_record] = line.getLineData();

  std::unique_ptr<VcfRecord> vcf_record_ptr(std::make_unique<VcfRecord>());

  vcf_record_ptr->line_num = line_num;

  std::vector<std::string_view> field_views = Utility::viewTokenizer(line_record, VCF_FIELD_DELIMITER_CHAR_);

  if (field_views.size() < MINIMUM_VCF_FIELDS_) {

    ExecEnv::log().error("ParseVCF::moveToVcfRecord; VCF file: {}, line: {}, record has less than the mandatory field count: {}"
                         , getFileName(), line_num, MINIMUM_VCF_FIELDS_);
    return nullptr;

  }

  vcf_record_ptr->contig_id = field_views[CONTIG_FIELD_IDX_];
  // Note all VCF file offsets are 1 based (obviously designed by geneticists).
  // However, all offsets are zero based in this parser.
  vcf_record_ptr->offset = std::stoull(std::string(field_views[OFFSET_FIELD_IDX_])) - 1;

  // The identifier field. Set field not present "." to the empty string.
  if (field_views[IDENT_FIELD_IDX_] == FIELD_NOT_PRESENT_) {

    vcf_record_ptr->id = "";

  } else {

    vcf_record_ptr->id = field_views[IDENT_FIELD_IDX_];
    vcf_record_ptr->id = Utility::trimEndWhiteSpace(vcf_record_ptr->id);

  }

  vcf_record_ptr->ref = field_views[REF_FIELD_IDX_];
  // A deletion variant can be signalled by a missing alt value.
  if (field_views[ALT_FIELD_IDX_] == FIELD_NOT_PRESENT_) {

    vcf_record_ptr->alt = "";

  } else {

    vcf_record_ptr->alt = field_views[ALT_FIELD_IDX_];

  }
  // The quality field can be omitted.
  if (field_views[QUALITY_FIELD_IDX_] == FIELD_NOT_PRESENT_) {

    vcf_record_ptr->qual = 0.0;

  } else {

    vcf_record_ptr->qual = std::stod(std::string(field_views[QUALITY_FIELD_IDX_]));

  }
  vcf_record_ptr->filter = field_views[FILTER_FIELD_IDX_];
  vcf_record_ptr->info = field_views[INFO_FIELD_IDX_];

  if (field_views.size() > MINIMUM_VCF_FIELDS_) {

    vcf_record_ptr->format = field_views[FORMAT_FIELD_IDX_];

    vcf_record_ptr->genotypeInfos.clear();
    for (size_t idx = (MINIMUM_VCF_FIELDS_ + 1); idx < field_views.size(); ++idx) {

      vcf_record_ptr->genotypeInfos.emplace_back(field_views[idx]);

    }

  }

  // Move the entire line record to the VCF field record.
  vcf_record_ptr->line_record_ptr = std::make_unique<std::string>(std::move(line_record));

  return vcf_record_ptr;

}

