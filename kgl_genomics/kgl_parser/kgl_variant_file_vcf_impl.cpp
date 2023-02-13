//
// Created by kellerberrin on 21/4/20.
//



#include <kel_exec_env.h>
#include "kel_utility.h"

#include "kgl_variant_file_vcf_impl.h"

#include <string>
#include <vector>
#include <memory>

// Implementation file classes need to be defined within namespaces.
namespace kellerberrin::genome {



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

RecordVCFIO::~RecordVCFIO() noexcept {

  // Ensure all worker threads are shutdown, suppress exceptions.
  try {

    launch_token_.wait();

  } catch(...) {}

}

bool RecordVCFIO::commenceVCFIO(const std::string& vcf_file_name) {

  if (not file_data_.commenceIO(vcf_file_name)) {

    // Enqueue an eof marker further up the pipeline.
    enqueueEOF();
    return false;

  }

  launch_token_ = detached_launch_.enqueueFuture(&RecordVCFIO::launchThreads, this);

  return true;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Launch the VCF record parser worker threads and await completion.

void RecordVCFIO::launchThreads() {

  WorkflowThreads vcf_record_threads{PARSER_THREADS_};
  std::vector<std::future<void>> thread_futures;
  for (size_t index = 0; index < vcf_record_threads.threadCount(); ++index) {

    thread_futures.push_back(vcf_record_threads.enqueueFuture(&RecordVCFIO::enqueueVCFRecord, this));

  }

  // Wait until processing is complete.
  for (auto const& future : thread_futures) {

    future.wait();

  }

  // Enqueue an eof marker further up the pipeline.
  enqueueEOF();

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parse the line record into a VCF record and enqueue the record.

void RecordVCFIO::enqueueVCFRecord() {

  while (true) {

    IOLineRecord line_record = file_data_.readIORecord();

    // Check for EOF condition.
    if (line_record.EOFRecord()) {

      // push the eof marker back on the queue and stop processing.
      file_data_.enqueueEOF();
      break;

    }

    // Get the line record.
    auto [line_count, line_string] = line_record.getLineData();

    // Check for zero length lines.
    if (line_string.empty()) {

      ExecEnv::log().warn("RecordVCFIO::enqueueVCFRecord; unexpected zero length line found at line: {}", line_count);
      continue;

    }

    // Skip header records by examining the first char for '#'
    if (line_string[0] == HEADER_CHAR_) {

      continue;

    }

    // Parse the line into basic VCF fields.
    std::optional<std::unique_ptr<VcfRecord>> vcf_record_opt = moveToVcfRecord(std::move(line_string));
    if (not vcf_record_opt) {

      ExecEnv::log().warn("FileVCFIO; Failed to parse VCF file: {} record line : {}", file_data_.fileName(), line_count);

    } else {

      // Queue the parsed VCF record for further processing.
      QueuedVCFRecord queue_record(std::pair<size_t, std::unique_ptr<VcfRecord>>(line_count, std::move(vcf_record_opt.value())));
      vcf_record_queue_.push(std::move(queue_record));

    }

  }

}

// Parse the VCF line into basic VCF fields.
std::optional<std::unique_ptr<VcfRecord>> RecordVCFIO::moveToVcfRecord(std::string line_record) {

  try {

    std::unique_ptr<VcfRecord> vcf_record_ptr(std::make_unique<VcfRecord>());

    std::vector<std::string_view> field_views = Utility::viewTokenizer(line_record, VCF_FIELD_DELIMITER_CHAR_);

    if (field_views.size() < MINIMUM_VCF_FIELDS_) {

      ExecEnv::log().error("FileVCFIO; VCF file: {}, record has less than the mandatory field count", file_data_.fileName());
      return std::nullopt;

    }

    vcf_record_ptr->contig_id = field_views[CONTIG_FIELD_IDX_];
    vcf_record_ptr->offset = std::stoull(std::string(field_views[OFFSET_FIELD_IDX_])) - 1; // all offsets are zero based.

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
  catch (const std::exception &e) {

    ExecEnv::log().error("FileVCFIO; Problem parsing record for VCF file: {}, Exception: {} thrown; VCF record ignored", file_data_.fileName(),  e.what());
    return std::nullopt;

  }

}

} // namespace