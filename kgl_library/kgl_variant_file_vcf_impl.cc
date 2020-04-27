//
// Created by kellerberrin on 21/4/20.
//



#include <kel_exec_env.h>

#include "kgl_variant_file_vcf_impl.h"

#include <string>
#include <vector>
#include <memory>
#include <thread>

// Implementation file classes need to be defined within namespaces.
namespace kellerberrin::genome {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool VCFParseHeader::parseHeader(const std::string& vcf_file_name, std::unique_ptr<BaseStreamIO>&& vcf_stream) {

  // Open input file.
  if (not vcf_stream->open(vcf_file_name)) {

    ExecEnv::log().critical("I/O error; could not open VCF file: {}", vcf_file_name);

  }

  try {

    long counter = 0;
    bool found_header = false;

    while (true) {

      IOLineRecord line_record = vcf_stream->readLine();
      if (not line_record) break;

      std::string& record_str = *line_record.value().second;

      size_t pos = record_str.find_first_of(KEY_SEPARATOR_);

      if (pos != std::string::npos) {

        std::string key = record_str.substr(0, pos);
        std::string value = record_str.substr(pos, std::string::npos);

        pos = value.find_first_of(KEY_SEPARATOR_);

        if (pos != std::string::npos) {

          value = value.erase(pos, pos + std::string(KEY_SEPARATOR_).length());

        }

        pos = key.find_first_of(KEY_PREFIX_);

        if (pos != std::string::npos) {

          key = key.erase(pos, pos + std::string(KEY_PREFIX_).length());

        }

        vcf_header_info_.emplace_back(key, value);

      }

      std::string line_prefix = record_str.substr(0, FIELD_NAME_FRAGMENT_LENGTH_);
      if (line_prefix == FIELD_NAME_FRAGMENT_) {

        found_header = true;
        std::vector<std::string> field_vector;
        if (not ParseVCFMiscImpl::tokenize(record_str, "\t",  field_vector)) {

          ExecEnv::log().error("Unable to parse VCF header line: {}", record_str);

        }
        size_t field_count = 0;
        for(auto const& field :field_vector) {

          if (field_count >= SKIP_FIELD_NAMES_) {

            vcf_genomes_.push_back(field);

          }

          ++field_count;

        }

        break; // #CHROM is the last field in the VCF header so stop processing.

      }

      ++counter;

    }

    if (not found_header) {

      ExecEnv::log().error("VCF Genome Names Not Found");

    } else {

      ExecEnv::log().info("{} Genomes in VCF Header {}, Header lines processed: {}", vcf_genomes_.size(), vcf_file_name, counter);

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCFParseHeader::parseHeader; VCF file: {}, unexpected I/O exception: {}", vcf_file_name, e.what());

  }

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

RecordVCFIO::~RecordVCFIO() {

  for (auto &thread : vcf_record_thread_vec_) {

    thread.join();

  }

}

void RecordVCFIO::commenceVCFIO(size_t reader_threads) {

  reader_threads_ = reader_threads;

  commenceIO(PARSER_THREADS_);

  for (size_t i = 0; i < PARSER_THREADS_; ++i) {

    vcf_record_thread_vec_.emplace_back(&RecordVCFIO::enqueueVCFRecord, this);

  }



}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the text or gzip stream as an rvalue to a pointer.

void RecordVCFIO::enqueueVCFRecord() {

  while (true) {

    IOLineRecord line_record = readIORecord();

    if (not line_record) { // check for EOF condition.

      // Only queue as many EOF (null optional) tokens as the number of reader threads.
      for (size_t i = 0; i < reader_threads_; ++i) {

        vcf_record_queue_.push(std::nullopt);

      }
      break;

    }

    // Skip header records by examining the first char for '#'
    std::string& line_string = *line_record.value().second;
    if (line_string[0] == HEADER_CHAR_) {

      continue;

    }

    std::unique_ptr<VcfRecord> vcf_record_ptr(std::make_unique<VcfRecord>());
    if (not parseVCFRecord(std::move(line_record.value().second), vcf_record_ptr)) {

      ExecEnv::log().warn("FileVCFIO; Failed to parse VCF file: {} record line : {}", fileName(), line_record.value().first);

    } else {

      QueuedVCFRecord queue_record(std::pair<size_t, std::unique_ptr<VcfRecord>>(line_record.value().first, std::move(vcf_record_ptr)));
      vcf_record_queue_.push(std::move(queue_record));

    }

  }

}


bool RecordVCFIO::parseVCFRecord( std::unique_ptr<const std::string> line_record_ptr,
                                  const std::unique_ptr<VcfRecord> &vcf_record_ptr) {

  std::vector<std::string> field_vector;

  if (not ParseVCFMiscImpl::tokenize(*line_record_ptr, VCF_FIELD_DELIMITER_, field_vector)) {

    ExecEnv::log().error("FileVCFIO; VCF file: {}, problem parsing record", fileName());
    return false;

  }

  if (field_vector.size() < MINIMUM_VCF_FIELDS_) {

    ExecEnv::log().error("FileVCFIO; VCF file: {}, record has less than the mandatory field count", fileName());
    return false;

  }

  if (not moveToVcfRecord(std::move(line_record_ptr),field_vector, *vcf_record_ptr)) {

    ExecEnv::log().error("FileVCFIO; VCF file: {}, cannot parse VCF record field", fileName());
    return false;

  }

  return true;

}


const VcfHeaderInfo& RecordVCFIO::VCFReadHeader() {

  parseheader_.parseHeader(fileName(), getSynchStream());
  return parseheader_.getHeaderInfo();

}


bool RecordVCFIO::moveToVcfRecord(std::unique_ptr<const std::string> line_record_ptr, std::vector<std::string> &fields, VcfRecord &vcf_record) {

  try {

    vcf_record.line_record_ptr = std::move(line_record_ptr);

    vcf_record.contig_id = fields[0];
    vcf_record.offset = std::stoull(fields[1]) - 1; // all offsets are zero based.
    vcf_record.id = std::move(fields[2]);
    vcf_record.ref = std::move(fields[3]);
    // A deletion variant can be signalled by a missing alt value.
    if (fields[4] == FIELD_NOT_PRESENT_) {

      vcf_record.alt = "";

    } else {

      vcf_record.alt = std::move(fields[4]);

    }
    // The quality field can be omitted.
    if (fields[5] == FIELD_NOT_PRESENT_) {

      vcf_record.qual = 0.0;

    } else {

      vcf_record.qual = std::stod(fields[5]);

    }
    vcf_record.filter = std::move(fields[6]);
    vcf_record.info = std::move(fields[7]);

    if (fields.size() > MINIMUM_VCF_FIELDS_) {

      vcf_record.format = std::move(fields[8]);

      vcf_record.genotypeInfos.clear();
      for (size_t idx = (MINIMUM_VCF_FIELDS_ + 1); idx < fields.size(); ++idx) {

        vcf_record.genotypeInfos.push_back(std::move(fields[idx]));

      }

    }

    return true;

  }
  catch (const std::exception &e) {

    ExecEnv::log().error("FileVCFIO; Problem parsing record for VCF file: {}, Exception: {} thrown; VCF record ignored", fileName(),  e.what());
    ExecEnv::log().error("FileVCFIO; VCF record line: {}", fileName(),  e.what());
    return false;

  }

}

} // namespace