//
// Created by kellerberrin on 18/4/20.
//


#include <kel_exec_env.h>
#include "kel_mt_queue.h"
#include "kgl_genome_types.h"
#include "kgl_variant_file_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <string>
#include <vector>
#include <memory>
#include <thread>

// Classes need to be defined within namespaces.
namespace kellerberrin::genome {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Text IO.
class TextStreamIO : public BaseStreamIO {

public:

  TextStreamIO() = default;

  TextStreamIO(const TextStreamIO &) = delete;

  ~TextStreamIO() override = default;

  bool open(const std::string &file_name) override;

  bool readLine(std::string &text_line) override { return not std::getline(file_, text_line).eof(); }

private:

  std::ifstream file_;

};


bool TextStreamIO::open(const std::string &file_name) {

  try {

    // Open input file.

    file_.open(file_name);
    if (not file_.good()) {

      ExecEnv::log().error("TextStreamIO; I/O error; could not open file: {}", file_name);
      return false;

    }
  }
  catch (std::exception const &e) {

    ExecEnv::log().error("TextStreamIO; Opening file: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gzip IO

class GZStreamIO : public BaseStreamIO {

public:

  GZStreamIO() = default;

  GZStreamIO(const GZStreamIO &) = delete;

  ~GZStreamIO() override = default;

  bool open(const std::string &file_name) override;

  bool readLine(std::string &text_line) override { return not std::getline(gz_file_, text_line).eof(); }

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream gz_file_;

};

bool GZStreamIO::open(const std::string &file_name) {

  try {

    // Open input file.

    file_.open(file_name, std::ios_base::in | std::ios_base::binary);

    if (not file_.good()) {

      ExecEnv::log().error("GZStreamIO; I/O error; could not open file: {}", file_name);
      return false;

    }

    gz_file_.push(bio::gzip_decompressor());
    gz_file_.push(file_);

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("GZStreamIO; Opening file: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}

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

      std::string record_str;
      if (not vcf_stream->readLine(record_str)) break;

      std::string contig_prefix = record_str.substr(0, CONTIG_NAME_FRAGMENT_LENGTH_);
      if (contig_prefix == CONTIG_NAME_FRAGMENT_) {

        std::string contig_string = Utility::trimAllWhiteSpace(record_str);

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

FileVCFIO::~FileVCFIO() {

  if (raw_io_thread_ptr_) raw_io_thread_ptr_->join();

  for (auto &thread : vcf_record_thread_vec_) {

    thread.join();

  }

}

void FileVCFIO::commenceIO() {

  std::string file_ext = Utility::toupper(Utility::fileExtension(fileName()));

  if (file_ext == VCF_FILE_EXTENSTION_) {

    parseheader_.parseHeader(fileName(), std::make_unique<TextStreamIO>());
    raw_io_thread_ptr_ = std::make_unique<std::thread>(&FileVCFIO::rawVCFIO, this, std::make_unique<TextStreamIO>());

  } else if (file_ext == GZ_FILE_EXTENSTION_) {

    parseheader_.parseHeader(fileName(), std::make_unique<GZStreamIO>());
    raw_io_thread_ptr_ = std::make_unique<std::thread>(&FileVCFIO::rawVCFIO, this, std::make_unique<GZStreamIO>());

  } else {

    ExecEnv::log().error("FileVCFIO; Invalid file name: {}", fileName());
    ExecEnv::log().critical(
    "FileVCFIO; Unsupported file type: '{}' for variant calling. Must be VCF or GZ ('.vcf' or '.gz')", file_ext);

  }

  for (size_t i = 0; i < PARSER_THREADS_; ++i) {

    vcf_record_thread_vec_.emplace_back(&FileVCFIO::enqueueVCFRecord, this);

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the text or gzip stream as a base pointer.
void FileVCFIO::rawVCFIO(std::unique_ptr<BaseStreamIO> &&vcf_stream) {

  try {

    if (not vcf_stream->open(fileName())) {

      ExecEnv::log().critical("I/O error; could not open VCF file: {}", fileName());

    }

    while (true) {

      std::unique_ptr<std::string> line_record_ptr = std::make_unique<std::string>();

      if (not vcf_stream->readLine(*line_record_ptr)) {

        // Enqueue the null eof indicator for each consumer thread.
        for (size_t i = 0; i < vcf_record_thread_vec_.size(); ++i) {

          raw_io_queue_.push(std::unique_ptr<std::string>(nullptr));

        }
        break;

      }

      // Check we have read a non-empty string.
      if (line_record_ptr->empty()) {

        ExecEnv::log().error("Empty VCF record string returned");
        continue;

      }

      // Skip header records.
      if ((*line_record_ptr)[0] == HEADER_CHAR_) {

        continue;

      }

      raw_io_queue_.push(std::move(line_record_ptr));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file unexpected Seqan I/O exception: {}", e.what());

  }

}


void FileVCFIO::enqueueVCFRecord() {

  size_t record_line = 0;
  std::unique_ptr<std::string> line_record_ptr;
  while (true) {

    raw_io_queue_.waitAndPop(line_record_ptr);

    if (not line_record_ptr) {

      vcf_record_queue_.push(VcfRecord::EOF_RECORD());
      break;

    }

    ++record_line;

    std::unique_ptr<VcfRecord> vcf_record_ptr(std::make_unique<VcfRecord>());
    if (not parseVCFRecord(line_record_ptr, vcf_record_ptr)) {

      ExecEnv::log().warn("Failed to parse VCF file: {} record line : {}", fileName(), record_line);

    } else {

      vcf_record_queue_.push(std::move(vcf_record_ptr));

    }

  }

}


bool FileVCFIO::parseVCFRecord(const std::unique_ptr<std::string> &line_record_ptr,
                               const std::unique_ptr<VcfRecord> &vcf_record_ptr) {

  std::vector<std::string> field_vector;
  if (not ParseVCFMiscImpl::tokenize(std::move(*line_record_ptr), VCF_FIELD_DELIMITER_, field_vector)) {

    ExecEnv::log().error("VCF file: {}, problem parsing record", fileName());
    return false;

  }

  if (field_vector.size() < MINIMUM_VCF_FIELDS_) {

    ExecEnv::log().error("VCF file: {}, record has less than the mandatory field count", fileName());
    return false;

  }

  if (not moveToVcfRecord(field_vector, *vcf_record_ptr)) {

    ExecEnv::log().error("VCF file: {}, cannot parse VCF record field", fileName());
    return false;

  }

  return true;

}


std::unique_ptr<VcfRecord> FileVCFIO::readVCFRecord() {

  std::unique_ptr<VcfRecord> vcf_record;
  vcf_record_queue_.waitAndPop(vcf_record);
  return vcf_record;

}


VcfHeaderInfo FileVCFIO::VCFReadHeader() {

  VcfHeaderInfo header_info;

  return header_info;

}


bool FileVCFIO::moveToVcfRecord(std::vector<std::string> &fields, VcfRecord &vcf_record) {

  try {

    vcf_record.contig_id = std::move(fields[0]);
    vcf_record.offset = std::stoull(fields[1]) - 1; // all offsets are zero based.
    vcf_record.id = std::move(fields[2]);
    vcf_record.ref = std::move(fields[3]);
    vcf_record.alt = std::move(fields[4]);
    vcf_record.qual = std::stod(fields[5]);
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

    ExecEnv::log().error("Problem parsing record for VCF file: {}, Exception: {} thrown; VCF record ignored", e.what(),
                         fileName());
    return false;

  }

}

} // namespace