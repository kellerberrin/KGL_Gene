//
// Created by kellerberrin on 19/4/20.
//

#ifndef KGL_VARIANT_FILE_IMPL_H
#define KGL_VARIANT_FILE_IMPL_H


#include <string>
#include <vector>
#include <fstream>
#include <thread>

#include "kel_exec_env.h"
#include "kel_mt_queue.h"
#include "kgl_genome_types.h"
#include "kgl_variant_file.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

namespace kellerberrin::genome {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plug one of the superclasses (defined in the implementation file) to read text or gzipped files.

class BaseStreamIO {

public:

  BaseStreamIO() = default;
  virtual ~BaseStreamIO() = default;

  virtual bool open(const std::string &file_name) = 0;
  virtual bool readLine(std::string &text_line) = 0;

private:

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parses the VCF header for Contig and Genome information.

class VCFParseHeader {

public:

  VCFParseHeader() = default;

  ~VCFParseHeader() = default;

  bool parseHeader(const std::string& file_name, std::unique_ptr<BaseStreamIO>&& vcf_stream);

  [[nodiscard]] const std::vector<std::string>& getGenomes() const { return vcf_genomes_; }

  [[nodiscard]] const VcfHeaderInfo& getHeaderInfo() const { return vcf_header_info_; }

private:

  std::vector<std::string> vcf_genomes_;                // Field (genome) names for each VCF record
  VcfHeaderInfo vcf_header_info_;
  // Parser constants.
  static constexpr const char* KEY_SEPARATOR_{"="};
  static constexpr const char* KEY_PREFIX_{"##"};
  static constexpr const char* FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.


};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IO object (multi-threaded) - uses standard File IO functions to parse VCF file.
// The code checks the file extension, if '.vcf' then it assumes a text file, if '.gz' then it assumes a gzipped file.

class FileVCFIO  {

public:

  explicit FileVCFIO(const std::string& vcf_file_name) : vcf_file_name_(vcf_file_name),
                                                         raw_io_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                         vcf_record_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  ~FileVCFIO();


  void commenceIO(); // Begin reading records, spawns threads.

  [[nodiscard]] std::unique_ptr<VcfRecord> readVCFRecord();

  [[nodiscard]] const VcfHeaderInfo& VCFReadHeader() const;

  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return parseheader_.getGenomes(); }

  [[nodiscard]] const std::string& fileName() const { return vcf_file_name_; }

private:

  std::string vcf_file_name_;

  BoundedMtQueue<std::unique_ptr<std::string>> raw_io_queue_; // The raw IO queue
  BoundedMtQueue<std::unique_ptr<VcfRecord>> vcf_record_queue_; // Parsed VCF record queue

  std::unique_ptr<std::thread> raw_io_thread_ptr_;
  std::vector<std::thread> vcf_record_thread_vec_;

  VCFParseHeader parseheader_;    // Get genome and contig information.

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records
  static constexpr const long PARSER_THREADS_{4};         // Threads parsing into vcf_records.
  static constexpr const char HEADER_CHAR_{'#'};          // If first char start with '#' then a header record.
  static constexpr const size_t MINIMUM_VCF_FIELDS_{8};   // At least 8 fields, any others are format and genotype fields (header specified).
  static constexpr const char* VCF_FIELD_DELIMITER_{"\t"};   // VCF Field separator.
  constexpr static const char* VCF_FILE_EXTENSTION_ = ".VCF"; // Valid file extensions (case insensitive)
  constexpr static const char* GZ_FILE_EXTENSTION_ = ".GZ"; // gzipped VCF file assumed.
  const std::string FIELD_NOT_PRESENT_{"."}; // no field value

  void rawVCFIO(std::unique_ptr<BaseStreamIO>&& vcf_stream); // read/decompress from disk.
  void enqueueVCFRecord(); // enqueue vcf_records.
  bool parseVCFRecord(const std::unique_ptr<std::string>& line_record_ptr, const std::unique_ptr<VcfRecord>& vcf_record_ptr);
  bool moveToVcfRecord(std::vector<std::string>& fields, VcfRecord& vcf_record);

};


} // end namespace

#endif //KGL_VARIANT_FILE_IMPL_H
