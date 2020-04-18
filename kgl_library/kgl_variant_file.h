//
// Created by kellerberrin on 15/4/20.
//

#ifndef KGL_VARIANT_FILE_H
#define KGL_VARIANT_FILE_H


#include <string>
#include <vector>

#include "kel_mt_queue.h"
#include "kgl_genome_types.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// Structure to return VCF header information as a vector of pair<key, value>
using VcfHeaderInfo = std::vector<std::pair<std::string, std::string>>;

// Parses the VCF header for Contig and Genome infromation.

class VCFParseHeader {

public:

  VCFParseHeader() = default;

  ~VCFParseHeader() = default;

  bool parseHeader(const std::string& file_name);

  [[nodiscard]] const std::vector<std::string>& getGenomes() const { return vcf_genomes_; }

private:

  std::vector<std::string> vcf_genomes_;                // Field (genome) names for each VCF record

  // Parser constants.
  static constexpr const char* CONTIG_NAME_FRAGMENT_{"##contig"};
  static constexpr const size_t CONTIG_NAME_FRAGMENT_LENGTH_{8};
  static constexpr const char* CONTIG_INFO_START_{"=<"};
  static constexpr const char* CONTIG_INFO_END_{">"};
  static constexpr const char* FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.


};

// Basic VCF Record LIne.
////////////////////////////////////////////////////////////////////////////////////////////////////////

class VcfRecord
{
public:

  // Default constructor.
  VcfRecord() : offset(INVALID_POS), qual(MISSING_QUAL) {}

  static std::unique_ptr<VcfRecord> EOF_RECORD() { return std::unique_ptr<VcfRecord>(); }

  // Numeric id of the reference sequence.
  ContigId_t contig_id;
  // Position on the reference.
  ContigOffset_t offset;
  // Textual identifier of the variant.
  std::string id;
  // Bases in the reference.
  std::string ref;
  // Bases in the alternatives, comma-separated.
  std::string alt;
  // Quality
  double qual;
  // Value of FILTER field.
  std::string filter;
  // Value of INFO field.
  std::string info;
  // Value of FORMAT field.
  std::string format;
  // The genotype infos.
  std::vector<std::string> genotypeInfos;


private:
  // Constant for invalid position.
  static constexpr const ContigOffset_t INVALID_POS = std::numeric_limits<ContigOffset_t>::max();
  // Undefined quality number.
  static constexpr const double MISSING_QUAL = std::numeric_limits<double>::lowest();

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual IO object (multi-threaded)


class BaseVCFIO {

public:

  BaseVCFIO(const BaseVCFIO&) = delete;
  explicit BaseVCFIO(std::string vcf_file_name) : vcf_file_name_(vcf_file_name) {

    parseFieldNames(vcf_file_name_);

  }
  virtual ~BaseVCFIO() = default;


  virtual void commenceIO() = 0;

  [[nodiscard]] virtual std::unique_ptr<VcfRecord> readVCFRecord() = 0;

  [[nodiscard]] virtual VcfHeaderInfo VCFReadHeader() = 0;

  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return parseheader_.getGenomes(); }

  [[nodiscard]] const std::string& fileName() const { return vcf_file_name_; }

private:

  std::string vcf_file_name_;
  VCFParseHeader parseheader_;    // Get genome and contig information.

  void parseFieldNames(const std::string& vcf_file_name) { parseheader_.parseHeader(vcf_file_name); }

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Seqan PIMPL facade IO object (multi-threaded)

class SeqanVCFIO : public BaseVCFIO {

public:

  explicit SeqanVCFIO(const std::string& vcf_file_name);
  ~SeqanVCFIO() override;


  void commenceIO() override; // Begin reading records, spawns threads.

  [[nodiscard]] std::unique_ptr<VcfRecord> readVCFRecord() override;

  [[nodiscard]] VcfHeaderInfo VCFReadHeader() override;


private:

  class SeqanVCFIOImpl;       // Forward declaration of seqan vcf read implementation class
  std::unique_ptr<SeqanVCFIOImpl> seqan_vcf_impl_ptr_;    // Seqan VCF parser PIMPL

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Uncompressed IO object (multi-threaded) - uses standard File IO functions to parse an uncompressed VCF file.

class FileVCFIO : public BaseVCFIO {

public:

  explicit FileVCFIO(const std::string& vcf_file_name) : BaseVCFIO(vcf_file_name),
                                                         raw_io_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                         vcf_record_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  ~FileVCFIO() override {

    raw_io_thread_ptr_->join();
    for (auto& thread : vcf_record_thread_vec_) {

      thread.join();

    }

  }


  void commenceIO() override; // Begin reading records, spawns threads.

  [[nodiscard]] std::unique_ptr<VcfRecord> readVCFRecord() override;

  [[nodiscard]] VcfHeaderInfo VCFReadHeader() override;


private:

  BoundedMtQueue<std::unique_ptr<std::string>> raw_io_queue_; // The raw IO queue
  BoundedMtQueue<std::unique_ptr<VcfRecord>> vcf_record_queue_; // Parsed VCF record queue

  std::unique_ptr<std::thread> raw_io_thread_ptr_;
  std::vector<std::thread> vcf_record_thread_vec_;

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records
  static constexpr const long PARSER_THREADS_{3};         // Threads parsing into vcf_records.
  static constexpr const char HEADER_CHAR_{'#'};          // If first char start with '*' then a header record.
  static constexpr const size_t MINIMUM_VCF_FIELDS_{8};   // At least 8 fields, any others are genotype fields (header specified).
  static constexpr const char* VCF_FIELD_DELIMITER_{"\t"};   // VCF Field separator.

  void rawVCFIO(); // enqueue seqan::vcf records.
  void enqueueVCFRecord(); // enqueue vcf_records.
  bool parseVCFRecord(const std::unique_ptr<std::string>& line_record_ptr, const std::unique_ptr<VcfRecord>& vcf_record_ptr);
  bool moveToVcfRecord(std::vector<std::string>& fields, VcfRecord& vcf_record);

};







} // end namespace


#endif //KGL_KGL_VARIANT_FILE_H
