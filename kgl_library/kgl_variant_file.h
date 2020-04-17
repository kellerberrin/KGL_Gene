//
// Created by kellerberrin on 15/4/20.
//

#ifndef KGL_VARIANT_FILE_H
#define KGL_VARIANT_FILE_H


#include <string>
#include <vector>

#include "kel_mt_queue.h"
#include "kgl_genome_types.h"

#include <seqan/vcf_io.h>


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
  VcfRecord(seqan::VcfRecord&& vcf_record, ContigId_t&& contig_id);

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
// Actual IO object )(multi-threaded)


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



class SeqanVCFIO : public BaseVCFIO {

public:

  SeqanVCFIO(const SeqanVCFIO&) = delete;
  explicit SeqanVCFIO(const std::string& vcf_file_name) : BaseVCFIO(vcf_file_name),
                                                          raw_io_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                          vcf_record_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                          vcfIn_ptr_(std::make_unique<seqan::VcfFileIn>(seqan::toCString(vcf_file_name)))
                                                          {}
  ~SeqanVCFIO() override {

    raw_io_thread_ptr_->join();
    vcf_record_thread_ptr_->join();

  }



  void commenceIO() override; // Begin reading records, spawns threads.

  [[nodiscard]] std::unique_ptr<VcfRecord> readVCFRecord() override;

  [[nodiscard]] VcfHeaderInfo VCFReadHeader() override;



private:

  BoundedMtQueue<std::unique_ptr<seqan::VcfRecord>> raw_io_queue_; // The raw IO queue
  BoundedMtQueue<std::unique_ptr<VcfRecord>> vcf_record_queue_; // Parsed VCF record queue
  std::unique_ptr<seqan::VcfFileIn> vcfIn_ptr_;
  VCFParseHeader parseheader_;    // Get genome and contig information.
  std::unique_ptr<std::thread> raw_io_thread_ptr_;
  std::unique_ptr<std::thread> vcf_record_thread_ptr_;

  static constexpr const long HIGH_TIDE_{100000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{10000};            // Low water mark to begin queueing VCF records

  [[nodiscard]] ContigId_t getContig(int32_t contig_idx) const;
  bool VCFRecordEOF();
  void rawVCFIO(); // enqueue seqan::vcf records.
  void enqueueVCFRecord(); // enqueue vcf_records.

  };




} // end namespace


#endif //KGL_KGL_VARIANT_FILE_H
