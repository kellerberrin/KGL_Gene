//
// Created by kellerberrin on 18/4/20.
//

#include <string>
#include <vector>

#include "kel_exec_env.h"

#include "kel_mt_queue.h"
#include "kgl_genome_types.h"
#include "kgl_variant_file.h"

#include <seqan/vcf_io.h>


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


class kgl::SeqanVCFIO::SeqanVCFIOImpl  {

public:

  SeqanVCFIOImpl(const SeqanVCFIOImpl&) = delete;
  explicit SeqanVCFIOImpl(const std::string& vcf_file_name) : raw_io_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                              vcf_record_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                              vcfIn_ptr_(std::make_unique<seqan::VcfFileIn>(seqan::toCString(vcf_file_name)))
  {}
  ~SeqanVCFIOImpl() {

    raw_io_thread_ptr_->join();
    vcf_record_thread_ptr_->join();

  }



  void commenceIO(); // Begin reading records, spawns threads.

  [[nodiscard]] std::unique_ptr<VcfRecord> readVCFRecord();

  [[nodiscard]] VcfHeaderInfo VCFReadHeader();



private:

  BoundedMtQueue<std::unique_ptr<seqan::VcfRecord>> raw_io_queue_; // The raw IO queue
  BoundedMtQueue<std::unique_ptr<VcfRecord>> vcf_record_queue_; // Parsed VCF record queue
  std::unique_ptr<seqan::VcfFileIn> vcfIn_ptr_;
  std::unique_ptr<std::thread> raw_io_thread_ptr_;
  std::unique_ptr<std::thread> vcf_record_thread_ptr_;

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records

  [[nodiscard]] ContigId_t getContig(int32_t contig_idx) const;
  bool VCFRecordEOF();
  void rawVCFIO(); // enqueue seqan::vcf records.
  void enqueueVCFRecord(); // enqueue vcf_records.
  static void moveToVcfRecord(VcfRecord& vcf_record, std::unique_ptr<seqan::VcfRecord>& seqan_vcf_record, ContigId_t&& contig_id);

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::SeqanVCFIO::SeqanVCFIOImpl::moveToVcfRecord(VcfRecord& vcf_record, std::unique_ptr<seqan::VcfRecord>& seqan_vcf_record_ptr, ContigId_t&& contig) {

  vcf_record.contig_id = contig;
  vcf_record.offset = seqan_vcf_record_ptr->beginPos;
  seqan::move(vcf_record.id, seqan_vcf_record_ptr->id);
  seqan::move(vcf_record.ref, seqan_vcf_record_ptr->ref);
  seqan::move(vcf_record.alt, seqan_vcf_record_ptr->alt);
  vcf_record.qual = seqan_vcf_record_ptr->qual;
  seqan::move(vcf_record.filter, seqan_vcf_record_ptr->filter);
  seqan::move(vcf_record.info, seqan_vcf_record_ptr->info);
  seqan::move(vcf_record.format, seqan_vcf_record_ptr->format);

  auto begin = seqan::begin(seqan_vcf_record_ptr->genotypeInfos);
  auto end = seqan::end(seqan_vcf_record_ptr->genotypeInfos);

  vcf_record.genotypeInfos.clear();
  for (auto it = begin; it !=end; ++it) {

    std::string genotype;
    seqan::move(genotype, *it);
    vcf_record.genotypeInfos.push_back(std::move(genotype));

  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::ContigId_t kgl::SeqanVCFIO::SeqanVCFIOImpl::getContig(int32_t contig_idx) const {

  std::string contig_id = seqan::toCString(contigNames(context(*vcfIn_ptr_))[contig_idx]);

  return contig_id;

}


void kgl::SeqanVCFIO::SeqanVCFIOImpl::commenceIO() {

  raw_io_thread_ptr_ = std::make_unique<std::thread>(&SeqanVCFIOImpl::rawVCFIO, this);
  vcf_record_thread_ptr_ = std::make_unique<std::thread>(&SeqanVCFIOImpl::enqueueVCFRecord, this);

}



void kgl::SeqanVCFIO::SeqanVCFIOImpl::rawVCFIO() {

  try {

    std::unique_ptr<seqan::VcfRecord> vcf_record_ptr;
    while (true) {

      if (VCFRecordEOF()) {

        raw_io_queue_.push(std::unique_ptr<seqan::VcfRecord>(nullptr));
        break;

      }

      vcf_record_ptr = std::make_unique<seqan::VcfRecord>();
      seqan::readRecord(*vcf_record_ptr, *vcfIn_ptr_);
      raw_io_queue_.push(std::move(vcf_record_ptr));

    }

  }
  catch (std::exception const &e) {

    kel::ExecEnv::log().critical("VCF file unexpected Seqan I/O exception: {}", e.what());

  }

}



void kgl::SeqanVCFIO::SeqanVCFIOImpl::enqueueVCFRecord() {

  std::unique_ptr<seqan::VcfRecord> seqan_vcf_record_ptr;
  while (true) {

    raw_io_queue_.waitAndPop(seqan_vcf_record_ptr);

    if (not seqan_vcf_record_ptr) {

      vcf_record_queue_.push(VcfRecord::EOF_RECORD());
      break;

    }

    std::unique_ptr<VcfRecord> vcf_record(std::make_unique<VcfRecord>());
    moveToVcfRecord(*vcf_record, seqan_vcf_record_ptr, getContig(seqan_vcf_record_ptr->rID));
    vcf_record_queue_.push(std::move(vcf_record));

  }

}


std::unique_ptr<kgl::VcfRecord> kgl::SeqanVCFIO::SeqanVCFIOImpl::readVCFRecord() {

  std::unique_ptr<kgl::VcfRecord> vcf_record;
  vcf_record_queue_.waitAndPop(vcf_record);
  return vcf_record;

}


bool kgl::SeqanVCFIO::SeqanVCFIOImpl::VCFRecordEOF() {

  return seqan::atEnd(*vcfIn_ptr_);

}


kgl::VcfHeaderInfo kgl::SeqanVCFIO::SeqanVCFIOImpl::VCFReadHeader() {

  seqan::VcfHeader vcf_header;
  seqan::readHeader(vcf_header, *vcfIn_ptr_);

  VcfHeaderInfo header_info;
  for (size_t idx = 0; idx != seqan::length(vcf_header); ++idx) {

    std::string key = seqan::toCString(vcf_header[idx].key);
    std::string value = seqan::toCString(vcf_header[idx].value);
    header_info.push_back(std::pair(key,value));

  }

  return header_info;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Seqan VCF Parser PIMPL


kgl::SeqanVCFIO::SeqanVCFIO(const std::string& vcf_file_name) : BaseVCFIO(vcf_file_name),
                                                                seqan_vcf_impl_ptr_(std::make_unique<kgl::SeqanVCFIO::SeqanVCFIOImpl>(vcf_file_name))  {}
kgl::SeqanVCFIO::~SeqanVCFIO() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.

// Functionality passed to the implementation object (defined above).

// Initialize and begin reading records, spawns threads
void kgl::SeqanVCFIO::commenceIO() {

  seqan_vcf_impl_ptr_->commenceIO();

}

std::unique_ptr<kgl::VcfRecord> kgl::SeqanVCFIO::readVCFRecord() {

  return seqan_vcf_impl_ptr_->readVCFRecord();

}

kgl::VcfHeaderInfo kgl::SeqanVCFIO::VCFReadHeader() {

  return seqan_vcf_impl_ptr_->VCFReadHeader();

}


