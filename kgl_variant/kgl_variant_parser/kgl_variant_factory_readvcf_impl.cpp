//
// Created by kellerberrin on 16/4/20.
//

#include "kgl_variant_factory_readvcf_impl.h"


namespace kgl = kellerberrin::genome;


void kgl::VCFReaderMT::readHeader(const std::string& file_name) {

  if (not parseheader_.parseHeader(file_name)) {

    ExecEnv::log().error("RecordVCFIO::VCFReadHeader, Problem parsing VCF header file: {}", file_name);

  }

  processVCFHeader(parseheader_.getHeaderInfo());

}


void kgl::VCFReaderMT::readVCFFile(const std::string& vcf_file_name) {


  ExecEnv::log().info("Parse Header VCF file: {}", vcf_file_name);

  readHeader(vcf_file_name);

  ExecEnv::log().info("Begin processing VCF file: {}", vcf_file_name);

  ExecEnv::log().info("Spawning: {} Consumer threads to process the VCF file", parser_threads_.threadCount());

  // Queue the worker thread tasks.
  std::vector<std::future<void>> thread_futures;
  for (size_t index = 0; index < parser_threads_.threadCount(); ++index) {

    thread_futures.push_back(parser_threads_.enqueueFuture(&VCFReaderMT::VCFConsumer, this));

  }

  // start reading records asynchronously.
  if (not vcf_parser_.open(vcf_file_name)) {

    ExecEnv::log().error("VCFReaderMT::readVCFFile; Problem opening VCF file: {}", vcf_file_name);
    return;

  }

  // Wait until processing is complete.
  for (auto const& future : thread_futures) {

    future.wait();

  }

}


void kgl::VCFReaderMT::VCFConsumer() {

  // Loop until EOF.
  size_t final_count = 0;
  while (true) {

    // Dequeue the vcf record.
    std::unique_ptr<const VCFRecord> vcf_record_ptr = vcf_parser_.readVCFRecord();
    // Terminate on EOF
    if (vcf_record_ptr->EOFRecord()) {

      vcf_parser_.enqueueEOF();
      break;  // Eof encountered, terminate processing.

    }

    if (vcf_record_ptr->contig_id.empty()) {

      ExecEnv::log().error("Empty VCF_record encountered; consumer thread terminates.");
      break;

    }

    // Call the consumer object with the dequeued record.
    ProcessVCFRecord(std::move(vcf_record_ptr));
    ++final_count;

  }

  ExecEnv::log().info("Final; Consumer thread processed: {} VCF records", final_count);


}


