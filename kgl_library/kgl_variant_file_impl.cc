//
// Created by kellerberrin on 18/4/20.
//

#include "kgl_variant_file_impl.h"
#include <vector>
#include <memory>
#include <thread>

namespace kgl = kellerberrin::genome;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::FileVCFIO::~FileVCFIO() {

  // Shutdown the IO thread.
  if (raw_io_thread_ptr_) raw_io_thread_ptr_->join();

}

void kgl::FileVCFIO::commenceIO(size_t reader_threads) {

  // The number of consumer threads reading the raw record queue; used for shutdown.
  reader_threads_ = reader_threads;
  // Commence reading with just one IO thread.
  raw_io_thread_ptr_ = std::make_unique<std::thread>(&FileVCFIO::rawVCFIO, this);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the text or gzip stream as an rvalue to a pointer.
void kgl::FileVCFIO::rawVCFIO() {

  try {

    std::optional<std::unique_ptr<BaseStreamIO>> vcf_stream_opt = BaseStreamIO::getReaderStream(fileName());

    if (not vcf_stream_opt) {

      ExecEnv::log().critical("FileVCFIO; I/O error; could not open VCF file: {}", fileName());

    }

    while (true) {

      IOLineRecord line_record = vcf_stream_opt.value()->readLine();
      if (not line_record) {

        // Enqueue the null eof indicator for each consumer thread.
        for (size_t i = 0; i < reader_threads_; ++i) {

          raw_io_queue_.push(std::nullopt);

        }

        // Terminate the read line loop.
        break;

      }

      // Check we have a valid line record.
      if (line_record.value().second->empty()) {

        ExecEnv::log().error("FileVCFIO; Empty VCF record on line: {} - ignored", line_record.value().first);
        continue;

      }

      raw_io_queue_.push(std::move(line_record));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("FileVCFIO; VCF file: {} unexpected I/O exception: {}", fileName(), e.what());

  }

}

