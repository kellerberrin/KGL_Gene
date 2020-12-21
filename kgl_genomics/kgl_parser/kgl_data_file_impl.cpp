//
// Created by kellerberrin on 18/4/20.
//

#include "kgl_data_file_impl.h"
#include <vector>
#include <memory>

namespace kgl = kellerberrin::genome;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::FileDataIO::~FileDataIO() noexcept {

  // Wait for completion, suppress any exceptions.
  try {

    launch_token_.wait();

  } catch(...) {}

}


bool kgl::FileDataIO::commenceIO(std::string read_file_name) {

  read_file_name_ = std::move(read_file_name);
  file_stream_opt_ = BaseStreamIO::getReaderStream(fileName());

  if (not file_stream_opt_) {

    ExecEnv::log().error("FileDataIO::rawDataIO; I/O error; could not open Data file: {}", fileName());
    enqueueEOF();
    return false;

  }

  // Activate the worker threads.
  launch_token_ = detached_launch_.enqueueTask(&FileDataIO::launchThreads, this);

  return true;

}



void kgl::FileDataIO::launchThreads() {

  ThreadPool vcf_record_threads{IO_THREAD_COUNT_};
  std::vector<std::future<void>> thread_futures;
  for (size_t index = 0; index < vcf_record_threads.threadCount(); ++index) {

    thread_futures.push_back(vcf_record_threads.enqueueTask(&FileDataIO::rawDataIO, this));

  }

  // Wait until processing is complete.
  for (auto const& future : thread_futures) {

    future.wait();

  }

  // Enqueue an eof marker further up the pipeline.
  enqueueEOF();

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the line records and just queue them.
void kgl::FileDataIO::rawDataIO() {

  try {

    while (true) {

      // ReadLine() should not block on eof.
      IOLineRecord line_record = file_stream_opt_.value()->readLine();
      if (not line_record) {

        // Terminate the read line loop.
        break;

      }

      // Check we have a valid line record.
      if (line_record.value().second->empty()) {

        ExecEnv::log().error("FileDataIO::rawDataIO; Empty Data record on line: {} - ignored", line_record.value().first);
        continue;

      }

      raw_io_queue_.push(std::move(line_record));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("FileDataIO::rawDataIO; Data file: {} unexpected I/O exception: {}", fileName(), e.what());

  }

}

