// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// Created by kellerberrin on 8/09/17.
//

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "kgl_logging.h"
#include "kgl_read_sam.h"
#include "kgl_consume_sam.h"



namespace py = pybind11;

namespace kgl = kellerberrin::genome;

// Simple class to bolt the Python producer and consumer together.
using PythonConsumer = kgl::ConsumeMTSAM<kgl::ConsumerNumpyRecord>;
using PythonProducer = kgl::ProduceMTSAM<PythonConsumer>;

class PythonProcessSam {

public:

  explicit PythonProcessSam(const std::string &log_file, int readQuality) : log(SAM_READ_MODULE_NAME_, log_file) {

    consumer_ptr_ = std::shared_ptr<PythonConsumer>(std::make_shared<PythonConsumer>(log));
    consumer_ptr_->readQuality(readQuality);
    producer_ptr_ = std::unique_ptr<PythonProducer>(std::make_unique<PythonProducer>(log, consumer_ptr_));

  }
  virtual ~PythonProcessSam() = default;

  inline void readSAMFile(const std::string& file_name) {

    producer_ptr_->readSamFile(file_name);
    consumer_ptr_->finalize();

    // Convert the inserted sequences to Python format.
    for(auto& contig : consumer_ptr_->contigDataMap().getMap()) {

      contig.second->getInsertArray().convertToQueue(contig.first, insert_queue_);

    }

  }

  void insertContig( const kgl::ContigId_t& contig_id,
                     kgl::NucleotideReadCount_t *data_ptr,
                     const kgl::ContigSize_t contig_size,
                     const kgl::ContigOffset_t num_nucleotides) {

    std::unique_ptr<kgl::ConsumerNumpyRecord> contig_matrix_ptr(std::make_unique<kgl::ConsumerNumpyRecord>(log,
                                                                                                 data_ptr,
                                                                                                 contig_size,
                                                                                                 num_nucleotides));
    contigDataMap().addContigBlock(contig_id, contig_matrix_ptr);

  }

  inline kgl::ContigDataMap<kgl::ConsumerNumpyRecord>& contigDataMap() { return consumer_ptr_->contigDataMap(); }
  inline const std::vector<std::string>& getQueueContigs() { return insert_queue_.getQueueContigs(); }
  inline const std::vector<std::size_t>& getQueueOffsets() { return insert_queue_.getQueueOffsets(); }
  inline const std::vector<std::string>& getQueueSequences() { return insert_queue_.getQueueSequences(); }

private:

  kgl::Logger log;                              // Must be declared First. Emit log messages to console and log file.
  static constexpr const char* SAM_READ_MODULE_NAME_{"SamRead"};  // Name of this module for the logger

  std::shared_ptr<PythonConsumer> consumer_ptr_; ;  // Thread safe SAM record consumer
  std::unique_ptr<PythonProducer> producer_ptr_;   //  Thread safe SAM record producer.
  kgl::InsertQueue insert_queue_;                // Inserted sequences converted for Python.

};


// Class to implement the Python bindings to the underlying C++ code.
class PythonProcessSamBind : public PythonProcessSam {

public:

  explicit PythonProcessSamBind(const std::string& log_file, int read_quality) : PythonProcessSam(log_file,
                                                                                                read_quality) {}
  ~PythonProcessSamBind() override = default;

  void registerContigNumpy( const kgl::ContigId_t &contig_name
                          , py::array_t<kgl::NucleotideReadCount_t> contig_numpy_data) {
    // Use the PyBind11 input variables to register the contig numpy with the underlying C++ object.
    auto numpy_dims = contig_numpy_data.mutable_unchecked<2>();
    auto numpy_data_ptr =  contig_numpy_data.mutable_data(0,0);
    auto contig_size = static_cast<kgl::ContigSize_t>(numpy_dims.shape(0));
    auto num_nucleotides = static_cast<kgl::ContigSize_t>(numpy_dims.shape(1));

    insertContig(contig_name, numpy_data_ptr, contig_size, num_nucleotides);

  }

private:

};


PYBIND11_MODULE(libread_sam, m) {

  m.doc() = "Python binding for 'libread_sam' using 'pybind11'"; // module docstring
   py::class_<PythonProcessSamBind>(m, "ProcessSamFile")
      .def(py::init<const std::string&, int>())
      .def("registerContigNumpy", &PythonProcessSamBind::registerContigNumpy)
      .def("readSamFile", &PythonProcessSamBind::readSAMFile)
      .def("getQueueContigs", &PythonProcessSamBind::getQueueContigs)
      .def("getQueueOffsets", &PythonProcessSamBind::getQueueOffsets)
      .def("getQueueSequences", &PythonProcessSamBind::getQueueSequences);

}


