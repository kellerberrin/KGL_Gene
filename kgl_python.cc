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
#include "kgl_exec_env.h"
#include "kgl_sam_read.h"



namespace py = pybind11;

namespace kgl = kellerberrin::genome;

// Simple class to bolt the Python producer and consumer together.
using NumpyDataBlock = kgl::ContigDataMap<kgl::ConsumerNumpyRecord>;
using PythonConsumer = kgl::ConsumeMTSAM<NumpyDataBlock>;
using PythonProducer = kgl::ProduceMTSAM<PythonConsumer>;

class PythonProcessSam {

public:

  explicit PythonProcessSam(const std::string &log_file) {

    kgl::ExecEnv::createLogger(SAM_READ_MODULE_NAME_, log_file); // Must be the first statement.
    numpy_data_ptr_ = std::shared_ptr<NumpyDataBlock >(std::make_shared<NumpyDataBlock>());
    consumer_ptr_ = std::shared_ptr<PythonConsumer>(std::make_shared<PythonConsumer>(kgl::ExecEnv::log(),
                                                                                     numpy_data_ptr_));
    producer_ptr_ = std::unique_ptr<PythonProducer>(std::make_unique<PythonProducer>(kgl::ExecEnv::log(),
                                                                                     consumer_ptr_));

  }
  virtual ~PythonProcessSam() = default;

  inline void readSAMFile(const std::string& file_name, int readQuality) {

    consumer_ptr_->readQuality(static_cast<unsigned char>(readQuality));
    producer_ptr_->readSamFile(file_name);
    consumer_ptr_->finalize();

    // Convert the inserted sequences to Python format.
    for(auto& contig : consumer_ptr_->contigDataMap().getMap()) {

      contig.second->getInsertArray().convertToQueue(contig.first, insert_queue_);

    }

  }

  bool insertContig( const kgl::ContigId_t& contig_id,
                     kgl::NucleotideReadCount_t *data_ptr,
                     const kgl::ContigSize_t contig_size,
                     const kgl::ContigOffset_t num_nucleotides) {

    return numpy_data_ptr_->insertContig(contig_id, contig_size, data_ptr, num_nucleotides);

  }

  inline const std::vector<std::string>& getQueueContigs() { return insert_queue_.getQueueContigs(); }
  inline const std::vector<std::size_t>& getQueueOffsets() { return insert_queue_.getQueueOffsets(); }
  inline const std::vector<std::string>& getQueueSequences() { return insert_queue_.getQueueSequences(); }

private:

  static constexpr const char* SAM_READ_MODULE_NAME_{"SamRead"};  // Name of this module for the logger

  std::shared_ptr<NumpyDataBlock> numpy_data_ptr_;  // Numpy contig data blocks passed in from python.
  std::shared_ptr<PythonConsumer> consumer_ptr_;   // Thread safe SAM record consumer
  std::unique_ptr<PythonProducer> producer_ptr_;   //  Thread safe SAM record producer.
  kgl::InsertQueue insert_queue_;                // Inserted sequences converted for Python.

};


// Class to implement the Python bindings to the underlying C++ code.
class PythonProcessSamBind : public PythonProcessSam {

public:

  explicit PythonProcessSamBind(const std::string& log_file) : PythonProcessSam(log_file) {}
  ~PythonProcessSamBind() override = default;

  bool registerContigNumpy( const kgl::ContigId_t &contig_name
                          , py::array_t<kgl::NucleotideReadCount_t> contig_numpy_data) {
    // Use the PyBind11 input variables to register the contig numpy with the underlying C++ object.
    auto numpy_dims = contig_numpy_data.mutable_unchecked<2>();
    auto numpy_data_ptr =  contig_numpy_data.mutable_data(0,0);
    auto contig_size = static_cast<kgl::ContigSize_t>(numpy_dims.shape(0));
    auto num_nucleotides = static_cast<kgl::ContigSize_t>(numpy_dims.shape(1));

    return insertContig(contig_name, numpy_data_ptr, contig_size, num_nucleotides);

  }

private:

};


PYBIND11_MODULE(libread_sam, m) {

  m.doc() = "Python binding for 'libread_sam' using 'pybind11'"; // module docstring
   py::class_<PythonProcessSamBind>(m, "ProcessSamFile")
      .def(py::init<const std::string&>())
      .def("registerContigNumpy", &PythonProcessSamBind::registerContigNumpy)
      .def("readSamFile", &PythonProcessSamBind::readSAMFile)
      .def("getQueueContigs", &PythonProcessSamBind::getQueueContigs)
      .def("getQueueOffsets", &PythonProcessSamBind::getQueueOffsets)
      .def("getQueueSequences", &PythonProcessSamBind::getQueueSequences);

}


