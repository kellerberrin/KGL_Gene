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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "kgl_process_sam.h"

namespace py = pybind11;
namespace kgl = kellerberrin::genome;


// Class to implement the Python bindings to the underlying C++ code.
class PythonProcessSamFile : public kgl::ProcessSAMFile {

public:

  explicit PythonProcessSamFile(const std::string& log_file, int read_quality) : ProcessSAMFile(log_file,
                                                                                                read_quality) {}
  ~PythonProcessSamFile() override = default;

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
   py::class_<PythonProcessSamFile>(m, "ProcessSamFile")
      .def(py::init<const std::string&, int>())
      .def("registerContigNumpy", &PythonProcessSamFile::registerContigNumpy)
      .def("readSamFile", &PythonProcessSamFile::readSAMFile)
      .def("getQueueContigs", &PythonProcessSamFile::getQueueContigs)
      .def("getQueueOffsets", &PythonProcessSamFile::getQueueOffsets)
      .def("getQueueSequences", &PythonProcessSamFile::getQueueSequences);

}


