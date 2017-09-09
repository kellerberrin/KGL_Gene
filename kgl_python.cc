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
//
//
// Created by kellerberrin on 8/09/17.
//

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "kgl_read_sam.h"

namespace py = pybind11;
namespace kgl = kellerberrin::genome;

void print_dict(const py::dict& dict) {
/* Easily interact with Python types */
  std::cout << "Dictionary Element" << std::endl;

  for (auto item : dict) {
    std::cout << "type=" << typeid(item.first).name() << "  key=" << std::string(py::str(item.first)) << std:: endl
              << "type=" << typeid(item.second).name() << "  value=" << std::string(py::str(item.second)) << std::endl;
  }
}

void contig_args(const std::vector<std::string> &contig_names, const std::vector<int64_t> &contig_sizes) {
/* Easily interact with Python types */
  std::cout << "List of Contigs" << std::endl;

  for (int i = 0; i < contig_sizes.size(); ++i) {
    std::cout << "  contig=" << contig_names[i] << " size=" << contig_sizes[i] << std:: endl;
  }
}


PYBIND11_MODULE(libread_sam, m) {
  m.doc() = "Python binding for 'libread_sam' using 'pybind11'"; // module docstring
  m.def("print_dict", &print_dict);
  m.def("contig_args", &contig_args);
  py::class_<kgl::ProcessSamFile>(m, "ProcessSamFile")
      .def(py::init<const std::string&>())
      .def("readSamFile", &kgl::ProcessSamFile::readSamFile);
}
