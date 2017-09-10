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
#include <pybind11/numpy.h>
#include "kgl_read_sam.h"

namespace py = pybind11;
namespace kgl = kellerberrin::genome;


using numpy_type = uint32_t;

py::array_t<numpy_type> contig_args( const std::vector<std::string> &contig_names
                                  , const std::vector<int64_t> &contig_sizes
                                  , int nucleotides) {

  std::cout << "List of Contigs, nucleotides:" << nucleotides << std::endl;

  for (int i = 0; i < contig_sizes.size(); ++i) {
    std::cout << "  contig=" << contig_names[i] << " size=" << contig_sizes[i] << std::endl;
  }

  // Allocate and initialize some data; make this big so
  // we can see the impact on the process memory use:
  constexpr size_t size = 100 * 1000 * 1000;
  auto foo = new numpy_type[size];
  for (size_t i = 0; i < size; i++) {
    foo[i] = (numpy_type) i;
  }


  // Create a Python object that will free the allocated
  // memory when destroyed:
  py::capsule free_when_done(foo, [](void *f) {
    auto boo = reinterpret_cast<numpy_type *>(f);
    std::cerr << "Element [3] = " << boo[3] << "\n";
    std::cerr << "freeing memory @ " << f << "\n";
    delete[] boo;
  });

  return py::array_t<numpy_type>(
      {100, 1000, 1000}, // shape
      {1000*1000*sizeof(numpy_type), 1000*sizeof(numpy_type), sizeof(numpy_type)}, // C-style contiguous strides
      foo, // the data pointer
      free_when_done); // numpy array references this parent

}



void numpy_args(py::array_t<numpy_type> a, numpy_type v, const std::string &contig_name) {

  auto r = a.mutable_unchecked<2>();
  auto ptr =  a.mutable_data(0,0);

  std::cerr << "ptr has type: " << typeid(ptr).name() << std::endl;

  std::cerr << "ptr " << ((typeid(ptr).name() == typeid(numpy_type*).name()) ? "matches uint32_t" : "does NOT match uint32_t") << std::endl;


  std::cerr << "dim1:" << r.shape(0) << " dim2:" << r.shape(1) <<  " contig:" << contig_name << std::endl;

  for (ssize_t i = 0; i < r.shape(0); i++) {

    for (ssize_t j = 0; j < r.shape(1); j++) {

      r(i, j) += v;

    }

  }

}


class NoLeak {

public:

  NoLeak(): pvec(std::make_unique<std::vector<std::string>>(10000000, "3.141592")) {}
  ~NoLeak() { std::cerr << "NoLeak out of scope" << std::endl; }

  std::vector<std::string>& getVec() { *pvec; }

private:

  std::unique_ptr<std::vector<std::string>> pvec;

};


PYBIND11_MODULE(libread_sam, m) {
  m.doc() = "Python binding for 'libread_sam' using 'pybind11'"; // module docstring
  m.def("contig_args", &contig_args);
  m.def("numpy_args", &numpy_args, py::arg().noconvert(), py::arg(), py::arg());
  py::class_<kgl::ProcessSamFile>(m, "ProcessSamFile")
      .def(py::init<const std::string&>())
      .def("readSamFile", &kgl::ProcessSamFile::readSamFile);
  py::class_<NoLeak>(m, "NoLeak")
      .def(py::init<>())
      .def("getVec", &NoLeak::getVec);
}


