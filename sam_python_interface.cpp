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
// Created by kellerberrin on 31/08/17.
//

#include <pybind11/pybind11.h>
#include <iostream>
#include <fstream>

#include "read_sam_file.h"
#include "execenv.h"


namespace py = pybind11;



void print_dict(py::dict dict) {
/* Easily interact with Python types */
    std::cout << "Dictionary Element" << std::endl;

    for (auto item : dict) {
        std::cout << "key=" << std::string(py::str(item.first)) << ", "
                  << "value=" << std::string(py::str(item.second)) << std::endl;
    }
}



PYBIND11_MODULE(libread_sam, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("log", &test_logging, "A function which tests C++ logging");
    m.def("print_dict", &print_dict);
    py::class_<ProcessSamFile>(m, "ProcessSamFile")
            .def(py::init<>())
            .def("readSamFile", &ProcessSamFile::readSamFile);
}
