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
// Created by kellerberrin on 9/09/17.
//

#ifndef KGL_PROCESS_SAM_H
#define KGL_PROCESS_SAM_H


#include <memory>
#include <string>
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class ProcessSamRecord {

public:

  explicit ProcessSamRecord(Logger& logger): log(logger) {}
  ~ProcessSamRecord() = default;

  void processSamRecord(std::unique_ptr<const std::string>& record_ptr);

private:

  Logger& log;      // Emit log messages to console and log file.

};

}   // namespace genome
}   // namespace kellerberrin


#endif // KGL_PROCESS_SAM_H
