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
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class ReadCountFilter : public VariantFilter {

public:

  explicit ReadCountFilter(NucleotideReadCount_t read_count) : read_count_(read_count) {}
  ~ReadCountFilter() override = default;

  bool applyFilter(ReadCountVariant& variant) const final { return variant.readCount() >= read_count_; }
  std::string filterName() const final;

private:

  NucleotideReadCount_t read_count_;

};


class MutantProportionFilter : public VariantFilter {

public:

  explicit MutantProportionFilter(double proportion) : mutant_proportion_(proportion) {}
  ~MutantProportionFilter() override = default;

  bool applyFilter(ReadCountVariant& variant) const final { return variant.proportion() >= mutant_proportion_; }
  std::string filterName() const final;

private:

  double mutant_proportion_;

};

}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
