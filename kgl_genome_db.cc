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
// Created by kellerberrin on 7/10/17.
//

#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


void kgl::GenomeSequences::addContigSequence(kgl::ContigId_t& contig_id, kgl::Sequence_t sequence) {

  using ContigPtr = std::unique_ptr<kgl::ContigRecord>;
  ContigPtr contig_ptr(std::make_unique<kgl::ContigRecord>(contig_id, std::move(sequence)));

  auto result = genome_sequence_map_.insert(std::make_pair(contig_id, std::move(contig_ptr)));
  if (not result.second) {

    log.error("addContigSequence(), Attempted to add duplicate contig; {}", contig_id);

  }

}
