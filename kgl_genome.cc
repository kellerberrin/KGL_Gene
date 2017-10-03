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
// Created by kellerberrin on 3/10/17.
//
#include "kgl_genome.h"

namespace kgl = kellerberrin::genome;

kgl::GenomeAnalysis::GenomeAnalysis(kgl::Logger& log, const kgl::ExecEnv::Args& args ) {

  // Parse the Fasta and GFF files.
  kgl::ParseGFFSFasta gff_fasta(log, args.gffFile, args.fastaFile);

  // Declare a SAM reader.
  kgl::LocalProcessSam process_sam(log, args.readQuality);

  // Register the Fasta contigs with the SAM reader.
  for (auto contig_pair : gff_fasta.getFasta().getFastaMap()) {

    process_sam.insertContig(contig_pair.first, contig_pair.second.length());

  }

  process_sam.readSAMFile(args.mutantFile);

  if (args.parentFile.length() > 0) {

    process_sam.readSAMFile(args.parentFile); // Optionally specified.

  }

}

