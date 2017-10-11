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

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeSequences> genome_db_ptr(std::make_shared<kgl::GenomeSequences>());

  { // Attach to a scoped reader to parse in Fasta and Gff files.

    kgl::ParseGffFasta(log).readFastaGffFile(args.fastaFile, args.gffFile, genome_db_ptr);
    // Wireup the genome database.
    genome_db_ptr->createVerifyGenomeDatabase();

  }

  // Create a data block to hold the read data.
  std::shared_ptr<kgl::ContigDataBlock> contig_data_ptr(std::make_shared<kgl::ContigDataBlock>());

  // Register with the genome database to setup the contig data blocks.
  genome_db_ptr->registerContigData(contig_data_ptr);

  { // Attach a scoped SAM reader to the data block and read the SAM file into the data block.

    kgl::LocalProcessSam(contig_data_ptr, log).readSAMFile(args.mutantFile, args.readQuality);
  }

}

