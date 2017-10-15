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
// Created by kellerberrin on 13/10/17.
//
#include "kgl_application.h"

namespace kgl = kellerberrin::genome;


kgl::GenomeApplication::GenomeApplication(kgl::Logger& log, const kgl::ExecEnv::Args& args ) {

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = kgl::ParseGffFasta(log).readFastaGffFile(args.fastaFile,
                                                                                                args.gffFile);
  // Wire-up the genome database.
  genome_db_ptr->createVerifyGenomeDatabase();

  // Create a data block to hold the read data.
  std::shared_ptr<kgl::ContigCountData> contig_data_ptr(std::make_shared<kgl::ContigCountData>());

  // Register with the genome database to setup the contig data blocks.
  genome_db_ptr->registerContigData(contig_data_ptr);

  // Attach a SAM reader to the contig data block and read in the SAM file.
  kgl::SamCountReader(contig_data_ptr, log).readSAMFile(args.mutantFile, args.readQuality);

  // Create a genome variant to hold the SNP variant data.
  std::shared_ptr<kgl::GenomeVariant> variant_ptr = kgl::GenomeAnalysis().simpleSNPVariants(contig_data_ptr,
                                                                                            genome_db_ptr);


}


