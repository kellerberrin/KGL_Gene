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

#ifndef KGL_GFF_FASTA_H
#define KGL_GFF_FASTA_H


#include <memory>
#include <string>
#include <map>
#include "kgl_genome_db.h"
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Creates an instance of a Genome database object.
// Parses the input gff(3) file and annotates it with a fasta sequence
class ParseGFFSFasta {

public:

  ParseGFFSFasta( Logger& logger) : log(logger) {}
  ~ParseGFFSFasta() = default;

  std::unique_ptr<GenomeSequences> readFastaFile(const std::string& fasta_file_name);
  std::unique_ptr<GenomeSequences> readFastaGffFile(const std::string& fasta_file_name,
                                                    const std::string& gff_file_name);


private:

  Logger& log; // Should declared first.

  void readGffFile(const std::string& gff_file_name, GenomeSequences& genome_sequences);

};


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_GFF_FASTA_H
