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
#include "kgl_genome_types.h"
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class GFFRecord {

public:

  explicit GFFRecord(Logger& logger) : log(logger) {};
  ~GFFRecord() = default;

  void readGFFFile(const std::string& gff_file_name);

private:

  Logger& log;

};



using FastaMap = std::map<const ContigId_t, Sequence_t>;
class FastaRecord {

public:

  explicit FastaRecord(Logger& logger): log(logger) {}
  ~FastaRecord() = default;

  void addSequence(const ContigId_t& contig, const Sequence_t& sequence);
  void readFastaFile(const std::string& fasta_file_name);
  const FastaMap& getFastaMap() const { return fasta_map_; }

private:

  Logger& log;

  FastaMap fasta_map_;

};

class ParseGFFSFasta { // parses the input gff(3) file and annotates it with a fasta sequence

public:

  ParseGFFSFasta( Logger& logger,
                  const std::string& gff_file_name,
                  const std::string& fasta_file_name): log(logger),
                                                       fasta_ptr_(std::make_unique<FastaRecord>(logger)),
                                                       gff_ptr_(std::make_unique<GFFRecord>(logger)) {

    fasta_ptr_->readFastaFile(fasta_file_name);
    gff_ptr_->readGFFFile(gff_file_name);

  }
  ~ParseGFFSFasta() = default;

  inline const GFFRecord& getGFF() { return *gff_ptr_; }
  inline const FastaRecord& getFasta() { return *fasta_ptr_; }

private:

  Logger& log; // Should declared first.

  std::unique_ptr<FastaRecord> fasta_ptr_;
  std::unique_ptr<GFFRecord> gff_ptr_;


};


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_GFF_FASTA_H
