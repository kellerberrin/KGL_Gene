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
// Created by kellerberrin on 3/10/17.
//

#include "kgl_gff_fasta.h"
#include <seqan/seq_io.h>

namespace kgl = kellerberrin::genome;


void kgl::GFFRecord::readGFFFile(const std::string& gff_file_name) {


}


void kgl::FastaRecord::addSequence(const kgl::ContigId_t& contig_id, const kgl::Sequence_t& sequence) {


  auto result = fasta_map_.insert(std::make_pair(contig_id, std::move(sequence)));

  if (not result.second) {

    log.error("addSequence(), Attempted to add duplicate contig; {}", contig_id);

  }


}

void kgl::FastaRecord::readFastaFile(const std::string& fasta_file_name) {

  seqan::SeqFileIn seqFileIn;
  if (!seqan::open(seqFileIn, fasta_file_name.c_str()))
  {
    log.critical("Could not open fasta file: {}", fasta_file_name);
  }

  log.info("Reading Fasta file: {}", fasta_file_name);

  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::Dna5String> seqs;

  try
  {
    readRecords(ids, seqs, seqFileIn);
  }
  catch (seqan::Exception const & e)
  {
    log.critical("Error: {} reading fasta file: {}", e.what(), fasta_file_name);
  }

  for (unsigned i = 0; i < length(ids); ++i) {

    kgl::ContigId_t contig_id = toCString(ids[i]);
    Sequence_t sequence;
    seqan::move(sequence, seqs[i]);
    addSequence(contig_id, sequence);
    log.info("Fasta Contig id: {}; Contig sequence length: {} ", contig_id, sequence.length());

  }

}




