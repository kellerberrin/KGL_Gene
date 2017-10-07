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


void kgl::ParseGFFSFasta::readGffFile(const std::string &gff_file_name, kgl::GenomeSequences& genome_sequences) {

  seqan::GffFileIn gff_file_in;
  if (!seqan::open(gff_file_in, gff_file_name.c_str())) {

    log.critical("Could not open Gff file: {}", gff_file_name);

  }

  log.info("Reading Gff file: {}", gff_file_name);

  seqan::GffFileOut gff_file_out(std::cout, seqan::Gff());

  seqan::GffRecord record;

  try {

    while(!seqan::atEnd(gff_file_in)) {

      seqan::readRecord(record, gff_file_in);
//      seqan::writeRecord(gff_file_out, record);

    }

  } catch(seqan::Exception const & e) {

    log.critical("Error: {} reading fasta file: {}", e.what(), gff_file_name);

  }

}



std::unique_ptr<kgl::GenomeSequences> kgl::ParseGFFSFasta::readFastaFile(const std::string& fasta_file_name) {

  seqan::SeqFileIn seq_file_in;
  if (!seqan::open(seq_file_in, fasta_file_name.c_str())) {

    log.critical("Could not open fasta file: {}", fasta_file_name);

  }

  log.info("Reading Fasta file: {}", fasta_file_name);

  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::Dna5String> seqs;

  try {

    readRecords(ids, seqs, seq_file_in);

  }
  catch (seqan::Exception const & e) {

    log.critical("Error: {} reading fasta file: {}", e.what(), fasta_file_name);

  }

  std::unique_ptr<kgl::GenomeSequences> genome_sequences_ptr(std::make_unique<kgl::GenomeSequences>(log));

  for (unsigned i = 0; i < length(ids); ++i) {

    kgl::ContigId_t contig_id = toCString(ids[i]);
    Sequence_t sequence;
    seqan::move(sequence, seqs[i]);
    genome_sequences_ptr->addContigSequence(contig_id, sequence);
    log.info("Fasta Contig id: {}; Contig sequence length: {} ", contig_id, sequence.length());

  }

  return genome_sequences_ptr;

}


std::unique_ptr<kgl::GenomeSequences> kgl::ParseGFFSFasta::readFastaGffFile(const std::string& fasta_file_name,
                                                                            const std::string& gff_file_name) {

  std::unique_ptr<kgl::GenomeSequences> genome_sequences_ptr = readFastaFile(fasta_file_name);
  readGffFile(gff_file_name, *genome_sequences_ptr);

  return genome_sequences_ptr;

}


