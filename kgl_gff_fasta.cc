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
#include <boost/tokenizer.hpp>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ParseGffFasta::GffFastaImpl does all the heavy lifting using 3rd a party library. In this case; Seqan.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class kgl::ParseGffFasta::GffFastaImpl {

public:

  GffFastaImpl(kgl::Logger& logger) : log(logger) {}
  ~GffFastaImpl() = default;

  void readGffFile(const std::string &gff_file_name, kgl::GenomeSequences& genome_sequences);
  std::unique_ptr<kgl::GenomeSequences> readFastaFile(const std::string& fasta_file_name);

private:

  Logger& log; // Should declared first.

  bool parseGffRecord(const kgl::GenomeSequences& genome_sequences,
                      seqan::GffRecord& record,
                      long gff_line_counter);


};


std::unique_ptr<kgl::GenomeSequences> kgl::ParseGffFasta::GffFastaImpl::readFastaFile(const std::string& fasta_file_name) {

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

  std::unique_ptr<kgl::GenomeSequences> genome_sequences_ptr(std::make_unique<kgl::GenomeSequences>());

  for (unsigned i = 0; i < length(ids); ++i) {

    std::string id_line = toCString(ids[i]);
    ContigId_t  contig_id = id_line.substr(0, id_line.find_first_of(" \t,")); // Only the identifier.
    Sequence_t sequence;
    seqan::move(sequence, seqs[i]);

    if (not genome_sequences_ptr->addContigSequence(contig_id, sequence)) {

      log.error("addContigSequence(), Attempted to add duplicate contig; {}", contig_id);

    }

    log.info("Fasta Contig id: {}; Contig sequence length: {} ", contig_id, sequence.length());

  }

  return genome_sequences_ptr;

}


void kgl::ParseGffFasta::GffFastaImpl::readGffFile(const std::string &gff_file_name,
                                                   kgl::GenomeSequences& genome_sequences) {

  seqan::GffFileIn gff_file_in;
  if (!seqan::open(gff_file_in, gff_file_name.c_str())) {

    log.critical("Could not open Gff file: {}", gff_file_name);

  }

  log.info("Reading Gff file: {}", gff_file_name);

  seqan::GffFileOut gff_file_out(std::cout, seqan::Gff());

  seqan::GffRecord record;

  try {

    long gff_line_counter = 0;

    while(!seqan::atEnd(gff_file_in)) {

      seqan::readRecord(record, gff_file_in);

      ++gff_line_counter;

      if (not parseGffRecord(genome_sequences, record, gff_line_counter)) {

        log.error("Parse error, unable to load feature at Gff line: {}", gff_line_counter);
        seqan::writeRecord(gff_file_out, record);

      }

    }

    log.info("Processed: {} Gff feature records", gff_line_counter);

  } catch(seqan::Exception const & e) {

    log.critical("Error: {} reading fasta file: {}", e.what(), gff_file_name);

  }

}



bool kgl::ParseGffFasta::GffFastaImpl::parseGffRecord(const kgl::GenomeSequences& genome_sequences,
                                                      seqan::GffRecord& record,
                                                      long gff_line_counter) {
  // Get the attributes.
  kgl::FeatureAttributes record_attributes;
  for (unsigned i = 0; i < length(record.tagNames); i++) {

    std::string key = toCString(record.tagNames[i]);
    std::string value = toCString(record.tagValues[i]);


    bt::char_separator<char> sep(",");
    bt::tokenizer<bt::char_separator<char>> tokenize(value, sep);
    for(auto iter = tokenize.begin(); iter != tokenize.end(); ++iter) {

      record_attributes.insertAttribute(key, *iter);

    }


  }
  // Create a sequence object.
  ContigOffset_t begin = record.beginPos;
  ContigOffset_t end = record.endPos;
  kgl::StrandSense strand;
  // Check validity of the strand character
  switch(record.strand) {

    case '+':
    case '-':
    case '.':
      strand = static_cast<kgl::StrandSense>(record.strand);
      break;

    default:
      log.error("Strand Sense character: {} is not one of ['+', '-', '.']", record.strand);
      return false;

  }
  FeatureSequence sequence (begin, end, strand);
  // Get the feature type and convert to upper case.
  kgl::FeatureType_t type = toCString(record.type);
  std::transform(type.begin(), type.end(), type.begin(), ::toupper);
  // Get (or construct) the feature ID.
  std::vector<kgl::FeatureIdent_t> feature_id_vec;
  kgl::FeatureIdent_t feature_id;
  if (not record_attributes.getIds(feature_id_vec)) {

    // Construct an id
    feature_id = type + std::to_string(begin);
    log.warn("Gff line: {}, 'ID' key not found; ID: {} generated", gff_line_counter, feature_id);

  } else if (feature_id_vec.size() > 1) {

    log.warn("Gff line: {}, Has {} 'ID' values, choosing first value"
        , gff_line_counter, feature_id_vec.size());

    feature_id = feature_id_vec[0];

  } else {

    feature_id = feature_id_vec[0];

  }
  // Get the contig id.
  kgl::ContigId_t contig_id = toCString(record.ref);
  // Get a pointer to the contig.
  std::shared_ptr<kgl::ContigRecord> contig_ptr;
  if (not genome_sequences.getContigSequence(contig_id, contig_ptr)) {

    log.error("Could not find contig: {}", contig_id);
    return false;

  }
  // Check for a valid phase field
  CDSPhaseType_t phase;
  bool valid_phase = false;
  switch(record.phase) {

    case '0':
      phase = 0;
      valid_phase = true;
      break;

    case '1':
      phase = 1;
      valid_phase = true;
      break;

    case '2':
      phase = 2;
      valid_phase = true;
      break;

    case '.':
      break;

    default:
      log.warn("Unexpected phase code: {} in Gff file should be one of '0', '1', '2', '.'", record.phase);
      break;


  }
  // Check that the type field contains "CDS"
  if (valid_phase) {

    if (type.find(CDSRecord::CDS_TYPE) == std::string::npos) {

      log.warn("Mis-match between valid phase: {} and record type: {}", phase, type);

    }

  }

  std::shared_ptr<kgl::FeatureRecord> feature_ptr;
  if (valid_phase) {
    // Create a CDS Feature.
    feature_ptr = std::make_shared<kgl::CDSRecord>(feature_id, phase, contig_ptr, sequence);
  }
  else if (type.find(GeneRecord::GENE_TYPE) != std::string::npos) {
    // Create a GENE feature
    feature_ptr = std::make_shared<kgl::GeneRecord>(feature_id, contig_ptr, sequence);

  } else {
    // Create a general feature
    feature_ptr = std::make_shared<kgl::FeatureRecord>(feature_id, type, contig_ptr, sequence);

  }
  // Add in the attributes.
  feature_ptr->setAttributes(record_attributes);
  // Annotate the contig.
  if (not contig_ptr->addFeature(feature_ptr)) {

    log.error("Could not add duplicate feature: {} to contig: {}", feature_id, contig_id);
    return false;

  }

  return true;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ParseGffFasta is a public facade class that passes the functionality onto ParseGffFasta::GffFastaImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::ParseGffFasta::ParseGffFasta(kgl::Logger& logger)
    : log(logger), gff_fasta_impl_ptr_(std::make_unique<kgl::ParseGffFasta::GffFastaImpl>(logger)) {}


kgl::ParseGffFasta::~ParseGffFasta() {}

// Functionality passed to the implmentation.
std::unique_ptr<kgl::GenomeSequences> kgl::ParseGffFasta::readFastaFile(const std::string& fasta_file_name) {

  return gff_fasta_impl_ptr_->readFastaFile(fasta_file_name);

}

// Functionality passed to the implmentation.
void kgl::ParseGffFasta::readGffFile(const std::string& gff_file_name, GenomeSequences& genome_sequences) {

  gff_fasta_impl_ptr_->readGffFile(gff_file_name, genome_sequences);

}


std::unique_ptr<kgl::GenomeSequences> kgl::ParseGffFasta::readFastaGffFile(const std::string& fasta_file_name,
                                                                            const std::string& gff_file_name) {

  std::unique_ptr<kgl::GenomeSequences> genome_sequences_ptr = readFastaFile(fasta_file_name);
  readGffFile(gff_file_name, *genome_sequences_ptr);

  return genome_sequences_ptr;

}


