//
// Created by kellerberrin on 3/10/17.
//

#include "kgl_exec_env.h"
#include "kgl_gff_fasta.h"
#include "kgl_sequence_base.h"

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>

#include <boost/tokenizer.hpp>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ParseGffFasta::GffFastaImpl does all the heavy lifting using 3rd a party library. In this case; Seqan.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class kgl::ParseGffFasta::GffFastaImpl {

public:

  GffFastaImpl() = default;
  ~GffFastaImpl() = default;

  void readGffFile(const std::string &gff_file_name, std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr);
  bool writeFastaFile(const std::string& fasta_file_name, const std::vector<WriteFastaSequence>& fasta_sequences);
  void readFastaFile(const std::string& fasta_file_name, std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr);
  bool readFastaFile(const std::string& fasta_file_name, std::vector<ReadFastaSequence>& fasta_sequences);

private:


  bool parseGffRecord(std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr,
                      seqan::GffRecord& record,
                      long gff_line_counter);


};



bool kgl::ParseGffFasta::GffFastaImpl::readFastaFile(const std::string& fasta_file_name,
                                                     std::vector<ReadFastaSequence>& fasta_sequences) {

  seqan::SeqFileIn seq_file_in;
  if (!seqan::open(seq_file_in, fasta_file_name.c_str())) {

    ExecEnv::log().error("Could not open fasta file: {}", fasta_file_name);
    return false;

  }

  ExecEnv::log().info("Reading Fasta file: {}", fasta_file_name);

  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> seqs;

  try {

    readRecords(ids, seqs, seq_file_in);

  }
  catch (seqan::Exception const & e) {

    ExecEnv::log().error("Error: {} reading fasta file: {}", e.what(), fasta_file_name);
    return false;

  }

  for (unsigned i = 0; i < length(ids); ++i) {

    std::string id_line;
    seqan::move(id_line, ids[i]);
    ReadFastaSequence fasta_sequence;
    fasta_sequence.first = id_line.substr(0, id_line.find_first_of(" \t,")); // Only the identifier.
    std::shared_ptr<AlphabetSequence_t> sequence_string_ptr(std::make_shared<AlphabetSequence_t>());
    seqan::move(*sequence_string_ptr, seqs[i]);  // convert from seqan
    fasta_sequence.second = sequence_string_ptr;
    fasta_sequences.push_back(fasta_sequence);

  }

  return true;

}


void kgl::ParseGffFasta::GffFastaImpl::readFastaFile(const std::string& fasta_file_name,
                                                     std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr) {

  std::vector<ReadFastaSequence> fasta_sequences;
  if (not readFastaFile(fasta_file_name, fasta_sequences)) {

    ExecEnv::log().critical("Could not read genome fasta file: {}", fasta_file_name);

  }

  for (auto sequence : fasta_sequences) {

    const std::string& contig_id = sequence.first;
    StringDNA5 DNA5sequence(*sequence.second); // convert to alphabet DNA5.
    std::shared_ptr<DNA5SequenceContig> sequence_ptr(std::make_shared<DNA5SequenceContig>(DNA5sequence));

    if (not genome_db_ptr->addContigSequence(contig_id, sequence_ptr)) {

      ExecEnv::log().error("addContigSequence(), Attempted to add duplicate contig; {}", contig_id);

    }

    ExecEnv::log().info("Fasta Contig id: {}; Contig sequence length: {} ", contig_id, sequence_ptr->length());

  }

}


void kgl::ParseGffFasta::GffFastaImpl::readGffFile(const std::string &gff_file_name,
                                                   std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr) {

  seqan::GffFileIn gff_file_in;
  if (!seqan::open(gff_file_in, gff_file_name.c_str())) {

    ExecEnv::log().critical("Could not open Gff file: {}", gff_file_name);

  }

  ExecEnv::log().info("Reading Gff file: {}", gff_file_name);

  seqan::GffFileOut gff_file_out(std::cout, seqan::Gff());

  seqan::GffRecord record;

  try {

    long gff_line_counter = 0;

    while(!seqan::atEnd(gff_file_in)) {

      seqan::readRecord(record, gff_file_in);

      ++gff_line_counter;

      if (not parseGffRecord(genome_db_ptr, record, gff_line_counter)) {

        ExecEnv::log().error("Parse error, unable to load feature at Gff line: {}", gff_line_counter);
        seqan::writeRecord(gff_file_out, record);

      }

    }

    ExecEnv::log().info("Processed: {} Gff feature records", gff_line_counter);

  } catch(seqan::Exception const & e) {

    ExecEnv::log().critical("Error: {} reading fasta file: {}", e.what(), gff_file_name);

  }

}



bool kgl::ParseGffFasta::GffFastaImpl::parseGffRecord(std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr,
                                                      seqan::GffRecord& record,
                                                      long gff_line_counter) {
  // Get the attributes.
  kgl::Attributes record_attributes;
  for (unsigned i = 0; i < length(record.tagNames); i++) {

    std::string key = seqan::toCString(record.tagNames[i]);
    std::string value = seqan::toCString(record.tagValues[i]);


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
      ExecEnv::log().error("Strand Sense character: {} is not one of ['+', '-', '.']", record.strand);
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
    ExecEnv::log().warn("Gff line: {}, 'ID' key not found; ID: {} generated", gff_line_counter, feature_id);

  } else if (feature_id_vec.size() > 1) {

    ExecEnv::log().warn("Gff line: {}, Has {} 'ID' values, choosing first value"
        , gff_line_counter, feature_id_vec.size());

    feature_id = feature_id_vec[0];

  } else {

    feature_id = feature_id_vec[0];

  }
  // Get the contig id.
  kgl::ContigId_t contig_id = toCString(record.ref);
  // Get a pointer to the contig.
  std::shared_ptr<const kgl::ContigFeatures> contig_ptr;
  if (not genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().error("Could not find contig: {}", contig_id);
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
      phase = 0;
      valid_phase = false;
      break;

    default:
      ExecEnv::log().warn("Unexpected phase code: {} in Gff file should be one of '0', '1', '2', '.'", record.phase);
      phase = 0;
      valid_phase = false;
      break;


  }
  // Check that the type field contains "CDS"
  if (valid_phase) {

    if (type.find(CDSFeature::CDS_TYPE) == std::string::npos) {

      ExecEnv::log().warn("Mis-match between valid phase: {} and record type: {}", phase, type);

    }

  }

  std::shared_ptr<kgl::Feature> feature_ptr;
  if (type.find(CDSFeature::CDS_TYPE) != std::string::npos) {
    // Create a CDS Feature.
    feature_ptr = std::make_shared<kgl::CDSFeature>(feature_id, phase, contig_ptr, sequence);
  }
  else if (type.find(mRNAFeature::MRNA_TYPE) != std::string::npos) {
    // Create a mRNA feature
    feature_ptr = std::make_shared<kgl::mRNAFeature>(feature_id, contig_ptr, sequence);
  }
  else if (type.find(EXONFeature::EXON_TYPE) != std::string::npos) {
    // Create a mRNA feature
    feature_ptr = std::make_shared<kgl::EXONFeature>(feature_id, contig_ptr, sequence);
  }
  else if (type.find(GeneFeature::GENE_TYPE) != std::string::npos) {
    // Create a GENE feature
    feature_ptr = std::make_shared<kgl::GeneFeature>(feature_id, contig_ptr, sequence);

  } else {
    // Create a general feature
    feature_ptr = std::make_shared<kgl::Feature>(feature_id, type, contig_ptr, sequence);

  }
  // Add in the attributes.
  feature_ptr->setAttributes(record_attributes);
  // Annotate the contig.
  std::shared_ptr<kgl::ContigFeatures> mutable_contig_ptr = std::const_pointer_cast<kgl::ContigFeatures>(contig_ptr);
  if (not mutable_contig_ptr->addFeature(feature_ptr)) {

    ExecEnv::log().error("Could not add duplicate feature: {} to contig: {}", feature_id, contig_id);
    return false;

  }

  return true;

}


bool kgl::ParseGffFasta::GffFastaImpl::writeFastaFile(const std::string& fasta_file_name,
                                                      const std::vector<WriteFastaSequence>& fasta_sequences) {

  if (fasta_sequences.empty()) {

    ExecEnv::log().warn("writeFastaFile(), No Fasta sequences to write to file: {}", fasta_file_name);
    return false;

  }

  seqan::SeqFileOut seq_file_out;
  if (!seqan::open(seq_file_out, fasta_file_name.c_str())) {

    ExecEnv::log().error("Could not open fasta file: {}", fasta_file_name);
    return false;

  }


  try {


    for (auto sequence : fasta_sequences) {

      seqan::CharString seqan_id = sequence.first;
      seqan::CharString seqan_seq = sequence.second->getSequenceAsString();
      seqan::writeRecord(seq_file_out, seqan_id, seqan_seq, seqan::Fasta());

    }

  }
  catch(seqan::Exception const & e) {

    ExecEnv::log().error("Error: {} writing to Fasta file: {}", e.what(), fasta_file_name);
    return false;

  }


  return true;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ParseGffFasta is a public facade class that passes the functionality onto ParseGffFasta::GffFastaImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::ParseGffFasta::ParseGffFasta() : gff_fasta_impl_ptr_(std::make_unique<kgl::ParseGffFasta::GffFastaImpl>()) {}
kgl::ParseGffFasta::~ParseGffFasta() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.

// Functionality passed to the implmentation.

std::shared_ptr<kgl::GenomeDatabase> kgl::ParseGffFasta::readFastaFile(const std::string& fasta_file_name) {

  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr(std::make_shared<kgl::GenomeDatabase>());
  gff_fasta_impl_ptr_->readFastaFile(fasta_file_name, genome_db_ptr);
  return genome_db_ptr;

}


bool kgl::ParseGffFasta::writeFastaFile(const std::string& fasta_file_name,
                                        const std::vector<WriteFastaSequence>& fasta_sequences) {

  return gff_fasta_impl_ptr_->writeFastaFile(fasta_file_name, fasta_sequences);

}


std::shared_ptr<kgl::GenomeDatabase> kgl::ParseGffFasta::readFastaGffFile(const std::string& fasta_file_name,
                                                                          const std::string& gff_file_name ) {

  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = readFastaFile(fasta_file_name);
  gff_fasta_impl_ptr_->readGffFile(gff_file_name, genome_db_ptr);
  return genome_db_ptr;

}


bool kgl::ParseGffFasta::readFastaFile(const std::string& fasta_file_name,
                                       std::vector<ReadFastaSequence>& fasta_sequences) {

  return gff_fasta_impl_ptr_->readFastaFile(fasta_file_name, fasta_sequences);

}



