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

  void readGffFile(const std::string &gff_file_name, std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr);
  void readFastaFile(const std::string& fasta_file_name, std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr);

private:

  Logger& log; // Should declared first.

  bool parseGffRecord(std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr,
                      seqan::GffRecord& record,
                      long gff_line_counter);


};


void kgl::ParseGffFasta::GffFastaImpl::readFastaFile(const std::string& fasta_file_name,
                                                     std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr) {

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

  for (unsigned i = 0; i < length(ids); ++i) {

    std::string id_line;
    seqan::move(id_line, ids[i]);
    ContigId_t  contig_id = id_line.substr(0, id_line.find_first_of(" \t,")); // Only the identifier.
    DNA5Sequence::SequenceString sequence;
    seqan::move(sequence, seqs[i]);
    std::shared_ptr<DNA5Sequence> sequence_ptr(std::make_shared<DNA5Sequence>(sequence));

    if (not genome_db_ptr->addContigSequence(contig_id, sequence_ptr)) {

      log.error("addContigSequence(), Attempted to add duplicate contig; {}", contig_id);

    }

    log.info("Fasta Contig id: {}; Contig sequence length: {} ", contig_id, sequence_ptr->length());

  }

}


void kgl::ParseGffFasta::GffFastaImpl::readGffFile(const std::string &gff_file_name,
                                                   std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr) {

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

      if (not parseGffRecord(genome_db_ptr, record, gff_line_counter)) {

        log.error("Parse error, unable to load feature at Gff line: {}", gff_line_counter);
        seqan::writeRecord(gff_file_out, record);

      }

    }

    log.info("Processed: {} Gff feature records", gff_line_counter);

  } catch(seqan::Exception const & e) {

    log.critical("Error: {} reading fasta file: {}", e.what(), gff_file_name);

  }

}



bool kgl::ParseGffFasta::GffFastaImpl::parseGffRecord(std::shared_ptr<kgl::GenomeDatabase>& genome_db_ptr,
                                                      seqan::GffRecord& record,
                                                      long gff_line_counter) {
  // Get the attributes.
  kgl::Attributes record_attributes;
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
  std::shared_ptr<kgl::ContigFeatures> contig_ptr;
  if (not genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

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

    if (type.find(CDSFeature::CDS_TYPE) == std::string::npos) {

      log.warn("Mis-match between valid phase: {} and record type: {}", phase, type);

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
std::shared_ptr<kgl::GenomeDatabase> kgl::ParseGffFasta::readFastaFile(const std::string& fasta_file_name) {

  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr(std::make_shared<kgl::GenomeDatabase>());
  gff_fasta_impl_ptr_->readFastaFile(fasta_file_name, genome_db_ptr);
  return genome_db_ptr;

}


std::shared_ptr<kgl::GenomeDatabase> kgl::ParseGffFasta::readFastaGffFile(const std::string& fasta_file_name,
                                                                          const std::string& gff_file_name ) {

  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = readFastaFile(fasta_file_name);
  gff_fasta_impl_ptr_->readGffFile(gff_file_name, genome_db_ptr);

  return genome_db_ptr;

}


