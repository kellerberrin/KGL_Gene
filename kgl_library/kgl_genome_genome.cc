//
// Created by kellerberrin on 12/11/17.
//


#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_gaf_parser.h"
#include "kgl_gff_fasta.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::GenomeDatabase> kgl::GenomeDatabase::createGenomeDatabase(const std::string& fasta_file,
                                                                               const std::string& gff_file,
                                                                               const std::string& gaf_file,
                                                                               const std::string& translation_table) {
  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = ParseGffFasta().readFastaGffFile(fasta_file, gff_file);

  // Optionally set the translation table (defaults to the standard table).
  if (not translation_table.empty()) {

    // Set the amino translation table
    genome_db_ptr->setTranslationTable(translation_table);

  }

  // Wire-up the genome database.
  genome_db_ptr->createVerifyGenomeDatabase();

  // Optionally read in gaf records to add into the genome database.
  if (not gaf_file.empty()) {

    genome_db_ptr->gene_ontology_.readGafFile(gaf_file);

  }

  // return a const pointer.
  return genome_db_ptr;

}



void kgl::GenomeDatabase::readAuxillary(std::shared_ptr<GenomeDatabase> genome_db_ptr, const std::string& tss_gff_file) {

  // Optionally read in tss_gff records into the genome database.
  if (not tss_gff_file.empty()) {

    ParseGffFasta().readTssGffFile(tss_gff_file, genome_db_ptr);

  }

  // Wire-up the genome auxillary database.
  genome_db_ptr->createVerifyAuxillary();

}


bool kgl::GenomeDatabase::addContigSequence(const kgl::ContigId_t& contig_id,
                                            std::shared_ptr<kgl::DNA5SequenceContig> sequence_ptr) {

  using ContigPtr = std::shared_ptr<kgl::ContigFeatures>;
  ContigPtr contig_ptr(std::make_shared<kgl::ContigFeatures>(contig_id, sequence_ptr));

  auto result = genome_sequence_map_.insert(std::make_pair(contig_id, std::move(contig_ptr)));

  return result.second;

}

bool kgl::GenomeDatabase::getContigSequence(const kgl::ContigId_t& contig_id,
                                            std::shared_ptr<const ContigFeatures>& contig_ptr) const {

  auto result_iter = genome_sequence_map_.find(contig_id);

  if (result_iter != genome_sequence_map_.end()) {

    contig_ptr = result_iter->second;
    return true;

  }

  return false;

}


void kgl::GenomeDatabase::createVerifyGenomeDatabase() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->verifyFeatureHierarchy();

  }

}


void kgl::GenomeDatabase::createVerifyAuxillary() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->verifyAuxillaryHierarchy();

  }

}


void kgl::GenomeDatabase::setTranslationTable(const std::string& table) {

  ExecEnv::log().info("GenomeDatabase::setTranslationTable; All contigs set to Amino translation table: {}", table);

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->setTranslationTable(table);

  }

}


// Given a sequence offset, returns a contig offset.
bool kgl::GenomeDatabase::contigOffset( const ContigId_t& contig_id,
                                        const FeatureIdent_t& gene_id,
                                        const FeatureIdent_t& sequence_id,
                                        ContigOffset_t sequence_offset,
                                        ContigOffset_t& contig_offset) const {

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("contigOffset(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("contigOffset(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }

  ContigSize_t coding_sequence_length;
  return SequenceOffset::refCodingSequenceContigOffset(coding_sequence_ptr, sequence_offset, contig_offset, coding_sequence_length);

}
