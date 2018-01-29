//
// Created by kellerberrin on 12/11/17.
//


#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_gaf_parser.h"
#include "kgl_gff_fasta.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<const kgl::GenomeDatabase> kgl::GenomeDatabase::createGenomeDatabase(const std::string& fasta_file,
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

  setupFeatureHierarchy();
  verifyFeatureHierarchy();

}

void kgl::GenomeDatabase::setupFeatureHierarchy() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->setupFeatureHierarchy();

  }

}


void kgl::GenomeDatabase::verifyFeatureHierarchy() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->verifyFeatureHierarchy();

  }

}

void kgl::GenomeDatabase::setTranslationTable(const std::string& table) {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->setTranslationTable(table);

  }

}


void kgl::GenomeDatabase::registerContigData(std::shared_ptr<kgl::ContigCountData>& contig_data_ptr) const {
// Create data blocks for each contig in the genome database
  for (const auto &contig_pair : genome_sequence_map_) {

    if (not contig_data_ptr->insertContig(contig_pair.first, contig_pair.second->sequence().length())) {

      kgl::ExecEnv::log().error("ContigCountData; attempted to add duplicate contig; {}", contig_pair.first);

    }

  }

}
