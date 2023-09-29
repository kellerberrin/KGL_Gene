//
// Created by kellerberrin on 12/11/17.
//


#include "kel_exec_env.h"
#include "kgl_gaf_parser.h"
#include "kgl_io_gff_fasta.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::shared_ptr<kgl::GenomeReference> kgl::GenomeReference::createGenomeDatabase(const std::string& organism,
                                                                                 const std::string& fasta_file,
                                                                                 const std::string& gff_file,
                                                                                 const std::string& gaf_file,
                                                                                 const std::string& translation_table) {

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeReference> genome_db_ptr = ParseGffFasta::readFastaGffFile(organism, fasta_file, gff_file);

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


bool kgl::GenomeReference::addContigSequence(const kgl::ContigId_t& contig_id,
                                             const std::string& description,
                                             std::shared_ptr<kgl::DNA5SequenceLinear> sequence_ptr) {

  using ContigPtr = std::shared_ptr<kgl::ContigReference>;
  ContigPtr contig_ptr(std::make_shared<kgl::ContigReference>(contig_id, sequence_ptr));
  contig_ptr->description(description);

  auto result = genome_sequence_map_.insert(std::make_pair(contig_id, std::move(contig_ptr)));

  return result.second;

}

std::optional<std::shared_ptr<const kgl::ContigReference>> kgl::GenomeReference::getContigSequence(const kgl::ContigId_t& contig_id) const {

  auto result_iter = genome_sequence_map_.find(contig_id);

  if (result_iter != genome_sequence_map_.end()) {

    return result_iter->second;

  }

  return std::nullopt;

}


void kgl::GenomeReference::createVerifyGenomeDatabase() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->verifyFeatureHierarchy();

  }

}



void kgl::GenomeReference::setTranslationTable(const std::string& table) {

  ExecEnv::log().info("GenomeReference::setTranslationTable; All contigs set to Amino translation table: {}", table);

  for (auto [contig_id, contig_ptr] : genome_sequence_map_) {

    if (not contig_ptr->setTranslationTable(table)) {

      ExecEnv::log().error("setTranslationTable(), Could not set translation table: {} for contig_ref_ptr: {}", table, contig_id);

    }

  }

}

bool kgl::GenomeReference::equivalent(const GenomeReference& lhs) const {

  bool equivalent_contigs{true};
  for (auto const& [contig_id, contig_reference] : genome_sequence_map_) {

    auto contig_opt = lhs.getContigSequence(contig_id);
    if (contig_opt) {

      if (not contig_opt.value()->equivalent(*contig_reference)) {

        ExecEnv::log().warn("GenomeReference::equivalent; contig id: {} has differing contig_ref_ptr/feature attributes", contig_id);
        equivalent_contigs = false;

      }

    } else {

      ExecEnv::log().warn("GenomeReference::equivalent; contig id: {} not found in comparison contig_ref_ptr", contig_id);
      equivalent_contigs = false;

    }

  }

  for (auto const& [contig_id, contig_reference] : lhs.genome_sequence_map_) {

    auto contig_opt = getContigSequence(contig_id);
    if (not contig_opt) {

      ExecEnv::log().warn("GenomeReference::equivalent; comparison contig_ref_ptr: {} not found", contig_id);
      equivalent_contigs = false;

    }

  }

  return equivalent_contigs;

}
