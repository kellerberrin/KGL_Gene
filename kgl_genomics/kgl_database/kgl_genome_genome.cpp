//
// Created by kellerberrin on 12/11/17.
//


#include "kel_exec_env.h"
#include "kel_patterns.h"
#include "kgl_gaf_parser.h"
#include "kgl_gff_fasta_seqan.h"
#include "kgl_sequence_offset.h"

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
  std::shared_ptr<kgl::GenomeReference> genome_db_ptr = ParseGffFastaSeqan().readFastaGffFile(organism, fasta_file, gff_file);

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeReference> new_genome_db_ptr = ParseGffFasta::readFastaGffFile(organism, fasta_file, gff_file);

  if (genome_db_ptr->compareReference(*new_genome_db_ptr)) {

    ExecEnv::log().info("GenomeReference::createGenomeDatabase; Equivalent Seqan and Native Genome Database: {}", organism);

  } else {

    ExecEnv::log().error("GenomeReference::createGenomeDatabase; mis-match between Seqan and Native Genome Database: {}", organism);

  }

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
                                             std::shared_ptr<kgl::DNA5SequenceContig> sequence_ptr) {

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

      ExecEnv::log().error("setTranslationTable(), Could not set translation table: {} for contig: {}", table, contig_id);

    }

  }

}


// Given a sequence offset, returns a contig offset.
bool kgl::GenomeReference::contigOffset(const ContigId_t& contig_id,
                                        const FeatureIdent_t& gene_id,
                                        const FeatureIdent_t& sequence_id,
                                        ContigOffset_t sequence_offset,
                                        ContigOffset_t& contig_offset) const {

  // Get the contig.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("contigOffset(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_opt.value()->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("contigOffset(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }

  ContigSize_t coding_sequence_length;
  return SequenceOffset::refCodingSequenceContigOffset(coding_sequence_ptr, sequence_offset, contig_offset, coding_sequence_length);

}

bool kgl::GenomeReference::compareReference(const GenomeReference& compare_genome) const {

  bool equivalent_contigs{true};
  for (auto const& [contig_id, contig_reference] : genome_sequence_map_) {

    auto contig_opt = compare_genome.getContigSequence(contig_id);
    if (contig_opt) {

      if (not contig_opt.value()->compareContig(*contig_reference)) {

        ExecEnv::log().warn("GenomeReference::compareReference; contig id: {} has differing contig sequences", contig_id);
        equivalent_contigs = false;

      }

    } else {

      ExecEnv::log().warn("GenomeReference::compareReference; contig id: {} not found in comparison contig", contig_id);
      equivalent_contigs = false;

    }

  }

  for (auto const& [contig_id, contig_reference] : compare_genome.genome_sequence_map_) {

    auto contig_opt = getContigSequence(contig_id);
    if (not contig_opt) {

      ExecEnv::log().warn("GenomeReference::compareReference; comparison contig: {} not found", contig_id);
      equivalent_contigs = false;

    }

  }

  return equivalent_contigs;

}
