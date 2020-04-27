//
// Created by kellerberrin on 12/11/17.
//


#include "kel_exec_env.h"
#include "kel_patterns.h"
#include "kgl_gaf_parser.h"
#include "kgl_gff_fasta.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::shared_ptr<kgl::GenomeDatabase> kgl::GenomeDatabase::createGenomeDatabase(const kgl::RuntimeProperties& runtime_options,
                                                                               const GenomeId_t& organism) {

  // Get the genome database runtime parameters.
  std::string fasta_file, gff_file, gaf_file, amino_translation_table;
  runtime_options.getGenomeDBFiles(organism , fasta_file, gff_file, gaf_file, amino_translation_table);



  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = createGenomeDatabase(organism,
                                                                            fasta_file,
                                                                            gff_file,
                                                                            gaf_file,
                                                                            amino_translation_table);

  if (not genome_db_ptr) {

    ExecEnv::log().critical("GenomeDatabase::createGenomeDatabase; severe error NULL genome object for organism: {}", organism);

  }

  if (not genome_db_ptr->readGenomeAuxiliary(runtime_options)) {

    ExecEnv::log().info("GenomeDatabase::createGenomeDatabase; Problem reading auxiliary genome data for organism: {}", organism);

  }

  ExecEnv::log().info("**** Successfully created a Genome Database for organism: {} ****", genome_db_ptr->genomeId());

  return genome_db_ptr;

}


std::shared_ptr<kgl::GenomeDatabase> kgl::GenomeDatabase::createGenomeDatabase(const std::string& organism,
                                                                               const std::string& fasta_file,
                                                                               const std::string& gff_file,
                                                                               const std::string& gaf_file,
                                                                               const std::string& translation_table) {
  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = ParseGffFasta().readFastaGffFile(organism,  fasta_file, gff_file);

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


bool kgl::GenomeDatabase::readGenomeAuxiliary(const RuntimeProperties& runtime_options) {

  std::vector<AuxFileInfo> aux_file_list;

  if (not runtime_options.getGenomeAuxFiles(genomeId(), aux_file_list)) {

    ExecEnv::log().error("GenomeDatabase::readGenomeAuxiliary; problem reading runtime genome auxiliary file list for organism: {}", genomeId());
    return false;

  }

  for (auto aux_file : aux_file_list) {

    if (aux_file.auxType() == AuxFileInfo::ADJALLEY_TSS_GFF_) {

      readAuxillary(aux_file.fileName());

    } else {

      ExecEnv::log().error("GenomeDatabase::readGenomeAuxiliary; genome aux file type: {} not supported for organism: {}", aux_file.auxType(), genomeId());

    }

  }

  return true;

}

void kgl::GenomeDatabase::readAuxillary(const std::string& tss_gff_file) {

  // Optionally read in tss_gff records into the genome database.
  if (not tss_gff_file.empty()) {

    ParseGffFasta().readTssGffFile(tss_gff_file, *this);

  }

  // Wire-up the genome auxillary database.
  createVerifyAuxillary();

}


bool kgl::GenomeDatabase::addContigSequence(const kgl::ContigId_t& contig_id,
                                            const std::string& description,
                                            std::shared_ptr<kgl::DNA5SequenceContig> sequence_ptr) {

  using ContigPtr = std::shared_ptr<kgl::ContigFeatures>;
  ContigPtr contig_ptr(std::make_shared<kgl::ContigFeatures>(contig_id, sequence_ptr));
  contig_ptr->description(description);

  auto result = genome_sequence_map_.insert(std::make_pair(contig_id, std::move(contig_ptr)));

  return result.second;

}

std::optional<std::shared_ptr<const kgl::ContigFeatures>> kgl::GenomeDatabase::getContigSequence(const kgl::ContigId_t& contig_id) const {

  auto result_iter = genome_sequence_map_.find(contig_id);

  if (result_iter != genome_sequence_map_.end()) {

    return result_iter->second;

  }

  return std::nullopt;

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

  for (auto [contig_id, contig_ptr] : genome_sequence_map_) {

    if (not contig_ptr->setTranslationTable(table)) {

      ExecEnv::log().error("setTranslationTable(), Could not set translation table: {} for contig: {}", table, contig_id);

    }

  }

}


// Given a sequence offset, returns a contig offset.
bool kgl::GenomeDatabase::contigOffset( const ContigId_t& contig_id,
                                        const FeatureIdent_t& gene_id,
                                        const FeatureIdent_t& sequence_id,
                                        ContigOffset_t sequence_offset,
                                        ContigOffset_t& contig_offset) const {

  // Get the contig.
  std::optional<std::shared_ptr<const ContigFeatures>> contig_opt = getContigSequence(contig_id);
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeCollection - A map of different organism genomes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::GenomeCollection::getGenome(const GenomeId_t& genome_id, std::shared_ptr<const GenomeDatabase>& genome_database) const {

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    genome_database = result->second;
    return true;

  } else {

    genome_database = nullptr;
    return false;

  }

}


bool kgl::GenomeCollection::addGenome(std::shared_ptr<const GenomeDatabase> genome_database) {

  auto result = genome_map_.insert(std::pair<GenomeId_t, std::shared_ptr<const GenomeDatabase>>(genome_database->genomeId(), genome_database));

  return result.second;

}


std::shared_ptr<const kgl::GenomeDatabase> kgl::GenomeCollection::getGenome(const std::string& GenomeID) const {

  std::shared_ptr<const GenomeDatabase> genome_ptr;

  if (not getGenome(GenomeID, genome_ptr)) {

    ExecEnv::log().critical("GenomeCollection::getGenome; genome: {} not found", GenomeID);

  }

  return genome_ptr;

}



std::shared_ptr<kgl::GenomeCollection> kgl::GenomeCollection::createGenomeCollection(const RuntimeProperties& runtime_options) {

  std::vector<std::string> genome_list;
  if (not runtime_options.getActiveGenomes(genome_list)) {

    ExecEnv::log().error("GenomeCollection::createGenomeCollection; Problem retrieving runtime genome list");

  }

  std::shared_ptr<kgl::GenomeCollection> genome_collection(std::make_shared<kgl::GenomeCollection>());

  for (auto genome : genome_list) {

    // Create the genome database.
    std::shared_ptr<GenomeDatabase> genome_ptr = GenomeDatabase::createGenomeDatabase(runtime_options, genome);

    if (not genome_collection->addGenome(genome_ptr)) {

      ExecEnv::log().error("GenomeCollection::createGenomeCollection; Unable to add Genome Database: {} (probable duplicate)", genome_ptr->genomeId());

    }

  }

  return genome_collection;

}

