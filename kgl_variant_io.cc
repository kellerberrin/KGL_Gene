//
// Created by kellerberrin on 12/11/17.
//


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_filter.h"
#include "kgl_gff_fasta.h"

namespace kgl = kellerberrin::genome;


bool kgl::GenomeVariant::outputCSV(const std::string& file_name, VariantOutputIndex output_index, bool detail) const {

  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << output(',', output_index, detail);

  return out_file.good();

}

bool kgl::GenomeVariant::mutantProtein( const ContigId_t& contig_id,
                                        const FeatureIdent_t& gene_id,
                                        const FeatureIdent_t& sequence_id,
                                        const std::shared_ptr<const GenomeDatabase>& genome_db,
                                        std::shared_ptr<AminoSequence>& amino_sequence) const {
  // Get the contig.
  std::shared_ptr<ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }


  // Get the DNA sequence
  std::shared_ptr<DNA5SequenceCoding> dna_sequence_ptr;
  if (not contig_ptr->getDNA5SequenceCoding(coding_sequence_ptr, dna_sequence_ptr))  {

    ExecEnv::log().warn("No valid DNA sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  // Extract the variants for processing.
  OffsetVariantMap variant_map;
  getCodingSortedVariants(contig_id, coding_sequence_ptr->start(), coding_sequence_ptr->end(), variant_map);

  if (not GenomeVariant::mutateDNA(variant_map, sequence_id, dna_sequence_ptr)) {

    ExecEnv::log().warn("Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                        contig_id, gene_id, sequence_id);
    return false;

  }

  amino_sequence = contig_ptr->getAminoSequence(dna_sequence_ptr);
  return true;

}

// Perform the actual mutation.
bool kgl::GenomeVariant::mutateDNA(const OffsetVariantMap& variant_map,
                                   const FeatureIdent_t& sequence_id,
                                   std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr) {

  // Mutate the base sequence.
  for (const auto& variant : variant_map) {

    if (not variant.second->mutateCodingSequence(sequence_id, dna_sequence_ptr)) {

      ExecEnv::log().warn("mutateDNA(), problem with variant contig: {}, offset: {}",
                          variant.second->contigId(), variant.second->offset());
      return false;

    }

  }

  return true;

}


std::ostream& operator<<(std::ostream &os, const kgl::GenomeVariant& genome_variant) {

  os << genome_variant.output(' ', kgl::VariantOutputIndex::START_1_BASED, false);
  os.flush();

  return os;

}

