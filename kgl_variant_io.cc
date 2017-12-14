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

bool kgl::GenomeVariant::mutantProtein( const std::string& sequence_name,
                                        const ContigId_t& contig_id,
                                        const FeatureIdent_t& gene_id,
                                        const FeatureIdent_t& sequence_id,
                                        const std::shared_ptr<const GenomeDatabase>& genome_db,
                                        std::shared_ptr<AminoSequence>& amino_sequence) const {

  // Filter the variants, assumes that gene is only unique to contig and sequence is only unique to gene.
  std::shared_ptr<GenomeVariant> variant_ptr = filterVariants(ContigFilter(contig_id));
  variant_ptr = variant_ptr->filterVariants(GeneFilter(gene_id));
  variant_ptr = variant_ptr->filterVariants(SequenceFilter(sequence_id));
  std::cout << *variant_ptr; // debug.
  // Extract the variants for processing.
  std::vector<std::shared_ptr<const Variant>> variant_vector;
  variant_ptr->getVariants(variant_vector);

  std::shared_ptr<DNA5SequenceCoding> sequence_ptr;
  if (genome_db->getDNA5SequenceCoding(contig_id, gene_id, sequence_id, sequence_ptr)) {

    // Mutate the base sequence.
    for (const auto& variant : variant_vector) {

      if (not variant->mutateCodingSequence(sequence_id, sequence_ptr)) {

        ExecEnv::log().warn("mutateCodingSequence(), problem with variant contig: {}, offset: {}",
                            variant->contigId(), variant->offset());
        return false;

      }

    }

    std::shared_ptr<ContigFeatures> contig_ptr;
    if (genome_db->getContigSequence(contig_id, contig_ptr)) {

      amino_sequence = contig_ptr->getAminoSequence(sequence_ptr);
      return true;

    } else {

      ExecEnv::log().warn("Could not find contig: {} in genome database", contig_id);
      return false;

    }

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, gene: {}, sequence id: {}",
                        contig_id, gene_id, sequence_id);
    return false;

  }

}



std::ostream& operator<<(std::ostream &os, const kgl::GenomeVariant& genome_variant) {

  os << genome_variant.output(' ', kgl::VariantOutputIndex::START_1_BASED, false);
  os.flush();

  return os;

}

