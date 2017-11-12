//
// Created by kellerberrin on 12/11/17.
//


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_gff_fasta.h"

namespace kgl = kellerberrin::genome;


bool kgl::GenomeVariant::outputCSV(const std::string& file_name, VariantOutputIndex output_index) const {

  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << output(',', output_index);

  return out_file.good();

}


bool kgl::GenomeVariant::writeMutantProtein(const std::string& fasta_file,
                                            const std::string& sequence_name,
                                            const ContigId_t& contig_id,
                                            const FeatureIdent_t& gene_id,
                                            const FeatureIdent_t& sequence_id,
                                            const std::shared_ptr<const GenomeDatabase>& genome_db) const {

  std::shared_ptr<DNA5Sequence> sequence_ptr;
  if (genome_db->getDNA5Sequence(contig_id, gene_id, sequence_id, sequence_ptr)) {

    std::shared_ptr<ContigFeatures> contig_ptr;
    if (genome_db->getContigSequence(contig_id, contig_ptr)) {

      std::shared_ptr<AminoSequence> amino_sequence = contig_ptr->getAminoSequence(sequence_ptr);

      FastaSequence fasta_sequence;
      fasta_sequence.first = sequence_name;
      fasta_sequence.second = amino_sequence;
      return ParseGffFasta().writeFastaFile(fasta_file, std::vector<FastaSequence>{fasta_sequence});

    } else {

      ExecEnv::log().warn("Could not find contig: {} in genome database", contig_id);
      return false;

    }

  } else {

    ExecEnv::log().warn("No valid sequence for contig: {}, gene: {}, sequence id: {} Fasta: {} not written",
                        contig_id, gene_id, sequence_id, fasta_file);
    return false;

  }

}


std::ostream& operator<<(std::ostream &os, const kgl::GenomeVariant& genome_variant) {

  os << genome_variant.output(' ', kgl::VariantOutputIndex::START_1_BASED);
  os.flush();

  return os;

}

