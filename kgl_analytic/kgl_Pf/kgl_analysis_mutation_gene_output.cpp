//
// Created by kellerberrin on 27/1/21.
//



#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_mutation.h"

#include <fstream>


namespace kgl = kellerberrin::genome;




// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::writeOutput(const std::shared_ptr<const GenomePEDData>& ped_data, const std::string& output_file_name, char output_delimiter) const {

  const double homozygous_bias{0.1};

  std::ofstream out_file(output_file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("GenomeMutation::writeOutput; could not open file: {} for output", output_file_name);
    return false;

  } else {

    ExecEnv::log().info("GenomeMutation writing output to file: {}", output_file_name);

  }

  if (not gene_vector_.empty()) {

    writeHeader(ped_data, out_file, output_delimiter, gene_vector_.front().gene_characteristic, gene_vector_.front().clinvar);

  }

  for (auto const& gene : gene_vector_) {


    gene.gene_characteristic.writeGene(out_file, output_delimiter);

    out_file << output_delimiter;
    out_file << gene.unique_variants << output_delimiter
        << gene.variant_count << output_delimiter
        << gene.span_variant_count << output_delimiter
        << gene.female_lof << output_delimiter
        << gene.male_lof << output_delimiter
        << gene.hom_lof << output_delimiter;

    out_file << gene.female_high_effect << output_delimiter
             << gene.male_high_effect << output_delimiter
             << gene.hom_high_effect << output_delimiter
             << gene.genome_count << output_delimiter
             << gene.genome_variant << output_delimiter
             << (static_cast<double>(gene.span_variant_count) / static_cast<double>(gene.gene_characteristic.geneSpan())) << output_delimiter;

    if (gene.gene_characteristic.nucleotides() > 0) {

      out_file << (static_cast<double>(gene.variant_count) / static_cast<double>(gene.gene_characteristic.nucleotides())) << output_delimiter;

    } else {

      out_file << 0 << output_delimiter;

    }

    out_file << gene.heterozygous << output_delimiter
             << gene.homozygous << output_delimiter;
    double ratio = static_cast<double>(gene.heterozygous) / (static_cast<double>(gene.homozygous) + homozygous_bias);
    out_file << ratio << output_delimiter
             << (gene.indel * 100) << output_delimiter
             << (gene.transition * 100) << output_delimiter
             << (gene.transversion * 100) << output_delimiter;

    gene.clinvar.writeOutput(ped_data, out_file, output_delimiter);

    out_file << '\n';

  } // Gene

  return true;

}


void kgl::GeneCharacteristic::writeGene(std::ostream& out_file, char output_delimiter) const {

  out_file << genome << output_delimiter
           << contig << output_delimiter
           << gene_id << output_delimiter
           << gene_name << output_delimiter
           << description << output_delimiter
           << biotype << output_delimiter
           << (valid_protein ? "Valid" : "Invalid") << output_delimiter
           << gaf_id << output_delimiter
           << gene_begin << output_delimiter
           << gene_end << output_delimiter
           << gene_span_ << output_delimiter
      << strand << output_delimiter
      << sequences << output_delimiter
      << seq_name << output_delimiter
      << nucleotides_ << output_delimiter
           << exons << output_delimiter
           << attributes;

}




void kgl::GeneCharacteristic::writeGeneHeader(std::ostream& out_file, char output_delimiter) const {

  out_file << "Genome" << output_delimiter
           << "Contig" << output_delimiter
           << "Gene" << output_delimiter
           << "Name" << output_delimiter
           << "Description" << output_delimiter
           << "BioType" << output_delimiter
           << "ValidProtein" << output_delimiter
           << "GafId" <<  output_delimiter
           << "Begin" << output_delimiter
           << "End" << output_delimiter
           << "Span" << output_delimiter
           << "Strand" << output_delimiter
           << "Sequences" << output_delimiter
           << "SeqName" << output_delimiter
           << "Nucleotides" << output_delimiter
           << "Exons" << output_delimiter
           << "Attributes";

}



void kgl::GenomeMutation::writeHeader(const std::shared_ptr<const GenomePEDData>& ped_data,
                                      std::ostream& out_file,
                                      char output_delimiter,
                                      const GeneCharacteristic& gene_characteristic,
                                      const GeneClinvar& clinvar) {

  gene_characteristic.writeGeneHeader(out_file, output_delimiter);

  out_file << output_delimiter
           << "UniqueVariants" << output_delimiter
           << "VariantCount" << output_delimiter
           << "SpanVariantCount" << output_delimiter
           << "FemalePhaseLoF" << output_delimiter
           << "MalePhaseLoF" << output_delimiter
           << "HomLoF" << output_delimiter
           << "FemaleHigh" << output_delimiter
           << "MaleHigh" << output_delimiter
           << "HomHigh" << output_delimiter
           << "GenomeCount" << output_delimiter
           << "GenomeVariant" << output_delimiter
           << "SpanDensity" << output_delimiter
           << "VariantDensity" << output_delimiter
           << "Heterozygous" << output_delimiter
           << "Homozygous" << output_delimiter
           << "Het/Hom" << output_delimiter
           << "Indel%" << output_delimiter
           << "Transition%" << output_delimiter
           << "Transversion%" << output_delimiter;

           clinvar.writeHeader(ped_data, out_file, output_delimiter);

  out_file << '\n';

}

