//
// Created by kellerberrin on 27/1/21.
//



#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_mutation.h"

#include <fstream>


namespace kgl = kellerberrin::genome;




// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::writeOutput(const std::shared_ptr<const GenomePEDData>& ped_data, const std::string& output_file_name, char output_delimiter) const {

  std::ofstream out_file(output_file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("GenomeMutation::writeOutput; could not open file: {} for output", output_file_name);
    return false;

  } else {

    ExecEnv::log().info("GenomeMutation writing output to file: {}", output_file_name);

  }

  if (not gene_vector_.empty()) {

    writeHeader(ped_data, out_file, output_delimiter, gene_vector_.front());

  }

  for (auto const& gene : gene_vector_) {

    gene.gene_characteristic.writeGene(out_file, output_delimiter);

    out_file << output_delimiter;

    gene.ontology.writeOntology(out_file, output_delimiter);

    out_file << output_delimiter;

    gene.gene_variants.writeVariantOutput(ped_data, out_file, output_delimiter);

    out_file << output_delimiter;

    gene.clinvar.writeOutput(ped_data, out_file, output_delimiter);

    out_file << '\n';

  } // Gene

  return true;

}


void kgl::GeneCharacteristic::writeGene(std::ostream& out_file, char output_delimiter) const {

  out_file << genome_ << output_delimiter
      << contig_ << output_delimiter
      << HGNC_id_ << output_delimiter
      << ensembl_id_ << output_delimiter
      << gene_name_ << output_delimiter
      << description_ << output_delimiter
      << biotype_ << output_delimiter
      << (valid_protein_ ? "Valid" : "Invalid") << output_delimiter
      << gaf_id_ << output_delimiter;

  out_file << gene_begin_ << output_delimiter
      << gene_end_ << output_delimiter
      << gene_span_ << output_delimiter
      << strand_ << output_delimiter
      << sequences_ << output_delimiter
      << seq_name_ << output_delimiter
      << nucleotides_ << output_delimiter
      << exons_ << output_delimiter
      << attributes_;

}




void kgl::GeneCharacteristic::writeGeneHeader(std::ostream& out_file, char output_delimiter) const {

  out_file << "Genome" << output_delimiter
           << "Contig" << output_delimiter
           << "HGNC_id" << output_delimiter
           << "Ensembl_Id" << output_delimiter
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


// Perform the genetic analysis per iteration.
bool kgl::GeneVariants::writeVariantOutput( const std::shared_ptr<const GenomePEDData>& ped_data,
                                            std::ostream& out_file,
                                            char output_delimiter) const {

  const double homozygous_bias{0.1};

  out_file << unique_variants_ << output_delimiter
           << variant_count_ << output_delimiter
           << span_variant_count_ << output_delimiter;

  for (auto const& [population, upper_tail] : upper_tail_) {

    out_file << upper_tail << output_delimiter;

  }

  for (auto const& [population, lower_tail] : lower_tail_) {

    out_file << lower_tail << output_delimiter;

  }

  out_file << all_lof_ << output_delimiter
           << female_lof_ << output_delimiter
           << male_lof_ << output_delimiter
           << hom_lof_ << output_delimiter;

  ethnic_lof_.writeOutput(ped_data, out_file, output_delimiter);

  out_file << all_high_effect_ << output_delimiter
           << hom_high_effect_ << output_delimiter;

  ethnic_high_.writeOutput(ped_data, out_file, output_delimiter);

  out_file << all_moderate_effect_ << output_delimiter
           << hom_moderate_effect_ << output_delimiter;

  ethnic_moderate_.writeOutput(ped_data, out_file, output_delimiter);

  out_file << genome_count_ << output_delimiter
           << genome_variant_ << output_delimiter;

  out_file << heterozygous_ << output_delimiter
           << homozygous_ << output_delimiter;

  double ratio = static_cast<double>(heterozygous_) / (static_cast<double>(homozygous_) + homozygous_bias);
  out_file << ratio;

  return true;

}





void kgl::GeneVariants::writeVariantHeader( const std::shared_ptr<const GenomePEDData>& ped_data,
                                            std::ostream& out_file,
                                            char output_delimiter) const {

  out_file << "UniqueVariants" << output_delimiter
           << "VariantCount" << output_delimiter
           << "SpanVariantCount" << output_delimiter;

  for (auto const& [upper_population, upper_tail] : upper_tail_) {

    out_file << upper_population << output_delimiter;

  }

  for (auto const& [lower_population, lower_tail] : lower_tail_) {

    out_file << lower_population << output_delimiter;

  }

  out_file << "AllLof" << output_delimiter
           << "FemalePhaseLoF" << output_delimiter
           << "MalePhaseLoF" << output_delimiter
           << "HomLoF" << output_delimiter;

  ethnic_lof_.writeHeader(ped_data, out_file, output_delimiter);

  out_file  << "AllHigh" << output_delimiter
            << "HomHigh" << output_delimiter;

  ethnic_high_.writeHeader(ped_data, out_file, output_delimiter);

  out_file  << "AllModerate" << output_delimiter
            << "HomModerate" << output_delimiter;

  ethnic_moderate_.writeHeader(ped_data, out_file, output_delimiter);

  out_file << "GenomeCount" << output_delimiter
           << "GenomeVariant" << output_delimiter
           << "Heterozygous" << output_delimiter
           << "Homozygous" << output_delimiter
           << "Het/Hom";

}



void kgl::GenomeMutation::writeHeader(const std::shared_ptr<const GenomePEDData>& ped_data,
                                      std::ostream& out_file,
                                      char output_delimiter,
                                      const GeneMutation& gene_mutation) {

  gene_mutation.gene_characteristic.writeGeneHeader(out_file, output_delimiter);

  out_file << output_delimiter;

  gene_mutation.ontology.writeOntologyHeader(out_file, output_delimiter);

  out_file << output_delimiter;

  gene_mutation.gene_variants.writeVariantHeader(ped_data, out_file, output_delimiter);

  out_file << output_delimiter;

  gene_mutation.clinvar.writeHeader(ped_data, out_file, output_delimiter);

  out_file << '\n';

}

