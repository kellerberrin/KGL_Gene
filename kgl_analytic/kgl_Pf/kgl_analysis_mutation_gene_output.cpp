//
// Created by kellerberrin on 27/1/21.
//



#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_mutation.h"

#include <fstream>


namespace kgl = kellerberrin::genome;




// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::writeOutput(const std::string& output_file_name, char output_delimiter) const {

  const double homozygous_bias{0.1};
  const double sex_bias{0.1};

  std::ofstream out_file(output_file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("GenomeMutation::writeOutput; could not open file: {} for output", output_file_name);
    return false;

  } else {

    ExecEnv::log().info("GenomeMutation writing output to file: {}", output_file_name);

  }

  writeHeader(out_file, output_delimiter);

  for (auto const& gene : gene_vector_) {

    out_file << gene.genome << output_delimiter
             << gene.contig << output_delimiter
             << gene.gene_id << output_delimiter
             << gene.gene_name << output_delimiter
             << gene.description << output_delimiter
             << gene.biotype << output_delimiter
             << (gene.valid_protein ? "Valid" : "Invalid") << output_delimiter
             << gene.gaf_id << output_delimiter
             << gene.gene_begin << output_delimiter
             << gene.gene_end << output_delimiter
             << gene.gene_span << output_delimiter
             << gene.strand << output_delimiter
             << gene.sequences << output_delimiter
        << gene.seq_name << output_delimiter
        << gene.nucleotides << output_delimiter
        << gene.exons << output_delimiter
        << gene.attribute_size << output_delimiter
        << gene.unique_variants << output_delimiter
        << gene.variant_count << output_delimiter
        << gene.span_variant_count << output_delimiter
        << gene.female_lof << output_delimiter
        << gene.male_lof << output_delimiter
        << gene.hom_lof << output_delimiter
        << gene.male_clinvar << output_delimiter
        << gene.female_clinvar << output_delimiter
        << gene.hom_clinvar << output_delimiter;

    std::string concat_desc{"\""};
    for (auto const& desc : gene.clinvar_desc) {

      if (desc != *gene.clinvar_desc.begin()) {

        concat_desc += CLINVAR_DESC_CONCAT;

      }

      concat_desc += desc;

    }

    concat_desc += "\"";

    out_file << concat_desc << output_delimiter
             << gene.female_high_effect << output_delimiter
             << gene.male_high_effect << output_delimiter
             << gene.hom_high_effect << output_delimiter
             << (static_cast<double>(gene.female_phase) / static_cast<double>(gene.male_phase + sex_bias)) << output_delimiter
             << gene.female_value << output_delimiter
             << gene.male_value  << output_delimiter
             << gene.genome_count << output_delimiter
             << gene.genome_variant << output_delimiter
             << (static_cast<double>(gene.span_variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter;

    if (gene.nucleotides > 0) {

      out_file << (static_cast<double>(gene.variant_count) / static_cast<double>(gene.nucleotides)) << output_delimiter;

    } else {

      out_file << 0 << output_delimiter;

    }

    out_file << gene.EAS << output_delimiter
             << gene.EUR << output_delimiter
             << gene.SAS << output_delimiter
             << gene.AMR << output_delimiter
             << gene.AFR << output_delimiter;

    out_file << gene.heterozygous << output_delimiter
             << gene.homozygous << output_delimiter;
    double ratio = static_cast<double>(gene.heterozygous) / (static_cast<double>(gene.homozygous) + homozygous_bias);
    out_file << ratio << output_delimiter
             << (gene.indel * 100) << output_delimiter
             << (gene.transition * 100) << output_delimiter
             << (gene.transversion * 100) << '\n';

  } // Gene

  return true;

}


void kgl::GenomeMutation::writeHeader(std::ostream& out_file, char output_delimiter) const {

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
           << "Attributes" << output_delimiter
           << "UniqueVariants" << output_delimiter
           << "VariantCount" << output_delimiter
           << "SpanVariantCount" << output_delimiter
           << "FemaleLoF" << output_delimiter
           << "MaleLoF" << output_delimiter
           << "HomLoF" << output_delimiter
           << "MaleClinvar" << output_delimiter
           << "FemaleClinvar" << output_delimiter
           << "HomClinvar" << output_delimiter
           << "ClinvarDesc" << output_delimiter
           << "FemaleHigh" << output_delimiter
           << "MaleHigh" << output_delimiter
           << "HomHigh" << output_delimiter
           << "F/MPhase" << output_delimiter
           << "Female" << output_delimiter
           << "Male" << output_delimiter
           << "GenomeCount" << output_delimiter
           << "GenomeVariant" << output_delimiter
           << "SpanDensity" << output_delimiter
           << "VariantDensity" << output_delimiter
           << "EAS" << output_delimiter
           << "EUR" << output_delimiter
           << "SAS" << output_delimiter
           << "AMR" << output_delimiter
           << "AFR" << output_delimiter
           << "Heterozygous" << output_delimiter
           << "Homozygous" << output_delimiter
           << "Het/Hom" << output_delimiter
           << "Indel%" << output_delimiter
           << "Transition%" << output_delimiter
           << "Transversion%" << '\n';

}

