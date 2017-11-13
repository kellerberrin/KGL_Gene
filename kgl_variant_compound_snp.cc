//
// Created by kellerberrin on 5/11/17.
//

#include "kgl_filter.h"
#include "kgl_variant_evidence.h"

namespace kgl = kellerberrin::genome;


std::string kgl::CompoundSNP::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index) << delimiter;
  ss << mutation(delimiter, output_index) << delimiter;
  ss << "Compound SNP >>>>>\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output(delimiter, output_index) << "\n";
  }
  ss << "<<<<< Compound SNP\n";
  return ss.str();

}

std::string kgl::CompoundSNP::mutation(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << "-" << "(" << variant_map_.size() << ")" << offsetOutput(contigOffset(), output_index) << delimiter;
  return ss.str();

}

bool kgl::CompoundSNP::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                            std::shared_ptr<DNA5Sequence>& mutated_sequence) const {

  ExecEnv::log().warn("mutateCodingSequence() not yet implemented for CompoundSNP");
  return false;

}

// Compound SNPs are created to hold 2 or 3 coding SNPs that modify an Amino Acid simultaneously.
// 1. The SNP can only be in the same coding sequence (same mRNA parent).
// 2. They must be within an offset of 2 of each other.
// 3. They must lie within the same codon boundary (a subset of point 2).

std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantAnalysis::compoundSNP(const std::shared_ptr<const GenomeVariant>& SNPs,
                                  const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) {

  std::shared_ptr<GenomeVariant>
  genome_compoundSNP_variants = kgl::GenomeVariant::emptyGenomeVariant(SNPs->genomeId(), genome_db_ptr);

  // Filter to only coding base mutation SNPs.
  std::shared_ptr<GenomeVariant> mutant_coding_SNPs = SNPs->filterVariants(MutantSNPFilter());
  mutant_coding_SNPs = mutant_coding_SNPs->filterVariants(InCDSFilter());

  CompoundVariantMap variant_map;
  for (auto contig_variant : mutant_coding_SNPs->contigMap()) {

    for (auto variant : contig_variant.second->getMap()) {
/*
      int64_t offset = variant.second->contigOffset();
      int64_t previous_offset = variant_map.rbegin()->second->contigOffset();
      int64_t offset_diff = ::llabs(offset - previous_offset);
      switch(variant_map.size()) {

        case 2:
          if (offset_diff == 1) {


          }

          break;

        case 1:
          break;

        case 0:
          break;

      } // switch map size
*/
    } // for variant

  } // for contig.

  return genome_compoundSNP_variants;

}



