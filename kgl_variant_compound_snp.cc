//
// Created by kellerberrin on 5/11/17.
//

#include "kgl_filter.h"
#include "kgl_variant_evidence.h"

namespace kgl = kellerberrin::genome;


std::string kgl::CompoundSNP::output() const {

  std::stringstream ss;
  ss << genomeOutput() << " ";
  ss << mutation() << "\n";
  ss << "Compound SNP >>>>>\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output() << "\n";
  }
  ss << "<<<<< Compound SNP\n";
  return ss.str();

}

std::string kgl::CompoundSNP::mutation() const {

  std::stringstream ss;
  ss << "-" << "(" << variant_map_.size() << ")" << contigOffset() << " ";
  return ss.str();

}

// Compound SNPs are created to hold 2 or 3 coding SNPs that modify an Amino Acid simultaneously.
// 1. The SNP can only be in a coding sequence.
// 2. They must be within an offset of 2 of each other.
// 3. They must lie within the same codon boundary (a subset of 2).

std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantAnalysis::compoundSNP(const std::shared_ptr<const GenomeVariant>& SNPs,
                                  const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) {

  std::shared_ptr<GenomeVariant>
  genome_compoundSNP_variants = kgl::GenomeVariant::emptyGenomeVariant("Compound SNP Variants",
                                                                       SNPs->genomeId(),
                                                                       genome_db_ptr);

  // Filter to only coding base mutation SNPs.
  std::shared_ptr<GenomeVariant> mutant_coding_SNPs = SNPs->filterVariants(MutantSNPFilter());
  mutant_coding_SNPs = mutant_coding_SNPs->filterVariants(InCDSFilter());

  CompoundVariantMap variant_map;
  for (auto contig_variant : mutant_coding_SNPs->contigMap()) {

    for (auto variant : contig_variant.second->getMap()) {

      if (variant_map.empty()) {


      }
      else {

        int64_t offset = variant.second->contigOffset();
        int64_t previous_offset = variant_map.rbegin()->second->contigOffset();
        int64_t offset_diff = ::llabs(offset - previous_offset);
        switch(variant_map.size()) {

          case 2:
            if (offset_diff == 1) {

              // check if this variant is in the same codon as the variants in the CompoundVariantMap.

              GeneVector variant_genes = variant.second->geneMembership();
              for (auto gene : variant_genes) {
                std::shared_ptr<const CodingSequenceArray> coding_seqs = kgl::GeneFeature::getCodingSequences(gene);
                ContigOffset_t codon_offset;
                ContigSize_t base_in_codon;
                for (const auto& sequence : coding_seqs->getMap()) {

                  if (kgl::CodingSequenceDNA5::codonOffset(sequence.second,
                                                           variant.second->contigOffset(),
                                                           codon_offset,
                                                           base_in_codon)) {


                  }

                }

              }

            }

            break;

          case 1:
            break;

          case 0:
            variant_map[variant.second->contigOffset()] = variant.second;
            break;

        }

      }

    } // for variant

  } // for contig.

  return genome_compoundSNP_variants;

}



