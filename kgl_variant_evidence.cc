//
// Created by kellerberrin on 29/10/17.
//

#include "kgl_variant_evidence.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


// Generate SNP variants.
std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantAnalysis::SNPVariants(const std::shared_ptr<const ContigCountData>& count_data,
                                  const std::shared_ptr<const GenomeDatabase>& genome_db) {

  std::shared_ptr<GenomeVariant> snp_variants = kgl::GenomeVariant::emptyGenomeVariant("simpleSNPVariants",
                                                                                       count_data->fileName(),
                                                                                       genome_db);
  size_t snp_count = 0;

  for (auto& contig_block : count_data->getMap()) {   // For each contig block.

    // Get the sequence.
    std::shared_ptr<ContigFeatures> contig_ptr;
    if (not genome_db->getContigSequence(contig_block.first, contig_ptr)) {

      ExecEnv::log().error("Contig: {} not found in SNPVariants()", contig_block.first);
      continue;

    } else {

      const DNA5Sequence& contig_sequence = contig_ptr->sequence();

      const auto &nucleotide_array = contig_block.second->getNucleotideArray();
      for (ContigOffset_t contig_offset = 0; contig_offset < nucleotide_array.contigSize(); ++contig_offset) {

        const NucleotideReadCount_t* nucleotide_count_ptr = nucleotide_array.readCount(contig_offset);

        ContigOffset_t max_count_offset = 0;
        NucleotideReadCount_t read_count = 0;
        NucleotideReadCount_t max_count = 0;
        for(ContigOffset_t idx = 0; idx <  NucleotideColumn_DNA5::NUCLEOTIDE_COLUMNS; ++idx) {

          if (max_count < nucleotide_count_ptr[idx]) {

            max_count_offset = idx;
            max_count = nucleotide_count_ptr[idx];

          }

          read_count += nucleotide_count_ptr[idx];

        }

        typename NucleotideColumn_DNA5::NucleotideType max_count_nucleotide
        = NucleotideColumn_DNA5::offsetToNucleotide(max_count_offset);
        typename NucleotideColumn_DNA5::NucleotideType reference_nucleotide = contig_sequence[contig_offset];
        NucleotideReadCount_t reference_count
        =  nucleotide_count_ptr[NucleotideColumn_DNA5::nucleotideToColumn(reference_nucleotide)];

        if (max_count > reference_count
            and max_count_nucleotide != contig_sequence[contig_offset]
            and read_count > 0) {

          std::shared_ptr<const Variant>
          snp_variant(std::make_shared<const SNPVariantDNA5>(contig_ptr,
                                                             contig_offset,
                                                             read_count,
                                                             max_count,
                                                             nucleotide_count_ptr,
                                                             NucleotideColumn_DNA5::NUCLEOTIDE_COLUMNS,
                                                             reference_nucleotide,
                                                             max_count_nucleotide));

          if (not snp_variants->addVariant(snp_variant)) { // Annotate with genome information

            ExecEnv::log().error("Unable to add SNP variant: {} to contig: {} - probable offset duplicate",
                                 snp_variant->output(), contig_ptr->contigId());

          }
          ++snp_count;

        }

      }  // for all sequence elements

      ExecEnv::log().info("Contig: {} has: {} raw SNPs", contig_ptr->contigId(), snp_count);
      snp_count = 0;

    } // found contig.

  }  // for all contigs.

  return snp_variants;

}


