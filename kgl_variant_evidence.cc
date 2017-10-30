//
// Created by kellerberrin on 29/10/17.
//

//
// Created by kellerberrin on 13/10/17.
//

#include "kgl_variant_evidence.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


// Generate SNP variants.
std::shared_ptr<kgl::GenomeVariant> kgl::VariantAnalysis::SNPVariants(std::shared_ptr<ContigCountData>& count_data,
                                                                     std::shared_ptr<GenomeDatabase>& genome_db) {

  std::shared_ptr<GenomeVariant> snp_variants(std::make_shared<GenomeVariant>("simpleSNPVariants",
                                                                              count_data->fileName()));

  for (auto& contig_block : count_data->getMap()) {   // For each contig block.

    // Get the sequence.
    std::shared_ptr<ContigFeatures> contig_ptr;
    if (not genome_db->getContigSequence(contig_block.first, contig_ptr)) {

      ExecEnv::log().error("Contig: {} not found in SNPVariants()", contig_block.first);
      continue;

    } else {

      const DNA5Sequence& contig_sequence = contig_ptr->sequence();
      std::shared_ptr<ContigVariant> contig_variant_ptr(std::make_shared<ContigVariant>(contig_ptr->contigId()));

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

          contig_variant_ptr->addVariant(contig_offset, snp_variant);

        }

      }  // for all sequence elements

      ExecEnv::log().info("Contig: {} has: {} raw SNPs",
                          contig_variant_ptr->contigId(),
                          contig_variant_ptr->variantCount());


      if (not snp_variants->addContigVariant(contig_variant_ptr)) {

        ExecEnv::log().info("Duplicate Contig: {} variant. SNP variants not added to genome variants",
                            contig_variant_ptr->contigId());

      }

    } // found contig.

  }  // for all contigs.

  return snp_variants;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Generate compound delete variants.
//  Given a list of raw delete variants:
// 1. find delete SNPs in coding sequences.
// 2. find 3 contiguous SNP deletions aligned on a Gene codon boundary that delete a single Amino Acid
// Then assemble these SNP deletions into a compound CodonDelete variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate SNP variants.
std::shared_ptr<kgl::GenomeVariant>
kgl::VariantAnalysis::codonDelete(std::shared_ptr<kgl::GenomeVariant> delete_SNPs,
                                  std::shared_ptr<ContigCountData>& count_data,
                                  std::shared_ptr<GenomeDatabase>& genome_db_ptr) {

  std::shared_ptr<kgl::GenomeVariant> genome_delete_variants;

  // Filter to only coding deletion SNPs.
  std::shared_ptr<kgl::GenomeVariant> coding_delete_SNPs = delete_SNPs->filterVariants(kgl::InCDSFilter(genome_db_ptr));
  coding_delete_SNPs =coding_delete_SNPs->filterVariants(kgl::DeleteSNPFilter());
  std::vector<CompoundVariantMap> contiguous_delete_vec;
  CompoundVariantMap contiguous_deletes;

  for (const auto& contig_variants : coding_delete_SNPs->contigMap()) {

    contiguous_deletes.clear();

    for (const auto& variants : contig_variants.second->getMap()) {

      if (contiguous_deletes.empty()) {

        auto result = contiguous_deletes.insert(variants);

        if (not result.second) {

          ExecEnv::log().error("generateCodonDelete(), unexpected - could not insert to empty map");

        }

      } else { // not empty

        // check if next delete contiguous.
        if ((contiguous_deletes.rbegin()->first + 1) == variants.first) {

          auto result = contiguous_deletes.insert(variants);

          if (not result.second) {

            ExecEnv::log().error("generateCodonDelete(), unexpected - could not insert to contiguous map");

          } // insert OK

        } else { // not contiguous

          if (contiguous_deletes.size() >= 2) {  // 2 or more added to the vector.

            contiguous_delete_vec.push_back(contiguous_deletes);

          }

          contiguous_deletes.clear();

        } // not contig.

      } // not empty.

    } // for all variants.

    if (contiguous_deletes.size() >= 2) {  // 2 or more added to the vector.

      contiguous_delete_vec.push_back(contiguous_deletes);

    }

    contiguous_deletes.clear();

  } // for all contigs.

  ExecEnv::log().info("Found: {} contiguous SNP deletions", contiguous_delete_vec.size());

  return genome_delete_variants;

}
