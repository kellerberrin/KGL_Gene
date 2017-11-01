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
std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantAnalysis::SNPVariants(std::shared_ptr<const ContigCountData> count_data,
                                  std::shared_ptr<const GenomeDatabase> genome_db) {

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
// 3. Then assemble these SNP deletions into a compound CodonDelete variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantAnalysis::codonDelete(std::shared_ptr<const GenomeVariant> delete_SNPs,
                                  std::shared_ptr<const ContigCountData> count_data,
                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  std::shared_ptr<const kgl::GenomeVariant> genome_delete_variants;
  std::vector<CompoundVariantMap> contiguous_delete_vec;

  // Aggregate the coding deletions.
  aggregateCodingDeletions(delete_SNPs, genome_db_ptr, contiguous_delete_vec);
  // Check the aggregates for the same Gene(s) membership.
  if (!membershipCodingDeletions(contiguous_delete_vec)) {

    ExecEnv::log().info("codonDelete(), problem verifying the gene structures of contiguous SNP deletes.");

  } else { // All good.

    generateCodonDeletes(genome_db_ptr, count_data, contiguous_delete_vec);

  }

  printCompoundVariant(contiguous_delete_vec);

  return genome_delete_variants;

}


// Generate SNP variants.
void kgl::VariantAnalysis::aggregateCodingDeletions(std::shared_ptr<const GenomeVariant>& delete_SNPs,
                                                    std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                                    std::vector<CompoundVariantMap>& contiguous_delete_vec) {

  // Filter to only coding deletion SNPs.
  delete_SNPs = delete_SNPs->filterVariants(kgl::DeleteSNPFilter());
  std::shared_ptr<kgl::GenomeVariant> coding_delete_SNPs = delete_SNPs->filterVariants(kgl::InCDSFilter(genome_db_ptr));
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
        if ((contiguous_deletes.rbegin()->second->contigOffset() + 1) == variants.second->contigOffset()) {

          auto result = contiguous_deletes.insert(variants);

          if (not result.second) {

            ExecEnv::log().error("generateCodonDelete(), unexpected - could not insert to contiguous map");

          } // insert OK

        } else { // not contiguous

          if (contiguous_deletes.size() >= 2) {  // 2 or more added to the vector.

            contiguous_delete_vec.push_back(contiguous_deletes);

          }

          contiguous_deletes.clear();

          auto result = contiguous_deletes.insert(variants);

          if (not result.second) {

            ExecEnv::log().error("generateCodonDelete(), unexpected - could not insert to contiguous map");

          } // insert OK

        } // not contig.

      } // not empty.

    } // for all variants.

    if (contiguous_deletes.size() >= 2) {  // 2 or more added to the vector.

      contiguous_delete_vec.push_back(contiguous_deletes);

    }

    contiguous_deletes.clear();

  } // for all contigs.

  ExecEnv::log().info("Found: {} contiguous SNP deletions", contiguous_delete_vec.size());

}

// Not strictly necessary. However, always best to check the contents of complex data structures.
// Checks contiguous SNPs for the same Gene membership.
bool
kgl::VariantAnalysis::membershipCodingDeletions(const std::vector<kgl::CompoundVariantMap>& contiguous_delete_vec) {

  bool result = true;

  if (contiguous_delete_vec.empty()) return true;


  for (auto variant_map : contiguous_delete_vec) {

    if (variant_map.size() < 2) {

      ExecEnv::log().error("membershipCodingDeletions(), contiguous map size: {} should be >= 2", variant_map.size());
      result = false;

    } else {  // correct size.

      auto it = variant_map.begin();
      GeneVector cmp_genes = it->second->geneMembership();

      while (it != variant_map.end()) {

        if (it->second->geneMembership().size() != cmp_genes.size()) {

          ExecEnv::log().error("membershipCodingDeletions(), mismatching gene vector size");
          result = false;

        } else { // correct size

          for (size_t idx = 0; idx < cmp_genes.size(); ++idx) {

            if (not cmp_genes[idx]->isGene() or cmp_genes[idx]->id() != it->second->geneMembership()[idx]->id()) {

              ExecEnv::log().error("membershipCodingDeletions(), mismatching genes id: {} and id: {}",
                                   cmp_genes[idx]->id(), it->second->geneMembership()[idx]->id());
              result = false;

            } // Check gene ids.

          } // for all genes

        } // Correct number of genes.

        ++it;

      } // while variant in variant_map

    } // correct size.

  } // for all variant_maps.

  return result;

}


void kgl::VariantAnalysis::generateCodonDeletes( std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                                 std::shared_ptr<const ContigCountData> count_data,
                                                 const std::vector<CompoundVariantMap>& contiguous_delete_vec) {

  size_t codon_lengthmod3_delete = 0;
  size_t codon_aligned_deletions = 0;
  size_t codon_offset_deletions = 0;
  size_t previous_delete_evidence = 0;
  size_t codon_offset_2_deletions = 0;
  size_t previous_codon_evidence = 0;
  size_t coding_deletions = 0;

  if (contiguous_delete_vec.empty()) return;

  for (auto variant_map : contiguous_delete_vec) {

    if ((variant_map.size() % 3) == 0) {

      ++codon_lengthmod3_delete;

    }

    auto it = variant_map.begin();
    GeneVector gene_vector = it->second->geneMembership();


    for (auto gene : gene_vector) {

      StrandSense gene_strand = gene->sequence().strand();
      SortedCDSVector coding_cds_vec;
      gene->getSortedCDS(coding_cds_vec);
      ContigOffset_t delete_offset;

      switch(gene_strand) {

        case StrandSense::UNKNOWN:
        case StrandSense::FORWARD:
          delete_offset = variant_map.begin()->second->contigOffset();
          break;

        case StrandSense::REVERSE:
          delete_offset = variant_map.rbegin()->second->contigOffset();
          break;
      }

      for (auto coding_cds : coding_cds_vec) {

        ContigOffset_t sequence_offset;
        ContigSize_t sequence_length;
        if (not DNA5Sequence::offsetWithinSequence(coding_cds, delete_offset, sequence_offset, sequence_length)) {

          ExecEnv::log().info("Variant compound delete contig: {}, offset:{}, size: {} not in gene: {} coding sequence",
                               it->second->contigId(), delete_offset, variant_map.size(), gene->id());
          gene->recusivelyPrintsubfeatures();
          gene->printCDSvector(coding_cds_vec);

        } else { // we have the delete sequence offset determine if

          ++coding_deletions;
          if ((sequence_offset % 3) == 0) {

            // create the compound delete.
            ++codon_aligned_deletions;

          }
          if ((variant_map.size() % 3) == 1) {

            ++codon_offset_deletions;

            // Get the gene strand
            ContigSize_t previous_offset;
            // see if the the previous offset read evidence has any deletions
            switch(gene_strand) {

              case StrandSense::UNKNOWN:
              case StrandSense::FORWARD:
                previous_offset = delete_offset - 1;
                break;

              case StrandSense::REVERSE:
                previous_offset = delete_offset + 1;
                break;
            }
            auto contig_block = count_data->findContigBlock(it->second->contigId());
            if (contig_block) {

              const auto &nucleotide_array = contig_block->getNucleotideArray();
              NucleotideReadCount_t delete_count;
              delete_count = nucleotide_array.readCount(previous_offset, NucleotideColumn_DNA5::DELETE_NUCLEOTIDE);
              if (delete_count > 1) {

                ++previous_delete_evidence;

              } else { // no delete evidence.


              }

            } else {

              ExecEnv::log().error("Read count data for contig: {} not found", it->second->contigId());

            }

          } // offset 2 from codon boundary.
          if ((variant_map.size() % 3) == 2) {

            ++codon_offset_2_deletions;
            ContigSize_t previous_offset;
            ContigSize_t codon_boundary_offset;

            switch(gene_strand) {

              case StrandSense::UNKNOWN:
              case StrandSense::FORWARD:
                previous_offset = delete_offset - 1;
                codon_boundary_offset = delete_offset - 2;
                break;

              case StrandSense::REVERSE:
                previous_offset = delete_offset + 1;
                codon_boundary_offset = delete_offset + 2;
                break;
            }

            // see if the previous 2 offset read evidence have deletions
            auto contig_block = count_data->findContigBlock(it->second->contigId());
            if (contig_block) {

              const auto &nucleotide_array = contig_block->getNucleotideArray();
              NucleotideReadCount_t previous_count, codon_count;
              previous_count = nucleotide_array.readCount(previous_offset, NucleotideColumn_DNA5::DELETE_NUCLEOTIDE);
              codon_count = nucleotide_array.readCount(codon_boundary_offset, NucleotideColumn_DNA5::DELETE_NUCLEOTIDE);
              if (previous_count > 1 and codon_count > 1) {

                ++previous_codon_evidence;

              } else { // no delete evidence.


              }

            } else {

              ExecEnv::log().error("Read count data for contig: {} not found", it->second->contigId());

            }

          } // offset 2 from codon boundary.

        } // Check coding sequence offset.

      } // coding sequences cds collection

    }  // all coding sequences per gene.

  } // for all genes

  ExecEnv::log().info("Total contiguous deletions: {}, mod(size,3) == 0 :{}, in coding sequence: {}",
                      contiguous_delete_vec.size(), codon_lengthmod3_delete, coding_deletions);
  ExecEnv::log().info("Contiguous deletions on codon boundary: {}", codon_aligned_deletions);
  ExecEnv::log().info("Contiguous deletions offset+1 from codon boundary: {} and offset-1 delete evidence : {}",
                      codon_offset_deletions, previous_delete_evidence);
  ExecEnv::log().info("Contiguous deletions offset+2 from codon boundary: {} and previous delete evidence : {}",
                      codon_offset_2_deletions, previous_codon_evidence);

}


void kgl::VariantAnalysis::printCompoundVariant(const std::vector<CompoundVariantMap>& map_vec) {

  for (auto map : map_vec) {

    ExecEnv::log().info("***************** Compound Variant size: {} ****************", map.size());

    if (map.empty()) continue;

    GeneVector gene_vec = map.begin()->second->geneMembership();

    if (gene_vec.empty()) continue;

    auto gene = *(gene_vec.begin());

    SortedCDSVector coding_cds_vec;
    gene->getSortedCDS(coding_cds_vec);

    if (coding_cds_vec.empty()) continue;

    StrandSense strand = gene->sequence().strand();
    ContigOffset_t contig_offset;

    switch(strand) {

      case StrandSense::UNKNOWN:
      case StrandSense::FORWARD:
        contig_offset = map.begin()->second->contigOffset();
        break;

      case StrandSense::REVERSE:
        contig_offset = map.rbegin()->second->contigOffset();
        break;

    }

    ContigOffset_t sequence_offset;
    ContigSize_t sequence_length;
    if (DNA5Sequence::offsetWithinSequence(coding_cds_vec.front(), contig_offset, sequence_offset, sequence_length)) {

      ExecEnv::log().info("Sequence length:{} ({}), offset: {} ({}), mod(3): {} strand: {}",
                          sequence_length, static_cast<int>(sequence_length/3), sequence_offset,
                          static_cast<int>(sequence_offset/3), (sequence_offset % 3), static_cast<char>(strand));

    }

    for (auto variant : map) {

      ExecEnv::log().info("Component variant: {}", variant.second->output());

    }

  }

}
