//
// Created by kellerberrin on 2/11/17.
//


#include "kgl_variant_evidence.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Generate compound delete variants.
//  Given a list of raw delete variants:
// 1. Find delete SNPs in coding sequences.
// 2. Find mod(size, 3) = 0 contiguous SNP deletions.
// 3. If aligned on a codon boundary delete size/3 Amino Acids.
// 4. Or if not aligned, modify 2 Amino Acids and delete size/3 -1 Amino Acids.
// 5. Assemble these SNP deletions into a compound CodonDelete variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::CompoundDelete::output() const {

  std::stringstream ss;
  ss << genomeOutput() << '\n';
  ss << "Compound Delete >>>>>\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output() << "\n";
  }
  ss << "<<<<< Compound Delete\n";
  return ss.str();

}


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantAnalysis::codonDelete(std::shared_ptr<const GenomeVariant> delete_SNPs,
                                  std::shared_ptr<const ContigCountData> count_data,
                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  std::shared_ptr<kgl::GenomeVariant>
  genome_delete_variants = kgl::GenomeVariant::emptyGenomeVariant("Compound Delete Variants",
                                                                  delete_SNPs->genomeId(),
                                                                  genome_db_ptr);
  std::vector<CompoundVariantMap> contiguous_delete_vec;

  // Aggregate the coding deletions.
  aggregateCodingDeletions(delete_SNPs, genome_db_ptr, contiguous_delete_vec);
  // Check the aggregates for the same Gene(s) membership.
  if (!membershipCodingDeletions(contiguous_delete_vec)) {

    ExecEnv::log().info("codonDelete(), problem verifying the gene structures of contiguous SNP deletes.");

  } else { // All good.

    generateCodonDeletes(genome_db_ptr, count_data, contiguous_delete_vec, genome_delete_variants);

  }

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
                                                 const std::vector<CompoundVariantMap>& contiguous_delete_vec,
                                                 std::shared_ptr<kgl::GenomeVariant> genome_variant_ptr) {

  if (contiguous_delete_vec.empty()) return;

  for (auto variant_map : contiguous_delete_vec) {

    ContigSize_t mod3_size = variant_map.size() % 3;
    std::shared_ptr<const Variant> compound_delete;

    switch(mod3_size) {

      case 0:
        genome_variant_ptr->addVariant(createCompoundDelete(variant_map));
        break;

      case 1:
        break;

      case 2:
        break;

    }

  } // for all genes

}


std::shared_ptr<const kgl::Variant> kgl::VariantAnalysis::createCompoundDelete(const CompoundVariantMap& variant_map) {

  auto it = variant_map.begin();
  GeneVector gene_vector = it->second->geneMembership();

  if (gene_vector.empty()) {

    ExecEnv::log().error("Variant has no associated gene: {}", it->second->output());
    return nullptr;
  }

  StrandSense gene_strand = gene_vector.front()->sequence().strand();
  std::shared_ptr<const ContigFeatures> contig_ptr;
  ContigOffset_t variant_offset;

  switch(gene_strand) {

    case StrandSense::UNKNOWN:
    case StrandSense::FORWARD:
      contig_ptr = variant_map.begin()->second->contig();
      variant_offset = variant_map.begin()->second->contigOffset();
      break;

    case StrandSense::REVERSE:
      contig_ptr = variant_map.rbegin()->second->contig();
      variant_offset = variant_map.rbegin()->second->contigOffset();
      break;
  }

  std::shared_ptr<const Variant> compound_delete(std::make_shared<const CompoundDelete>(contig_ptr,
                                                                                        variant_offset,
                                                                                        variant_map));

  return compound_delete;

}


/*
void kgl::VariantAnalysis::generateCodonDeletes( std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                                 std::shared_ptr<const ContigCountData> count_data,
                                                 const std::vector<CompoundVariantMap>& contiguous_delete_vec) {

  size_t codon_aligned_deletions = 0;
  size_t codon_offset_deletions = 0;
  size_t previous_delete_evidence = 0;
  size_t codon_offset_2_deletions = 0;
  size_t previous_codon_evidence = 0;
  size_t coding_deletions = 0;

  if (contiguous_delete_vec.empty()) return;

  for (auto variant_map : contiguous_delete_vec) {

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

//          ExecEnv::log().info("Variant compound delete contig: {}, offset:{}, size: {} not in gene: {} coding sequence",
//                              it->second->contigId(), delete_offset, variant_map.size(), gene->id());
//          gene->recusivelyPrintsubfeatures();
//          gene->printCDSvector(coding_cds_vec);

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
*/

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
