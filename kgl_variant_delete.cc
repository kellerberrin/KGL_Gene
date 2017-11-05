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
  ss << genomeOutput() << " ";
  ss << mutation() << "\n";
  ss << "Compound Delete >>>>>\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output() << "\n";
  }
  ss << "<<<<< Compound Delete\n";
  return ss.str();

}

std::string kgl::CompoundDelete::mutation() const {

  std::stringstream ss;
  ss << "-" << "(" << variant_map_.size() << ")" << contigOffset() << " ";
  return ss.str();

}


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantAnalysis::codonDelete(const std::shared_ptr<const GenomeVariant>& delete_SNPs,
                                  const std::shared_ptr<const ContigCountData>& count_data,
                                  const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) {

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
void kgl::VariantAnalysis::aggregateCodingDeletions(const std::shared_ptr<const GenomeVariant>& delete_SNPs,
                                                    const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                                    std::vector<CompoundVariantMap>& contiguous_delete_vec) {

  // Filter to only coding deletion SNPs.
  std::shared_ptr<kgl::GenomeVariant> filtered_delete_SNPs = delete_SNPs->filterVariants(kgl::DeleteSNPFilter());
  std::shared_ptr<kgl::GenomeVariant> coding_delete_SNPs = filtered_delete_SNPs->filterVariants(kgl::InCDSFilter());
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


void kgl::VariantAnalysis::generateCodonDeletes( const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                                 const std::shared_ptr<const ContigCountData>& count_data,
                                                 const std::vector<CompoundVariantMap>& contiguous_delete_vec,
                                                 std::shared_ptr<kgl::GenomeVariant> genome_variant_ptr) {

  if (contiguous_delete_vec.empty()) return;

  for (auto variant_map : contiguous_delete_vec) {

    ContigSize_t mod3_size = variant_map.size() % 3;
    std::shared_ptr<const Variant> compound_delete;

    switch (mod3_size) {

      case 0: {

          std::shared_ptr<const Variant> delete_variant = createCompoundDelete(variant_map);
          if (not genome_variant_ptr->addVariant(delete_variant)) {

            ExecEnv::log().error("Unable to add compound delete variant: {} - probable offset duplicate",
                                 delete_variant->output());

          }

        }
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
