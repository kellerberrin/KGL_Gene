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


std::string kgl::CompoundDelete::output(char delimiter) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter) << delimiter;
  ss << mutation() << "\n";
  ss << "Compound Delete >>>>>\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output(delimiter);
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
  // Generate the Codon_deletes
  generateCodonDeletes(genome_db_ptr, count_data, contiguous_delete_vec, genome_delete_variants);

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
          if (delete_variant != nullptr) {

            if (not genome_variant_ptr->addVariant(delete_variant)) {

              ExecEnv::log().error("Unable to add compound delete variant: {} - probable offset duplicate",
                                   delete_variant->output(' '));

            }

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


std::shared_ptr<const kgl::Variant>
kgl::VariantAnalysis::createCompoundDelete(const CompoundVariantMap& variant_map) {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  if (variant_map.begin()->second->geneMembership().empty()) {

    ExecEnv::log().error("Variant has no associated gene");
    return nullptr;

  }

  std::shared_ptr<const GeneFeature> gene_ptr = variant_map.begin()->second->geneMembership().front();

  StrandSense gene_strand = gene_ptr->sequence().strand();
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

