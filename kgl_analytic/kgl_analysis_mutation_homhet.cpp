//
// Created by kellerberrin on 14/8/20.
//



#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Calculate the HetHom Ratio
bool kgl::InbreedingAnalysis::hetHomRatio(std::shared_ptr<const DiploidPopulation> population) const {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38_->getMap()) {

    ThreadPool thread_pool;
    std::vector<std::future<future_ret_tuple>> future_vector;
    std::map<GenomeId_t, std::tuple<size_t, size_t, double, double>> genome_variant_count;
    for (auto const& [genome, genome_ptr] : population->getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        std::future<future_ret_tuple> future = thread_pool.enqueueTask(&InbreedingAnalysis::processContig, this, genome_contig_id, genome_ptr);
        future_vector.push_back(std::move(future));

      }

    }

    for (auto& future : future_vector) {

      auto [genome, contig_flag, hetero_count, homo_count, expected_hetero, expected_homo] = future.get();

      if (contig_flag) {

        genome_variant_count[genome] = std::tuple<size_t, size_t, double, double>(hetero_count, homo_count, expected_hetero, expected_homo);

      }

    }

    if (not genome_variant_count.empty()) {

      // Append the results.
      std::ofstream outfile;
      outfile.open(output_file_name_, std::ofstream::out |  std::ofstream::app);

      outfile << genome_contig_ptr->contigId() << DELIMITER_ << "Population" << DELIMITER_
              << "Description" << DELIMITER_ << "SuperPopulation" << DELIMITER_ << "Description" << DELIMITER_
              << "HetCount" << DELIMITER_ << "HomCount" << DELIMITER_ << "ObsHet/Hom" << DELIMITER_
              << "ExpectedHet" << DELIMITER_ << "ExpectedHom" << DELIMITER_ << "ExHet/Hom" << DELIMITER_ << "InBreeding\n";

      for (auto const& [genome_id, het_hom_count] : genome_variant_count) {

        auto const& [het_count, hom_count, expected_het, expected_hom] = het_hom_count;

        auto result = ped_data_->getMap().find(genome_id);

        if (result == ped_data_->getMap().end()) {

          ExecEnv::log().error("InbreedingAnalysis::hetHomRatio, Genome sample: {} does not have a PED record", genome_id);

        }

        auto const& [sample_id, ped_record] = *result;

        outfile << genome_id << DELIMITER_;
        outfile << ped_record.population() << DELIMITER_;
        outfile << ped_record.populationDescription() << DELIMITER_;
        outfile << ped_record.superPopulation() << DELIMITER_;
        outfile << ped_record.superDescription() << DELIMITER_;
        outfile << het_count << DELIMITER_;
        outfile << hom_count << DELIMITER_;
        outfile << static_cast<double>(het_count)/static_cast<double>(hom_count) << DELIMITER_;
        outfile << expected_het << DELIMITER_;
        outfile << expected_hom << DELIMITER_;
        outfile << static_cast<double>(expected_het)/static_cast<double>(expected_hom) << DELIMITER_;
        outfile << 1.0 - (static_cast<double>(expected_het)/static_cast<double>(het_count));
        outfile << '\n';

      }

      outfile.flush();

    }

  }

  return true;

}


kgl::InbreedingAnalysis::future_ret_tuple
kgl::InbreedingAnalysis::processContig(ContigId_t contig_id, std::shared_ptr<const DiploidGenome> genome_ptr) const {

  auto contig_opt = genome_ptr->getContig(contig_id);

  if (contig_opt) {

    // variants are restricted to SNPs
    auto contig_ptr = contig_opt.value()->filterVariants(SNPFilter());


    double expected_homozygous{0};
    double expected_heterozygous{0};
    size_t heterozygous_count{0};
    size_t homozygous_count{0};
    size_t multiple_snp_count{0};
    ContigOffset_t variant_next_offset{0};
    for (auto const&[offset, offset_ptr] : contig_ptr->getMap()) {

      if (offset <= variant_next_offset) {

        // Insufficient spacing, Skip this variant
        continue;

      } else {

        // Accept the variant and set the spacing for the next variant.
        variant_next_offset = offset + VARIANT_SPACING_;

      }

      OffsetVariantArray variants = offset_ptr->getVariantArray();

      if (variants.size() == 2 or  variants.size() == 1) {

        // Get the allele frequency.
        //        auto [result, allele_frequency] = alleleFrequency_1000Genome(genome_ptr->genomeId(), variants[0]);
        auto [result, allele_frequency] = alleleFrequency_Gnomad(genome_ptr->genomeId(), variants[0]);
        //        auto [result, allele_frequency] = alleleFrequency_SNPdb(genome_ptr->genomeId(), variants[0]);

        if (result and allele_frequency >= MIN_MAF) {

          expected_heterozygous += 2.0 * allele_frequency * (1.0 - allele_frequency);
          expected_homozygous += allele_frequency * allele_frequency;

          if (variants.size() == 1) {

            ++heterozygous_count;

          } else if (variants.size() == 2) {

            if (variants[0]->homozygous(*variants[1])) {

              ++homozygous_count;

            }

          }

        }

      } else {

        multiple_snp_count += variants.size();

      }

    }

    ExecEnv::log().info("InbreedingAnalysis::hetHomRatio, Processed Genome: {}, Contig: {}, All Variants: {}, Filtered Variants: {}, multiple snp: {}",
                        genome_ptr->genomeId(), contig_id, contig_opt.value()->variantCount(), contig_ptr->variantCount(), multiple_snp_count);

    return { genome_ptr->genomeId(), true, heterozygous_count, homozygous_count, expected_heterozygous, expected_homozygous};

  }

  return { genome_ptr->genomeId(), false, 0, 0, 0.0, 0.0}; // contig not present.

}



std::tuple<bool, double> kgl::InbreedingAnalysis::alleleFrequency_Gnomad(GenomeId_t genome_id,
                                                                       std::shared_ptr<const Variant> variant_ptr) const {

  // Lookup the corresponding variant in the Gnomad database
  auto joined_variant_opt = lookupUnphasedVariant(variant_ptr);

  if (not joined_variant_opt) {

    return {false, 0.0};  // variant not found.

  }

  auto result = ped_data_->getMap().find(genome_id);

  if (result == ped_data_->getMap().end()) {

    ExecEnv::log().error("InbreedingAnalysis::alleleFrequency_Gnomad, Genome sample: {} does not have a PED record", genome_id);
    return {false, 0.0};

  } else {

    auto [sample_id, ped_record] = *result;

    return processFloatField(joined_variant_opt.value(), lookupSuperPopulationField(ped_record.superPopulation()));

  }

}



std::tuple<bool, double> kgl::InbreedingAnalysis::alleleFrequency_SNPdb(GenomeId_t,
                                                                      std::shared_ptr<const Variant> variant_ptr) const {

  // Lookup the AF frequency.
  auto joined_variant_opt = lookupUnphasedVariant(variant_ptr);

  if (joined_variant_opt) {

    return processStringField(joined_variant_opt.value(), SNP_DB_FREQ_);

  } else {

    return {false, 0.0};

  }

}


std::tuple<bool, double>  kgl::InbreedingAnalysis::alleleFrequency_1000Genome(GenomeId_t genome_id,
                                                                            std::shared_ptr<const Variant> variant_ptr) const {

  auto result = ped_data_->getMap().find(genome_id);

  if (result == ped_data_->getMap().end()) {

    ExecEnv::log().error("InbreedingAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);
    return {false, 0.0};

  } else {

    // The Info AF field for this the Genome super population.
    std::string super_population_AF_field = result->second.superPopulation() + GENOME_1000_FREQ_SUFFIX_;
    // Lookup the AF frequency.
    return processFloatField(variant_ptr, super_population_AF_field);

  }

}



std::optional<std::shared_ptr<const kgl::Variant>>
kgl::InbreedingAnalysis::lookupUnphasedVariant(std::shared_ptr<const Variant> variant_ptr) const {

  if (unphased_population_->getMap().size() == 1) {

    auto [genome_id, genome_ptr] = *(unphased_population_->getMap().begin());

    auto contig_opt = genome_ptr->getContig(variant_ptr->contigId());

    if (contig_opt) {

      auto variant_opt = contig_opt.value()->findVariant(*variant_ptr);

      return variant_opt;

    } else {

      return std::nullopt;

    }


  } else {

    ExecEnv::log().error("InbreedingAnalysis::lookupUnphasedVariant, Joining Population has {} Genomes, Expected 1",
                         unphased_population_->getMap().size());
    return std::nullopt;

  }

}

