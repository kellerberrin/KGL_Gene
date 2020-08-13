//
// Created by kellerberrin on 3/7/20.
//

#include "kel_thread_pool.h"
#include "kgl_analysis_mutation.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const RuntimeParameterMap& named_parameters,
                                               std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(REFERENCE_GENOME_);

  if (ref_genome_opt) {

    genome_GRCh38_ = ref_genome_opt.value();

  } else {

    ExecEnv::log().error("MutationAnalysis::initializeAnalysis, Could not find Genome: {} Analysis: {} disabled.", REFERENCE_GENOME_, ident());
    return false;

  }

  if (not getParameters(work_directory, named_parameters)) {

    return false;

  }

// Clear the data file a
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out | std::ofstream::trunc);

  return true;

}

// This function superclasses the data objects and stores them for further use.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const DataObjectBase> data_object_ptr) {

  ExecEnv::log().info("Analysis: {}, begin processing data file", ident(), data_object_ptr->Id());


  if (data_object_ptr->dataType() == DataTypeEnum::DiploidPopulation) {

    std::shared_ptr<const DiploidPopulation> diploid_population = std::dynamic_pointer_cast<const DiploidPopulation>(data_object_ptr);

    if (diploid_population) {

      ExecEnv::log().info("Analysis: {}, Generate Hetreozygous/Homozygous ratio statistics for file: {}", ident(), data_object_ptr->Id());

      diploid_population_ = diploid_population;

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a Diploid Population", ident(), data_object_ptr->Id());
      return false;

    }

  }

  if (data_object_ptr->dataType() == DataTypeEnum::UnphasedPopulation) {

    std::shared_ptr<const UnphasedPopulation> unphased_population = std::dynamic_pointer_cast<const UnphasedPopulation>(data_object_ptr);

    if (unphased_population) {

      unphased_population_ = unphased_population;

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not an Unphased Population", ident(), data_object_ptr->Id());
      return false;

    }

  }

  if (data_object_ptr->dataType() == DataTypeEnum::PedAncestor) {

    std::shared_ptr<const GenomePEDData> ped_data = std::dynamic_pointer_cast<const GenomePEDData>(data_object_ptr);

    if (ped_data) {

      ped_data_ = ped_data;
      ExecEnv::log().info("Analysis: {}, ped file: {} contains: {} PED records", ident(), ped_data->Id(), ped_data->getMap().size());

    } else {

      ExecEnv::log().error("MutationAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a PED Ancestor Object", ident(), data_object_ptr->Id());
      return false;

    }

  }

  ExecEnv::log().info("Analysis: {}, completed data file: {}", ident(), data_object_ptr->Id());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  if (diploid_population_ and unphased_population_ and ped_data_) {

    ExecEnv::log().info("Filtered Population : {} Count: {} and Joined Population: {} Count: {} both active",
                        diploid_population_->variantCount(), diploid_population_->populationId(),
                        unphased_population_->populationId(), unphased_population_->variantCount());

    if (not hetHomRatio(diploid_population_)) {

      ExecEnv::log().error("Analysis: {},  problem creating Het/Hom ratio", ident());
      return false;

    }
    joinPopulations();
    checkPED();



  } else {

    ExecEnv::log().info("Failed to create Filtered Population and Joined Populations Id: {}", ident());

  }

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}


bool kgl::MutationAnalysis::getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters) {

  // Get the output filename
  auto result = named_parameters.find(OUTPUT_FILE_);
  if (result == named_parameters.end()) {
    ExecEnv::log().error("Analytic: {}; Expected Parameter: {} to be defined. {} is deactivated. Available named Parameters:", ident(), OUTPUT_FILE_, ident());
    for (auto const& [parameter_ident, parameter_value] : named_parameters) {

      ExecEnv::log().info("Analysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

    }
    return false;
  }
  output_file_name_ = Utility::filePath(result->second, work_directory);

  ExecEnv::log().info("Analysis: {}, initialized with output file: {}", ident(), output_file_name_);

  return true;

}


// Calculate the HetHom Ratio
bool kgl::MutationAnalysis::hetHomRatio(std::shared_ptr<const DiploidPopulation> population) {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38_->getMap()) {

    ThreadPool thread_pool;
    std::vector<std::future<future_ret_tuple>> future_vector;
    std::map<GenomeId_t, std::tuple<size_t, size_t, double, double>> genome_variant_count;
    for (auto const& [genome, genome_ptr] : population->getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        std::future<future_ret_tuple> future = thread_pool.enqueueTask(&MutationAnalysis::processContig, this, genome_contig_id, genome_ptr);
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

          ExecEnv::log().error("MutationAnalysis::hetHomRatio, Genome sample: {} does not have a PED record", genome_id);

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


kgl::MutationAnalysis::future_ret_tuple
kgl::MutationAnalysis::processContig(ContigId_t contig_id, std::shared_ptr<const DiploidGenome> genome_ptr) {

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

        if (result and allele_frequency >= MAX_MAF_) {

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

    ExecEnv::log().info("MutationAnalysis::hetHomRatio, Processed Genome: {}, Contig: {}, All Variants: {}, Filtered Variants: {}, multiple snp: {}",
                        genome_ptr->genomeId(), contig_id, contig_opt.value()->variantCount(), contig_ptr->variantCount(), multiple_snp_count);

    return { genome_ptr->genomeId(), true, heterozygous_count, homozygous_count, expected_heterozygous, expected_homozygous};

  }

  return { genome_ptr->genomeId(), false, 0, 0, 0.0, 0.0}; // contig not present.

}



std::tuple<bool, double> kgl::MutationAnalysis::alleleFrequency_Gnomad(GenomeId_t genome_id, std::shared_ptr<const Variant> variant_ptr) {

  // Lookup the corresponding variant in the Gnomad database
  auto joined_variant_opt = lookupUnphasedVariant(variant_ptr);

  if (not joined_variant_opt) {

    return {false, 0.0};  // variant not found.

  }

  auto result = ped_data_->getMap().find(genome_id);

  if (result == ped_data_->getMap().end()) {

    ExecEnv::log().error("MutationAnalysis::alleleFrequency_Gnomad, Genome sample: {} does not have a PED record", genome_id);
    return {false, 0.0};

  } else {

    auto [sample_id, ped_record] = *result;

    if (ped_record.superPopulation() == SUPER_POP_AFR_GNOMAD_.first) {

      return processFloatField(joined_variant_opt.value(), SUPER_POP_AFR_GNOMAD_.second);

    } else if (ped_record.superPopulation() == SUPER_POP_AMR_GNOMAD_.first) {

      return processFloatField(joined_variant_opt.value(), SUPER_POP_AMR_GNOMAD_.second);

    } else if (ped_record.superPopulation() == SUPER_POP_EAS_GNOMAD_.first) {

      return processFloatField(joined_variant_opt.value(), SUPER_POP_EAS_GNOMAD_.second);

    } else if (ped_record.superPopulation() == SUPER_POP_EUR_GNOMAD_.first) {

      return processFloatField(joined_variant_opt.value(), SUPER_POP_EUR_GNOMAD_.second);

    } else if (ped_record.superPopulation() == SUPER_POP_SAS_GNOMAD_.first) {

      return processFloatField(joined_variant_opt.value(), SUPER_POP_SAS_GNOMAD_.second);

    } else  {

      ExecEnv::log().error("MutationAnalysis::alleleFrequency_Gnomad, Sample Id: {} PED Record Unknown Super Population: {}"
                           , sample_id, ped_record.superPopulation());
      return {false, 0.0};

    }

  }

}



std::tuple<bool, double> kgl::MutationAnalysis::alleleFrequency_SNPdb(GenomeId_t, std::shared_ptr<const Variant> variant_ptr) {

    // Lookup the AF frequency.
  auto joined_variant_opt = lookupUnphasedVariant(variant_ptr);

  if (joined_variant_opt) {

    return processStringField(joined_variant_opt.value(), SNP_DB_FREQ_);

  } else {

    return {false, 0.0};

  }

}


std::tuple<bool, double>  kgl::MutationAnalysis::alleleFrequency_1000Genome(GenomeId_t genome_id, std::shared_ptr<const Variant> variant_ptr) {

  auto result = ped_data_->getMap().find(genome_id);

  if (result == ped_data_->getMap().end()) {

    ExecEnv::log().error("MutationAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);
    return {false, 0.0};

  } else {

    // The Info AF field for this the Genome super population.
    std::string super_population_AF_field = result->second.superPopulation() + GENOME_1000_FREQ_SUFFIX_;
    // Lookup the AF frequency.
    return processFloatField(variant_ptr, super_population_AF_field);

  }

}


std::tuple<bool, double> kgl::MutationAnalysis::processFloatField( const std::shared_ptr<const Variant>& variant_ptr
                                                                  , const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, field_name);

  if (field_opt) {

    std::vector<double> field_vec = InfoEvidenceAnalysis::varianttoFloats(field_opt.value());

    if (field_vec.size() == 1) {

      return {true, field_vec.front() };

    } else if (field_vec.size() == 0) {

      // Missing value
      return {false, 0.0};

    } else {

      std::string vector_str;
      for (auto const& str : field_vec) {

        vector_str += str;
        vector_str += ";";

      }

      ExecEnv::log().error("MutationAnalysis::processField, Field: {} expected vector size 1, get vector size: {}, vector: {}",
                           field_name, field_vec.size(), vector_str);
      return {false, 0.0};

    }

  } else {

    ExecEnv::log().error("MutationAnalysis::processField, Field: {} not found for Variant: {}",
                         field_name, variant_ptr->output(',',VariantOutputIndex::START_0_BASED, false));

    return {false, 0.0};

  }

}



std::tuple<bool, double> kgl::MutationAnalysis::processStringField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, field_name);

  if (field_opt) {

    std::vector<std::string> field_vec = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());

    if (field_vec.size() != 1) {

      ExecEnv::log().error("MutationAnalysis::processStringField, Field: {} expected vector size 1, get vector size: {}",
                           field_name, field_vec.size());
      return {false, 0.0};

    } else {

      try {

        double frequency = std::stod(field_vec.front());
        return {true, frequency };

      }
      catch(std::exception& e) {

        ExecEnv::log().error("MutationAnalysis::processStringField, Field: {} problem converting: {} to double, exception: {}",
                             field_name, field_vec.front(), e.what());
        return {false, 0.0};

      }

    }

  } else {

    ExecEnv::log().error("MutationAnalysis::processField, Field: {} not found for Variant: {}",
                         field_name, variant_ptr->output(',',VariantOutputIndex::START_0_BASED, false));
    return {false, 0.0};
  }

}




std::optional<std::shared_ptr<const kgl::Variant>>
kgl::MutationAnalysis::lookupUnphasedVariant(std::shared_ptr<const Variant> variant_ptr) {

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

    ExecEnv::log().error("utationAnalysis::lookupUnphasedVariant, Joining Population has {} Genomes, Expected 1",
                         unphased_population_->getMap().size());
    return std::nullopt;

  }

}



// Join Diploid and a single genome population such as Gnomad or Clinvar.
void kgl::MutationAnalysis::joinPopulations() const {


  JoinSingleGenome join_pop(diploid_population_, unphased_population_, ped_data_);

  join_pop.joinPopulations();

  ExecEnv::log().info("Joined: {} Genomes in Population: {}", diploid_population_->getMap().size(), diploid_population_->populationId());

}


// Calculate the HetHom Ratio
bool kgl::MutationAnalysis::hetHomRatioLocus(const std::shared_ptr<const DiploidPopulation>& population) const {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38_->getMap()) {

    // Generate the super population allele locus lists
    std::map<std::string, std::shared_ptr<const ContigVariant>> locus_map;
    locus_map[SUPER_POP_AFR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_AFR_GNOMAD_.second, MAX_MAF_);
    locus_map[SUPER_POP_AMR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_AMR_GNOMAD_.second, MAX_MAF_);
    locus_map[SUPER_POP_EAS_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_EAS_GNOMAD_.second, MAX_MAF_);
    locus_map[SUPER_POP_EUR_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_EUR_GNOMAD_.second, MAX_MAF_);
    locus_map[SUPER_POP_SAS_GNOMAD_.first] = getLocusList(genome_contig_id, VARIANT_SPACING_,
                                                          SUPER_POP_SAS_GNOMAD_.second, MAX_MAF_);

    // Use a thread pool to calculate inbreeding and relatedness.
    ThreadPool thread_pool;
    std::vector<std::future<locus_ret_tuple>> future_vector;
    for (auto const&[genome, genome_ptr] : population->getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        auto result = ped_data_->getMap().find(genome);
        if (result == ped_data_->getMap().end()) {

          ExecEnv::log().error("MutationAnalysis::hetHomRatioLocu, Genome sample: {} does not have a PED record", genome);
          continue;

        }
        auto const& [sample_id, ped_record] = *result;
        auto locus_result = locus_map.find(ped_record.superPopulation());
        if (locus_result == locus_map.end()) {

          ExecEnv::log().error("MutationAnalysis::hetHomRatioLocu, Locus set not found for super population: {}", ped_record.superPopulation());
          continue;

        }
        auto const& [super_pop_id, locus_list] = *locus_result;

        std::future<locus_ret_tuple> future = thread_pool.enqueueTask(&MutationAnalysis::processLocusContig, this,
                                                                      genome, contig_opt.value(), super_pop_id, locus_list);
        future_vector.push_back(std::move(future));

      }

    }

    // Retrieve the thread results into a map.
    std::map<GenomeId_t, std::tuple<size_t, size_t, double, double>> genome_results_map;
    for (auto &future : future_vector) {

      auto[genome, hetero_count, homo_count, expected_hetero, expected_homo] = future.get();
      genome_results_map[genome] = std::tuple<size_t, size_t, double, double>(hetero_count, homo_count, expected_hetero,
                                                                              expected_homo);

    }

    writeResults(genome_contig_id, genome_results_map);

  }

  return true;

}

bool kgl::MutationAnalysis::writeResults( const ContigId_t& contig_id,
                                          const std::map<GenomeId_t, std::tuple<size_t, size_t, double, double>>& genome_results_map) const {

  // Append the results.
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out |  std::ofstream::app);

  outfile << contig_id << DELIMITER_ << "Population" << DELIMITER_
          << "Description" << DELIMITER_ << "SuperPopulation" << DELIMITER_ << "Description" << DELIMITER_
          << "HetCount" << DELIMITER_ << "HomCount" << DELIMITER_ << "ObsHet/Hom" << DELIMITER_
          << "ExpectedHet" << DELIMITER_ << "ExpectedHom" << DELIMITER_ << "ExHet/Hom" << DELIMITER_ << "InBreeding\n";

  for (auto const& [genome_id, het_hom_count] : genome_results_map) {

    auto const& [het_count, hom_count, expected_het, expected_hom] = het_hom_count;

    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().error("MutationAnalysis::hetHomRatio, Genome sample: {} does not have a PED record", genome_id);

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

  return true;

}


kgl::MutationAnalysis::locus_ret_tuple
kgl::MutationAnalysis::processLocusContig( const GenomeId_t& genome_id,
                                           const std::shared_ptr<const DiploidContig>& contig_ptr,
                                           const std::string& super_population,
                                           const std::shared_ptr<const ContigVariant>& locus_list) const {


  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {
  // Join on the diploid contig.

    auto diploid_variant_opt = contig_ptr->findOffsetArray(offset);

    if (diploid_variant_opt) {
    // Determine if the sample alternate allele is Hom/Het or Mixed.
      auto const& diploid_offset = diploid_variant_opt.value();
      auto locus_variant_array = offset_ptr->getVariantArray();
      if (diploid_offset.size() == 1) {
      // The sample is alt allele heterozygous
      // Find the matching locus allele
        for (auto const& locus_variant : locus_variant_array) {

          if (diploid_offset[0]->homozygous(*locus_variant)) {
          // Found the matching locus allele.


          }

        }

      } else if (diploid_offset.size() == 2) {

        if (diploid_offset[0]->homozygous(*diploid_offset[1])) {
        // The sample is alt allele homozygous


        } else {
        // The sample has different alt alleles.

        }

      } else {

        ExecEnv::log().error("MutationAnalysis::processLocusContig; Diploid genome: {} has: {} SNPs at offset: {} contig: {}",
                             genome_id, diploid_offset.size(), offset, contig_ptr->contigId());
        continue;
      }

    } else {
    // The sample has a homozygous reference allele at this location.


    }

  }

  return { genome_id, 0, 0, 0.0, 0.0}; // contig not present.

}


// Get a list of hom/het SNPs with a spcified spacing to minimise linkage dis-equilibrium
// and at a specified frequency for the super population. Used as a template for calculating
// the inbreeding coefficient and sample relatedness
std::shared_ptr<const kgl::ContigVariant> kgl::MutationAnalysis::getLocusList( const ContigId_t& contig_id,
                                                                               ContigOffset_t spacing,
                                                                               const std::string& super_pop_freq,
                                                                               double min_frequency) const {

  // Annotate the variant list with the super population frequency identifier
  std::shared_ptr<ContigVariant> locus_list(std::make_shared<ContigVariant>(super_pop_freq));

  if (unphased_population_->getMap().size() == 1) {

    auto [genome_id, genome_ptr] = *(unphased_population_->getMap().begin());

    auto contig_opt = genome_ptr->getContig(contig_id);

    if (contig_opt) {

      // Filter for SNP (is this necessary?).
      auto snp_contig_ptr = contig_opt.value()->filterVariants(SNPFilter());
      // Filter for minimum AF frequency
      snp_contig_ptr = snp_contig_ptr->filterVariants(InfoGEQFloatFilter(super_pop_freq, min_frequency));

      ContigOffset_t previous_offset{0};
      for (auto const& [offset, offset_ptr] : snp_contig_ptr->getMap()) {

        if (offset >= previous_offset + spacing) {

          OffsetVariantArray variant_array = offset_ptr->getVariantArray();

          for (auto const& variant_ptr : variant_array) {

            if (not locus_list->addVariant(variant_ptr)) {

              ExecEnv::log().error("MutationAnalysis::getLocusList, Could not add variant: {}",
                                   variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          } // for variant array

        } // if spacing

      } // for offset

    } // if contig_opt

  } else {

    ExecEnv::log().error("MutationAnalysis::getLocusList, Unphased Population has {} Genomes, Expected 1",
                         unphased_population_->getMap().size());

  }

  return locus_list;

}


bool kgl::JoinSingleGenome::joinPopulations() {

  // The unphased population only has 1 genome. We find the variant by offset.
  if (joining_population_->getMap().size() != 1) {

    ExecEnv::log().error("JoinSingleGenome::lookupJoinedPop, unphased population : {} expected 1 genome, actually contains: {}",
                         joining_population_->populationId(), joining_population_->getMap().size());

    return false;

  }

  auto [joining_genome_id, joining_genome_ptr] = *(joining_population_->getMap().begin());

  ThreadPool thread_pool;
  std::vector<std::future<std::tuple<std::string, size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [joined_genome, joined_genome_ptr] : joined_population_->getMap()) {

    // function, object_ptr, arg1, ..., argn
    std::future<std::tuple<std::string, size_t, size_t>> future = thread_pool.enqueueTask(&JoinSingleGenome::processGenome, this, joined_genome_ptr, joining_genome_ptr);
    future_vector.push_back(std::move(future));

  }

  for (auto& future : future_vector) {

    auto [genome_id, total_variants, joined_variants] = future.get();

    double percent = (static_cast<double>(joined_variants) / static_cast<double>(total_variants)) * 100.0;

    std::string population_desc;
    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().warn("MutationAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);
      population_desc = "";

    } else {

      population_desc = result->second.population();

    }

    ExecEnv::log().info("Population: {}, Genome: {}, Total Variants: {}, Joined Variants: {} ({}%)",
                        population_desc, genome_id, total_variants, joined_variants, percent);

  }

  return true;

}


std::tuple<std::string, size_t, size_t>
kgl::JoinSingleGenome::processGenome(std::shared_ptr<const DiploidGenome> diploid_genome, std::shared_ptr<const GenomeVariant> unphased_genome_ptr) {

  LocalGenomeJoin genome_join(unphased_genome_ptr);

  diploid_genome->processAll(genome_join, &LocalGenomeJoin::lookupJoinedPop);

  return {diploid_genome->genomeId(), genome_join.variantsProcessed(), genome_join.joinedVariantsFound()};

}

// Check that all samples (genomes) have a corresponding PED record.
void kgl::MutationAnalysis::checkPED() const {

  size_t PED_record_count = 0;
  for (auto const& [genome_id, genome_ptr] : diploid_population_->getMap()) {

    auto result = ped_data_->getMap().find(genome_id);

    if (result == ped_data_->getMap().end()) {

      ExecEnv::log().warn("MutationAnalysis::checkPED, Genome sample: {} does not have a PED record", genome_id);

    } else {

      ++PED_record_count;

    }

  }

  ExecEnv::log().info("Genome samples with PED records: {}", PED_record_count);

}


// Joins a single genome population (Gnomad, Clinvar) to another (generally phased Diploid) population.
bool kgl::LocalGenomeJoin::lookupJoinedPop(std::shared_ptr<const Variant> variant_ptr) {

  ++variants_processed_;

  auto contig_opt = unphased_genome_ptr_->getContig(variant_ptr->contigId());

  if (contig_opt) {

    auto variant_opt = contig_opt.value()->findVariant(*variant_ptr);

    if (variant_opt) {

      ++joined_variants_found_;

    }

  }

  return true;

}



