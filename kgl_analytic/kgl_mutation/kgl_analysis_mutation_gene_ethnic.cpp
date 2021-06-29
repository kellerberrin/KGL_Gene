//
// Created by kellerberrin on 31/1/21.
//

#include "kgl_analysis_mutation_gene_ethnic.h"


namespace kgl = kellerberrin::genome;




bool kgl::GeneEthnicitySex::pedAnalysis(const GenomeId_t& genome_id,
                                      size_t count,
                                      const std::shared_ptr<const HsGenomeAux>& genome_aux_data) {

  if (count == 0) {

    return true;

  }

  if (not genome_aux_data) {

    ExecEnv::log().critical("GenomeMutation::pedAnalysis; PED file pointer net defined");
    return false;

  }

  auto record_opt = genome_aux_data->getGenome(genome_id);

  if (not record_opt) {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; Genome sample: {} does not have a PED record", genome_id);
    return false;

  }

  auto const ped_record = record_opt.value();

  if (ped_record.sexType() == AuxSexType::MALE) {

    male_ += count;

  } else {

    female_ += count;

  }

  std::string super_pop = ped_record.superPopulation();
  auto find_result = super_population_.find(super_pop);
  if (find_result != super_population_.end()) {

    auto& [super_pop, super_count] = *find_result;
    super_count += count;

  } else {

    ExecEnv::log().error("GenomeMutation::pedAnalysis; Unable to find super population: {}", super_pop);

  }

  std::string pop = ped_record.population();
  auto pop_result = population_.find(pop);
  if (pop_result != population_.end()) {

    auto& [pop, pop_count] = *pop_result;
    pop_count += count;

  } else {

    ExecEnv::log().error("GenomeMutation::pedAnalysis; Unable to find population: {}", pop);

  }

  return true;

}


void kgl::GeneEthnicitySex::updatePopulations(const std::shared_ptr<const HsGenomeAux>& genome_aux_data) {

  // Generate the template populations.
  std::map<std::string, size_t> population_map;
  for (auto const& [population, description] : genome_aux_data->populationList()) {

    population_map[population] = 0;

  }

  // Generate the template super populations.
  std::map<std::string, size_t> super_population_map;
  for (auto const& [super_population, description] : genome_aux_data->superPopulationList()) {

    super_population_map[super_population] = 0;

  }


  super_population_ = super_population_map;
  population_ = population_map;

}


void kgl::GeneEthnicitySex::writeHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                         std::ostream& out_file,
                                         char output_delimiter) const {

  if (display_flags_ & DISPLAY_SEX_FLAG) {

    writeSexHeader(out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_SUPER_POP_FLAG) {

    writeSuperPopHeader(genome_aux_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_POPULATION_FLAG) {

    writePopHeader(genome_aux_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }

}



void kgl::GeneEthnicitySex::writeSexHeader( std::ostream& out_file, char output_delimiter) const {

  out_file << header_prefix_ << "Male" << output_delimiter
           << header_prefix_ << "Female";

}

void kgl::GeneEthnicitySex::writeSuperPopHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                                 std::ostream& out_file,
                                                 char output_delimiter) const {

  for (auto const& [population, description] : genome_aux_data->superPopulationList()) {

    if (population != genome_aux_data->superPopulationList().begin()->first) {

      out_file << output_delimiter;

    }

    out_file << header_prefix_ << population;


  }

}

void kgl::GeneEthnicitySex::writePopHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                            std::ostream& out_file,
                                            char output_delimiter) const {

  for (auto const& [population, description] : genome_aux_data->populationList()) {

    if (population != genome_aux_data->populationList().begin()->first) {

      out_file << output_delimiter;

    }

    out_file << header_prefix_ << population;


  }

}



void kgl::GeneEthnicitySex::writeOutput(const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                        std::ostream& out_file,
                                        char output_delimiter) const {

  if (display_flags_ & DISPLAY_SEX_FLAG) {

    kgl::GeneEthnicitySex::writeSex(out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_SUPER_POP_FLAG) {

    writeSuperPop(genome_aux_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_POPULATION_FLAG) {

    writePop(genome_aux_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }


}


void kgl::GeneEthnicitySex::writeSex( std::ostream& out_file,
                                      char output_delimiter) const {

  out_file << male_ << output_delimiter
           << female_;

}


void kgl::GeneEthnicitySex::writeSuperPop( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                           std::ostream& out_file,
                                           char output_delimiter) const {

  if (super_population_.size() != genome_aux_data->superPopulationList().size()) {

    ExecEnv::log().error("GeneEthnicitySex::writeSuperPop; Mismatch between data super population size: {}, and Ped size: {}",
                         super_population_.size(), genome_aux_data->superPopulationList().size());

  }

  for (auto const& [super_pop, super_pop_count] : super_population_) {

    if (super_pop != super_population_.begin()->first) {

      out_file << output_delimiter;

    }

    out_file << super_pop_count;

  }

}


void kgl::GeneEthnicitySex::writePop( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                      std::ostream& out_file,
                                      char output_delimiter) const {

  if (population_.size() != genome_aux_data->populationList().size()) {

    ExecEnv::log().error("GeneEthnicitySex::writePop; Mismatch between data super population size: {}, and Ped size: {}",
                         population_.size(), genome_aux_data->populationList().size());

  }

  for (auto const& [pop, pop_count] : population_) {

    if (pop != population_.begin()->first) {

      out_file << output_delimiter;

    }

    out_file << pop_count;

  }

}


size_t kgl::GeneEthnicitySex::superPopulationCount(const std::string& super_population) const {

  auto result = super_population_.find(super_population);
  if (result == super_population_.end()) {

    ExecEnv::log().error("GeneEthnicitySex::superPopulationCount; could not find super population: {}", super_population);
    return 0;

  } else {

    auto const& [pop, count] = *result;
    return count;

  }

}


size_t kgl::GeneEthnicitySex::populationCount(const std::string& population) const {

  auto result = population_.find(population);
  if (result == population_.end()) {

    ExecEnv::log().error("GeneEthnicitySex::populationCount; could not find population: {}", population);
    return 0;

  } else {

    auto const& [pop, count] = *result;
    return count;

  }

}

size_t kgl::GeneEthnicitySex::superPopulationTotal() const {

  size_t total{0};
  for (auto const& [population, count] : super_population_) {

    total += count;

  }

  return total;

}

size_t kgl::GeneEthnicitySex::populationTotal() const {


  size_t total{0};
  for (auto const& [population, count] : population_) {

    total += count;

  }

  return total;

}

bool kgl::GeneEthnicitySex::auditTotals() const {

  if (superPopulationTotal() != populationTotal()) {

    ExecEnv::log().error("GeneEthnicitySex::auditTotals; super population total: {} not equal to population total: {}", superPopulationTotal(), populationTotal());
    return false;

  }

  return true;

}